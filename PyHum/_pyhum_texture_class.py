'''
pyhum_texture.py
Part of PyHum software 

INFO:

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.1      Revision: Sept, 2014

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 series instruments. 

SYNTAX:
python pyhum_texture.py -i datfile -s sonpath
where datfile is the .DAT file associated with the survey, and sonpath is the (absolute or relative) path to where the associated .SON files are

Optional arguments:

EXAMPLES:
1) show help
python pyhum_texture.py -h

2) run the provided test case with all defaults
python pyhum_texture.py -i ./test.DAT -s ./test_data/ (linux)
python pyhum_texture.py -i test.DAT -s \test_data\ (windows)

OUTPUTS:
Files are created which contain the raw and parsed meta data. They are prefixed by the root of the input file (*) followed by:

These are python/matlab/octave .mat data format. To read use, for example:
data = loadmat('test.mat')

If doplot =1 (see above) the program will also create some rudimentary plots of the data (mainly to check everything is ok). These are stored in the same directory as the .son files and are hopefully self explanatory

Installation:

PYTHON LIBRARIES YOU MAY NEED TO INSTALL TO USE PyHum:
1) Pyproj: http://code.google.com/p/pyproj/
2) SciPy: http://www.scipy.org/scipylib/download.html
3) Numpy: http://www.scipy.org/scipylib/download.html
4) Matplotlib: http://matplotlib.org/downloads.html
5) Scikit-learn: http://scikit-learn.org/stable/
6) Python Image LIbrary (PIL) http://www.pythonware.com/products/pil/

All of the above are available through pip (https://pypi.python.org/pypi/pip) and easy_install (https://pythonhosted.org/setuptools/easy_install.html)

OTHER LIBRARIES (CYTHON) NEED TO BE COMPILED FOR SPEED:
1) pyread.pyx
2) cwt.pyx
3) replace_nans.pyx
4) ppdrc.pyx
5) spec_noise.pyx

'''

# =========================================================
# ====================== libraries ======================
# =========================================================

# operational
import sys, getopt, os
from time import clock, time
from scipy.io import loadmat, savemat
from PIL.Image import open as imopen
#from joblib import Parallel, delayed, cpu_count
from sklearn.externals.joblib import Parallel, delayed, cpu_count
from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory

import warnings
warnings.simplefilter('ignore', RuntimeWarning)

# numerical
import numpy as np
import cwt
import ppdrc
import spec_noise
import replace_nans
import PyHum.utils as humutils
from scipy.ndimage import binary_dilation, binary_erosion, binary_fill_holes
from scipy.ndimage.filters import median_filter
from fractions import gcd
   
# plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

__all__ = [
    'humtexture',
    'custom_save',
    'parallel_me',
    ]

#################################################
def humtexture(humfile, sonpath, c, t, f, win, shift, doplot, density, numclasses, maxscale, notes, shorepick, do_two):
                                  
      # prompt user to supply file if no input file given
      if not humfile:
         print 'An input file is required!!!!!!'
         Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
         inputfile = askopenfilename(filetypes=[("DAT files","*.DAT")]) 

      # prompt user to supply directory if no input sonpath is given
      if not sonpath:
         print 'A *.SON directory is required!!!!!!'
         Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
         sonpath = askdirectory() 

      # print given arguments to screen and convert data type where necessary
      if humfile:
         print 'Input file is %s' % (humfile)
      if sonpath:
         print 'Sonar file path is %s' % (sonpath)
      if c:
         c = np.asarray(c,float)
         print 'Sound speed is %s m/s' % (str(c))
      if t:
         t = np.asarray(t,float)
         print 'Transducer length is %s m' % (str(t))
      if f:
         f = np.asarray(f,int)
         print 'Frequency is %s kHz' % (str(f))
      if win:
         win = np.asarray(win,int)
         print 'Window is %s square pixels' % (str(win))
      if shift:
         shift = np.asarray(shift,int)
         print 'Min shift is %s pixels' % (str(shift))
      if density:
         density = np.asarray(density,int)
         print 'Image will be sampled every %s pixels' % (str(density))
      if numclasses:
         numclasses = np.asarray(numclasses,int)
         print 'Number of sediment classes: %s' % (str(numclasses))
      if maxscale:
         maxscale = np.asarray(maxscale,int)
         print 'Max scale as inverse fraction of data length: %s' % (str(maxscale))
      if notes:
         notes = np.asarray(notes,int)
         print 'Notes per octave: %s' % (str(notes))
      if doplot:
         doplot = int(doplot)
         if doplot==0:
            print "Plots will not be made"
      if shorepick:
         shorepick = np.asarray(shorepick,int)
         if shorepick==1:
            print 'Shore picking is manual'
         else: 
            print 'Shore picking is auto'
      if do_two:
         do_two = np.asarray(do_two,int)
         if do_two==1:
            print 'Two sweeps of data (takes twice as long)'
         else: 
            print 'One sweep of data'      
      
      
      print '[Default] Number of processors is %s' % (str(cpu_count()))

      if not c:
         c = 1450
         print '[Default] Sound speed is %s m/s' % (str(c))

      if not t:
         t = 0.108
         print '[Default] Transducer length is %s m' % (str(t))

      if not f:
         f = 455
         print '[Default] Frequency is %s kHz' % (str(f))

      if not win:
         win = 100
         print '[Default] Window is %s square pixels' % (str(win))

      if not shift:
         shift = 2
         print '[Default] Min shift is %s pixels' % (str(shift))

      if not density:
         density = win/2
         print '[Default] Image will be sampled every %s pixels' % (str(density))

      if not numclasses:
         numclasses = 4
         print '[Default] Number of sediment classes: %s' % (str(numclasses))

      if not maxscale:
         maxscale = 20
         print '[Default] Max scale as inverse fraction of data length: %s ' % (str(maxscale))

      if not notes:
         notes = 4
         print '[Default] Notes per octave: %s ' % (str(notes))

      if not doplot:
         if doplot != 0:
            doplot = 1
            print "[Default] Plots will be made"

      if not shorepick:
         shorepick = 0
         print '[Default] Shore picking is auto'

      if not do_two:
         do_two = 0
         print '[Default] One sweep of data'


      ########################################################
      ########################################################

      # start timer
      if os.name=='posix': # true if linux/mac or cygwin on windows
         start1 = time()
      else: # windows
         start1 = clock()

      theta3dB = np.arcsin(c/(t*(f*1000))) # *(180/pi) # to see in degs
      ft = (np.pi/2)*(1/theta3dB)

      # if son path name supplied has no separator at end, put one on
      if sonpath[-1]!=os.sep:
         sonpath = sonpath + os.sep

      base = humfile.split('.DAT') # get base of file name for output
      base = base[0].split('/')[-1]

      # convert to float16 from float64 to conserve memory
      try:
         port_mg = np.asarray(loadmat(sonpath+base+'port_la.mat')['port_mg_la'],'float16')
      except:
         sonpath = os.getcwd()+os.sep
         port_mg = np.asarray(loadmat(sonpath+base+'port_la.mat')['port_mg_la'],'float16')
         
      star_mg = np.asarray(loadmat(sonpath+base+'star_la.mat')['star_mg_la'],'float16')

      if shorepick==1:
         raw_input("Shore picking (starboard), are you ready? 30 seconds. Press Enter to continue...")
         shoreline_star={}
         fig = plt.figure()
         ax = plt.gca()
         im = ax.imshow(star_mg, cmap = 'gray', origin = 'upper')
         plt.axis('normal'); plt.axis('tight')
         pts1 = plt.ginput(n=300, timeout=30) # it will wait for 200 clicks or 30 seconds
         x1=map(lambda x: x[0],pts1) # map applies the function passed as 
         y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
         shoreline_star = np.interp(np.r_[:np.shape(star_mg)[1]],x1,y1)
         plt.close()
         del fig

         raw_input("Shore picking (port), are you ready? 30 seconds. Press Enter to continue...")
         shoreline_port={}
         fig = plt.figure()
         ax = plt.gca()
         im = ax.imshow(port_mg, cmap = 'gray', origin = 'upper')
         plt.axis('normal'); plt.axis('tight')
         pts1 = plt.ginput(n=300, timeout=30) # it will wait for 200 clicks or 30 seconds
         x1=map(lambda x: x[0],pts1) # map applies the function passed as 
         y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
         shoreline_port = np.interp(np.r_[:np.shape(port_mg)[1]],x1,y1)
         plt.close()
         del fig

         shoreline_port = np.asarray(shoreline_port,'int')
         shoreline_star = np.asarray(shoreline_star,'int')
         # shift proportionally depending on where the bed is
         for k in xrange(np.shape(star_mg)[1]):
            star_mg[shoreline_star[k]:,k] = np.nan
            port_mg[shoreline_port[k]:,k] = np.nan

         del shoreline_port, shoreline_star

      merge = np.vstack((np.flipud(port_mg),star_mg))
      del port_mg, star_mg
      merge = np.asarray(merge, 'float64')

      # convert to float16 from float64 to conserve memory
      port_mg = loadmat(sonpath+base+'port_la.mat')['port_mg_la']
      star_mg = loadmat(sonpath+base+'star_la.mat')['star_mg_la']
      merge_mask = np.vstack((np.flipud(port_mg),star_mg))
      del port_mg, star_mg
      #merge_mask = np.asarray(merge_mask, 'float64')

      merge[merge_mask==0] = 0
      del merge_mask

      print "Pre-processing data ... "
      # do some pre-processing to try to remove acoustic shadows
      # first, make a mask of merged scans
      mask = np.asarray(merge!=0,'int8') # only 8bit precision needed

      if shorepick==0:
         # get a 5-value segmentation   
         wc, values = humutils.cut_kmeans(merge,5)
   
         # the lowest value is zero, get the next lowest
         bw = wc==np.sort(values)[-4]
         del wc, values

         # erode and dilate to remove splotches of no data
         bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((3,3))), structure=np.ones((13,13)))
   
         # fill holes
         bw2 = binary_fill_holes(bw2, structure=np.ones((3,3))).astype(int)
         del bw
         bw2 = np.asarray(bw2!=0,'int8') # we only need 8 bit precision

         merge[bw2==1] = 0 #blank out bad data
      else:
         bw2 = np.asarray(np.isnan(merge),'int8')

      merge[np.isnan(merge)] = 0
      Ny, Nx = np.shape(merge)

      print "padding with spectral noise ..."
      # get some noise to plug the gaps for processing
      sn = spec_noise.Noise(humutils.im_resize(merge,int(Nx/10),int(Ny/10)),1)
      sn = humutils.im_resize(sn.getres(),Nx,Ny)
   
      # plug gaps that are not real data
      sn[((bw2==1) - (mask==1)) ==1] = 0
   
      sn = np.asarray(sn,'float16')
      #sn[mask==1] = 0

      merge2 = sn+merge
      #merge2[merge2<=0] = 100*np.random.rand()
   
      del sn, merge

      dat = ppdrc.ppdrc(merge2, 768)
      merge3 = dat.getdata()
      del dat, merge2

      Ny, Nx = np.shape(merge3)

      # end of pre-processing
   
      # get optimal number of slices
      if Nx%2==0:
         H = []
         for k in xrange(3,30):
            H.append(gcd(Nx,Nx/k))
         hslice = np.max(H)
         print "processing into %s slices" % str(Nx/hslice) 
      else:
         merge3 = np.hstack( (merge3,np.ones((Ny,1))) )
         Ny, Nx = np.shape(merge3)
         H = []
         for k in xrange(3,30):
            H.append(gcd(Nx,Nx/k))
         hslice = np.max(H)
         print "processing into %s slices" % str(Nx/hslice)      
   
      # get windowed data
      Zt,ind = humutils.sliding_window(merge3,(Ny,hslice))
    
      del merge3
      numsegs = len(Zt)

      # start timer
      if os.name=='posix': # true if linux/mac or cygwin on windows
         start = time()
      else: # windows
         start = clock()

      SRT = []
      for kk in xrange(numsegs):

         print "Chunk %s out of %s" % (str(kk+1), str(numsegs))

         # get windowed data
         try:
            print "Sliding window with shift= %s" % (str(shift))
            if numsegs==1:
               Z,ind = humutils.sliding_window(np.asarray(Zt,'int8'),(win,win),(shift,shift)) #float16
            else:
               Z,ind = humutils.sliding_window(np.asarray(Zt[kk],'int8'),(win,win),(shift,shift)) #float16
         except:
            print "Sliding window failed with shift= %s, memory error" % (str(shift))
         
            shift = shift+1       
            Z = None
            while Z is None:         
               try:
                  print "trying shift= %s" % (str(shift))
                  if numsegs==1:
                     Z,ind = humutils.sliding_window(np.asarray(Zt,'int8'),(win,win),(shift,shift)) #float16
                  else:
                     Z,ind = humutils.sliding_window(np.asarray(Zt[kk],'int8'),(win,win),(shift,shift)) #float16
               except:
                  shift = shift+1    
             
         filename = sonpath+'tmpfile'+humutils.id_generator()+'.dat'    
         print "Memory mapping the data"
         # create memory mapped file for Z
         fp = np.memmap(filename, dtype='int8', mode='w+', shape=np.shape(Z))
         fp[:] = Z[:]
         del fp
         shapeZ = np.shape(Z)
         del Z
         #we are only going to access the portion of memory required
         newfp = np.memmap(filename, dtype='int8', mode='r', shape=shapeZ)
	  
         try:
            print "Carrying out wavelet calculations ... round 1 ..."
            print "%s windows to process with a density of %s" % (str(len(newfp)), str(density)) #% (str(len(Z)), str(density))
            # do the wavelet clacs and get the stats
            d = Parallel(n_jobs = -1, verbose=1)(delayed(parallel_me)(newfp[k], maxscale, notes, win, density) for k in xrange(len(newfp)))
         except:
            print "memory error: trying a lower pre-dispatch"
            print "Carrying out wavelet calculations ... port ... takes a long time"
            d = Parallel(n_jobs = 1, verbose=1)(delayed(parallel_me)(newfp[k], maxscale, notes, win, density) for k in xrange(len(newfp)))

         # reshape arrays
         #mnsz = reshape(mnsz , ( ind[0], ind[1] ) )
         srt = np.reshape(d , ( ind[0], ind[1] ) )
         del d

         if do_two==1: # analyze transpose as well
            try:
               print "Carrying out wavelet calculations ... round 2 ..."
               print "%s windows to process with a density of %s" % (str(len(newfp)), str(density))
               # do the wavelet clacs and get the stats
               d = Parallel(n_jobs = -1, verbose=1)(delayed(parallel_me)(newfp[k].T, maxscale, notes, win, density) for k in xrange(len(newfp)))
            except:
               print "memory error: trying a lower pre-dispatch"
               print "Carrying out wavelet calculations ... port ... takes a long time"
               d = Parallel(n_jobs = 1, verbose=1)(delayed(parallel_me)(newfp[k].T, maxscale, notes, win, density) for k in xrange(len(newfp)))

            srt2 = np.reshape(d , ( ind[0], ind[1] ) )
            del d
            #del newfp #Z

            SRT.append(srt+srt2)
            del srt, srt2
         
         else: # just one
      
            SRT.append(srt)
            del srt #, newfp
      
         os.remove(filename)

      #del newfp
      if os.name=='posix': # true if linux/mac
         elapsed = (time() - start)
      else: # windows
         elapsed = (clock() - start)
      print "Processing took ", elapsed/60 , "minutes to analyse"

      SRT = np.hstack(SRT)

      Snn = shift*(SRT.copy()*1/(ft**0.5))
      del SRT
      Snn = np.asarray(Snn,'float64')

      # replace nans using infilling algorithm
      rn = replace_nans.RN(Snn,1000,0.01,2,'localmean')
      Snn = rn.getdata()
      del rn   

      Snn = median_filter(Snn,(int(Nx/100),int(Ny/100)))
   
      S = humutils.im_resize(Snn,Nx,Ny)
      del Snn

      S[mask==0] = np.nan
      S[bw2==1] = np.nan
   
      dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])

      star_mg_la = np.asarray(loadmat(sonpath+base+'star_la.mat')['star_mg_la'],'float64')
      port_mg_la = np.asarray(loadmat(sonpath+base+'port_la.mat')['port_mg_la'],'float64')

      merge = np.vstack((np.flipud(port_mg_la),star_mg_la))

      extent = np.shape(merge)[0]   
      del star_mg_la, port_mg_la
   
      merge[mask==0] = np.nan
      merge[bw2==1] = np.nan   
      del bw2, mask   
      merge = np.asarray(merge,'float32')

      ########################################################
      ########################################################
      if doplot==1:
         print "Plotting ... "
         # create fig 1
         fig = plt.figure()
         fig.subplots_adjust(wspace = 0, hspace=0.075)
         plt.subplot(2,1,1)
         ax = plt.gca()
         im = ax.imshow(-merge,cmap='gray',extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)],origin='upper')
         plt.ylabel('Horizontal distance (m)'); 
         plt.axis('tight')

         plt.tick_params(\
          axis='x',          # changes apply to the x-axis
          which='both',      # both major and minor ticks are affected
          bottom='off',      # ticks along the bottom edge are off
          labelbottom='off') # labels along the bottom edge are off

         if humfile=='test.DAT':
            plt.ylim(-20, 20)

         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.05)
         plt.colorbar(im, cax=cax)

         plt.subplot(2,1,2)
         ax = plt.gca()
         plt.imshow(merge,cmap='gray',extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)],origin='upper')
         im = ax.imshow(S, alpha=0.5,extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap='hot')
         plt.ylabel('Horizontal distance (m)'); 
         plt.xlabel('Distance along track (m)')
         plt.axis('tight')

         if humfile=='test.DAT':
            plt.ylim(-20, 20)

         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.05)
         plt.colorbar(im, cax=cax)

         custom_save(sonpath,'class')
         del fig

         ########################################################
         levels = np.linspace(np.nanmin(S),np.nanmax(S),numclasses+1)
         X, Y = np.meshgrid(dist_m, np.linspace(-extent*(1/ft),extent*(1/ft),extent))

         fig = plt.figure()
         plt.subplot(2,1,1)
         ax = plt.gca()
         CS = plt.contourf(X, Y, S, levels, extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)], cmap='hot',origin='upper')
         plt.ylabel('Horizontal distance (m)'); plt.xlabel('Distance along track (m)')
         plt.axis('tight')

         if humfile=='test.DAT':
            plt.ylim(-20, 20)
            plt.gca().invert_yaxis()

         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.05)
         plt.colorbar(CS, cax=cax)

         custom_save(sonpath,'class_contours')
         del fig
         del X, Y, CS, levels

      #######################################################
      # k-means 
      S4k = S.copy()
      S4k[np.isnan(S4k)] = 0
      wc, values = humutils.cut_kmeans(S4k,numclasses+1)
      wc[S4k==0] = np.nan
      del S4k

      if doplot==1:
         if humfile=='test.DAT':
            cols = ['w','y',[.5,.5,.5],'r'] # for size fractions
            gscmap = colors.ListedColormap(cols)

         fig = plt.figure()
         plt.subplot(2,1,1)
         ax = plt.gca()
         plt.imshow(-merge,cmap='gray',extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)],origin='upper')
         if humfile=='test.DAT':
            im = ax.imshow(wc, alpha=0.4, extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap=gscmap)
         else:
            im = ax.imshow(wc, alpha=0.4, extent=[min(dist_m), max(dist_m), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap='hot')   
         plt.ylabel('Horizontal distance (m)')
         plt.xlabel('Distance along track (m)')
         plt.axis('tight')

         if humfile=='test.DAT':
            plt.ylim(-20, 20)

         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.05)
         plt.colorbar(im, cax=cax)

         custom_save(sonpath,'class_kmeans')
         del fig

         ########################################################
         ########################################################
         # save away data
         savemat(sonpath+base+'class.mat', mdict={'numclasses': numclasses, 'S': S, 'kclass': wc, 'ft': ft, 'dist_m': dist_m, 'merge': merge },oned_as='row')

         if os.name=='posix': # true if linux/mac
            elapsed = (time() - start1)
         else: # windows
            elapsed = (clock() - start1)
         print "Processing took ", elapsed/60 , "minutes to analyse"

         print "Done!"

# =========================================================
def custom_save(figdirec,root):
    try:
       plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
    except:
       plt.savefig(os.getcwd()+os.sep+root,bbox_inches='tight',dpi=400)      

# =========================================================
def parallel_me(x, maxscale, notes, win, density):
   dat = cwt.Cwt(x, maxscale, notes, win, density)
   return dat.getvar()

