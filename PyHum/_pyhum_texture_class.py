'''
pyhum_texture.py
Part of PyHum software 

INFO:

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.2.3      Revision: Apr, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 and 1198 series instruments. 

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
import sys, getopt, os, time
from scipy.io import loadmat, savemat
from joblib import Parallel, delayed, cpu_count
from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory

import warnings
warnings.simplefilter('ignore', RuntimeWarning)

# numerical
import numpy as np
import cwt
import replace_nans
import PyHum.utils as humutils
from scipy.ndimage.filters import median_filter

# plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

__all__ = [
    'humtexture',
    'custom_save',
    'parallel_me',
    ]

#################################################
def humtexture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes):
                                  
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
      
      
      print '[Default] Number of processors is %s' % (str(cpu_count()))

      if not win:
         win = 100
         print '[Default] Window is %s square pixels' % (str(win))

      if not shift:
         shift = 10
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

      ########################################################
      ########################################################
      
      # start timer
      if os.name=='posix': # true if linux/mac or cygwin on windows
         start = time.time()
      else: # windows
         start = time.clock()

      # if son path name supplied has no separator at end, put one on
      if sonpath[-1]!=os.sep:
         sonpath = sonpath + os.sep

      base = humfile.split('.DAT') # get base of file name for output
      base = base[0].split(os.sep)[-1]

      ft = 1/loadmat(sonpath+base+'meta.mat')['pix_m']
      pix_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['pix_m'])
      dep_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dep_m'])
      dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])

      ### port
      print "processing port side ..."
      # load memory mapped scan ... port
      shape_port = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_port'])
      if shape_port!='':
         port_fp = np.memmap(sonpath+base+'_data_port_la.dat', dtype='float32', mode='r', shape=tuple(shape_port))
         port_fp2 = np.memmap(sonpath+base+'_data_port_l.dat', dtype='float32', mode='r', shape=tuple(shape_port))

      ### star
      print "processing starboard side ..."
      # load memory mapped scan ... port
      shape_star = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_star'])
      if shape_star!='':
         star_fp = np.memmap(sonpath+base+'_data_star_la.dat', dtype='float32', mode='r', shape=tuple(shape_star))
         star_fp2 = np.memmap(sonpath+base+'_data_star_l.dat', dtype='float32', mode='r', shape=tuple(shape_star))

      shape = shape_port.copy()
      shape[1] = shape_port[1] + shape_star[1]

      # create memory mapped file for Sp
      fp = np.memmap(sonpath+base+'_data_class.dat', dtype='float32', mode='w+', shape=tuple(shape))

      #SRT = []
      for p in xrange(len(port_fp)):

         Z,ind = humutils.sliding_window(np.vstack((np.flipud(port_fp[p]), star_fp[p])),(win,win),(shift,shift))

         try:
            print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
            # do the wavelet clacs and get the stats
            d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win, density) for k in xrange(len(Z)))
         except:
            print "memory error: trying serial"
            d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win, density) for k in xrange(len(Z)))

         srt = np.reshape(d , ( ind[0], ind[1] ) )
         del d

         try:
            print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
            # do the wavelet clacs and get the stats
            d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win, density) for k in xrange(len(Z)))
         except:
            print "memory error: trying serial"
            d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win, density) for k in xrange(len(Z)))

         srt2 = np.reshape(d , ( ind[0], ind[1] ) )
         del d
         Z = None

         SRT = srt+srt2
         del srt, srt2

         Snn = SRT.copy() 
         del SRT

         # replace nans using infilling algorithm
         rn = replace_nans.RN(Snn.astype('float64'),1000,0.01,2,'localmean')
         Snn = rn.getdata()
         del rn   

         Ny, Nx = np.shape( np.vstack((np.flipud(port_fp[p]), star_fp[p])) )
         Snn = median_filter(Snn,(int(Nx/100),int(Ny/100)))
   
         Sp = humutils.im_resize(Snn,Nx,Ny)
         del Snn

         Sp[np.isnan(np.vstack((np.flipud(port_fp[p]), star_fp[p])))] = np.nan
         Sp[np.isnan(np.vstack((np.flipud(port_fp2[p]), star_fp2[p])))] = np.nan

         extent = shape_port[1]
         Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
         yvec = np.linspace(pix_m,extent*pix_m,extent)
         d = dep_m[shape_port[-1]*p:shape_port[-1]*(p+1)]

         R_fp = np.memmap(sonpath+base+'_data_range.dat', dtype='float32', mode='r', shape=tuple(shape_star))

         #R = np.ones(np.shape(Sp))
         #for k in range(len(d)): 
         #   R[:,k] = np.hstack((np.flipud(d[k]/yvec), d[k]/yvec))

         #if len(d)<np.shape(port_fp[p])[1]:
         #   d = np.append(d,d[-1])
         #Zbed = np.squeeze(d*ft)

         #R1 = R[extent:,:]
         #R2 = np.flipud(R[:extent,:])

         ## shift proportionally depending on where the bed is
         #for k in xrange(np.shape(R1)[1]):
         #   R1[:,k] = np.r_[R1[Zbed[k]:,k], np.zeros( (np.shape(R1)[0] -  np.shape(R1[Zbed[k]:,k])[0] ,) )]

         #for k in xrange(np.shape(R2)[1]):
         #   R2[:,k] = np.r_[R2[Zbed[k]:,k], np.zeros( (np.shape(R2)[0] -  np.shape(R2[Zbed[k]:,k])[0] ,) )]

         #R = np.vstack((np.flipud(R2),R1))
         #del R1, R2

         R = np.vstack((np.flipud(R_fp[0]),R_fp[0]))
         
         R[R>0.8] = np.nan

         rn = replace_nans.RN(R.astype('float64'),1000,0.01,2,'localmean')
         R = rn.getdata()
         del rn   

         Sp = (Sp**2) * np.cos(R) / shift**2

         fp[p] = Sp.astype('float32')
         del Sp

      del fp # flush data to file

      class_fp = np.memmap(sonpath+base+'_data_class.dat', dtype='float32', mode='r', shape=tuple(shape))

      dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])

      ########################################################
      ########################################################
      if doplot==1:

         for p in xrange(len(star_fp)):
            plot_class(dist_m, shape_port, port_fp[p], star_fp[p], class_fp[p], ft, humfile, sonpath, base, p)

         for p in xrange(len(star_fp)):
            plot_contours(dist_m, shape_port, class_fp[p], ft, humfile, sonpath, base, numclasses, p)


      #######################################################
      # k-means 
      fp = np.memmap(sonpath+base+'_data_kclass.dat', dtype='float32', mode='w+', shape=tuple(shape))

      for p in xrange(len(port_fp)):
         Sk = class_fp[p].copy()
         Sk[np.isnan(Sk)] = 0
         wc, values = humutils.cut_kmeans(Sk,numclasses+1)
         wc[Sk==0] = np.nan
         del Sk
         fp[p] = wc.astype('float32')
         del wc

      del fp

      kclass_fp = np.memmap(sonpath+base+'_data_kclass.dat', dtype='float32', mode='r', shape=tuple(shape))

      ########################################################
      if doplot==1:

         for p in xrange(len(star_fp)):
            plot_kmeans(dist_m, shape_port, port_fp[p], star_fp[p], kclass_fp[p], ft, humfile, sonpath, base, p)

      if os.name=='posix': # true if linux/mac
         elapsed = (time.time() - start)
      else: # windows
         elapsed = (time.clock() - start)
      print "Processing took ", elapsed , "seconds to analyse"

      print "Done!"
      
# =========================================================
def custom_save(figdirec,root):
    '''
    save with no bounding box
    '''
    plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)

# =========================================================
def parallel_me(x, maxscale, notes, win, density):
   dat = cwt.Cwt(x, maxscale, notes, win, density)
   return dat.getvar()

# =========================================================
def plot_class(dist_m, shape_port, dat_port, dat_star, dat_class, ft, humfile, sonpath, base, p):

   Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
   extent = shape_port[1]

   print "Plotting ... "
   # create fig 1
   fig = plt.figure()
   fig.subplots_adjust(wspace = 0, hspace=0.075)
   plt.subplot(2,1,1)
   ax = plt.gca()
   im = ax.imshow(np.vstack((np.flipud(dat_port),dat_star)) ,cmap='gray',extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper')
   plt.ylabel('Horizontal distance (m)'); 
   plt.axis('tight')

   plt.tick_params(\
   axis='x',          # changes apply to the x-axis
   which='both',      # both major and minor ticks are affected
   bottom='off',      # ticks along the bottom edge are off
   labelbottom='off') # labels along the bottom edge are off

   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   plt.colorbar(im, cax=cax)

   plt.subplot(2,1,2)
   ax = plt.gca()
   plt.imshow(np.vstack((np.flipud(dat_port), dat_star)),cmap='gray',extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper')
   im = ax.imshow(dat_class, alpha=0.5,extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap='YlOrRd', vmin=0.5, vmax=3)
   plt.ylabel('Horizontal distance (m)'); 
   plt.xlabel('Distance along track (m)')
   plt.axis('tight')

   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   plt.colorbar(im, cax=cax)

   custom_save(sonpath,base+'class'+str(p))
   del fig

# =========================================================
def plot_contours(dist_m, shape_port, dat_class, ft, humfile, sonpath, base, numclasses, p):

   Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
   extent = shape_port[1] 

   levels = np.linspace(np.nanmin(dat_class),np.nanmax(dat_class),numclasses+1)

   fig = plt.figure()
   plt.subplot(2,1,1)
   ax = plt.gca()
   CS = plt.contourf(dat_class.astype('float64'), levels, extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)], cmap='YlOrRd',origin='upper', vmin=0.5, vmax=3)
   plt.ylabel('Horizontal distance (m)'); plt.xlabel('Distance along track (m)')
   plt.axis('tight')

   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   plt.colorbar(CS, cax=cax)

   custom_save(sonpath,base+'class_contours'+str(p))
   del fig
   del CS, levels

# =========================================================
def plot_kmeans(dist_m, shape_port, dat_port, dat_star, dat_kclass, ft, humfile, sonpath, base, p):

   Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
   extent = shape_port[1] 

   levels = [0.5,0.75,1.25,1.5,1.75,2,3]

   fig = plt.figure()
   plt.subplot(2,1,1)
   ax = plt.gca()
   plt.imshow(np.vstack((np.flipud(dat_port), dat_star)), cmap='gray',extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper')
   
   CS = plt.contourf(np.flipud(dat_kclass), levels, alpha=0.4, extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap='YlOrRd', vmin=0.5, vmax=3)   
   plt.ylabel('Horizontal distance (m)')
   plt.xlabel('Distance along track (m)')
   plt.axis('tight')

   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   plt.colorbar(CS, cax=cax)

   custom_save(sonpath,base+'class_kmeans'+str(p))
   del fig


