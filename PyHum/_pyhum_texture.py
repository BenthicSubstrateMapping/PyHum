## PyHum (Python program for Humminbird(R) data processing) 
## has been developed at the Grand Canyon Monitoring & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
## Project homepage: <https://github.com/dbuscombe-usgs/PyHum>
##
##This software is in the public domain because it contains materials that originally came from 
##the United States Geological Survey, an agency of the United States Department of Interior. 
##For more information, see the official USGS copyright policy at 
##http://www.usgs.gov/visual-id/credit_usgs.html#copyright
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

#"""
# ____        _   _                         
#|  _ \ _   _| | | |_   _ _ __ ___    _   _ 
#| |_) | | | | |_| | | | | '_ ` _ \  (_) (_)
#|  __/| |_| |  _  | |_| | | | | | |  _   _ 
#|_|    \__, |_| |_|\__,_|_| |_| |_| (_) (_)
#       |___/                               
#
#   __            __                
#  / /____  _  __/ /___  __________ 
# / __/ _ \| |/_/ __/ / / / ___/ _ \
#/ /_/  __/>  </ /_/ /_/ / /  /  __/
#\__/\___/_/|_|\__/\__,_/_/   \___/ 
#                                   
#
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
#+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#|d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
#|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+

#"""

# =========================================================
# ====================== libraries ======================
# =========================================================

# operational
import os, time #sys, getopt, 
from scipy.io import loadmat #, savemat
from joblib import Parallel, delayed, cpu_count
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory
except:
   pass
import warnings
warnings.simplefilter('ignore', RuntimeWarning)
import PyHum.io as io

# numerical
import numpy as np
import PyHum.utils as humutils
from scipy.ndimage.filters import median_filter

import PyHum.cwt as cwt
import PyHum.replace_nans as replace_nans

# plotting
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
try:
   from mpl_toolkits.axes_grid1 import make_axes_locatable
except:
   pass
   
#import dask.array as da
   
# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

import warnings
warnings.filterwarnings("ignore")

#################################################
def texture(humfile, sonpath, win=100, shift=10, doplot=1, density=50, numclasses=4, maxscale=20, notes=4):
          
      '''
      Create a texture lengthscale map using the algorithm detailed by Buscombe et al. (2015)
      This textural lengthscale is not a direct measure of grain size. Rather, it is a statistical 
      representation that integrates over many attributes of bed texture, of which grain size is the most important. 
      The technique is a physically based means to identify regions of texture within a sidescan echogram, 
      and could provide a basis for objective, automated riverbed sediment classification.

      Syntax
      ----------
      [] = PyHum.texture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

      Parameters
      ----------
      humfile : str
       path to the .DAT file
      sonpath : str
       path where the *.SON files are
      win : int, *optional* [Default=100]
       pixel in pixels of the moving window
      shift : int, *optional* [Default=10]
       shift in pixels for moving window operation
      doplot : int, *optional* [Default=1]
       if 1, make plots, otherwise do not make plots
      density : int, *optional* [Default=win/2]
       echogram will be sampled every 'density' pixels
      numclasses : int, *optional* [Default=4]
       number of 'k means' that the texture lengthscale will be segmented into
      maxscale : int, *optional* [Default=20]
       Max scale as inverse fraction of data length for wavelet analysis
      notes : int, *optional* [Default=100]
       notes per octave for wavelet analysis

      Returns
      -------
      sonpath+base+'_data_class.dat': memory-mapped file
        contains the texture lengthscale map

      sonpath+base+'_data_kclass.dat': memory-mapped file
        contains the k-means segmented texture lengthscale map

      References
      ----------
      .. [1] Buscombe, D., Grams, P.E., and Smith, S.M.C., 2015, Automated riverbed sediment
       classification using low-cost sidescan sonar. Journal of Hydraulic Engineering 10.1061/(ASCE)HY.1943-7900.0001079, 06015019.
      '''

      # prompt user to supply file if no input file given
      if not humfile:
         print 'An input file is required!!!!!!'
         Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
         humfile = askopenfilename(filetypes=[("DAT files","*.DAT")]) 

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

      # remove underscores, negatives and spaces from basename
      base = humutils.strip_base(base)   

      meta = loadmat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')))

      ft = 1/loadmat(sonpath+base+'meta.mat')['pix_m']
      #pix_m = np.squeeze(meta['pix_m'])
      #dep_m = np.squeeze(meta['dep_m'])
      dist_m = np.squeeze(meta['dist_m'])

      ### port
      print "processing port side ..."
      # load memory mapped scan ... port
      shape_port = np.squeeze(meta['shape_port'])
      if shape_port!='':

         if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat'))):
            port_fp = io.get_mmap_data(sonpath, base, '_data_port_lar.dat', 'float32', tuple(shape_port))         
         else:
            port_fp = io.get_mmap_data(sonpath, base, '_data_port_la.dat', 'float32', tuple(shape_port))

         port_fp2 = io.get_mmap_data(sonpath, base, '_data_port_l.dat', 'float32', tuple(shape_port))

      ### star
      print "processing starboard side ..."
      # load memory mapped scan ... port
      shape_star = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_star'])
      if shape_star!='':
         if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat'))):
            star_fp = io.get_mmap_data(sonpath, base, '_data_star_lar.dat', 'float32', tuple(shape_star))
         else:
            star_fp = io.get_mmap_data(sonpath, base, '_data_star_la.dat', 'float32', tuple(shape_star))

         star_fp2 = io.get_mmap_data(sonpath, base, '_data_star_l.dat', 'float32', tuple(shape_star))

      if len(shape_star)>2:
         shape = shape_port.copy()
         shape[1] = shape_port[1] + shape_star[1]
      else:
         shape = []
         shape.append(1)
         shape.append(shape_port[0])
         shape.append(shape_port[1])
         shape[1] = shape_port[0] + shape_star[0]

      # create memory mapped file for Sp
      with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'w+') as ff:
         fp = np.memmap(ff, dtype='float32', mode='w+', shape=tuple(shape))

      if len(shape_star)>2:
         #SRT = []
         for p in xrange(len(port_fp)):
            
            Z,ind = humutils.sliding_window_sliced(np.vstack((np.flipud(port_fp[p]), star_fp[p])), density, (win,win),(shift,shift))
            
            Snn = get_srt(Z,ind,maxscale, notes, win) #, density)
 
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

            #extent = shape_port[1]
            #Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
            #yvec = np.linspace(pix_m,extent*pix_m,extent)
            #d = dep_m[shape_port[-1]*p:shape_port[-1]*(p+1)]

            R_fp = io.get_mmap_data(sonpath, base, '_data_range.dat', 'float32', tuple(shape_star))

            R = np.vstack((np.flipud(R_fp[0]),R_fp[0]))
            
            R = R/np.max(R)

            #R[R>0.8] = np.nan

            rn = replace_nans.RN(R.astype('float64'),1000,0.01,2,'localmean')
            R = rn.getdata()
            del rn   

            Sp = (Sp**2) * np.cos(R) / shift**2

            fp[p] = Sp.astype('float32')
            del Sp

         del fp # flush data to file

         class_fp = io.get_mmap_data(sonpath, base, '_data_class.dat', 'float32', tuple(shape))

      else: 

            Z,ind = humutils.sliding_window_sliced(np.vstack((np.flipud(port_fp), star_fp)), density, (win,win),(shift,shift))  
            
            Snn = get_srt(Z,ind,maxscale, notes, win) #, density)
            
            # replace nans using infilling algorithm
            rn = replace_nans.RN(Snn.astype('float64'),1000,0.01,2,'localmean')
            Snn = rn.getdata()
            del rn   

            Ny, Nx = np.shape( np.vstack((np.flipud(port_fp), star_fp)) )

            Snn = median_filter(Snn,(int(Nx/100),int(Ny/100)))
   
            Sp = humutils.im_resize(Snn,Nx,Ny)
            del Snn

            Sp[np.isnan(np.vstack((np.flipud(port_fp), star_fp)))] = np.nan
            Sp[np.isnan(np.vstack((np.flipud(port_fp2), star_fp2)))] = np.nan

            #extent = shape_port[0]
            #Zdist = dist_m
            #yvec = np.linspace(pix_m,extent*pix_m,extent)
            #d = dep_m

            R_fp = io.get_mmap_data(sonpath, base, '_data_range.dat', 'float32', tuple(shape_star))

            R = np.vstack((np.flipud(R_fp),R_fp))
            R = R/np.max(R)
            
            #R[R>0.8] = np.nan

            rn = replace_nans.RN(R.astype('float64'),1000,0.01,2,'localmean')
            R = rn.getdata()
            del rn   

            Sp = (Sp**2) * np.cos(R) / shift**2

            with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'w+') as ff:
               np.save(ff, np.squeeze(Sp).astype('float32'))

            with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'r') as ff:
               class_fp = np.load(ff)

            del Sp

      dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])

      ########################################################
      if doplot==1:

         if len(shape_star)>2:
            for p in xrange(len(star_fp)):
               plot_class(dist_m, shape_port, port_fp[p], star_fp[p], class_fp[p], ft, humfile, sonpath, base, p)
         else:
            plot_class(dist_m, shape_port, port_fp, star_fp, class_fp, ft, humfile, sonpath, base, 0)

         if len(shape_star)>2:
            for p in xrange(len(star_fp)):
               plot_contours(dist_m, shape_port, class_fp[p], ft, humfile, sonpath, base, numclasses, p)
         else:
            plot_contours(dist_m, shape_port, class_fp, ft, humfile, sonpath, base, numclasses, 0)
        

      #######################################################
      # k-means 
      
      if len(shape_star)>2:
         with open(os.path.normpath(os.path.join(sonpath,base+'_data_kclass.dat')), 'w+') as ff:
            fp = np.memmap(ff, dtype='float32', mode='w+', shape=tuple(shape))

         for p in xrange(len(port_fp)):
            wc = get_kclass(class_fp[p].copy(), numclasses)
            fp[p] = wc.astype('float32')
            del wc

         del fp

         kclass_fp = io.get_mmap_data(sonpath, base, '_data_kclass.dat', 'float32', tuple(shape))
            
      else:
         wc = get_kclass(class_fp.copy(), numclasses)

         with open(os.path.normpath(os.path.join(sonpath,base+'_data_kclass.dat')), 'w+') as ff:
            np.save(ff, np.squeeze(wc).astype('float32'))

         del wc
         
         with open(os.path.normpath(os.path.join(sonpath,base+'_data_kclass.dat')), 'r') as ff:
            kclass_fp = np.load(ff)
            
      ########################################################
      if doplot==1:

         if len(shape_star)>2:
            for p in xrange(len(star_fp)):
               plot_kmeans(dist_m, shape_port, port_fp[p], star_fp[p], kclass_fp[p], ft, humfile, sonpath, base, p)
         else:
            plot_kmeans(dist_m, shape_port, port_fp, star_fp, kclass_fp, ft, humfile, sonpath, base, 0)         

      if os.name=='posix': # true if linux/mac
         elapsed = (time.time() - start)
      else: # windows
         elapsed = (time.clock() - start)
      print "Processing took ", elapsed , "seconds to analyse"

      print "Done!"
    
    
# =========================================================
def get_srt(Z,ind,maxscale, notes, win): #, density):    
    try:
       #print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
       print "%s windows to process" % (str(len(Z)))                              
       # do the wavelet clacs and get the stats
       d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win) for k in xrange(len(Z))) #density
    except:
       print "memory error: trying serial"
       d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win) for k in xrange(len(Z))) #density

    srt = np.reshape(d , ( ind[0], ind[1] ) )
    del d

    return srt #srt+srt2
            
# =========================================================
def get_kclass(Sk, numclasses):   
    Sk[np.isnan(Sk)] = 0
    wc, values = humutils.cut_kmeans(Sk,numclasses+1)
    wc[Sk==0] = np.nan
    return wc  
   
# =========================================================
def custom_save(figdirec,root):
    '''
    save with no bounding box
    '''
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=400)

# =========================================================
def parallel_me(x, maxscale, notes, win): #, density):
   dat = cwt.Cwt(x, maxscale, notes, win) #, density)
   return dat.getvar()

# =========================================================
def plot_class(dist_m, shape_port, dat_port, dat_star, dat_class, ft, humfile, sonpath, base, p):

   if len(shape_port)>2:
      Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
      extent = shape_port[1]
   else:
      Zdist = dist_m
      extent = shape_port[0]

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

   try:
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(im, cax=cax)
   except:
      plt.colorbar()

   plt.subplot(2,1,2)
   ax = plt.gca()
   plt.imshow(np.vstack((np.flipud(dat_port), dat_star)),cmap='gray',extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper')
   im = ax.imshow(dat_class, alpha=0.5,extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap='YlOrRd')#, vmin=0.5, vmax=3)
   plt.ylabel('Horizontal distance (m)'); 
   plt.xlabel('Distance along track (m)')
   plt.axis('tight')

   try:
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(im, cax=cax)
   except:
      plt.colorbar()

   if len(shape_port)>2:
      custom_save(sonpath,base+'class'+str(p))
   else:
      custom_save(sonpath,base+'class')
   del fig

# =========================================================
def plot_contours(dist_m, shape_port, dat_class, ft, humfile, sonpath, base, numclasses, p):

   if len(shape_port)>2:
      Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
      extent = shape_port[1]
   else:
      Zdist = dist_m
      extent = shape_port[0]

   levels = np.linspace(np.nanmin(dat_class),np.nanmax(dat_class),numclasses+1)

   fig = plt.figure()
   plt.subplot(2,1,1)
   ax = plt.gca()
   CS = plt.contourf(dat_class.astype('float64'), levels, extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)], cmap='YlOrRd',origin='upper')#vmin=0.5, vmax=3)
   plt.ylabel('Horizontal distance (m)'); plt.xlabel('Distance along track (m)')
   plt.axis('tight')

   try:
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(CS, cax=cax)
   except:
      plt.colorbar()

   custom_save(sonpath,base+'class_contours'+str(p))
   del fig
   del CS, levels

# =========================================================
def plot_kmeans(dist_m, shape_port, dat_port, dat_star, dat_kclass, ft, humfile, sonpath, base, p):

   if len(shape_port)>2:
      Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
      extent = shape_port[1]
   else:
      Zdist = dist_m
      extent = shape_port[0]

   #levels = [0.5,0.75,1.25,1.5,1.75,2,3]

   fig = plt.figure()
   plt.subplot(2,1,1)
   ax = plt.gca()
   plt.imshow(np.vstack((np.flipud(dat_port), dat_star)), cmap='gray',extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper')
   
   CS = plt.contourf(np.flipud(dat_kclass), alpha=0.4, extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)],origin='upper', cmap='YlOrRd')#, levels, vmin=0.5, vmax=3)   
   plt.ylabel('Horizontal distance (m)')
   plt.xlabel('Distance along track (m)')
   plt.axis('tight')

   try:
      divider = make_axes_locatable(ax)
      cax = divider.append_axes("right", size="5%", pad=0.05)
      plt.colorbar(CS, cax=cax)
   except:
      plt.colorbar()

   custom_save(sonpath,base+'class_kmeans'+str(p))
   del fig

# =========================================================
# =========================================================
if __name__ == '__main__':

   texture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

            #tmp = np.vstack((np.flipud(port_fp[p]), star_fp[p]))
            #tmpshape = np.shape(tmp)
            ## create memory mapped file for Sp
            #with open(os.path.normpath(os.path.join(sonpath,base+'_tmp.dat')), 'w+') as ff:
            #   tmpfp = np.memmap(ff, dtype='float32', mode='w+', shape=tmpshape)

            #tmpfp = tmp.astype('float32')
            #del tmp

            #del tmpfp # flush data to file
            #tmpfp = io.get_mmap_data(sonpath, base, '_tmp.dat', 'float32', tmpshape)

            #Z,ind = humutils.sliding_window(tmpfp,(win,win),(shift,shift))        
            #os.remove(os.path.normpath(os.path.join(sonpath,base+'_tmp.dat')))
                                
            #try:
            #   Z,ind = humutils.sliding_window(np.vstack((np.flipud(port_fp[p]), star_fp[p])),(win,win),(shift,shift)) 
            #except:
            #   print "memory error ... using dask array on scan"
            #   dat = da.from_array(np.vstack((np.flipud(port_fp[p]), star_fp[p])), chunks=1000)
            #   Z,ind = humutils.sliding_window(dat,(win,win),(shift,shift))
            #   del dat

         #Sk = class_fp.copy()
         #Sk[np.isnan(Sk)] = 0
         #wc, values = humutils.cut_kmeans(Sk,numclasses+1)
         #wc[Sk==0] = np.nan
         #del Sk
         
            #Sk = class_fp[p].copy()
            #Sk[np.isnan(Sk)] = 0
            #wc, values = humutils.cut_kmeans(Sk,numclasses+1)
            #wc[Sk==0] = np.nan
            #del Sk
            
#            try:
#               print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
#            # do the wavelet clacs and get the stats
#               d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win, density) for k in xrange(len(Z)))
#            except:
#               print "memory error: trying serial"
#               d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win, density) for k in xrange(len(Z)))

#            srt = np.reshape(d , ( ind[0], ind[1] ) )
#            del d

#            try:
#               print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
#            # do the wavelet clacs and get the stats
#               d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win, density) for k in xrange(len(Z)))
#            except:
#               print "memory error: trying serial"
#               d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win, density) for k in xrange(len(Z)))

#            srt2 = np.reshape(d , ( ind[0], ind[1] ) )
#            del d
#            Z = None

#            SRT = srt+srt2
#            del srt, srt2

#            Snn = SRT.copy() 
#            del SRT

#            try:
#               print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
#            # do the wavelet clacs and get the stats
#               d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win, density) for k in xrange(len(Z)))
#            except:
#               print "memory error: trying serial"
#               d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k], maxscale, notes, win, density) for k in xrange(len(Z)))

#            srt = np.reshape(d , ( ind[0], ind[1] ) )
#            del d

#            try:
#               print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
#            # do the wavelet clacs and get the stats
#               d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win, density) for k in xrange(len(Z)))
#            except:
#               print "memory error: trying serial"
#               d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win, density) for k in xrange(len(Z)))

#            srt2 = np.reshape(d , ( ind[0], ind[1] ) )
#            del d
#            Z = None

#            SRT = srt+srt2
#            del srt, srt2

#            Snn = SRT.copy() 
#            del SRT

            #try:
            #   Z,ind = humutils.sliding_window(np.vstack((np.flipud(port_fp[p]), star_fp[p])),(win,win),(shift,shift))              
            #except:
            #   print "memory-mapping failed in sliding window - trying memory intensive version"
            #   Z,ind = humutils.sliding_window_nomm(np.vstack((np.flipud(port_fp[p]), star_fp[p])),(win,win),(shift,shift))

           #try:
            #   #print "%s windows to process with a density of %s" % (str(len(Z)), str(density)) #% (str(len(Z)), str(density))
            #   print "%s windows to process" % (str(len(Z)))               
            ## do the wavelet clacs and get the stats
            #   d = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win) for k in xrange(len(Z))) #density
            #except:
            #   print "memory error: trying serial"
            #   d = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k].T, maxscale, notes, win) for k in xrange(len(Z))) #density

            #srt2 = np.reshape(d , ( ind[0], ind[1] ) )
            #del d
            #Z = None
            
            #try:
            #   Z,ind = humutils.sliding_window(np.vstack((np.flipud(port_fp[p]), star_fp[p])),(win,win),(shift,shift))   
            #except:
            #   print "memory-mapping failed in sliding window - trying memory intensive version"
            #   Z,ind = humutils.sliding_window_nomm(np.vstack((np.flipud(port_fp[p]), star_fp[p])),(win,win),(shift,shift))             
            
