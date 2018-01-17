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
#          __            __                             ___     
#         / /____  _  __/ /___  __________        _____/ (_)____
#        / __/ _ \| |/_/ __/ / / / ___/ _ \______/ ___/ / / ___/
#       / /_/  __/>  </ /_/ /_/ / /  /  __/_____(__  ) / / /__  
#       \__/\___/_/|_|\__/\__,_/_/   \___/     /____/_/_/\___/ 
#                                   
#
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
#+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#|d|a|n|i|e|l|.|b|u|s|c|o|m|b|e|@|n|a|u|.|e|d|u|
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#"""

# =========================================================
# ====================== libraries ======================
# =========================================================

# operational
from __future__ import print_function
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

import PyHum.cwt as cwt

from skimage.segmentation import slic

# plotting
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
try:
   from mpl_toolkits.axes_grid1 import make_axes_locatable
except:
   pass
   
# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

import warnings
warnings.filterwarnings("ignore")


#################################################
def texture_slic(humfile, sonpath, doplot=1, numclasses=4, maxscale=20, notes=4):
          
      '''
      Create a texture lengthscale map using the algorithm detailed by Buscombe et al. (2015)
      This textural lengthscale is not a direct measure of grain size. Rather, it is a statistical 
      representation that integrates over many attributes of bed texture, of which grain size is the most important. 
      The technique is a physically based means to identify regions of texture within a sidescan echogram, 
      and could provide a basis for objective, automated riverbed sediment classification.

      Syntax
      ----------
      [] = PyHum.texture(humfile, sonpath, doplot, numclasses, maxscale, notes)

      Parameters
      ----------
      humfile : str
       path to the .DAT file
      sonpath : str
       path where the *.SON files are
      doplot : int, *optional* [Default=1]
       if 1, make plots, otherwise do not make plots
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
         print('An input file is required!!!!!!')
         Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
         humfile = askopenfilename(filetypes=[("DAT files","*.DAT")]) 

      # prompt user to supply directory if no input sonpath is given
      if not sonpath:
         print('A *.SON directory is required!!!!!!')
         Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
         sonpath = askdirectory() 

      # print given arguments to screen and convert data type where necessary
      if humfile:
         print('Input file is %s' % (humfile))
         
      if sonpath:
         print('Sonar file path is %s' % (sonpath))

      if numclasses:
         numclasses = np.asarray(numclasses,int)
         print('Number of sediment classes: %s' % (str(numclasses)))
         
      if maxscale:
         maxscale = np.asarray(maxscale,int)
         print('Max scale as inverse fraction of data length: %s' % (str(maxscale)))
         
      if notes:
         notes = np.asarray(notes,int)
         print('Notes per octave: %s' % (str(notes)))
         
      if doplot:
         doplot = int(doplot)
         if doplot==0:
            print("Plots will not be made")   
      
      
      print('[Default] Number of processors is %s' % (str(cpu_count())))
                        
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
      print("processing port side ...")
      # load memory mapped scan ... port
      shape_port = np.squeeze(meta['shape_port'])
      if shape_port!='':

         if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat'))):
            port_fp = io.get_mmap_data(sonpath, base, '_data_port_lar.dat', 'float32', tuple(shape_port))         
         else:
            port_fp = io.get_mmap_data(sonpath, base, '_data_port_la.dat', 'float32', tuple(shape_port))
          
         #port_fp2 = io.get_mmap_data(sonpath, base, '_data_port_l.dat', 'float32', tuple(shape_port))

      ### star
      print("processing starboard side ...")
      # load memory mapped scan ... port
      shape_star = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_star'])
      if shape_star!='':
         if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat'))):
            star_fp = io.get_mmap_data(sonpath, base, '_data_star_lar.dat', 'float32', tuple(shape_star))
         else:
            star_fp = io.get_mmap_data(sonpath, base, '_data_star_la.dat', 'float32', tuple(shape_star))

         #star_fp2 = io.get_mmap_data(sonpath, base, '_data_star_l.dat', 'float32', tuple(shape_star))

      if len(shape_star)>2:
         shape = shape_port.copy()
         shape[1] = shape_port[1] + shape_star[1]
      else:
         shape = []
         shape.append(1)
         shape.append(shape_port[0])
         shape.append(shape_port[1])
         shape[1] = shape_port[0] + shape_star[0]

      #work on the entire scan
      #im = humutils.rescale(np.vstack((np.flipud(np.hstack(port_fp)), np.hstack(star_fp))),0,1)
      im = np.vstack((np.flipud(np.hstack(port_fp)), np.hstack(star_fp)))
      im[np.isnan(im)] = 0
      im = humutils.rescale(im,0,1)

      #get SLIC superpixels
      segments_slic = slic(im, n_segments=int(im.shape[0]/10), compactness=.1)

      #pre-allocate texture lengthscale array
      tl = np.zeros(im.shape, dtype = "float64")

      #cycle through each segment and compute tl
      for k in np.unique(segments_slic):
         mask = np.zeros(im.shape[:2], dtype = "uint8")
         mask[segments_slic == k] = 255
         cmask, cim = crop_toseg(mask, im)
         tl[segments_slic == k] = parallel_me(cim, maxscale, notes, np.shape(cim)[0])

      R_fp = io.get_mmap_data(sonpath, base, '_data_range.dat', 'float32', tuple(shape_star))
      R = np.vstack((np.flipud(np.hstack(R_fp)), np.hstack(R_fp)))
      R = R/np.max(R)

      #correct for range and scale
      tl = tl * np.cos(R) * (1/ft)
      tl[im==0] = np.nan 
      tl[np.isnan(im)] = np.nan 

      # create memory mapped file for Sp
      with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'w+') as ff:
         fp = np.memmap(ff, dtype='float32', mode='w+', shape=tuple(shape))

      counter = 0
      if len(shape_star)>2:
         for p in range(len(port_fp)):
            if p==0:
               n,m = np.shape(np.vstack((np.flipud(port_fp[p]), star_fp[p])))
            else:
               n,m = np.shape(np.vstack((np.flipud(port_fp[p]), star_fp[p])))
            Sp = tl[:n, counter:counter+m]
            counter = counter+m
            fp[p] = Sp.astype('float32')
            del Sp
         del fp # flush data to file

         class_fp = io.get_mmap_data(sonpath, base, '_data_class.dat', 'float32', tuple(shape))

      else:

            with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'w+') as ff:
               np.save(ff, np.squeeze(Sp).astype('float32'))

            with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'r') as ff:
               class_fp = np.load(ff)

      dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])

      ########################################################
      if doplot==1:

         if len(shape_star)>2:
            for p in range(len(star_fp)):
               plot_class(dist_m, shape_port, port_fp[p], star_fp[p], class_fp[p], ft, humfile, sonpath, base, p)
         else:
            plot_class(dist_m, shape_port, port_fp, star_fp, class_fp, ft, humfile, sonpath, base, 0)

         if len(shape_star)>2:
            for p in range(len(star_fp)):
               plot_contours(dist_m, shape_port, class_fp[p], ft, humfile, sonpath, base, numclasses, p)
         else:
            plot_contours(dist_m, shape_port, class_fp, ft, humfile, sonpath, base, numclasses, 0)
        

      #######################################################
      # k-means 
      
      if len(shape_star)>2:
         with open(os.path.normpath(os.path.join(sonpath,base+'_data_kclass.dat')), 'w+') as ff:
            fp = np.memmap(ff, dtype='float32', mode='w+', shape=tuple(shape))

         for p in range(len(port_fp)):
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
            for p in range(len(star_fp)):
               plot_kmeans(dist_m, shape_port, port_fp[p], star_fp[p], kclass_fp[p], ft, humfile, sonpath, base, p)
         else:
            plot_kmeans(dist_m, shape_port, port_fp, star_fp, kclass_fp, ft, humfile, sonpath, base, 0)         

      if os.name=='posix': # true if linux/mac
         elapsed = (time.time() - start)
      else: # windows
         elapsed = (time.clock() - start)
      print("Processing took "+str(elapsed)+"seconds to analyse")

      print("Done!")
    

# =========================================================
def crop_toseg(mask, im):
   true_points = np.argwhere(mask)
   top_left = true_points.min(axis=0)
   # take the largest points and use them as the bottom right of your crop
   bottom_right = true_points.max(axis=0)
   return mask[top_left[0]:bottom_right[0]+1, top_left[1]:bottom_right[1]+1], im[top_left[0]:bottom_right[0]+1, top_left[1]:bottom_right[1]+1]  

            
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
def parallel_me(x, maxscale, notes, win): 
   dat = cwt.Cwt(x, maxscale, notes, win) 
   return dat.getvar()

# =========================================================
def plot_class(dist_m, shape_port, dat_port, dat_star, dat_class, ft, humfile, sonpath, base, p):

   if len(shape_port)>2:
      Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
      extent = shape_port[1]
   else:
      Zdist = dist_m
      extent = shape_port[0]

   print("Plotting ... ")
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
      custom_save(sonpath,base+'slic_class'+str(p))
   else:
      custom_save(sonpath,base+'slic_class')
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

   custom_save(sonpath,base+'class_slic_contours'+str(p))
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

   custom_save(sonpath,base+'class_slic_kmeans'+str(p))
   del fig

# =========================================================
# =========================================================
if __name__ == '__main__':

   texture_slic(humfile, sonpath, doplot, numclasses, maxscale, notes)




