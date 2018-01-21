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
#                        __              __                  
#   _________ ___  _____/ /_  ____ _____/ /___ _      _______
#  / ___/ __ `__ \/ ___/ __ \/ __ `/ __  / __ \ | /| / / ___/
# / /  / / / / / (__  ) / / / /_/ / /_/ / /_/ / |/ |/ (__  ) 
#/_/  /_/ /_/ /_/____/_/ /_/\__,_/\__,_/\____/|__/|__/____/  
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

#operational
from __future__ import print_function
from __future__ import division
from scipy.io import loadmat #savemat, 
import os, time #, sys, getopt
#import shutil
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory
except:
   pass


#numerical
import numpy as np
import PyHum.utils as humutils
#from pyhum_utils import sliding_window, im_resize, cut_kmeans
from joblib import Parallel, delayed #, cpu_count
from scipy.ndimage import binary_dilation, binary_erosion, binary_fill_holes, grey_erosion

#import stdev
from skimage.feature import greycomatrix, greycoprops

#plotting
import matplotlib.pyplot as plt
#import matplotlib.colors as colors

# suppress divide and invalid warnings
#np.seterr(divide='ignore')
#np.seterr(invalid='ignore')
np.seterr(all='ignore')

import warnings
warnings.filterwarnings("ignore")

# =========================================================
# =============== begin program ======================
# ========================================================

#################################################
def rmshadows(humfile, sonpath, win=31, shadowmask=0, doplot=1, dissim=3, correl=0.2, contrast=6, energy=0.15, mn=4):
    '''
    Remove dark shadows in scans caused by shallows, shorelines, and attenuation of acoustics with distance
    Manual or automated processing options available
    Works on the radiometrically corrected outputs of the correct module

    Syntax
    ----------
    [] = PyHum.rmshadows(humfile, sonpath, win, shadowmask, doplot)

    Parameters
    ----------
    humfile : str
       path to the .DAT file
    sonpath : str
       path where the *.SON files are
    win : int, *optional* [Default=100]
       window size (pixels) for the automated shadow removal algorithm
    shadowmask : int, *optional* [Default=0]
       1 = do manual shadow masking, otherwise do automatic shadow masking
    doplot : int, *optional* [Default=1]
       1 = make plots, otherwise do not

    Returns
    -------
    sonpath+base+'_data_star_la.dat': memory-mapped file
        contains the starboard scan with water column removed and 
        radiometrically corrected, and shadows removed

    sonpath+base+'_data_port_la.dat': memory-mapped file
        contains the portside scan with water column removed and
        radiometrically corrected, and shadows removed

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

    if win:
       win = np.asarray(win,int)
       print('Window is %s square pixels' % (str(win)))
       
    if shadowmask:
       shadowmask = np.asarray(shadowmask,int)
       if shadowmask==1:
          print('Shadow masking is manual')
       else: 
          print('Shadow masking is auto')
          
    if doplot:
       doplot = int(doplot)
       if doplot==0:
          print("Plots will not be made")

    if dissim:
       dissim = np.asarray(dissim,int)
       print('Threshold dissimilarity (shadow is <) is %s' % (str(dissim)))

    if correl:
       correl = np.asarray(correl,int)
       print('Threshold correlation (shadow is <) is %s' % (str(correl)))

    if contrast:
       contrast = np.asarray(contrast,int)
       print('Threshold contrast (shadow is <) is %s' % (str(contrast)))

    if energy:
       energy = np.asarray(energy,int)
       print('Threshold energy (shadow is >) is %s' % (str(energy)))

    if mn:
       mn = np.asarray(mn,int)
       print('Threshold mean intensity (shadow is <) is %s' % (str(mn)))

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

    base = humutils.strip_base(base)

    meta = loadmat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')))

    # load memory mapped scans
    shape_port = np.squeeze(meta['shape_port'])
    if shape_port!='':
       #port_fp = np.memmap(sonpath+base+'_data_port_la.dat', dtype='float32', mode='r', shape=tuple(shape_port))
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')), 'r') as ff:
          port_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_port))

    shape_star = np.squeeze(meta['shape_star'])
    if shape_star!='':
       #star_fp = np.memmap(sonpath+base+'_data_star_la.dat', dtype='float32', mode='r', shape=tuple(shape_star))
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')), 'r') as ff:
          star_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_star))

    dist_m = np.squeeze(meta['dist_m'])
    ft = 1/(meta['pix_m'])
    extent = shape_star[1] 

    if shadowmask==1: #manual

       Zt = []
       if len(np.shape(star_fp))>2:
          for p in range(len(star_fp)):
             raw_input("Shore picking "+str(p+1)+" of "+str(len(star_fp))+" (starboard), are you ready? 60 seconds. Press Enter to continue...")
             shoreline_star={}
             fig = plt.figure()
             ax = plt.gca()
             ax.imshow(star_fp[p], cmap = 'gray') #, origin = 'upper') #im = 
             plt.axis('normal'); plt.axis('tight')
             pts1 = plt.ginput(n=300, timeout=75) # it will wait for 200 clicks or 75 seconds
             x1=map(lambda x: x[0],pts1) # map applies the function passed as 
             y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
             shoreline_star = np.interp(np.r_[:np.shape(star_fp[p])[1]],x1,y1)
             plt.close()
             del fig

             star_mg = star_fp[p].copy()
             
             shoreline_star = np.asarray(shoreline_star,'int')
             # shift proportionally depending on where the bed is
             for k in range(np.shape(star_mg)[1]):
                star_mg[shoreline_star[k]:,k] = np.nan

             del shoreline_star

             Zt.append(star_mg)
             
       else:

          raw_input("Shore picking "+str(len(star_fp))+" of "+str(len(star_fp))+" (starboard), are you ready? 60 seconds. Press Enter to continue...")
          shoreline_star={}
          fig = plt.figure()
          ax = plt.gca()
          ax.imshow(star_fp, cmap = 'gray') #, origin = 'upper') #im = 
          plt.axis('normal'); plt.axis('tight')
          pts1 = plt.ginput(n=300, timeout=75) # it will wait for 200 clicks or 75 seconds
          x1=map(lambda x: x[0],pts1) # map applies the function passed as 
          y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
          shoreline_star = np.interp(np.r_[:np.shape(star_fp)[1]],x1,y1)
          plt.close()
          del fig

          star_mg = star_fp.copy()

          shoreline_star = np.asarray(shoreline_star,'int')
          # shift proportionally depending on where the bed is
          for k in range(np.shape(star_mg)[1]):
             star_mg[shoreline_star[k]:,k] = np.nan

          del shoreline_star

          Zt.append(star_mg)

       ## create memory mapped file for Z
       #p = np.memmap(sonpath+base+'_data_star_la.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
       #fp[:] = Zt[:]
       #del fp

       Zt = np.squeeze(Zt)
             
       # create memory mapped file for Zs
       #fp = np.memmap(sonpath+base+'_data_star_lar.dat', dtype='float32', mode='w+', shape=np.shape(Zs))
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       del Zt

       #shutil.move(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat')), os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')))


       Zt = []
       if len(np.shape(star_fp))>2:
          for p in range(len(port_fp)):

             raw_input("Shore picking "+str(p+1)+" of "+str(len(port_fp))+" (port), are you ready? 60 seconds. Press Enter to continue...")
             shoreline_port={}
             fig = plt.figure()
             ax = plt.gca()
             ax.imshow(port_fp[p], cmap = 'gray') #, origin = 'upper') #im = 
             plt.axis('normal'); plt.axis('tight')
             pts1 = plt.ginput(n=300, timeout=75) # it will wait for 200 clicks or 75 seconds
             x1=map(lambda x: x[0],pts1) # map applies the function passed as 
             y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
             shoreline_port = np.interp(np.r_[:np.shape(port_fp[p])[1]],x1,y1)
             plt.close()
             del fig

             port_mg = port_fp[p].copy()

             shoreline_port = np.asarray(shoreline_port,'int')
             # shift proportionally depending on where the bed is
             for k in range(np.shape(port_mg)[1]):
                port_mg[shoreline_port[k]:,k] = np.nan

             del shoreline_port

             Zt.append(port_mg)
          
       else:

          raw_input("Shore picking "+str(len(port_fp))+" of "+str(len(port_fp))+" (port), are you ready? 60 seconds. Press Enter to continue...")
          shoreline_port={}
          fig = plt.figure()
          ax = plt.gca()
          ax.imshow(port_fp, cmap = 'gray') #, origin = 'upper') #im = 
          plt.axis('normal'); plt.axis('tight')
          pts1 = plt.ginput(n=300, timeout=75) # it will wait for 200 clicks or 75 seconds
          x1=map(lambda x: x[0],pts1) # map applies the function passed as 
          y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
          shoreline_port = np.interp(np.r_[:np.shape(port_fp)[1]],x1,y1)
          plt.close()
          del fig

          port_mg = port_fp.copy()

          shoreline_port = np.asarray(shoreline_port,'int')
          # shift proportionally depending on where the bed is
          for k in range(np.shape(port_mg)[1]):
             port_mg[shoreline_port[k]:,k] = np.nan

          del shoreline_port

          Zt.append(port_mg)

       Zt = np.squeeze(Zt)
       ## create memory mapped file for Z
       #fp = np.memmap(sonpath+base+'_data_port_la.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
       #fp[:] = Zt[:]
       #del fp

       # create memory mapped file for Zp
       #fp = np.memmap(sonpath+base+'_data_port_lar.dat', dtype='float32', mode='w+', shape=np.shape(Zp))
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       del Zt    

       #shutil.move(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat')), os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')))

    else: #auto

       Zs = []; Zp = []
       if len(np.shape(star_fp))>2:
          for p in range(len(star_fp)):
             merge = np.vstack((np.flipud(port_fp[p]),star_fp[p]))
             merge = np.asarray(merge, 'float64')

             merge_mask = np.vstack((np.flipud(port_fp[p]),star_fp[p]))

             merge[merge_mask==0] = 0
             del merge_mask

             mask = np.asarray(merge!=0,'int8') # only 8bit precision needed

             merge[np.isnan(merge)] = 0

             #Z,ind = humutils.sliding_window(merge,(win,win),(win/2,win/2))
             Z,ind = humutils.sliding_window(merge,(win,win),(win,win))

             #zmean = np.reshape(zmean, ( ind[0], ind[1] ) )
             Ny, Nx = np.shape(merge)
             #zmean[np.isnan(zmean)] = 0
          
             try: #parallel processing with all available cores     
                w = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k], dissim, correl, contrast, energy, mn) for k in range(len(Z)))
             except: #fall back to serial
                w = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k], dissim, correl, contrast, energy, mn) for k in range(len(Z)))          
          
             zmean = np.reshape(w , ( ind[0], ind[1] ) )
             del w
        
             M = humutils.im_resize(zmean,Nx,Ny)
             M[mask==0] = 0
             del zmean

             bw = M>0.5  
             del M

             # erode and dilate to remove splotches of no data
             bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((3,3))), structure=np.ones((13,13)))             
             #bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((win/4,win/4))), structure=np.ones((win/4,win/4)))
             ##bw2 = binary_erosion(bw,structure=np.ones((win*2,win*2)))
                         
             ## fill holes
             bw2 = binary_fill_holes(bw2, structure=np.ones((win,win))).astype(int)
             merge2 = grey_erosion(merge,structure=np.ones((win,win)))
                
             #del bw
             #bw2 = np.asarray(bw2!=0,'int8') # we only need 8 bit precision

             bw2 = np.asarray(bw!=0,'int8') # we only need 8 bit precision
             del bw

             merge[bw2==1] = 0 #blank out bad data
             merge[merge2==np.min(merge2)] = 0 #blank out bad data
             del merge2
         
             ## do plots of merged scans
             if doplot==1:

                Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]

                fig = plt.figure()
                plt.imshow(merge, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
                plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                plt.axis('normal'); plt.axis('tight')
                custom_save(sonpath,'merge_corrected_rmshadow_scan'+str(p))
                del fig

             Zp.append(np.flipud(merge[:shape_port[1],:]))
             Zs.append(merge[shape_port[1]:,:])
             del merge, bw2

       else:

          merge = np.vstack((np.flipud(port_fp),star_fp))
          merge = np.asarray(merge, 'float64')

          merge_mask = np.vstack((np.flipud(port_fp),star_fp))

          merge[merge_mask==0] = 0
          del merge_mask

          mask = np.asarray(merge!=0,'int8') # only 8bit precision needed

          merge[np.isnan(merge)] = 0

          #Z,ind = humutils.sliding_window(merge,(win,win),(win/2,win/2))
          Z,ind = humutils.sliding_window(merge,(win,win),(win,win))

          #zmean = np.reshape(zmean, ( ind[0], ind[1] ) )
          Ny, Nx = np.shape(merge)
          #zmean[np.isnan(zmean)] = 0
          
          try: #parallel processing with all available cores     
             w = Parallel(n_jobs = -1, verbose=0)(delayed(parallel_me)(Z[k], dissim, correl, contrast, energy, mn) for k in range(len(Z)))
          except: #fall back to serial
             w = Parallel(n_jobs = 1, verbose=0)(delayed(parallel_me)(Z[k], dissim, correl, contrast, energy, mn) for k in range(len(Z)))          
          
          zmean = np.reshape(w , ( ind[0], ind[1] ) )
          del w
        
          M = humutils.im_resize(zmean,Nx,Ny)
          M[mask==0] = 0
          del zmean

          bw = M>0.5  
          del M

          # erode and dilate to remove splotches of no data
          bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((3,3))), structure=np.ones((13,13)))             
          #bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((win/4,win/4))), structure=np.ones((win/4,win/4)))
          ##bw2 = binary_erosion(bw,structure=np.ones((win*2,win*2)))
                         
          ## fill holes
          bw2 = binary_fill_holes(bw2, structure=np.ones((win,win))).astype(int)
          merge2 = grey_erosion(merge,structure=np.ones((win,win)))
                
          #del bw
          #bw2 = np.asarray(bw2!=0,'int8') # we only need 8 bit precision

          bw2 = np.asarray(bw!=0,'int8') # we only need 8 bit precision
          del bw

          merge[bw2==1] = 0 #blank out bad data
          merge[merge2==np.min(merge2)] = 0 #blank out bad data
          del merge2

          # erode and dilate to remove splotches of no data
          #bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((3,3))), structure=np.ones((13,13)))
          #bw2 = binary_dilation(binary_erosion(bw,structure=np.ones((win,win))), structure=np.ones((win*2,win*2)))
          #bw2 = binary_erosion(bw,structure=np.ones((win,win)))
                       
          # fill holes
          #bw2 = binary_fill_holes(bw2, structure=np.ones((3,3))).astype(int)
          #del bw
          #bw2 = np.asarray(bw2!=0,'int8') # we only need 8 bit precision

          #merge[bw2==1] = 0 #blank out bad data
         
          ## do plots of merged scans
          if doplot==1:

             Zdist = dist_m
             fig = plt.figure()
             plt.imshow(merge, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
             plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

             plt.axis('normal'); plt.axis('tight')
             custom_save(sonpath,'merge_corrected_rmshadow_scan'+str(0))
             del fig

          Zp.append(np.flipud(merge[:shape_port[0],:]))
          Zs.append(merge[shape_port[0]:,:])
          del merge, bw2

       Zp = np.squeeze(Zp)
       Zs = np.squeeze(Zs)
       # create memory mapped file for Zp
       #fp = np.memmap(sonpath+base+'_data_port_lar.dat', dtype='float32', mode='w+', shape=np.shape(Zp))
       #with open(sonpath+base+'_data_port_lar.dat', 'w+') as f:
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zp))
       fp[:] = Zp[:]
       del fp
       del Zp    

       #shutil.move(sonpath+base+'_data_port_lar.dat', sonpath+base+'_data_port_la.dat')
       #shutil.move(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat')), os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')))

       # create memory mapped file for Zs
       #fp = np.memmap(sonpath+base+'_data_star_lar.dat', dtype='float32', mode='w+', shape=np.shape(Zs))
       #with open(sonpath+base+'_data_star_lar.dat', 'w+') as f:
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zs))
       fp[:] = Zs[:]
       del fp
       del Zs

       #shutil.move(sonpath+base+'_data_star_lar.dat', sonpath+base+'_data_star_la.dat')
       #shutil.move(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat')), os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')))

    if os.name=='posix': # true if linux/mac
       elapsed = (time.time() - start)
    else: # windows
       elapsed = (time.clock() - start)
    print("Processing took "+str(elapsed)+ "seconds to analyse")

    print("Done!")
    print("===================================================")

# =========================================================
def custom_save(figdirec,root):
    '''
    save with no bounding box
    '''
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=400)

# =========================================================
def parallel_me(Z, dissim, correl, contrast, energy, mn):
    try:
       glcm = greycomatrix(Z, [5], [0], 256, symmetric=True, normed=True)
       if (greycoprops(glcm, 'dissimilarity')[0, 0] < dissim) and (greycoprops(glcm, 'correlation')[0, 0] < correl) and (greycoprops(glcm, 'contrast')[0, 0] < contrast) and (greycoprops(glcm, 'energy')[0, 0] > energy) and (np.mean(Z)<mn):
          return 1
       else:
          return 0
    except:
       return 0
       
# =========================================================
# =========================================================
if __name__ == '__main__':

   rmshadows(humfile, sonpath, win, shadowmask, doplot, dissim, correl, contrast, energy, mn)

