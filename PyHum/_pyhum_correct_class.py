'''
pyhum_correct.py
Part of PyHum software 

INFO:


Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.1.7      Revision: Mar, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 and 1198 series instruments. 

SYNTAX:
python pyhum_correct.py -i datfile -s sonpath
where datfile is the .DAT file associated with the survey, and sonpath is the (absolute or relative) path to where the associated .SON files are

Optional arguments:


EXAMPLES:
1) show help
python pyhum_correct.py -h

2) run the provided test case with all defaults
python pyhum_correct.py -i ./test.DAT -s ./test_data/ (linux)
python pyhum_correct.py -i test.DAT -s \test_data\ (windows)


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
- use the shell script "compile_pyhum.sh" on linux/mac

'''

# =========================================================
# ====================== libraries ======================
# =========================================================

#operational
from __future__ import division
from scipy.io import savemat, loadmat
import os, time, sys, getopt
from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory
from joblib import Parallel, delayed, cpu_count

#numerical
import numpy as np
#from pyhum_utils import rm_spikes, sliding_window, runningMeanFast, dpboundary, rescale
import PyHum.utils as humutils
from scipy.stats import nanmean, nanmedian
import ppdrc

#plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

# =========================================================
# =============== begin program ======================
# ========================================================

__all__ = [
    'humcorrect',
    'custom_save',
    'remove_water',
    'correct_scans',
    ]

#################################################
def humcorrect(humfile, sonpath, maxW, doplot):

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
   if maxW:
      maxW = np.asarray(maxW,float)
      print 'Max. transducer power is %s W' % (str(maxW))
   if doplot:
      doplot = int(doplot)
      if doplot==0:
         print "Plots will not be made"

   if not maxW:
      maxW = 1000
      print '[Default] Max. transducr power is %s W' % (str(maxW))

   if not doplot:
      if doplot != 0:
         doplot = 1
         print "[Default] Plots will be made"

   ##show_stages = 1

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

   # add wattage to metadata dict 
   meta = loadmat(sonpath+base+'meta.mat')
   meta['maxW'] = maxW
   savemat(sonpath+base+'meta.mat', meta ,oned_as='row')
   del meta

   bed = np.squeeze(loadmat(sonpath+base+'meta.mat')['bed'])
   ft = 1/loadmat(sonpath+base+'meta.mat')['pix_m'] #np.squeeze(loadmat(sonpath+base+'meta.mat')['ft'])
   dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])

   # load memory mapped scans
   shape_port = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_port'])
   if shape_port!='':
      port_fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='r', shape=tuple(shape_port))

   shape_star = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_star'])
   if shape_star!='':
      star_fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='r', shape=tuple(shape_star))

   extent = shape_star[1] #np.shape(data_port)[0]

   bed = np.asarray(bed,'int')+int(0.25*ft)

   # calculate in dB
   ######### star
   Zt = remove_water(star_fp, bed, shape_star, maxW)

   # create memory mapped file for Z
   fp = np.memmap(sonpath+base+'_data_star_l.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
   fp[:] = Zt[:]
   del fp
   shape_star = np.shape(Zt)
   del Zt
   #we are only going to access the portion of memory required
   star_fp = np.memmap(sonpath+base+'_data_star_l.dat', dtype='float32', mode='r', shape=shape_star)

   Zt = correct_scans(star_fp)
   for p in xrange(len(Zt)):

      dat1 = Zt[p].astype('float64')
      dat1[np.isnan(dat1)] = 0
      dat1 = ppdrc.ppdrc(dat1, shape_star[-1]/4)
      dat1 = humutils.rescale(dat1.getdata(),0,255)
      dat1[np.isnan(Zt[p])] = np.nan
      Zt[p] = dat1.astype('float32')
      del dat1

   # create memory mapped file for Z
   fp = np.memmap(sonpath+base+'_data_star_la.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
   fp[:] = Zt[:]
   del fp
   shape_star = np.shape(Zt)
   del Zt
   #we are only going to access the portion of memory required
   star_fp = np.memmap(sonpath+base+'_data_star_la.dat', dtype='float32', mode='r', shape=shape_star)

   ######### port
   Zt = remove_water(port_fp, bed, shape_port, maxW)

   # create memory mapped file for Z
   fp = np.memmap(sonpath+base+'_data_port_l.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
   fp[:] = Zt[:]
   del fp
   shape_port = np.shape(Zt)
   del Zt
   #we are only going to access the portion of memory required
   port_fp = np.memmap(sonpath+base+'_data_port_l.dat', dtype='float32', mode='r', shape=shape_port)

   Zt = correct_scans(port_fp)
   for p in xrange(len(Zt)):

      dat1 = Zt[p].astype('float64')
      dat1[np.isnan(dat1)] = 0
      dat1 = ppdrc.ppdrc(dat1, shape_port[-1]/4)
      dat1 = humutils.rescale(dat1.getdata(),0,255)
      dat1[np.isnan(Zt[p])] = np.nan
      Zt[p] = dat1.astype('float32')
      del dat1

   # create memory mapped file for Z
   fp = np.memmap(sonpath+base+'_data_port_la.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
   fp[:] = Zt[:]
   del fp
   shape_port = np.shape(Zt)
   del Zt
   #we are only going to access the portion of memory required
   port_fp = np.memmap(sonpath+base+'_data_port_la.dat', dtype='float32', mode='r', shape=shape_port)

   ## do plots of merged scans
   if doplot==1:

      for p in xrange(len(port_fp)):
         merge = np.vstack((np.flipud(port_fp[p]), star_fp[p]))

         Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]

         fig = plt.figure()
         plt.imshow(merge, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
         plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

         plt.axis('normal'); plt.axis('tight')
         custom_save(sonpath,'merge_corrected_scan'+str(p))
         del fig

#      # treats each chunk in parallel for speed
#      try:
#         d = Parallel(n_jobs = min(cpu_count(),len(star_fp)), verbose=0)(delayed(plot_merged_scans)(port_fp[p], star_fp[p], dist_m, shape_port, ft, sonpath, p) for p in xrange(len(star_fp)))
#      except:
#         print "memory error: trying serial"
#         d = Parallel(n_jobs = 1, verbose=0)(delayed(plot_merged_scans)(port_fp[p], star_fp[p], dist_m, shape_port, ft, sonpath, p) for p in xrange(len(star_fp)))

   # load memory mapped scans
   shape_low = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_low'])
   if shape_low!='':
      low_fp = np.memmap(sonpath+base+'_data_dwnlow.dat', dtype='int16', mode='r', shape=tuple(shape_low))

   shape_hi = np.squeeze(loadmat(sonpath+base+'meta.mat')['shape_hi'])
   if shape_hi!='':
      hi_fp = np.memmap(sonpath+base+'_data_dwnhi.dat', dtype='int16', mode='r', shape=tuple(shape_hi))

   if 'low_fp' in locals():
      ######### low
      Zt = remove_water(low_fp, bed, shape_low, maxW)

      # create memory mapped file for Z
      fp = np.memmap(sonpath+base+'_data_dwnlow_l.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
      fp[:] = Zt[:]
      del fp
      shape_low = np.shape(Zt)
      del Zt
      #we are only going to access the portion of memory required
      low_fp = np.memmap(sonpath+base+'_data_dwnlow_l.dat', dtype='float32', mode='r', shape=shape_low)

      Zt = correct_scans(low_fp)
      for p in xrange(len(Zt)):

         dat1 = Zt[p].astype('float64')
         dat1[np.isnan(dat1)] = 0
         dat1 = ppdrc.ppdrc(dat1, shape_low[-1]/4)
         dat1 = humutils.rescale(dat1.getdata(),0,255)
         dat1[np.isnan(Zt[p])] = np.nan
         Zt[p] = dat1.astype('float32')
         del dat1

      # create memory mapped file for Z
      fp = np.memmap(sonpath+base+'_data_dwnlow_la.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
      fp[:] = Zt[:]
      del fp
      shape_low = np.shape(Zt)
      del Zt
      #we are only going to access the lowion of memory required
      low_fp = np.memmap(sonpath+base+'_data_dwnlow_la.dat', dtype='float32', mode='r', shape=shape_low)

      if doplot==1:

         for p in xrange(len(low_fp)):

            Zdist = dist_m[shape_low[-1]*p:shape_low[-1]*(p+1)]

            fig = plt.figure()
            plt.imshow(low_fp[p], cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0])
            plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

            plt.axis('normal'); plt.axis('tight')
            custom_save(sonpath,'dwnlow_corrected_scan'+str(p))
            del fig

#         # treats each chunk in parallel for speed
#         try:
#            d = Parallel(n_jobs = min(cpu_count(),len(low_fp)), verbose=0)(delayed(plot_dwnlow_scans)(low_fp[p], dist_m[shape_low[-1]*p:shape_low[-1]*(p+1)], shape_low, ft, sonpath, p) for p in xrange(len(low_fp)))
#         except:
#            print "memory error: trying serial"
#            d = Parallel(n_jobs = 1, verbose=0)(delayed(plot_dwnlow_scans)(low_fp[p], dist_m[shape_low[-1]*p:shape_low[-1]*(p+1)], shape_low, ft, sonpath, p) for p in xrange(len(low_fp)))


   if 'hi_fp' in locals():
      ######### hi
      Zt = remove_water(hi_fp, bed, shape_hi, maxW)

      # create memory mapped file for Z
      fp = np.memmap(sonpath+base+'_data_dwnhi_l.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
      fp[:] = Zt[:]
      del fp
      shape_hi = np.shape(Zt)
      del Zt
      #we are only going to access the portion of memory required
      hi_fp = np.memmap(sonpath+base+'_data_dwnhi_l.dat', dtype='float32', mode='r', shape=shape_hi)

      Zt = correct_scans(hi_fp)
      for p in xrange(len(Zt)):
         dat1 = Zt[p].astype('float64')
         dat1[np.isnan(dat1)] = 0
         dat1 = ppdrc.ppdrc(dat1, shape_hi[-1]/4)
         dat1 = humutils.rescale(dat1.getdata(),0,255)
         dat1[np.isnan(Zt[p])] = np.nan
         Zt[p] = dat1.astype('float32')
         del dat1

      # create memory mapped file for Z
      fp = np.memmap(sonpath+base+'_data_dwnhi_la.dat', dtype='float32', mode='w+', shape=np.shape(Zt))
      fp[:] = Zt[:]
      del fp
      shape_hi = np.shape(Zt)
      del Zt
      #we are only going to access the hiion of memory required
      hi_fp = np.memmap(sonpath+base+'_data_dwnhi_la.dat', dtype='float32', mode='r', shape=shape_hi)

      if doplot==1:

         for p in xrange(len(hi_fp)):
  
            Zdist = dist_m[shape_hi[-1]*p:shape_hi[-1]*(p+1)]

            fig = plt.figure()
            plt.imshow(hi_fp[p], cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0])
            plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

            plt.axis('normal'); plt.axis('tight')
            custom_save(sonpath,'dwnhi_corrected_scan'+str(p))
            del fig

#         # treats each chunk in parallel for speed
#         try:
#            d = Parallel(n_jobs = min(cpu_count(),len(hi_fp)), verbose=0)(delayed(plot_dwnhi_scans)(hi_fp[p], dist_m[shape_hi[-1]*p:shape_hi[-1]*(p+1)], shape_hi, ft, sonpath, p) for p in xrange(len(hi_fp)))
#         except:
#            print "memory error: trying serial"
#            d = Parallel(n_jobs = 1, verbose=0)(delayed(plot_dwnhi_scans)(hi_fp[p], dist_m[shape_hi[-1]*p:shape_hi[-1]*(p+1)], shape_hi, ft, sonpath, p) for p in xrange(len(hi_fp)))

   if os.name=='posix': # true if linux/mac
      elapsed = (time.time() - start)
   else: # windows
      elapsed = (time.clock() - start)
   print "Processing took ", elapsed , "seconds to analyse"
   print "Done!"


# =========================================================
def custom_save(figdirec,root):
   try:
      plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
   except:
      plt.savefig(os.path.expanduser("~")+os.sep+root,bbox_inches='tight',dpi=400)      


# =========================================================
def remove_water(fp, bed, shape, maxW):
   Zt = []

   for p in xrange(len(fp)):
      data_dB = fp[p]*(10*np.log10(maxW)/255)

      Zbed = np.squeeze(bed[shape[-1]*p:shape[-1]*(p+1)])

      # shift proportionally depending on where the bed is
      for k in xrange(np.shape(data_dB)[1]):
         try:
            data_dB[:,k] = np.r_[data_dB[Zbed[k]:,k], np.zeros( (np.shape(data_dB)[0] -  np.shape(data_dB[Zbed[k]:,k])[0] ,) )]
         except:
            data_dB[:,k] = np.ones(np.shape(data_dB)[0])

      Zt.append(data_dB)    

   return Zt

# =========================================================
def correct_scans(fp):
   Zt = []

   for p in xrange(len(fp)):

      mg = 10**np.log10(np.asarray(fp[p],'float32')+0.001)
      mg[fp[p]==0] = np.nan
      extent = np.shape(mg)[0]
         
      beam_av = nanmedian(mg, axis=1)+1 #nanmean(mg, axis=1)+1
      beam_av[np.isinf(beam_av)] = np.nan
      index = np.where( 1-np.sum(np.asarray(np.isnan(mg),'int'),axis=1)/extent <0)
      beam_av[index] = np.nan
      ref = nanmedian(beam_av) #nanmean(beam_av)

      # correct for directivity index
      mg_la = np.ones(np.shape(mg))
      # apply attenuation corrections
      for k in xrange(len(beam_av)):
         mg_la[k,:] = (mg[k,:]/beam_av[k])*ref

      beam_av = nanmedian(mg_la, axis=1)+1 #nanmean(mg, axis=1)+1
      beam_av[np.isinf(beam_av)] = np.nan
      index = np.where( 1-np.sum(np.asarray(np.isnan(mg_la),'int'),axis=1)/extent <0)
      beam_av[index] = np.nan
      ref = nanmedian(beam_av) #nanmean(beam_av)

      # apply attenuation corrections
      for k in xrange(len(beam_av)):
         mg_la[k,:] = (mg_la[k,:]/beam_av[k])*ref

      del mg, ref, beam_av, index
      Zt.append(mg_la)
      del mg_la

   return Zt

## =========================================================
#def plot_merged_scans(dat_port, dat_star, dist_m, shape_port, ft, sonpath, p):

#   Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
#   extent = shape_port[1] #np.shape(merge)[0]

#   fig = plt.figure()
#   plt.imshow(np.vstack((np.flipud(dat_port), dat_star)), cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
#   plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

#   plt.axis('normal'); plt.axis('tight')
#   custom_save(sonpath,'merge_corrected_scan'+str(p))
#   del fig

## =========================================================
#def plot_dwnlow_scans(dat_dwnlow, Zdist, shape_low, ft, sonpath, p):

#   #Zdist = dist_m[shape_low[-1]*p:shape_low[-1]*(p+1)]
#   extent = shape_low[1] #np.shape(merge)[0]

#   fig = plt.figure()
#   plt.imshow(dat_dwnlow, cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0])

#   plt.axis('normal'); plt.axis('tight')
#   custom_save(sonpath,'dwnlow_corrected_scan'+str(p))
#   del fig

## =========================================================
#def plot_dwnhi_scans(dat_dwnhi, Zdist, shape_hi, ft, sonpath, p):

#   #Zdist = dist_m[shape_hi[-1]*p:shape_hi[-1]*(p+1)]
#   extent = shape_hi[1] #np.shape(merge)[0]

#   fig = plt.figure()
#   plt.imshow(dat_dwnhi, cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0])
#   plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

#   plt.axis('normal'); plt.axis('tight')
#   custom_save(sonpath,'dwnhi_corrected_scan'+str(p))
#   del fig




