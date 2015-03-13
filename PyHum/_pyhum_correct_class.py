'''
pyhum_correct.py
Part of PyHum software 

INFO:


Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.1.4      Revision: Mar, 2015

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

#numerical
import numpy as np
#from pyhum_utils import rm_spikes, sliding_window, runningMeanFast, dpboundary, rescale
import PyHum.utils as humutils
from scipy.stats import nanmean
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

#################################################
class humcorrect:

   def __init__(self, humfile, sonpath, maxW, doplot):

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

      show_stages = 1

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
      try:
         meta = loadmat(sonpath+base+'meta.mat')
      except:
         meta = loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')
      meta['maxW'] = maxW

      try:
         savemat(sonpath+base+'meta.mat', meta ,oned_as='row')
      except:
         savemat(os.path.expanduser("~")+os.sep+base+'meta.mat', meta ,oned_as='row')

      del meta

      try:
         bed = np.squeeze(loadmat(sonpath+base+'meta.mat')['bed'])
         ft = 1/loadmat(sonpath+base+'meta.mat')['pix_m'] #np.squeeze(loadmat(sonpath+base+'meta.mat')['ft'])
         dist_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['dist_m'])
         data_port = np.asarray(np.squeeze(loadmat(sonpath+base+'.mat')['data_port']),'float16')
         data_star = np.asarray(np.squeeze(loadmat(sonpath+base+'.mat')['data_star']),'float16')
      except:
         bed = np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['bed'])
         ft = 1/loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['pix_m'] #np.squeeze(loadmat(sonpath+base+'meta.mat')['ft'])
         dist_m = np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['dist_m'])
         data_port = np.asarray(np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'.mat')['data_port']),'float16')
         data_star = np.asarray(np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'.mat')['data_star']),'float16')

      extent = np.shape(data_port)[0]

      # calculate in dB
      # port and star
      data_star_dB = data_star*(10*np.log10(maxW)/255)
      del data_star

      data_port_dB = data_port*(10*np.log10(maxW)/255)
      del data_port

      if show_stages==1:
         nx, ny = np.shape(data_star_dB)
         # make a plot of the starboard side with the bed pick
         fig = plt.figure()
         plt.subplot(2,2,1)
         plt.imshow(data_star_dB[:,:min(ny,10000)],cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
         plt.colorbar()
         plt.axis('normal'); plt.axis('tight')
         plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')

         plt.subplot(2,2,3)
         plt.imshow(data_port_dB[:,:min(ny,10000)],cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
         plt.colorbar()
         plt.axis('normal'); plt.axis('tight')
         plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')

         self._custom_save(sonpath,'dB')
         del fig

      bed = np.asarray(bed,'int')+int(0.25*ft)
      # shift proportionally depending on where the bed is
      for k in xrange(np.shape(data_star_dB)[1]):
         data_star_dB[:,k] = np.r_[data_star_dB[bed[k]:,k], np.zeros( (np.shape(data_star_dB)[0] -  np.shape(data_star_dB[bed[k]:,k])[0] ,) )]
         data_port_dB[:,k] = np.r_[data_port_dB[bed[k]:,k], np.zeros( (np.shape(data_port_dB)[0] -  np.shape(data_port_dB[bed[k]:,k])[0] ,) )]


      star_mg = 10**np.log10(np.asarray(data_star_dB,'float64')+0.001)
      port_mg = 10**np.log10(np.asarray(data_port_dB,'float64')+0.001)

      #star_mg = 20*np.log10(np.asarray(data_star_dB,'float64')+0.001)
      #port_mg = 20*np.log10(np.asarray(data_port_dB,'float64')+0.001)

      star_mg[data_star_dB==0] = np.nan
      port_mg[data_port_dB==0] = np.nan

      del data_star_dB, data_port_dB

      if show_stages==1:
         nx, ny = np.shape(star_mg)
         # make a plot of the starboard side with the bed pick
         fig = plt.figure()
         plt.subplot(2,2,1)
         plt.imshow(star_mg[:,:min(ny,10000)],cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
         plt.colorbar()
         plt.axis('normal'); plt.axis('tight')
         plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')

         plt.subplot(2,2,3)
         plt.imshow(port_mg[:,:min(ny,10000)],cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
         plt.colorbar()
         plt.axis('normal'); plt.axis('tight')
         plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')

         self._custom_save(sonpath,'dB_directivity')
         del fig

      ## get average dB per 'beam' and the reference
      # which is simply the mean
      beam_av_port = nanmean(port_mg, axis=1)+1
      beam_av_port[np.isinf(beam_av_port)] = np.nan
      index = np.where( 1-np.sum(np.asarray(np.isnan(port_mg),'int'),axis=1)/extent <0)
      beam_av_port[index] = np.nan
      ref_port = nanmean(beam_av_port)

      beam_av_star = nanmean(star_mg, axis=1)+1
      beam_av_star[np.isinf(beam_av_star)] = np.nan
      index = np.where( 1-np.sum(np.asarray(np.isnan(star_mg),'int'),axis=1)/extent <0)
      beam_av_star[index] = np.nan
      ref_star = nanmean(beam_av_star)

      # correct for directivity index
      star_mg_la = np.ones(np.shape(star_mg))
      port_mg_la = np.ones(np.shape(port_mg))
      # apply attenuation corrections
      for k in xrange(len(beam_av_star)):
         star_mg_la[k,:] = (star_mg[k,:]/beam_av_star[k])*ref_star
         port_mg_la[k,:] = (port_mg[k,:]/beam_av_port[k])*ref_port

      del star_mg, port_mg

      # scale star by port intensities
      cf = beam_av_port/beam_av_star
      for k in xrange(len(beam_av_star)):
         star_mg_la[k,:] = star_mg_la[k,:]/cf[k]

      if show_stages==1:
         nx, ny = np.shape(star_mg_la)
         # make a plot of the starboard side with the bed pick
         fig = plt.figure()
         plt.subplot(2,2,1)
         plt.imshow(star_mg_la[:,:min(ny,10000)],cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
         plt.colorbar()
         plt.axis('normal'); plt.axis('tight')
         plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')

         plt.subplot(2,2,3)
         plt.imshow(port_mg_la[:,:min(ny,10000)],cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
         plt.colorbar()
         plt.axis('normal'); plt.axis('tight')
         plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')

         self._custom_save(sonpath,'dB_directivity_attenuation')
         del fig

      del beam_av_port, beam_av_star, ref_star, ref_port, cf

      if doplot==1:
         merge = np.vstack((np.flipud(port_mg_la),star_mg_la))
         merge[np.isnan(merge)] = 0

         # phase preserving dynamic range compression
         dat = ppdrc.ppdrc(merge, nx/2)
         merge2 = humutils.rescale(dat.getdata(),0,255)
         del dat
         del port_mg_la, star_mg_la

         merge[merge==0] = np.nan
         merge2[np.isnan(merge)] = np.nan

         ################################################## ppdrc-corrected
         nx, ny = np.shape(merge)
         if ny>10000:
            Z,inds = humutils.sliding_window(merge2,(nx,10000))
            del ny, inds
            Zdist,inds = humutils.sliding_window(np.squeeze(dist_m),(10000))
            del inds 
            flag=1
         else:
            Z = merge2.copy()  
            Zdist = np.squeeze(dist_m.copy())
            flag=0
         del merge2

         if flag==1:
            if len(Z)==nx:
               # show the merged corrected sidescan as greyscale
               fig = plt.figure()
               plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
               plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

               if humfile=='test.DAT':
                  plt.ylim(-20, 20)

               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'merge_corrected_scan')
               del fig

            else:
               for k in xrange(len(Z)):
                  # show the merged corrected sidescan as greyscale
                  fig = plt.figure()
                  plt.imshow(Z[k], cmap='gray', extent=[min(Zdist[k]), max(Zdist[k]), -extent*(1/ft), extent*(1/ft)])
                  plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                  if humfile=='test.DAT':
                     plt.ylim(-20, 20)

                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'merge_corrected_scan'+str(k))
                  del fig

               fig = plt.figure()
               plt.imshow(merge[:,-(np.shape(merge)[1]-len(Z)*10000):], cmap='gray', extent=[min(dist_m[-(np.shape(merge)[1]-len(Z)*10000):]), max(dist_m[-(np.shape(merge)[1]-len(Z)*10000):]), -extent*(1/ft), extent*(1/ft)])
               plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')
               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'merge_corrected_scan'+str(k+1))
               del fig

         else:
            fig = plt.figure()
            plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
            plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

            if humfile=='test.DAT':
               plt.ylim(-20, 20)

            plt.axis('normal'); plt.axis('tight')
            self._custom_save(sonpath,'merge_corrected_scan_ppdrc')
            del fig   
         del Z, Zdist

       ################################################## non-ppdrc-corrected
         nx, ny = np.shape(merge)
         if ny>10000:
            Z,inds = humutils.sliding_window(merge,(nx,10000))
            del ny, inds
            Zdist,inds = humutils.sliding_window(np.squeeze(dist_m),(10000))
            del inds 
            flag=1
         else:
            Z = merge.copy()  
            Zdist = np.squeeze(dist_m.copy())
            flag=0

         if flag==1:
            if len(Z)==nx:
               # show the merged corrected sidescan as greyscale
               fig = plt.figure()
               plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
               plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

               if humfile=='test.DAT':
                  plt.ylim(-20, 20)

               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'merge_corrected_scan')
               del fig

            else:
               for k in xrange(len(Z)):
                  # show the merged corrected sidescan as greyscale
                  fig = plt.figure()
                  plt.imshow(Z[k], cmap='gray', extent=[min(Zdist[k]), max(Zdist[k]), -extent*(1/ft), extent*(1/ft)])
                  plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                  if humfile=='test.DAT':
                     plt.ylim(-20, 20)

                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'merge_corrected_scan'+str(k))
                  del fig

               fig = plt.figure()
               plt.imshow(merge[:,-(np.shape(merge)[1]-len(Z)*10000):], cmap='gray', extent=[min(dist_m[-(np.shape(merge)[1]-len(Z)*10000):]), max(dist_m[-(np.shape(merge)[1]-len(Z)*10000):]), -extent*(1/ft), extent*(1/ft)])
               plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')
               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'merge_corrected_scan'+str(k+1))
               del fig

         else:
            fig = plt.figure()
            plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
            plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

            if humfile=='test.DAT':
               plt.ylim(-20, 20)

            plt.axis('normal'); plt.axis('tight')
            self._custom_save(sonpath,'merge_corrected_scan')
            del fig   
         del Z, Zdist

         # retreive port and starboard from merge
         port_mg_la = np.flipud(merge[:extent,:])
         star_mg_la = merge[extent:,:]
         del merge

      port_mg_la[np.isnan(port_mg_la)] = 0
      star_mg_la[np.isnan(star_mg_la)] = 0

      try:
         savemat(sonpath+base+'port_la.mat', mdict={'port_mg_la': np.asarray(port_mg_la, 'float16')},oned_as='row')
      except:
         savemat(os.path.expanduser("~")+os.sep+base+'port_la.mat', mdict={'port_mg_la': np.asarray(port_mg_la, 'float16')},oned_as='row')      

      del port_mg_la

      try:
         savemat(sonpath+base+'star_la.mat', mdict={'star_mg_la': np.asarray(star_mg_la, 'float16')},oned_as='row')
      except:
         savemat(os.path.expanduser("~")+os.sep+base+'star_la.mat', mdict={'star_mg_la': np.asarray(star_mg_la, 'float16')},oned_as='row')


      del star_mg_la

      ###########################################################
      ###########################################################
      # load in high freq echogram
      try:
         data_dwnhi = np.asarray(np.squeeze(loadmat(sonpath+base+'.mat')['data_dwnhi']),'float16')
      except:
         data_dwnhi = np.asarray(np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'.mat')['data_dwnhi']),'float16')

      if np.size(data_dwnhi)>0:
         # calculate in dB
         data_dwnhi_dB = data_dwnhi*(10*np.log10(maxW)/255)
         del data_dwnhi

         i = np.linspace(1,np.shape(data_dwnhi_dB)[1],len(bed))
         bedi = np.interp(np.linspace(1,np.shape(data_dwnhi_dB)[1],np.shape(data_dwnhi_dB)[1]), i, bed) 
         del i
         bedi = np.asarray(bedi,'int')

         # shift proportionally depending on where the bed is
         for k in xrange(np.shape(data_dwnhi_dB)[1]):
            data_dwnhi_dB[:,k] = np.r_[data_dwnhi_dB[bed[k]:,k], np.zeros( (np.shape(data_dwnhi_dB)[0] -  np.shape(data_dwnhi_dB[bed[k]:,k])[0] ,) )]

         dwnhi_mg = 10**np.log10(np.asarray(data_dwnhi_dB,'float64')+0.001)
         #dwnhi_mg = 20*np.log10(np.asarray(data_dwnhi_dB,'float64')+0.001)
         dwnhi_mg[data_dwnhi_dB==0] = np.nan
         del data_dwnhi_dB

         ## get average dB per 'beam' and the reference
         # which is simply the mean
         beam_av_dwnhi = nanmean(dwnhi_mg, axis=1)+1
         beam_av_dwnhi[np.isinf(beam_av_dwnhi)] = np.nan
         index = np.where( 1-np.sum(np.asarray(np.isnan(dwnhi_mg),'int'),axis=1)/extent <0)
         beam_av_dwnhi[index] = np.nan
         ref_dwnhi = nanmean(beam_av_dwnhi)

         # correct for directivity index
         dwnhi_mg_la = np.ones(np.shape(dwnhi_mg))
         # apply attenuation corrections
         for k in xrange(len(beam_av_dwnhi)):
            dwnhi_mg_la[k,:] = (dwnhi_mg[k,:]/beam_av_dwnhi[k])*ref_dwnhi

         del dwnhi_mg

         if doplot==1:
            nx, ny = np.shape(dwnhi_mg_la)
            if ny>10000:
               Z,inds = humutils.sliding_window(dwnhi_mg_la,(nx,10000))
               del ny, inds
               Zdist,inds = humutils.sliding_window(np.squeeze(dist_m),(10000))
               del inds 
               if len(Z) != len(dwnhi_mg_la):
                  flag = 1
               else:
                  flag=0
            else:
               Z = dwnhi_mg_la.copy()  
               Zdist = np.squeeze(dist_m.copy())
               flag=0

            if flag==1:
               if len(Z)==nx:
                  fig = plt.figure()
                  plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), 0])
                  plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                  if humfile=='test.DAT':
                     plt.ylim(-20, 20)

                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'dwnhi_corrected_scan')
                  del fig

               else:
                  for k in xrange(len(Z)):
                     fig = plt.figure()
                     plt.imshow(Z[k], cmap='gray', extent=[min(Zdist[k]), max(Zdist[k]), -extent*(1/ft), 0])
                     plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                     if humfile=='test.DAT':
                        plt.ylim(-20, 20)

                     plt.axis('normal'); plt.axis('tight')
                     self._custom_save(sonpath,'dwnhi_corrected_scan'+str(k))
                     del fig

                  fig = plt.figure()
                  plt.imshow(dwnhi_mg_la[:,-(np.shape(dwnhi_mg_la)[1]-len(Z)*10000):], cmap='gray', extent=[min(dist_m[-(np.shape(dwnhi_mg_la)[1]-len(Z)*10000):]), max(dist_m[-(np.shape(dwnhi_mg_la)[1]-len(Z)*10000):]), -extent*(1/ft), 0])
                  plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')
                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'dwnhi_corrected_scan'+str(k+1))
                  del fig

            else:
               fig = plt.figure()
               try:
                  plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), 0])
               except:
                  plt.imshow(Z, cmap='gray', extent=[min(Zdist[0]), max(Zdist[0]), -extent*(1/ft), 0])
               plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

               if humfile=='test.DAT':
                  plt.ylim(-20, 20)

               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'dwnhi_corrected_scan')
               del fig   
            del Z, Zdist

         dwnhi_mg_la[np.isnan(dwnhi_mg_la)] = 0
         try:
            savemat(sonpath+base+'dwnhi_la.mat', mdict={'dwnhi_mg_la': np.asarray(dwnhi_mg_la, 'float16')},oned_as='row')
         except:
            savemat(os.path.expanduser("~")+os.sep+base+'dwnhi_la.mat', mdict={'dwnhi_mg_la': np.asarray(dwnhi_mg_la, 'float16')},oned_as='row')

         del dwnhi_mg_la

           #################################################################################
      # load in low freq echogram
      try:
         data_dwnlow = np.asarray(np.squeeze(loadmat(sonpath+base+'.mat')['data_dwnlow']),'float16')
      except:
         data_dwnlow = np.asarray(np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'.mat')['data_dwnlow']),'float16')
     
      if np.size(data_dwnlow)>0:
         data_dwnlow_dB = data_dwnlow*(10*np.log10(maxW)/255)
         del data_dwnlow

         i = np.linspace(1,np.shape(data_dwnlow_dB)[1],len(bed))
         bedi = np.interp(np.linspace(1,np.shape(data_dwnlow_dB)[1],np.shape(data_dwnlow_dB)[1]), i, bed) 
         del i
         bedi = np.asarray(bedi,'int')

         # shift proportionally depending on where the bed is
         for k in xrange(np.shape(data_dwnlow_dB)[1]):
            data_dwnlow_dB[:,k] = np.r_[data_dwnlow_dB[bed[k]:,k], np.zeros( (np.shape(data_dwnlow_dB)[0] -  np.shape(data_dwnlow_dB[bed[k]:,k])[0] ,) )]

         dwnlow_mg = 10**np.log10(np.asarray(data_dwnlow_dB,'float64')+0.001)
         #dwnlow_mg = 20*np.log10(np.asarray(data_dwnlow_dB,'float64')+0.001)
         dwnlow_mg[data_dwnlow_dB==0] = np.nan
         del data_dwnlow_dB

         ## get average dB per 'beam' and the reference
         # which is simply the mean
         beam_av_dwnlow = nanmean(dwnlow_mg, axis=1)+1
         beam_av_dwnlow[np.isinf(beam_av_dwnlow)] = np.nan
         index = np.where( 1-np.sum(np.asarray(np.isnan(dwnlow_mg),'int'),axis=1)/extent <0)
         beam_av_dwnlow[index] = np.nan
         ref_dwnlow = nanmean(beam_av_dwnlow)

         # correct for directivity index
         dwnlow_mg_la = np.ones(np.shape(dwnlow_mg))
         # apply attenuation corrections
         for k in xrange(len(beam_av_dwnlow)):
            dwnlow_mg_la[k,:] = (dwnlow_mg[k,:]/beam_av_dwnlow[k])*ref_dwnlow

         del dwnlow_mg

         # scale dwnlow by dwnhi intensities
         cf = beam_av_dwnhi/beam_av_dwnlow
         for k in xrange(len(beam_av_dwnlow)):
            dwnlow_mg_la[k,:] = dwnlow_mg_la[k,:]/cf[k]

         del beam_av_dwnlow, ref_dwnlow, cf

         if doplot==1:

            nx, ny = np.shape(dwnlow_mg_la)
            if ny>10000:
               Z,inds = humutils.sliding_window(dwnlow_mg_la,(nx,10000))
               del ny, inds
               Zdist,inds = humutils.sliding_window(np.squeeze(dist_m),(10000))
               del inds 
               if len(Z) != len(dwnlow_mg_la):
                  flag = 1
               else:
                  flag=0
            else:
               Z = dwnlow_mg_la.copy()  
               Zdist = np.squeeze(dist_m.copy())
               flag=0

            if flag==1:
               if len(Z)==nx:
                  fig = plt.figure()
                  plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), 0])
                  plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                  if humfile=='test.DAT':
                     plt.ylim(-20, 20)

                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'dwnlow_corrected_scan')
                  del fig

               else:
                  for k in xrange(len(Z)):
                     fig = plt.figure()
                     plt.imshow(Z[k], cmap='gray', extent=[min(Zdist[k]), max(Zdist[k]), -extent*(1/ft), 0])
                     plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

                     if humfile=='test.DAT':
                        plt.ylim(-20, 20)

                     plt.axis('normal'); plt.axis('tight')
                     self._custom_save(sonpath,'dwnlow_corrected_scan'+str(k))
                     del fig

                  fig = plt.figure()
                  plt.imshow(dwnlow_mg_la[:,-(np.shape(dwnlow_mg_la)[1]-len(Z)*10000):], cmap='gray', extent=[min(dist_m[-(np.shape(dwnlow_mg_la)[1]-len(Z)*10000):]), max(dist_m[-(np.shape(dwnlow_mg_la)[1]-len(Z)*10000):]), -extent*(1/ft), 0])
                  plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')
                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'dwnlow_corrected_scan'+str(k+1))
                  del fig

            else:
               fig = plt.figure()
               try:
                  plt.imshow(Z, cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), 0])
               except:
                  plt.imshow(Z, cmap='gray', extent=[min(Zdist[0]), max(Zdist[0]), -extent*(1/ft), 0])
               plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

               if humfile=='test.DAT':
                  plt.ylim(-20, 20)

               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'dwnlow_corrected_scan')
               del fig   
            del Z, Zdist

         dwnlow_mg_la[np.isnan(dwnlow_mg_la)] = 0

         try:
            savemat(sonpath+base+'dwnlow_la.mat', mdict={'dwnlow_mg_la': np.asarray(dwnlow_mg_la, 'float16')},oned_as='row')
         except:
            savemat(os.path.expanduser("~")+os.sep+base+'dwnlow_la.mat', mdict={'dwnlow_mg_la': np.asarray(dwnlow_mg_la, 'float16')},oned_as='row')


         del dwnlow_mg_la

      if os.name=='posix': # true if linux/mac
         elapsed = (time.time() - start)
      else: # windows
         elapsed = (time.clock() - start)
      print "Processing took ", elapsed , "seconds to analyse"

      print "Done!"

   # =========================================================
   def _custom_save(self, figdirec,root):
      try:
         plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
      except:
         plt.savefig(os.path.expanduser("~")+os.sep+root,bbox_inches='tight',dpi=400)      

