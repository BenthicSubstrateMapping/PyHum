'''
pyhum_read.py
Part of PyHum software 

INFO:
Python script to read Humminbird DAT and associated SON files, and export data

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.0      Revision: June, 2014

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 series instruments. 

SYNTAX:
python pyhum_read.py -i datfile -s sonpath
where datfile is the .DAT file associated with the survey, and sonpath is the (absolute or relative) path to where the associated .SON files are

Optional arguments:
-d measured draft of instrument in metres [Default = 0]
-c coordinate transformation. By default coordinates are output in WGS84. If you would like an additional coordinate transformation, specify the EPSG ID number for use in cs2cs (pyproj). See: http://cs2cs.mygeodata.eu/; http://www.eye4software.com/resources/stateplane/ [Default is Arizona State Plane: -c "epsg:26949"]
-p make simple plots of data [Default = 1 (yes); 0 = no]

EXAMPLES:
1) show help
python pyhum_read.py -h

2) run the provided test case with all defaults
python pyhum_read.py -i ./test.DAT -s ./test_data/ (linux)
python pyhum_read.py -i test.DAT -s \test_data\ (windows)

3) run a file and output eastings/northings in NAD Colorado Central (26954) with a draft of 0.4m
python pyhum_read.py -i ./test.DAT -s ./test_data/ -c "epsg:26954" -d 0.4

4) run a file and output eastings/northings in OSGB36 British National Grid (27700) with a draft of 1m 
python pyhum_read.py -i ./test.DAT -s ./test_data/ -c "epsg:27700" -d 1 

5) run the provided test case with all defaults except no plots
python pyhum_read.py -i ./test.DAT -s ./test_data/ -p 0

OUTPUTS:
Files are created which contain the raw and parsed meta data. They are prefixed by the root of the input file (*) followed by:
1) *.mat = all raw data (port side scan, starboard side scan, downward low freq, downward high freq)
2) *meta.mat = longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)

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
import glob, sys, getopt
from scipy.io import savemat, loadmat
import os, time
from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory

#numerical
import pyread
from numpy import shape
#from pyhum_utils import sliding_window
import PyHum.utils as humutils

#plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

#################################################
class humread:

   def __init__(self, humfile, sonpath, cs2cs_args, draft, doplot):

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
         print 'Son files are in %s' % (sonpath)
      if cs2cs_args:
         print 'cs2cs arguments are %s' % (cs2cs_args)
      if draft:
         draft = float(draft)
         print 'Draft: %s' % (str(draft))
      if doplot:
         doplot = int(doplot)
         if doplot==0:
            print "Plots will not be made"

      if not draft:
         draft = 0
         print '[Default] Draft = %s metres' % (str(draft))
      if not cs2cs_args:
         # arguments to pass to cs2cs for coordinate transforms
         cs2cs_args = "epsg:26949"
         print '[Default] cs2cs arguments are %s' % (cs2cs_args)
      if not doplot:
         if doplot != 0:
            doplot = 1
            print "[Default] Plots will be made"

      ## for debugging
      #humfile = r"test.DAT"; sonpath = "test_data"
      #cs2cs_args = "epsg:26949"; doplot = 1; draft = 0

      # start timer
      if os.name=='posix': # true if linux/mac or cygwin on windows
         start = time.time()
      else: # windows
         start = time.clock()

      # number of bytes in a header packet in SON file
      headbytes = 67

      # if son path name supplied has no separator at end, put one on
      if sonpath[-1]!=os.sep:
         sonpath = sonpath + os.sep

      base = humfile.split('.DAT') # get base of file name for output
      base = base[0].split(os.sep)[-1]

      # get the SON files from this directory
      sonfiles = glob.glob(sonpath+'*.SON')
      if not sonfiles:
          sonfiles = glob.glob(os.getcwd()+os.sep+sonpath+'*.SON')

      print "WARNING: Because files have to be read in byte by byte,"
      print "this could take a very long time ..."

      data = pyread.pyread(sonfiles, humfile, headbytes, cs2cs_args)

      data_port = data.getportscans()
      data_star = data.getstarscans()
      data_dwnlow = data.getlowscans()
      data_dwnhi = data.gethiscans()

      dat = data.gethumdat() 
      metadat = data.getmetadata()
      del data

      #imshow(data_port)
      try:
         savemat(sonpath+base+'.mat', mdict={'dat': dat, 'data_port': data_port, 'data_star': data_star, 'data_dwnlow': data_dwnlow, 'data_dwnhi': data_dwnhi},oned_as='row')
         savemat(sonpath+base+'meta.mat', metadat ,oned_as='row')
      except:
         sonpath = os.getcwd()+os.sep   
         savemat(sonpath+base+'.mat', mdict={'dat': dat, 'data_port': data_port, 'data_star': data_star, 'data_dwnlow': data_dwnlow, 'data_dwnhi': data_dwnhi},oned_as='row')
         savemat(sonpath+base+'meta.mat', metadat ,oned_as='row')
         
      if doplot==1:

         import ppdrc
         # phase preserving dynamic range compression
         dat = ppdrc.ppdrc(data_port.astype('float64'), 768)
         data_port = dat.getdata()
         del dat

         # port
         nx, ny = shape(data_port)
         if ny>10000:
            Z,inds = humutils.sliding_window(data_port,(nx,10000))
            del nx, ny, inds
            if len(Z) != len(data_port):
               flag = 1
            else:
               flag=0
            del data_port
         else:
            Z = data_port.copy()
            del data_port
            flag = 0
   
         if flag==1:
            for k in xrange(len(Z)):
               fig = plt.figure()
               plt.imshow(Z[k],cmap='gray', origin = 'upper'); plt.colorbar(); 
               plt.title('Portside Raw Scan'+str(k))
               plt.xlabel('Ping Number (Time)')
               plt.ylabel('Range (Distance)')
               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'raw_port'+str(k))
               del fig
         else:
            fig = plt.figure()
            plt.imshow(Z,cmap='gray', origin = 'upper'); plt.colorbar(); 
            plt.title('Portside Raw Scan')
            plt.xlabel('Ping Number (Time)')
            plt.ylabel('Range (Distance)')
            plt.axis('normal'); plt.axis('tight')
            self._custom_save(sonpath,'raw_port')
            del fig

         del Z

         dat = ppdrc.ppdrc(data_star.astype('float64'), 768)
         data_star = dat.getdata()
         del dat
         
         # starboard
         nx, ny = shape(data_star)
         if ny>10000:
            Z,inds = humutils.sliding_window(data_star,(nx,10000))
            del nx, ny, inds
            if len(Z) != len(data_star):
               flag = 1
            else:
               flag=0
            del data_star
         else:
            Z = data_star.copy()
            del data_star
            flag = 0
   
         if flag==1:
            for k in xrange(len(Z)):
               fig = plt.figure()
               plt.imshow(Z[k],cmap='gray', origin = 'upper'); plt.colorbar(); 
               plt.title('Starboardside Raw Scan'+str(k))
               plt.xlabel('Ping Number (Time)')
               plt.ylabel('Range (Distance)')
               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'raw_star'+str(k))
               del fig
         else:
            fig = plt.figure()
            plt.imshow(Z,cmap='gray', origin = 'upper'); plt.colorbar(); 
            plt.title('Starboardside Raw Scan')
            plt.xlabel('Ping Number (Time)')
            plt.ylabel('Range (Distance)')
            plt.axis('normal'); plt.axis('tight')
            self._custom_save(sonpath,'raw_star')
            del fig

         del Z

         if 'data_dwnlow' in locals():

            dat = ppdrc.ppdrc(data_dwnlow.astype('float64'), 768)
            data_dwnlow = dat.getdata()
            del dat

            # dwnlow
            nx, ny = shape(data_dwnlow)
            if ny>10000:
               Z,inds = humutils.sliding_window(data_dwnlow,(nx,10000))
               del nx, ny, inds
               if len(Z) != len(data_dwnlow):
                  flag = 1
               else:
                  flag=0
               del data_dwnlow
            else:
               Z = data_dwnlow.copy()
               del data_dwnlow
               flag = 0
   
            if flag==1:
               for k in xrange(len(Z)):
                  fig = plt.figure()
                  plt.imshow(Z[k],cmap='gray', origin = 'upper'); plt.colorbar(); 
                  plt.title('Downward Low Raw Scan'+str(k))
                  plt.xlabel('Ping Number (Time)')
                  plt.ylabel('Range (Distance)')
                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'raw_dwnlow'+str(k))
                  del fig
            else:
               fig = plt.figure()
               plt.imshow(Z,cmap='gray', origin = 'upper'); plt.colorbar(); 
               plt.title('Downward Low Raw Scan')
               plt.xlabel('Ping Number (Time)')
               plt.ylabel('Range (Distance)')
               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'raw_dwnlow')
               del fig

         if 'data_dwnhi' in locals():

            dat = ppdrc.ppdrc(data_dwnhi.astype('float64'), 768)
            data_dwnhi = dat.getdata()
            del dat

            # dwnhigh
            nx, ny = shape(data_dwnhi)
            if ny>10000:
               Z,inds = humutils.sliding_window(data_dwnhi,(nx,10000))
               del nx, ny, inds
               if len(Z) != len(data_dwnhi):
                  flag = 1
               else:
                  flag=0
               del data_dwnhi
            else:
               Z = data_dwnhi.copy()
               del data_dwnhi
               flag = 0
   
            if flag==1:
               for k in xrange(len(Z)):
                  fig = plt.figure()
                  plt.imshow(Z[k],cmap='gray', origin = 'upper'); plt.colorbar(); 
                  plt.title('Downward High Raw Scan'+str(k))
                  plt.xlabel('Ping Number (Time)')
                  plt.ylabel('Range (Distance)')
                  plt.axis('normal'); plt.axis('tight')
                  self._custom_save(sonpath,'raw_dwnhi'+str(k))
                  del fig
            else:
               fig = plt.figure()
               plt.imshow(Z,cmap='gray', origin = 'upper'); plt.colorbar(); 
               plt.title('Downward High Raw Scan')
               plt.xlabel('Ping Number (Time)')
               plt.ylabel('Range (Distance)')
               plt.axis('normal'); plt.axis('tight')
               self._custom_save(sonpath,'raw_dwnhi')
               del fig

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
         plt.savefig(os.getcwd()+os.sep+root,bbox_inches='tight',dpi=400)



