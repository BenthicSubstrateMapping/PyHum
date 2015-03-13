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
Version: 1.0.9      Revision: Mar, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 and 1198 series instruments. 

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
import csv

#numerical
import pyread
import PyHum.utils as humutils
#from skimage.measure import LineModel, ransac
import numpy as np
import pyproj

#plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

#################################################
class humread:

   def __init__(self, humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr):

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
      if c:
         c = float(c)
         print 'Celerity of sound: %s m/s' % (str(c))
      if doplot:
         doplot = int(doplot)
         if doplot==0:
            print "Plots will not be made"
      if flip_lr:
         flip_lr = int(flip_lr)
         if flip_lr==1:
            print "Port and starboard will be flipped"
      if t:
         t = np.asarray(t,float)
         print 'Transducer length is %s m' % (str(t))
      if f:
         f = np.asarray(f,int)
         print 'Frequency is %s kHz' % (str(f))
      if bedpick:
         bedpick = np.asarray(bedpick,int)
         if bedpick==1:
            print 'Bed picking is auto'
         else:
            print 'Bed picking is manual'

      if not t:
         t = 0.108
         print '[Default] Transducer length is %s m' % (str(t))
      if not f:
         f = 455
         print '[Default] Frequency is %s kHz' % (str(f))
      if not c:
         c = 1450.0
         print '[Default] Celerity of sound = %s m/s' % (str(c))
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
      if not flip_lr:
         if flip_lr != 1:
            flip_lr = 0
            print "[Default] No port/starboard flipping"
      if not bedpick:
         bedpick = 1
         print '[Default] Bed picking is auto'

      ## for debugging
      #humfile = r"test.DAT"; sonpath = "test_data"
      #cs2cs_args = "epsg:26949"; doplot = 1; draft = 0

      # start timer
      if os.name=='posix': # true if linux/mac or cygwin on windows
         start = time.time()
      else: # windows
         start = time.clock()

      tonemap=1

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

      data = pyread.pyread(sonfiles, humfile, c, headbytes, cs2cs_args)


      try:
         if flip_lr==0:
            data_port = data.getportscans()
         else:
            data_port = data.getstarscans()
      except:
         data_port = ''
         print "portside scan not available"
      try:
         if flip_lr==0:
            data_star = data.getstarscans()
         else:
            data_star = data.getportscans()
      except:
         data_star = ''
         print "starboardside scan not available"
      try:
         data_dwnlow = data.getlowscans()
      except:
         data_dwnlow = ''
         print "low-freq. scan not available"
      try:
         data_dwnhi = data.gethiscans()
      except:
         data_dwnhi = ''
         print "high-freq. scan not available"

      dat = data.gethumdat() 
      metadat = data.getmetadata()
      del data

      savemat(sonpath+base+'.mat', mdict={'dat': dat, 'data_port': data_port, 'data_star': data_star, 'data_dwnlow': data_dwnlow, 'data_dwnhi': data_dwnhi},oned_as='row')

      try:
         es = humutils.runningMeanFast(metadat['e'],len(metadat['e'])/100)
         ns = humutils.runningMeanFast(metadat['n'],len(metadat['n'])/100)
      except:
         es = metadat['e']
         ns = metadat['n']
    
      metadat['es'] = es
      metadat['ns'] = ns

      trans =  pyproj.Proj(init=cs2cs_args)
      lon, lat = trans(es, ns, inverse=True)
      metadat['lon'] = lon
      metadat['lat'] = lat

      try:
         import simplekml
         # create kml for loading path into google earth
         kml = simplekml.Kml()
         ls = kml.newlinestring(name='trackline')
         ls.coords = zip(lon,lat)
         ls.extrude = 1
         ls.altitudemode = simplekml.AltitudeMode.relativetoground
         ls.style.linestyle.width = 5
         ls.style.linestyle.color = simplekml.Color.red
         kml.save(sonpath+base+"trackline.kml")
      except:
         print "install simplekml for kml plots"

      #e = np.squeeze(loadmat(sonpath+base+'meta.mat')['es'])
      #n = np.squeeze(loadmat(sonpath+base+'meta.mat')['ns'])

      #dist_m = np.cumsum(np.sqrt(np.gradient(e)**2 + np.gradient(n)**2))

      dist = np.zeros(len(lat))
      for k in xrange(len(lat)-1):
         dist[k] = self._distBetweenPoints(lat[k], lat[k+1], lon[k], lon[k+1])

      dist_m = np.cumsum(dist)

      # theta at 3dB in the horizontal
      theta3dB = np.arcsin(c/(t*(f*1000)))
      #resolution of 1 sidescan pixel to nadir
      ft = (np.pi/2)*(1/theta3dB)

      dep_m = np.squeeze(metadat['dep_m']) #loadmat(sonpath+base+'meta.mat')['dep_m'])

      dep_m = humutils.rm_spikes(dep_m,2)

      metadat['dist_m'] = dist_m

      data_port = np.asarray(np.squeeze(loadmat(sonpath+base+'.mat')['data_port']),'float16')
      data_star = np.asarray(np.squeeze(loadmat(sonpath+base+'.mat')['data_star']),'float16')

      extent = np.shape(data_port)[0]

      if bedpick == 1: # auto

         # get bed from depth trace
         bed = ft*dep_m
         # narrow image to within range of estimated bed
         imu = data_port[int(np.min(bed)):int(np.max(bed)),:]
       # use dynamic boundary tracing to get 2nd estimate of bed  
         x = np.squeeze(int(np.min(bed))+humutils.dpboundary(-imu.T))
         del imu 

         if doplot==1:

            dist_mi = np.linspace(np.min(dist_m),np.max(dist_m),len(dist_m))

            nx, ny = np.shape(data_port)
            if ny>10000:
               Z,inds = humutils.sliding_window(data_port,(nx,10000))
               del inds

            nx, ny = np.shape(data_star)
            if ny>10000:
               Zstar,inds = humutils.sliding_window(data_star,(nx,10000))
               del inds

            if ny>10000:
               Zx,inds = humutils.sliding_window(np.squeeze(x),(10000))
               del inds

            ny = np.shape(bed)[0]
            if ny>10000:   
               Zbed,inds = humutils.sliding_window(np.squeeze(bed),(10000))
               del inds

            if ny>10000:   
               Zdist,inds = humutils.sliding_window(np.squeeze(dist_mi),(10000))
               del inds   
            #del ny

            if 'Z' in locals():
               if len(Z)==nx:
                  fig = plt.figure()
                  fig.subplots_adjust(wspace = 0.1, hspace=0.1)
                  plt.subplot(2,2,1)
                  ax = plt.gca()
                  im = ax.imshow(np.flipud(Z),cmap='gray',extent=[min(Zdist), max(Zdist), 0, extent*(1/ft)],origin='upper')
                  plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')  
                  plt.axis('normal'); plt.axis('tight')

                  plt.subplot(2,2,3)
                  ax = plt.gca()
                  im = ax.imshow(Zstar,cmap='gray',extent=[min(Zdist), max(Zdist), extent*(1/ft), 0],origin='upper')
                  plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
                  plt.axis('normal'); plt.axis('tight')

                  axR=plt.subplot(1,2,2); 
                  axR.yaxis.tick_right()
                  axR.yaxis.set_label_position("right")
                  axR.imshow(Zstar,cmap='gray',extent=[min(Zdist), max(Zdist), extent*(1/ft), 0],origin='upper')
                  plt.plot(Zdist,Zbed/ft,'k')
                  plt.plot(Zdist,Zx/ft,'r')
                  plt.axis('normal'); plt.axis('tight')
                  plt.ylim(10,0)
                  plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
                  self._custom_save(sonpath,'bed_2picks')
                  del fig

               else:
                  for k in xrange(len(Z)):
                     fig = plt.figure()
                     fig.subplots_adjust(wspace = 0.1, hspace=0.1)
                     plt.subplot(2,2,1)
                     ax = plt.gca()
                     im = ax.imshow(np.flipud(Z[k]),cmap='gray',extent=[min(Zdist[k]), max(Zdist[k]), 0, extent*(1/ft)],origin='upper')
                     plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')  
                     plt.axis('normal'); plt.axis('tight')

                     plt.subplot(2,2,3)
                     ax = plt.gca()
                     im = ax.imshow(Zstar[k],cmap='gray',extent=[min(Zdist[k]), max(Zdist[k]), extent*(1/ft), 0],origin='upper')
                     plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
                     plt.axis('normal'); plt.axis('tight')

                     axR=plt.subplot(1,2,2); 
                     axR.yaxis.tick_right()
                     axR.yaxis.set_label_position("right")
                     axR.imshow(Zstar[k],cmap='gray',extent=[min(Zdist[k]), max(Zdist[k]), extent*(1/ft), 0],origin='upper')
                     plt.plot(Zdist[k],Zbed[k]/ft,'k')
                     plt.plot(Zdist[k],Zx[k]/ft,'r')
                     plt.axis('normal'); plt.axis('tight')
                     plt.ylim(10,0)
                     plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
                     self._custom_save(sonpath,'bed_2picks'+str(k))
                     del fig

            else:
               fig = plt.figure()
               fig.subplots_adjust(wspace = 0.1, hspace=0.1)
               plt.subplot(2,2,1)
               ax = plt.gca()
               im = ax.imshow(np.flipud(data_port),cmap='gray',extent=[min(dist_m), max(dist_m), 0, extent*(1/ft)],origin='upper')
               plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')  
               plt.axis('normal'); plt.axis('tight')

               plt.subplot(2,2,3)
               ax = plt.gca()
               im = ax.imshow(data_star,cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
               plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
               plt.axis('normal'); plt.axis('tight')

               axR=plt.subplot(1,2,2); 
               axR.yaxis.tick_right()
               axR.yaxis.set_label_position("right")
               axR.imshow(data_star,cmap='gray',extent=[min(dist_m), max(dist_m), extent*(1/ft), 0],origin='upper')
               plt.plot(dist_mi,bed/ft,'k')
               plt.plot(dist_mi,x/ft,'r')
               plt.axis('normal'); plt.axis('tight')
               plt.ylim(np.max((np.max(bed/ft),np.max(x/ft)))+1,0)
               plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
               self._custom_save(sonpath,'bed_2picks')
               del fig

            if 'Zstar' and 'Z' and 'Zx' and 'Zbed' and 'Zdist' in locals():
               del Zstar, Z, Zx, Zbed, Zdist

            del dist_mi

         # 'real' bed is estimated to be the minimum of the two
         #bed = np.max(np.vstack((bed,np.squeeze(x))),axis=0) 
         bed = np.min(np.vstack((bed,np.squeeze(x))),axis=0) 
         del x


      else: #manual

         nx, ny = np.shape(data_port)
         if ny>10000:
            Z,inds = humutils.sliding_window(data_port,(nx,10000))
            del inds

         if len(Z) != len(data_port):
            beds=[]
            for k in xrange(len(Z)):
               raw_input("Bed picking "+str(k)+" of "+str(len(Z))+", are you ready? 60 seconds. Press Enter to continue...")
               bed={}
               fig = plt.figure()
               ax = plt.gca()
               im = ax.imshow(Z[k], cmap = 'gray', origin = 'upper')
               pts1 = plt.ginput(n=300, timeout=60) # it will wait for 200 clicks or 60 seconds
               x1=map(lambda x: x[0],pts1) # map applies the function passed as 
               y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
               bed = np.interp(np.r_[:10000],x1,y1)
               plt.close()
               del fig
               beds.append(bed)
            # last chunk
            raw_input("Bed picking "+str(len(Z))+" of "+str(len(Z))+", are you ready? 60 seconds. Press Enter to continue...")
            fig = plt.figure()
            ax = plt.gca()
            im = ax.imshow(data_port[:,-(np.shape(data_port)[1]-len(Z)*10000):], cmap = 'gray', origin = 'upper')
            pts1 = plt.ginput(n=300, timeout=60) # it will wait for 200 clicks or 60 seconds
            x1=map(lambda x: x[0],pts1) # map applies the function passed as 
            y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
            bed = np.interp(np.r_[:np.shape(data_port)[1]-len(Z)*10000],x1,y1)
            plt.close()
            del fig
            beds.append(bed)
            del Z

            bed = np.hstack(beds)

         else:
            raw_input("Bed picking, are you ready? 60 seconds. Press Enter to continue...")
            bed={}
            fig = plt.figure()
            ax = plt.gca()
            im = ax.imshow(data_port, cmap = 'gray', origin = 'upper')
            pts1 = plt.ginput(n=300, timeout=60) # it will wait for 200 clicks or 60 seconds
            x1=map(lambda x: x[0],pts1) # map applies the function passed as 
            y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
            bed = np.interp(np.r_[:np.shape(data_port)[1]],x1,y1)
            plt.close()
            del fig

      # now revise the depth in metres
      dep_m = (1/ft)*bed

      if doplot==1:
         nx, ny = np.shape(data_star)
         # make a plot of the starboard side with the bed pick
         fig = plt.figure()
         plt.subplot(2,2,1)
         plt.imshow(data_star[:,:min(ny,10000)],cmap='gray')
         plt.plot(bed[:min(ny,10000)],'r')
         plt.axis('normal'); plt.axis('tight')
         #plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')

         plt.subplot(2,2,3)
         plt.imshow(data_port[:,:min(ny,10000)],cmap='gray')
         plt.plot(bed[:min(ny,10000)],'r')
         plt.axis('normal'); plt.axis('tight')
         #plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')

         self._custom_save(sonpath,'bed_pick')
         del fig

      metadat['dist_m'] = dist_m
      metadat['dep_m'] = dep_m
      metadat['pix_m'] = 1/ft
      metadat['bed'] = bed

      metadat['c'] = c
      metadat['t'] = t
      metadat['f'] = f

      savemat(sonpath+base+'meta.mat', metadat ,oned_as='row')

      # read meta-data back in and append new variables
      #metadat = loadmat(sonpath+base+'meta.mat')
      #metadat['es'] = es; del es
      #metadat['ns'] = ns; del ns

      #savemat(sonpath+base+'meta.mat', metadat ,oned_as='row')
      #del metadat

      heading = np.squeeze(loadmat(sonpath+base+'meta.mat')['heading'])

      f = open(sonpath+base+'rawdat.csv', 'wt')
      writer = csv.writer(f)
      writer.writerow( ('longitude', 'latitude', 'easting', 'northing', 'depth (m)', 'distance (m)', 'heading (deg.)' ) )
      for i in range(0, len(lon)):
         writer.writerow(( float(lon[i]),float(lat[i]),float(es[i]),float(ns[i]),float(dep_m[i]),float(dist_m[i]), float(heading[i]) ))
      f.close()

      del heading, lat, lon, dep_m, dist_m

      if doplot==1:

         fig = plt.figure()
         fig.subplots_adjust(wspace = 0.5, hspace=0.5)
         plt.subplot(221)
         plt.plot(metadat['e'],metadat['n'],'k')
         plt.plot(es,ns,'r.')
         #plt.plot(esi,nsi,'b.')
         plt.xlabel('Easting (m)')
         plt.ylabel('Northing (m)')
         plt.setp(plt.xticks()[1], rotation=30)
         plt.axis('normal'); plt.axis('tight')
         #custom_save(sonpath,'raw_filt_pos_en')
         #del fig

         plt.subplot(222)
         plt.plot(metadat['lon'],metadat['lat'],'k')
         plt.xlabel('Longitude')
         plt.ylabel('Latitude')
         plt.axis('normal'); plt.axis('tight')
         plt.setp(plt.xticks()[1], rotation=30)
         self._custom_save(sonpath,'raw_filt_pos')
         del fig

         if data_port!='':
            if tonemap==1:
               import ppdrc
               # phase preserving dynamic range compression
               dat = ppdrc.ppdrc(data_port.astype('float64'), 768)
               data_port = dat.getdata()
               del dat

            # port
            nx, ny = np.shape(data_port)
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

         if data_star!='':
            if tonemap==1:
               import ppdrc
               dat = ppdrc.ppdrc(data_star.astype('float64'), 768)
               data_star = dat.getdata()
               del dat

            # starboard
            nx, ny = np.shape(data_star)
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

         if data_dwnlow!='':
  
            if tonemap==1:
               import ppdrc
               dat = ppdrc.ppdrc(data_dwnlow.astype('float64'), 768)
               data_dwnlow = dat.getdata()
               del dat

            # dwnlow
            nx, ny = np.shape(data_dwnlow)
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

         if data_dwnhi!='':
  
            if tonemap==1:
               import ppdrc
               dat = ppdrc.ppdrc(data_dwnhi.astype('float64'), 768)
               data_dwnhi = dat.getdata()
               del dat

            # dwnhigh
            nx, ny = np.shape(data_dwnhi)
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
   def _distBetweenPoints(self, pos1_lat, pos2_lat, pos1_lon, pos2_lon):
      return 6378137.0 * 2.0 * np.arcsin(np.sqrt(np.power(np.sin((np.deg2rad(pos1_lat) - np.deg2rad(pos2_lat)) / 2.0), 2.0) + np.cos(np.deg2rad(pos1_lat)) * np.cos(np.deg2rad(pos2_lat)) * np.power(np.sin((np.deg2rad(pos1_lon) - np.deg2rad(pos2_lon)) / 2.0), 2.0)))

   # =========================================================
   def _custom_save(self, figdirec,root):
      try:
         plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
      except:
         plt.savefig(os.getcwd()+os.sep+root,bbox_inches='tight',dpi=400)



