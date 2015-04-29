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
#                        __
#   ________  ____ _____/ /
#  / ___/ _ \/ __ `/ __  / 
# / /  /  __/ /_/ / /_/ /  
#/_/   \___/\__,_/\__,_/   
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

#operational
import glob, sys, getopt
from scipy.io import savemat, loadmat
import os, time
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory
except:
   pass
import csv
from fractions import gcd
from joblib import Parallel, delayed, cpu_count

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
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
import simplekml

__all__ = [
    'read',
    'custom_save',
    'distBetweenPoints',
    'makechunks',
    'plot_2bedpicks',
    'plot_bedpick',
    ]

#################################################
def read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunksize):

    '''
    Read a .DAT and associated set of .SON files recorded by a Humminbird(R)
    instrument. 
    
    Parse the data into a set of memory mapped files that will
    subsequently be used by the other functions of the PyHum module. 
    
    Export time-series data and metadata in other formats. 
    
    Create a kml file for visualising boat track

    Syntax
    ----------
    [] = PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunksize)

    Parameters
    ------------
    humfile : str
       path to the .DAT file
    sonpath : str
       path where the *.SON files are
    cs2cs_args : int, *optional* [Default="epsg:26949"]
       arguments to create coordinates in a projected coordinate system
       this argument gets given to pyproj to turn wgs84 (lat/lon) coordinates
       into any projection supported by the proj.4 libraries
    c : float, *optional* [Default=1450.0]
       speed of sound in water (m/s). Defaults to a value of freshwater
    draft : float, *optional* [Default=0]
       draft from water surface to transducer face (m)
    doplot : float, *optional* [Default=1]
       if 1, plots will be made
    t : float, *optional* [Default=0.108]
       length of transducer array (m).
       Default value is that of the 998 series Humminbird(R)
    f : float, *optional* [Default=455]
       frequency of sidescan transducer in kHz
    bedpick : int, *optional* [Default=1]
       if 1, bedpicking with be carried out automatically
       if 0, user will be prompted to pick the bed location on screen
    flip_lr : int, *optional* [Default=0]
       if 1, port and starboard scans will be flipped
       (for situations where the transducer is flipped 180 degrees)
    chunksize : int, *optional* [Default=0]
       if not 0, the data will be parsed into 'chunks' of data which
       are 'chunksize' scans long. A scan is a ping, or the simultaneous
       acquisition of a port and starboard scan. A typical value to keep 
       data chunks a manageable (small) size is 10,000 - 50,000
     
    Returns
    ---------
    sonpath+base+'_data_port.dat': memory-mapped file
        contains the raw echogram from the port side
        sidescan sonar (where present)

    sonpath+base+'_data_port.dat': memory-mapped file
        contains the raw echogram from the starboard side
        sidescan sonar (where present)

    sonpath+base+'_data_dwnhi.dat': memory-mapped file
        contains the raw echogram from the high-frequency
        echosounder (where present)

    sonpath+base+'_data_dwnlow.dat': memory-mapped file
        contains the raw echogram from the low-frequency
        echosounder (where present)
        
    sonpath+base+"trackline.kml": google-earth kml file
        contains the trackline of the vessel during data
        acquisition
     
    sonpath+base+'rawdat.csv': comma separated value file
        contains time-series data. columns corresponding to
        longitude
        latitude
        easting (m)
        northing (m)
        depth to bed (m)
        alongtrack cumulative distance (m)
        vessel heading (deg.)
     
    sonpath+base+'meta.mat': .mat file
        matlab format file containing a dictionary object
        holding metadata information. Fields are:
        e : ndarray, easting (m)
        n : ndarray, northing (m)
        es : ndarray, low-pass filtered easting (m)
        ns : ndarray, low-pass filtered northing (m)
        lat : ndarray, latitude
        lon : ndarray, longitude
        shape_port : tuple, shape of port scans in memory mapped file
        shape_star : tuple, shape of starboard scans in memory mapped file
        shape_hi : tuple, shape of high-freq. scans in memory mapped file
        shape_low : tuple, shape of low-freq. scans in memory mapped file
        dep_m : ndarray, depth to bed (m)
        dist_m : ndarray, distance along track (m)
        heading : ndarray, heading of vessel (deg. N)
        pix_m: float, size of 1 pixel in across-track dimension (m)
        bed : ndarray, depth to bed (m)
        c : float, speed of sound in water (m/s)
        t : length of sidescan transducer array (m)
        f : frequency of sidescan sound (kHz)
        spd : ndarray, vessel speed (m/s)
        time_s : ndarray, time elapsed (s)
        caltime : ndarray, unix epoch time (s)
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

    if chunksize:
      chunksize = int(chunksize)
      if chunksize==0:
         print "Chunk size will be determined automatically"
      else:
         print 'Chunk size: %s' % (str(chunksize))


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
      
    if not chunksize:
      chunksize = 0
      print '[Default] Chunk size will be determined automatically'
      

    ## for debugging
    #humfile = r"test.DAT"; sonpath = "test_data"
    #cs2cs_args = "epsg:26949"; doplot = 1; draft = 0
    #chunksize=0; c=1450; bedpick=1; fliplr=1

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

    data = pyread.pyread(sonfiles, humfile, c, headbytes, cs2cs_args)

    dat = data.gethumdat() 
    metadat = data.getmetadata()

    try:
       if flip_lr==0:
          data_port = data.getportscans().astype('int16')
       else:
          data_port = data.getstarscans().astype('int16')

       Zt, ind_port = makechunks(data_port, chunksize)
       del data_port

       # create memory mapped file for Z
       fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_port = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       port_fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='r', shape=shape_port)

    except:
       data_port = ''
       print "portside scan not available"

    try:
       if flip_lr==0:
          data_star = data.getstarscans().astype('int16')
       else:
          data_star = data.getportscans().astype('int16')

       Zt, ind_star = makechunks(data_star, chunksize)
       del data_star

       # create memory mapped file for Z
       fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_star = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       star_fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='r', shape=shape_star)

    except:
       data_star = ''
       print "starboardside scan not available"

    try:
       data_dwnlow = data.getlowscans().astype('int16')


       if chunksize != 0:
          Zt, ind_low = makechunks(data_dwnlow, chunksize/2)
       else:
          Zt, ind_low = makechunks(data_dwnlow, chunksize)
       del data_dwnlow

       # create memory mapped file for Z
       fp = np.memmap(sonpath+base+'_data_dwnlow.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_low = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       dwnlow_fp = np.memmap(sonpath+base+'_data_dwnlow.dat', dtype='int16', mode='r', shape=shape_low)

    except:
       data_dwnlow = ''
       print "low-freq. scan not available"

    try:
       data_dwnhi = data.gethiscans().astype('int16')

       if chunksize != 0:
          Zt, ind_hi = makechunks(data_dwnhi, chunksize/2)
       else:
          Zt, ind_hi = makechunks(data_dwnhi, chunksize)
       del data_dwnhi

       # create memory mapped file for Z
       fp = np.memmap(sonpath+base+'_data_dwnhi.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_hi = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       dwnhi_fp = np.memmap(sonpath+base+'_data_dwnhi.dat', dtype='int16', mode='r', shape=shape_hi)

    except:
       data_dwnhi = ''
       print "high-freq. scan not available"

    del data

    if chunksize!=0:
       nrec = len(dwnhi_fp)*chunksize
    else:
       nrec = len(metadat['n'])

    try:
       es = humutils.runningMeanFast(metadat['e'][:nrec],len(metadat['e'][:nrec])/100)
       ns = humutils.runningMeanFast(metadat['n'][:nrec],len(metadat['n'][:nrec])/100)
    except:
       es = metadat['e'][:nrec]
       ns = metadat['n'][:nrec]
    
    metadat['es'] = es
    metadat['ns'] = ns

    trans =  pyproj.Proj(init=cs2cs_args)
    lon, lat = trans(es, ns, inverse=True)
    metadat['lon'] = lon
    metadat['lat'] = lat

    if 'shape_port' in locals():
       metadat['shape_port'] = shape_port
    else:
       metadat['shape_port'] = ''   

    if 'shape_star' in locals():
       metadat['shape_star'] = shape_star
    else:
       metadat['shape_star'] = ''   

    if 'shape_hi' in locals():
       metadat['shape_hi'] = shape_hi
    else:
       metadat['shape_hi'] = ''   

    if 'shape_low' in locals():
       metadat['shape_low'] = shape_low
    else:
       metadat['shape_low'] = ''   

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

    dist = np.zeros(len(lat))
    for k in xrange(len(lat)-1):
       dist[k] = distBetweenPoints(lat[k], lat[k+1], lon[k], lon[k+1])

    dist_m = np.cumsum(dist)

    # theta at 3dB in the horizontal
    theta3dB = np.arcsin(c/(t*(f*1000)))
    #resolution of 1 sidescan pixel to nadir
    ft = (np.pi/2)*(1/theta3dB)

    dep_m = np.squeeze(metadat['dep_m'][:nrec]) #loadmat(sonpath+base+'meta.mat')['dep_m'])
    dep_m = humutils.rm_spikes(dep_m,2)

    metadat['dist_m'] = dist_m

    if bedpick == 1: # auto

       # get bed from depth trace
       bed = ft*dep_m

       imu = []
       for k in xrange(len(port_fp)):
          imu.append(port_fp[k][int(np.min(bed)):int(np.max(bed)),:])
       imu = np.hstack(imu)

       ## narrow image to within range of estimated bed
       #imu = data_port[int(np.min(bed)):int(np.max(bed)),:]
       # use dynamic boundary tracing to get 2nd estimate of bed  
       x = np.squeeze(int(np.min(bed))+humutils.dpboundary(-imu.T))
       del imu 

       if doplot==1:
          for k in xrange(len(star_fp)):
             plot_2bedpicks(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k)

          # treats each chunk in parallel for speed
#          try:
#             d = Parallel(n_jobs = min(cpu_count(),len(star_fp)), verbose=0)(delayed(plot_2bedpicks)(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k) for k in xrange(len(star_fp)))
#          except:
#             print "memory error: trying serial"
#             d = Parallel(n_jobs = 1, verbose=0)(delayed(plot_2bedpicks)(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k) for k in xrange(len(star_fp)))

       # 'real' bed is estimated to be the minimum of the two
       #bed = np.max(np.vstack((bed,np.squeeze(x))),axis=0) 
       bed = np.min(np.vstack((bed,np.squeeze(x[:len(bed)]))),axis=0) 
       del x

    else: #manual

       beds=[]
       for k in xrange(len(port_fp)):
          raw_input("Bed picking "+str(k+1)+" of "+str(len(port_fp))+", are you ready? 30 seconds. Press Enter to continue...")
          bed={}
          fig = plt.figure()
          ax = plt.gca()
          im = ax.imshow(port_fp[k], cmap = 'gray', origin = 'upper')
          pts1 = plt.ginput(n=300, timeout=30) # it will wait for 200 clicks or 60 seconds
          x1=map(lambda x: x[0],pts1) # map applies the function passed as 
          y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
          bed = np.interp(np.r_[:ind_port[-1]],x1,y1)
          plt.close()
          del fig
          beds.append(bed)
          extent = np.shape(port_fp[k])[0]
       bed = np.asarray(np.hstack(beds),'float')

    # now revise the depth in metres
    dep_m = (1/ft)*bed

    if doplot==1:

       for k in xrange(len(star_fp)):
          plot_bedpick(port_fp[k], star_fp[k], (1/ft)*bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k)

#       # treats each chunk in parallel for speed
#       try:
#          d = Parallel(n_jobs = min(cpu_count(),len(star_fp)), verbose=0)(delayed(plot_bedpick)(port_fp[k], star_fp[k], (1/ft)*bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k) for k in xrange(len(star_fp)))
#       except:
#          print "memory error: trying serial"
#          d = Parallel(n_jobs = 1, verbose=0)(delayed(plot_bedpick)(port_fp[k], star_fp[k], (1/ft)*bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k) for k in xrange(len(star_fp)))

    #heading = np.squeeze(loadmat(sonpath+base+'meta.mat')['heading'])[:nrec]
    metadat['heading'] = metadat['heading'][:nrec]
    metadat['dist_m'] = dist_m[:nrec]
    metadat['dep_m'] = dep_m[:nrec]
    metadat['pix_m'] = 1/ft
    metadat['bed'] = bed[:nrec]
    metadat['c'] = c
    metadat['t'] = t
    metadat['f'] = f

    metadat['spd'] = metadat['spd'][:nrec]
    metadat['time_s'] = metadat['time_s'][:nrec]
    metadat['e'] = metadat['e'][:nrec]
    metadat['n'] = metadat['n'][:nrec]
    metadat['caltime'] = metadat['caltime'][:nrec]

    savemat(sonpath+base+'meta.mat', metadat ,oned_as='row')

    f = open(sonpath+base+'rawdat.csv', 'wt')
    writer = csv.writer(f)
    writer.writerow( ('longitude', 'latitude', 'easting', 'northing', 'depth (m)', 'distance (m)', 'heading (deg.)' ) )
    for i in range(0, len(lon)):
       writer.writerow(( float(lon[i]),float(lat[i]),float(es[i]),float(ns[i]),float(dep_m[i]),float(dist_m[i]), float(metadat['heading'][i]) ))
    f.close()

    del lat, lon, dep_m #, dist_m

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
       custom_save(sonpath,'raw_filt_pos')
       plt.close(); del fig

       if 'dwnlow_fp' in locals():

          for k in xrange(len(dwnlow_fp)):
             fig = plt.figure()
             plt.imshow(dwnlow_fp[k],cmap='gray')
             plt.axis('normal'); plt.axis('tight')
             plt.xlabel('Ping Number (Time)')
             plt.ylabel('Range (Distance)')

             custom_save(sonpath,'raw_dwnlow'+str(k))
             plt.close(); del fig

       if 'dwnhi_fp' in locals():

          for k in xrange(len(dwnhi_fp)):
             fig = plt.figure()
             plt.imshow(dwnhi_fp[k],cmap='gray')
             plt.axis('normal'); plt.axis('tight')
             plt.xlabel('Ping Number (Time)')
             plt.ylabel('Range (Distance)')

             custom_save(sonpath,'raw_dwnhi'+str(k))
             plt.close(); del fig

    if os.name=='posix': # true if linux/mac
       elapsed = (time.time() - start)
    else: # windows
       elapsed = (time.clock() - start)
    print "Processing took ", elapsed , "seconds to analyse"

    print "Done!"


# =========================================================
def custom_save(figdirec,root):
    plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)

# =========================================================
def distBetweenPoints(pos1_lat, pos2_lat, pos1_lon, pos2_lon):
   return 6378137.0 * 2.0 * np.arcsin(np.sqrt(np.power(np.sin((np.deg2rad(pos1_lat) - np.deg2rad(pos2_lat)) / 2.0), 2.0) + np.cos(np.deg2rad(pos1_lat)) * np.cos(np.deg2rad(pos2_lat)) * np.power(np.sin((np.deg2rad(pos1_lon) - np.deg2rad(pos2_lon)) / 2.0), 2.0)))

# =========================================================
def makechunks(dat, chunksize=0):
   Ny, Nx = np.shape(dat)
   
   if chunksize==0:
      # get optimal number of slices
      if Nx%2==0:
         H = []
         for k in xrange(2,50):
            H.append(gcd(Nx,Nx/k))
         hslice = np.max(H)
      else:
         dat = np.hstack( (dat,np.ones((Ny,1))) )
         Ny, Nx = np.shape(dat)
         H = []
         for k in xrange(2,50):
            H.append(gcd(Nx,Nx/k))
         hslice = np.max(H)
   
      # get windowed data
      Zt,ind = humutils.sliding_window(dat,(Ny,hslice))

   else:
      if chunksize>np.shape(data_port)[1]:
         Zt,ind = humutils.sliding_window(dat,(Ny,chunksize))
      else:
         print "Error: chunk size is larger than number of scan lines. Please choose smaller chunk size ... exiting"

   return Zt, ind


# =========================================================
def plot_2bedpicks(dat_port, dat_star, Zbed, Zdist, Zx, ft, shape_port, sonpath, k):

   extent = shape_port[1] #np.shape(merge)[0]

   fig = plt.figure()
   fig.subplots_adjust(wspace = 0.1, hspace=0.1)
   plt.subplot(2,2,1)
   ax = plt.gca()
   im = ax.imshow(np.flipud(dat_port),cmap='gray',extent=[min(Zdist), max(Zdist), 0, extent*(1/ft)],origin='upper')
   plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')  
   plt.axis('normal'); plt.axis('tight')

   plt.subplot(2,2,3)
   ax = plt.gca()
   im = ax.imshow(dat_star,cmap='gray',extent=[min(Zdist), max(Zdist), extent*(1/ft), 0],origin='upper')
   plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
   plt.axis('normal'); plt.axis('tight')

   axR=plt.subplot(1,2,2); 
   axR.yaxis.tick_right()
   axR.yaxis.set_label_position("right")
   axR.imshow(dat_star,cmap='gray',extent=[min(Zdist), max(Zdist), extent*(1/ft), 0],origin='upper')
   plt.plot(Zdist,Zbed/ft,'k')
   plt.plot(Zdist,Zx[:len(Zdist)]/ft,'r')
   plt.axis('normal'); plt.axis('tight')
   plt.ylim(10,0)
   plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
   custom_save(sonpath,'bed_2picks'+str(k))
   plt.close(); del fig

# =========================================================
def plot_bedpick(dat_port, dat_star, Zbed, Zdist, ft, shape_port, sonpath, k):

   extent = shape_port[1] #np.shape(merge)[0]

   fig = plt.figure()
   plt.subplot(2,2,1)
   plt.imshow(np.flipud(dat_star),cmap='gray', extent=[min(Zdist), max(Zdist), 0, extent*(1/ft)], origin='upper')
   plt.plot(np.linspace(min(Zdist), max(Zdist),len(Zbed)), Zbed,'r')
   plt.axis('normal'); plt.axis('tight')
   plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')

   plt.subplot(2,2,3)
   plt.imshow(dat_port,cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0], origin='upper')
   plt.plot(np.linspace(min(Zdist), max(Zdist),len(Zbed)), Zbed,'r')
   plt.axis('normal'); plt.axis('tight')
   plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')

   custom_save(sonpath,'bed_pick'+str(k))
   plt.close(); del fig



