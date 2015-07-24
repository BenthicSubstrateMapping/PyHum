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
#from joblib import Parallel, delayed, cpu_count

#numerical
import pyread
import PyHum.utils as humutils
#from skimage.measure import LineModel, ransac
import numpy as np
import pyproj

#import ppdrc
from scipy.ndimage.filters import median_filter

#plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)
import simplekml

import warnings
warnings.filterwarnings("ignore")


__all__ = [
    'read',
    'custom_save',
    'distBetweenPoints',
    'makechunks',
    'plot_2bedpicks',
    'plot_bedpick',
    ]

#################################################
def read(humfile, sonpath, cs2cs_args="epsg:26949", c=1450.0, draft=0.3, doplot=1, t=0.108, f=455, bedpick=1, flip_lr=0, chunksize=0, model=998, calc_bearing = 0, filt_bearing = 0, cog = 1):

    '''
    Read a .DAT and associated set of .SON files recorded by a Humminbird(R)
    instrument. 
    
    Parse the data into a set of memory mapped files that will
    subsequently be used by the other functions of the PyHum module. 
    
    Export time-series data and metadata in other formats. 
    
    Create a kml file for visualising boat track

    Syntax
    ----------
    [] = PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunksize, model, calc_bearing, filt_bearing, cog)

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
    draft : float, *optional* [Default=0.3]
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
    model: int, *optional* [Default=998]
       A 3 or 4 number code indicating the model number 
       Examples: 998, 997, 1198, 1199
    cog : int, *optional* [Default=1]
       if 1, heading calculated assuming GPS course-over-ground rather than
       using a compass
    calc_bearing : float, *optional* [Default=0]
       if 1, bearing will be calculated from coordinates
    filt_bearing : float, *optional* [Default=0]
       if 1, bearing will be filtered
            
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

    if model:
       model = int(model)
       print "Data is from the %s series"  % (str(model))

    if cog:
       cog = int(cog)
       if cog==1:
          print "Heading based on course-over-ground" 

    if calc_bearing:
       calc_bearing = int(calc_bearing)
       if calc_bearing==1:
          print "Bearing will be calculated from coordinates"     
 
    if filt_bearing:
       filt_bearing = int(filt_bearing)
       if filt_bearing==1:
          print "Bearing will be filtered"    

    ## for debugging
    #humfile = r"test.DAT"; sonpath = "test_data"
    #cs2cs_args = "epsg:26949"; doplot = 1; draft = 0
    #chunksize=0; c=1450; bedpick=1; fliplr=1


    try:
       from mpl_toolkits.basemap import Basemap
       m = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], 
          resolution = 'i', llcrnrlon=10, llcrnrlat=10, urcrnrlon=30, urcrnrlat=30)
       del m
    except:
       print "Error: the epsg code you have chosen is not compatible with Basemap"
       print "please choose a different epsg code (http://spatialreference.org/)"
       print "program will now close"
       sys.exit()


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
    if base.find('_')>-1:
       base = base[:base.find('_')]

    if base.find('-')>-1:
       base = base[:base.find('-')]

    if base.find(' ')>-1:
       base = base[:base.find(' ')]

    # get the SON files from this directory
    sonfiles = glob.glob(sonpath+'*.SON')
    if not sonfiles:
        sonfiles = glob.glob(os.getcwd()+os.sep+sonpath+'*.SON')

    print "WARNING: Because files have to be read in byte by byte,"
    print "this could take a very long time ..."

    data = pyread.pyread(sonfiles, humfile, c, model, cs2cs_args)

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
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_port = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), 'r') as ff:
          port_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_port)

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
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_star = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), 'r') as ff:
          star_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_star)

    except:
       data_star = ''
       print "starboardside scan not available"


    if 'star_fp' in locals() and 'port_fp' in locals():
       # check that port and starboard are same size
       # and trim if not
       if np.shape(star_fp)!=np.shape(port_fp):
          if np.shape(port_fp[0])[1] > np.shape(star_fp[0])[1]:
             tmp = port_fp.copy()
             tmp2 = np.empty_like(star_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(star_fp[k])[1]]
             del tmp

             #port_fp.flush()
             #del port_fp
             #if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat'))):
             #   os.remove(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')))
                
             # create memory mapped file for Z
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_port2.dat')), 'w+') as ff:
                fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(tmp2))
             fp[:] = tmp2[:]
             del fp
             shape_port = np.shape(tmp2)
             del tmp2
             #we are only going to access the portion of memory required
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_port2.dat')), 'r') as ff:
                port_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_port)
             ind_port = list(ind_port)
             ind_port[-1] = np.shape(star_fp[0])[1]
             ind_port = tuple(ind_port)

          elif np.shape(port_fp[0])[1] < np.shape(star_fp[0])[1]:
             tmp = star_fp.copy()
             tmp2 = np.empty_like(port_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(port_fp[k])[1]]
             del tmp

             #star_fp.flush()
             #del star_fp
             #if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat'))):
             #   os.remove(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')))
                
             # create memory mapped file for Z
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_star2.dat')), 'w+') as ff:
                fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(tmp2))
             fp[:] = tmp2[:]
             del fp
             shape_star = np.shape(tmp2)
             del tmp2
             #we are only going to access the portion of memory required
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_star2.dat')), 'r') as ff:
                port_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_star)
             ind_star = list(ind_star)
             ind_star[-1] = np.shape(port_fp[0])[1]
             ind_star = tuple(ind_star)

    try:
       data_dwnlow = data.getlowscans().astype('int16')


       if chunksize != 0:
          Zt, ind_low = makechunks(data_dwnlow, chunksize/2)
       else:
          Zt, ind_low = makechunks(data_dwnlow, chunksize)
       del data_dwnlow

       # create memory mapped file for Z
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_low = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required      
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')), 'r') as ff:
          dwnlow_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_low)

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
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(Zt))

       fp[:] = Zt[:]
       del fp
       shape_hi = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
          dwnhi_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_hi)

    except:
       data_dwnhi = ''
       print "high-freq. scan not available"


    if 'dwnhi_fp' in locals() and 'dwnlow_fp' in locals():
       # check that low and high are same size
       # and trim if not
       if np.shape(dwnhi_fp)!=np.shape(dwnlow_fp):
          if np.shape(dwnhi_fp[0])[1] > np.shape(dwnlow_fp[0])[1]:
             tmp = dwnhi_fp.copy()
             tmp2 = np.empty_like(dwnlow_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(dwnlow_fp[k])[1]]
             del tmp

             #dwnhi_fp.flush()
             #del dwnhi_fp
             #if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat'))):
             #   os.remove(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')))
                
             # create memory mapped file for Z
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi2.dat')), 'w+') as ff:
                fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(tmp2))
             fp[:] = tmp2[:]
             del fp
             shape_dwnhi = np.shape(tmp2)
             del tmp2
             #we are only going to access the portion of memory required
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi2.dat')), 'r') as ff:
                dwnhi_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_dwnhi)
             ind_hi = list(ind_hi)
             ind_hi[-1] = np.shape(dwnlow_fp[0])[1]
             ind_hi = tuple(ind_hi)

          elif np.shape(dwnhi_fp[0])[1] < np.shape(dwnlow_fp[0])[1]:
             tmp = dwnlow_fp.copy()
             tmp2 = np.empty_like(dwnhi_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(dwnhi_fp[k])[1]]
             del tmp

             #dwnlow_fp.flush()
             #del dwnlow_fp
             #if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat'))):
             #   os.remove(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')))
                
             # create memory mapped file for Z
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow2.dat')), 'w+') as ff:
                fp = np.memmap(ff, dtype='int16', mode='w+', shape=np.shape(tmp2))
             fp[:] = tmp2[:]
             del fp
             shape_dwnlow = np.shape(tmp2)
             del tmp2
             #we are only going to access the portion of memory required
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow2.dat')), 'r') as ff:
                dwnlow_fp = np.memmap(ff, dtype='int16', mode='r', shape=shape_dwnlow)
             ind_low = list(ind_low)
             ind_low[-1] = np.shape(dwnhi_fp[0])[1]
             ind_low = tuple(ind_low)

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

    try:
       trans =  pyproj.Proj(init=cs2cs_args)
    except:
       trans =  pyproj.Proj(cs2cs_args.lstrip(), inverse=True)       

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
       #kml.save(sonpath+base+"trackline.kml")
       kml.save(os.path.normpath(os.path.join(sonpath,base+'trackline.kml')))
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
    dep_m = humutils.runningMeanFast(dep_m, 3)

    metadat['dist_m'] = dist_m

    if 'port_fp' in locals() and 'star_fp' in locals():

       if bedpick == 1: # auto

          buff = 10

          # get bed from depth trace
          bed = ft*dep_m

          imu = []

          for k in xrange(len(port_fp)):
             #imu.append(port_fp[k][int(np.min(bed)):int(np.max(bed)),:])
             imu.append(port_fp[k][np.max([0,int(np.min(bed))-buff]):int(np.max(bed))+buff,:])
          imu = np.hstack(imu)

          imu = np.asarray(imu, 'float64')

          #imu = ppdrc.ppdrc(imu, np.shape(imu)[1]/2).getdata()

          imu = median_filter(imu,(20,20))

          ## narrow image to within range of estimated bed
          #imu = data_port[int(np.min(bed)):int(np.max(bed)),:]
          # use dynamic boundary tracing to get 2nd estimate of bed  
          x = np.squeeze(int(np.min(bed))+humutils.dpboundary(-imu.T)) - buff
          #x = np.squeeze(humutils.dpboundary(-imu.T))
          del imu 

          if len(x)<len(bed):
             x = np.append(x,x[-1]*np.ones(len(bed)-len(x)))
          elif len(x)>len(bed):
             bed = np.append(bed,bed[-1]*np.ones(len(x)-len(bed)))

          # if standard deviation of auto bed pick is too small, then use acoustic bed pick
          if np.std(x)<5:
             print "stdev of auto bed pick is low, using acoustic pick"
             x = bed.copy()

          if len(dist_m)<len(bed):
             dist_m = np.append(dist_m,dist_m[-1]*np.ones(len(bed)-len(dist_m)))

          if doplot==1:
             # treats each chunk in parallel for speed
             #try:
             #   d = Parallel(n_jobs = -1, verbose=0)(delayed(plot_2bedpicks)(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k) for k in xrange(len(star_fp)))
             #except:
             for k in xrange(len(star_fp)):
                plot_2bedpicks(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k)

          # 'real' bed is estimated to be the minimum of the two
          #bed = np.max(np.vstack((bed,np.squeeze(x))),axis=0) 
          bed = np.min(np.vstack((bed[:nrec],np.squeeze(x[:nrec]))),axis=0) 
          bed = humutils.runningMeanFast(bed, 3)

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
          # treats each chunk in parallel for speed
          #try:
          #   d = Parallel(n_jobs = -1, verbose=0)(delayed(plot_bedpick)(port_fp[k], star_fp[k], (1/ft)*bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k) for k in xrange(len(star_fp)))
          #except:
          for k in xrange(len(star_fp)):
             plot_bedpick(port_fp[k], star_fp[k], (1/ft)*bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k)


       metadat['bed'] = bed[:nrec]

    else:
       metadat['bed'] = dep_m[:nrec]*ft

    #heading = np.squeeze(loadmat(sonpath+base+'meta.mat')['heading'])[:nrec]
    metadat['heading'] = metadat['heading'][:nrec]

    # over-ride measured bearing and calc from positions
    if calc_bearing==1:
       lat = np.squeeze(metadat['lat'])
       lon = np.squeeze(metadat['lon']) 

       #point-to-point bearing
       bearing = np.zeros(len(lat))
       for k in xrange(len(lat)-1):
          bearing[k] = bearingBetweenPoints(lat[k], lat[k+1], lon[k], lon[k+1])
       del lat, lon

    else:
       # reported bearing by instrument (Kalman filtered?)
       bearing = np.squeeze(metadat['heading'])

    # if stdev in heading is large, there's probably noise that needs to be filtered out
    if np.std(bearing)>180:
       print "WARNING: large heading stdev - attempting filtering"
       from sklearn.cluster import MiniBatchKMeans
       # can have two modes
       data = np.column_stack([bearing, bearing])
       k_means = MiniBatchKMeans(2)
       # fit the model
       k_means.fit(data) 
       values = k_means.cluster_centers_.squeeze()
       labels = k_means.labels_

       if np.sum(labels==0) > np.sum(labels==1):
          bearing[labels==1] = np.nan
       else:
          bearing[labels==0] = np.nan

       nans, y= humutils.nan_helper(bearing)
       bearing[nans]= np.interp(y(nans), y(~nans), bearing[~nans]) 

    if filt_bearing ==1:
       bearing = humutils.runningMeanFast(bearing, len(bearing)/100)

    if cog==1:
       theta = np.asarray(bearing, 'float')/(180/np.pi)
       #course over ground is given as a compass heading (ENU) from True north, or Magnetic north.
       #To get this into NED (North-East-Down) coordinates, you need to rotate the ENU 
       # (East-North-Up) coordinate frame. 
       #Subtract pi/2 from your heading
       theta = theta - np.pi/2
       # (re-wrap to Pi to -Pi)
       theta = np.unwrap(-theta)
       metadat['heading'] = theta * (180/np.pi)
    else:
       metadat['heading'] = bearing


    metadat['dist_m'] = dist_m[:nrec]
    metadat['dep_m'] = dep_m[:nrec]
    metadat['pix_m'] = 1/ft
    metadat['bed'] = metadat['bed'][:nrec]
    metadat['c'] = c
    metadat['t'] = t
    metadat['f'] = f

    metadat['spd'] = metadat['spd'][:nrec]
    metadat['time_s'] = metadat['time_s'][:nrec]
    metadat['e'] = metadat['e'][:nrec]
    metadat['n'] = metadat['n'][:nrec]
    metadat['caltime'] = metadat['caltime'][:nrec]

    #savemat(sonpath+base+'meta.mat', metadat ,oned_as='row')
    savemat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')), metadat ,oned_as='row')

    #f = open(sonpath+base+'rawdat.csv', 'wt')
    f = open(os.path.normpath(os.path.join(sonpath,base+'rawdat.csv')), 'wt')
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
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=400)

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
      if chunksize<np.shape(dat)[1]:
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

# =========================================================
def bearingBetweenPoints(pos1_lat, pos2_lat, pos1_lon, pos2_lon):
   lat1 = np.deg2rad(pos1_lat)
   lon1 = np.deg2rad(pos1_lon)
   lat2 = np.deg2rad(pos2_lat)
   lon2 = np.deg2rad(pos2_lon)

   bearing = np.arctan2(np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1), np.sin(lon2 - lon1) * np.cos(lat2))

   db = np.rad2deg(bearing)
   return (90.0 - db + 360.0) % 360.0
   
# =========================================================
# =========================================================
if __name__ == '__main__':

   read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunksize, model, calc_bearing, filt_bearing, cog)


       #dwnlow_fp = np.memmap(sonpath+base+'_data_dwnlow.dat', dtype='int16', mode='r', shape=shape_low)
       #dwnlow_fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')), dtype='int16', mode='r', shape=shape_low)
       #fp = np.memmap(sonpath+base+'_data_dwnlow.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       #fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')), dtype='int16', mode='w+', shape=np.shape(Zt))
             #star_fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='r', shape=shape_star)
             #star_fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), dtype='int16', mode='r', shape=shape_star)
             #fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='w+', shape=np.shape(tmp2))
             #fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), dtype='int16', mode='w+', shape=np.shape(tmp2))
             #port_fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='r', shape=shape_port)
             #port_fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), dtype='int16', mode='r', shape=shape_port)
             #fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='w+', shape=np.shape(tmp2))
             #fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), dtype='int16', mode='w+', shape=np.shape(tmp2))
       #star_fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='r', shape=shape_star)
       #star_fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), dtype='int16', mode='r', shape=shape_star)
       #fp = np.memmap(sonpath+base+'_data_star.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       #fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), dtype='int16', mode='w+', shape=np.shape(Zt))
       #port_fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='r', shape=shape_port)
       #port_fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), dtype='int16', mode='r', shape=shape_port)
       #fp = np.memmap(sonpath+base+'_data_port.dat', dtype='int16', mode='w+', shape=np.shape(Zt))
       #fp = np.memmap(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), dtype='int16', mode='w+', shape=np.shape(Zt)) 
