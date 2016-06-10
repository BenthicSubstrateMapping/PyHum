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
import glob, sys #, getopt
from scipy.io import savemat
import os, time
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory

except:
   pass
import csv
import PyHum.io as io

#numerical
#import pyread
import PyHum.pyread_single as pyread_single
from joblib import Parallel, delayed, cpu_count

import PyHum.utils as humutils
import numpy as np
import pyproj
from numpy.lib.stride_tricks import as_strided as ast


try:
   from pykdtree.kdtree import KDTree
   pykdtree=1
except:
   #print "install pykdtree for faster kd-tree operations: https://github.com/storpipfugl/pykdtree"
   from scipy.spatial import cKDTree as KDTree
   pykdtree=0


#plotting
import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

#################################################
def read(humfile, sonpath, cs2cs_args="epsg:26949", c=1450.0, draft=0.3, doplot=1, t=0.108, bedpick=1, flip_lr=0, model=998, calc_bearing = 0, filt_bearing = 0, chunk='d100'): #cog = 1,

    '''
    Read a .DAT and associated set of .SON files recorded by a Humminbird(R)
    instrument.

    Parse the data into a set of memory mapped files that will
    subsequently be used by the other functions of the PyHum module.

    Export time-series data and metadata in other formats.

    Create a kml file for visualising boat track

    Syntax
    ----------
    [] = PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, bedpick, flip_lr, chunksize, model, calc_bearing, filt_bearing, chunk)

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
    bedpick : int, *optional* [Default=1]
       if 1, bedpicking with be carried out automatically
       if 0, user will be prompted to pick the bed location on screen
    flip_lr : int, *optional* [Default=0]
       if 1, port and starboard scans will be flipped
       (for situations where the transducer is flipped 180 degrees)
    model: int, *optional* [Default=998]
       A 3 or 4 number code indicating the model number
       Examples: 998, 997, 1198, 1199
    calc_bearing : float, *optional* [Default=0]
       if 1, bearing will be calculated from coordinates
    filt_bearing : float, *optional* [Default=0]
       if 1, bearing will be filtered
    chunk : str, *optional* [Default='d100' (distance, 100 m)]
       letter, followed by a number.
       There are the following letter options:
       'd' - parse chunks based on distance, then number which is distance in m
       'p' - parse chunks based on number of pings, then number which is number of pings
       'h' - parse chunks based on change in heading, then number which is the change in heading in degrees
       '1' - process just 1 chunk

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

    if bedpick:
      bedpick = np.asarray(bedpick,int)
      if bedpick==1:
         print 'Bed picking is auto'
      elif bedpick==0:
         print 'Bed picking is manual'
      else:
         print 'User will be prompted per chunk about bed picking method'

    if chunk:
       chunk = str(chunk)
       if chunk[0]=='d':
          chunkmode=1
          chunkval = int(chunk[1:])
          print 'Chunks based on distance of %s m' % (str(chunkval))
       elif chunk[0]=='p':
          chunkmode=2
          chunkval = int(chunk[1:])
          print 'Chunks based on %s pings' % (str(chunkval))
       elif chunk[0]=='h':
          chunkmode=3
          chunkval = int(chunk[1:])
          print 'Chunks based on heading devation of %s degrees' % (str(chunkval))
       elif chunk[0]=='1':
          chunkmode=4
          chunkval = 1
          print 'Only 1 chunk will be produced'
       else:
          print "Chunk mode not understood - should be 'd', 'p', or 'h' - using defaults"
          chunkmode=1
          chunkval = 100
          print 'Chunks based on distance of %s m' % (str(chunkval))

    if model:
       model = int(model)
       print "Data is from the %s series"  % (str(model))

#    if cog:
#       cog = int(cog)
#       if cog==1:
#          print "Heading based on course-over-ground"

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
    #c=1450; bedpick=1; fliplr=1; chunk = 'd100'
    #model=998; cog=1; calc_bearing=0; filt_bearing=0

    f = 455

    try:
       print "Checking the epsg code you have chosen for compatibility with Basemap ... "
       from mpl_toolkits.basemap import Basemap
       m = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1],
          resolution = 'i', llcrnrlon=10, llcrnrlat=10, urcrnrlon=30, urcrnrlat=30)
       del m
       print "... epsg code compatible"
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

    # get the SON files from this directory
    sonfiles = glob.glob(sonpath+'*.SON')
    if not sonfiles:
        sonfiles = glob.glob(os.getcwd()+os.sep+sonpath+'*.SON')

    base = humfile.split('.DAT') # get base of file name for output
    base = base[0].split(os.sep)[-1]

    # remove underscores, negatives and spaces from basename
    base = humutils.strip_base(base)

    print "WARNING: Because files have to be read in byte by byte,"
    print "this could take a very long time ..."

    #reading each sonfile in parallel should be faster ...
    try:
       o = Parallel(n_jobs = np.min([len(sonfiles), cpu_count()]), verbose=0)(delayed(getscans)(sonfiles[k], humfile, c, model, cs2cs_args) for k in xrange(len(sonfiles)))
       X, Y, A, B = zip(*o)

       for k in xrange(len(Y)):
          if Y[k] == 'sidescan_port':
             dat = A[k] #data.gethumdat()
             metadat = B[k] #data.getmetadata()
             if flip_lr==0:
                data_port = X[k].astype('int16')
             else:
                data_star = X[k].astype('int16')

          elif Y[k] == 'sidescan_starboard':
             if flip_lr==0:
                data_star = X[k].astype('int16')
             else:
                data_port = X[k].astype('int16')

          elif Y[k] == 'down_lowfreq':
             data_dwnlow = X[k].astype('int16')

          elif Y[k] == 'down_highfreq':
             data_dwnhi = X[k].astype('int16')

       del X, Y, A, B, o
       old_pyread = 0

       if 'data_port' not in locals():
          data_port = ''
          print "portside scan not available"

       if 'data_star' not in locals():
          data_star = ''
          print "starboardside scan not available"

       if 'data_dwnhi' not in locals():
          data_dwnlow = ''
          print "high-frq. downward scan not available"

       if 'data_dwnlow' not in locals():
          data_dwnlow = ''
          print "low-frq. downward scan not available"

    except: # revert back to older version if paralleleised version fails

       print "something went wrong with the parallelised version of pyread ..."

       import pyread
       data = pyread.pyread(sonfiles, humfile, c, model, cs2cs_args)

       dat = data.gethumdat()

       metadat = data.getmetadata()

       old_pyread = 1

    nrec = len(metadat['n'])

    metadat['instr_heading'] = metadat['heading'][:nrec]

    #metadat['heading'] = humutils.get_bearing(calc_bearing, filt_bearing, cog, metadat['lat'], metadat['lon'], metadat['instr_heading'])

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

    metadat['heading'] = humutils.get_bearing(calc_bearing, filt_bearing, metadat['lat'], metadat['lon'], metadat['instr_heading']) #cog

    dist_m = humutils.get_dist(lat, lon)
    metadat['dist_m'] = dist_m

    if calc_bearing==1: # recalculate speed, m/s
       ds=np.gradient(np.squeeze(metadat['time_s']))
       dx=np.gradient(np.squeeze(metadat['dist_m']))
       metadat['spd'] = dx[:nrec]/ds[:nrec]

    # theta at 3dB in the horizontal
    theta3dB = np.arcsin(c/(t*(f*1000)))
    #resolution of 1 sidescan pixel to nadir
    ft = (np.pi/2)*(1/theta3dB) #/ (f/455)

    dep_m = humutils.get_depth(metadat['dep_m'][:nrec])

    if old_pyread == 1: #older pyread version

       # port scan
       try:
          if flip_lr==0:
             data_port = data.getportscans().astype('int16')
          else:
             data_port = data.getstarscans().astype('int16')
       except:
          data_port = ''
          print "portside scan not available"

    if data_port!='':

       Zt, ind_port = makechunks_scan(chunkmode, chunkval, metadat, data_port, 0)
       del data_port

       ## create memory mapped file for Z
       shape_port = io.set_mmap_data(sonpath, base, '_data_port.dat', 'int16', Zt)

       ##we are only going to access the portion of memory required
       port_fp = io.get_mmap_data(sonpath, base, '_data_port.dat', 'int16', shape_port)

    if old_pyread == 1: #older pyread version
       # starboard scan
       try:
          if flip_lr==0:
             data_star = data.getstarscans().astype('int16')
          else:
             data_star = data.getportscans().astype('int16')
       except:
          data_star = ''
          print "starboardside scan not available"

    if data_star!='':

       Zt, ind_star = makechunks_scan(chunkmode, chunkval, metadat, data_star, 1)
       del data_star

       # create memory mapped file for Z
       shape_star = io.set_mmap_data(sonpath, base, '_data_star.dat', 'int16', Zt)

       star_fp = io.get_mmap_data(sonpath, base, '_data_star.dat', 'int16', shape_star)

    if 'star_fp' in locals() and 'port_fp' in locals():
       # check that port and starboard are same size
       # and trim if not
       if np.shape(star_fp)!=np.shape(port_fp):
          print "port and starboard scans are different sizes ... rectifying"
          if np.shape(port_fp[0])[1] > np.shape(star_fp[0])[1]:
             tmp = port_fp.copy()
             tmp2 = np.empty_like(star_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(star_fp[k])[1]]
             del tmp

             # create memory mapped file for Z
             shape_port = io.set_mmap_data(sonpath, base, '_data_port2.dat', 'int16', tmp2)
             #shape_star = shape_port.copy()
             shape_star = tuple(np.asarray(shape_port).copy())

             ##we are only going to access the portion of memory required
             port_fp = io.get_mmap_data(sonpath, base, '_data_port2.dat', 'int16', shape_port)

             ind_port = list(ind_port)
             ind_port[-1] = np.shape(star_fp[0])[1]
             ind_port = tuple(ind_port)

          elif np.shape(port_fp[0])[1] < np.shape(star_fp[0])[1]:
             tmp = star_fp.copy()
             tmp2 = np.empty_like(port_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(port_fp[k])[1]]
             del tmp

             # create memory mapped file for Z
             shape_port = io.set_mmap_data(sonpath, base, '_data_star2.dat', 'int16', tmp2)
             #shape_star = shape_port.copy()
             shape_star = tuple(np.asarray(shape_port).copy())

             #we are only going to access the portion of memory required
             star_fp = io.get_mmap_data(sonpath, base, '_data_star2.dat', 'int16', shape_star)

             ind_star = list(ind_star)
             ind_star[-1] = np.shape(port_fp[0])[1]
             ind_star = tuple(ind_star)

    if old_pyread == 1: #older pyread version
       # low-freq. sonar
       try:
          data_dwnlow = data.getlowscans().astype('int16')
       except:
          data_dwnlow = ''
          print "low-freq. scan not available"

    if data_dwnlow!='':

       Zt, ind_low = makechunks_scan(chunkmode, chunkval, metadat, data_dwnlow, 2)
       del data_dwnlow

       # create memory mapped file for Z
       shape_low = io.set_mmap_data(sonpath, base, '_data_dwnlow.dat', 'int16', Zt)

       ##we are only going to access the portion of memory required
       dwnlow_fp = io.get_mmap_data(sonpath, base, '_data_dwnlow.dat', 'int16', shape_low)

    if old_pyread == 1: #older pyread version
       # hi-freq. sonar
       try:
          data_dwnhi = data.gethiscans().astype('int16')
       except:
          data_dwnhi = ''
          print "high-freq. scan not available"

    if data_dwnhi!='':

       Zt, ind_hi = makechunks_scan(chunkmode, chunkval, metadat, data_dwnhi, 3)
       del data_dwnhi

       # create memory mapped file for Z
       shape_hi = io.set_mmap_data(sonpath, base, '_data_dwnhi.dat', 'int16', Zt)

       dwnhi_fp = io.get_mmap_data(sonpath, base, '_data_dwnhi.dat', 'int16', shape_hi)

    if 'dwnhi_fp' in locals() and 'dwnlow_fp' in locals():
       # check that low and high are same size
       # and trim if not
       if (np.shape(dwnhi_fp)!=np.shape(dwnlow_fp)) and (chunkmode!=4):
          print "dwnhi and dwnlow are different sizes ... rectifying"
          if np.shape(dwnhi_fp[0])[1] > np.shape(dwnlow_fp[0])[1]:
             tmp = dwnhi_fp.copy()
             tmp2 = np.empty_like(dwnlow_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(dwnlow_fp[k])[1]]
             del tmp

             # create memory mapped file for Z
             shape_low = io.set_mmap_data(sonpath, base, '_data_dwnhi2.dat', 'int16', tmp2)
             #shape_hi = shape_low.copy()
             shape_hi = tuple(np.asarray(shape_low).copy())

             ##we are only going to access the portion of memory required
             dwnhi_fp = io.get_mmap_data(sonpath, base, '_data_dwnhi2.dat', 'int16', shape_hi)

             ind_hi = list(ind_hi)
             ind_hi[-1] = np.shape(dwnlow_fp[0])[1]
             ind_hi = tuple(ind_hi)

          elif np.shape(dwnhi_fp[0])[1] < np.shape(dwnlow_fp[0])[1]:
             tmp = dwnlow_fp.copy()
             tmp2 = np.empty_like(dwnhi_fp)
             for k in xrange(len(tmp)):
                 tmp2[k] = tmp[k][:,:np.shape(dwnhi_fp[k])[1]]
             del tmp

             # create memory mapped file for Z
             shape_low = io.set_mmap_data(sonpath, base, '_data_dwnlow2.dat', 'int16', tmp2)
             #shape_hi = shape_low.copy()
             shape_hi = tuple(np.asarray(shape_low).copy())

             ##we are only going to access the portion of memory required
             dwnlow_fp = io.get_mmap_data(sonpath, base, '_data_dwnlow2.dat', 'int16', shape_low)

             ind_low = list(ind_low)
             ind_low[-1] = np.shape(dwnhi_fp[0])[1]
             ind_low = tuple(ind_low)

    if old_pyread == 1: #older pyread version
       del data

    if ('shape_port' in locals()) and (chunkmode!=4):
       metadat['shape_port'] = shape_port
       nrec = metadat['shape_port'][0] * metadat['shape_port'][2]
    elif ('shape_port' in locals()) and (chunkmode==4):
       metadat['shape_port'] = shape_port
       nrec = metadat['shape_port'][1]
    else:
       metadat['shape_port'] = ''

    if ('shape_star' in locals()) and (chunkmode!=4):
       metadat['shape_star'] = shape_star
       nrec = metadat['shape_star'][0] * metadat['shape_star'][2]
    elif ('shape_star' in locals()) and (chunkmode==4):
       metadat['shape_star'] = shape_star
       nrec = metadat['shape_star'][1]
    else:
       metadat['shape_star'] = ''

    if ('shape_hi' in locals()) and (chunkmode!=4):
       metadat['shape_hi'] = shape_hi
       #nrec = metadat['shape_hi'][0] * metadat['shape_hi'][2] * 2
    elif ('shape_hi' in locals()) and (chunkmode==4):
       metadat['shape_hi'] = shape_hi
    else:
       metadat['shape_hi'] = ''

    if ('shape_low' in locals()) and (chunkmode!=4):
       metadat['shape_low'] = shape_low
       #nrec = metadat['shape_low'][0] * metadat['shape_low'][2] * 2
    elif ('shape_low' in locals()) and (chunkmode==4):
       metadat['shape_low'] = shape_low
    else:
       metadat['shape_low'] = ''

    #make kml boat trackline
    humutils.make_trackline(lon,lat, sonpath, base)

    if 'port_fp' in locals() and 'star_fp' in locals():

       #if not os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'meta.mat'))):
       if 2>1:
          if bedpick == 1: # auto

             x, bed = humutils.auto_bedpick(ft, dep_m, chunkmode, port_fp, c)

             if len(dist_m)<len(bed):
                dist_m = np.append(dist_m,dist_m[-1]*np.ones(len(bed)-len(dist_m)))

             if doplot==1:
                if chunkmode!=4:
                   for k in xrange(len(star_fp)):
                      plot_2bedpicks(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k, chunkmode)
                else:
                   plot_2bedpicks(port_fp, star_fp, bed, dist_m, x, ft, shape_port, sonpath, 0, chunkmode)

             # 'real' bed is estimated to be the minimum of the two
             bed = np.min(np.vstack((bed[:nrec],np.squeeze(x[:nrec]))),axis=0)
             bed = humutils.runningMeanFast(bed, 3)

          elif bedpick>1: # user prompt

             x, bed = humutils.auto_bedpick(ft, dep_m, chunkmode, port_fp, c)

             if len(dist_m)<len(bed):
                dist_m = np.append(dist_m,dist_m[-1]*np.ones(len(bed)-len(dist_m)))

             # 'real' bed is estimated to be the minimum of the two
             bed = np.min(np.vstack((bed[:nrec],np.squeeze(x[:nrec]))),axis=0)
             bed = humutils.runningMeanFast(bed, 3)

             # manually intervene
             fig = plt.figure()
             ax = plt.gca()
             if chunkmode !=4:
                im = ax.imshow(np.hstack(port_fp), cmap = 'gray', origin = 'upper')
             else:
                im = ax.imshow(port_fp, cmap = 'gray', origin = 'upper')
             plt.plot(bed,'r')
             plt.axis('normal'); plt.axis('tight')

             pts1 = plt.ginput(n=300, timeout=30) # it will wait for 200 clicks or 60 seconds
             x1=map(lambda x: x[0],pts1) # map applies the function passed as
             y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
             plt.close()
             del fig

             if x1 != []: # if x1 is not empty
                tree = KDTree(zip(np.arange(1,len(bed)), bed))
                try:
                   dist, inds = tree.query(zip(x1, y1), k = 100, eps=5, n_jobs=-1)
                except:
                   dist, inds = tree.query(zip(x1, y1), k = 100, eps=5)

                b = np.interp(inds,x1,y1)
                bed2 = bed.copy()
                bed2[inds] = b
                bed = bed2

             if doplot==1:
                if chunkmode!=4:
                   for k in xrange(len(star_fp)):
                      plot_2bedpicks(port_fp[k], star_fp[k], bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], x[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k, chunkmode)
                else:
                   plot_2bedpicks(port_fp, star_fp, bed, dist_m, x, ft, shape_port, sonpath, 0, chunkmode)

          else: #manual

             beds=[]

             if chunkmode!=4:
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
             else:
                raw_input("Bed picking - are you ready? 30 seconds. Press Enter to continue...")
                bed={}
                fig = plt.figure()
                ax = plt.gca()
                im = ax.imshow(port_fp, cmap = 'gray', origin = 'upper')
                pts1 = plt.ginput(n=300, timeout=30) # it will wait for 200 clicks or 60 seconds
                x1=map(lambda x: x[0],pts1) # map applies the function passed as
                y1=map(lambda x: x[1],pts1) # first parameter to each element of pts
                bed = np.interp(np.r_[:ind_port[-1]],x1,y1)
                plt.close()
                del fig
                beds.append(bed)
                extent = np.shape(port_fp)[1]
                bed = np.asarray(np.hstack(beds),'float')

          # now revise the depth in metres
          dep_m = (1/ft)*bed

          if doplot==1:
             if chunkmode!=4:
                for k in xrange(len(star_fp)):
                   plot_bedpick(port_fp[k], star_fp[k], (1/ft)*bed[ind_port[-1]*k:ind_port[-1]*(k+1)], dist_m[ind_port[-1]*k:ind_port[-1]*(k+1)], ft, shape_port, sonpath, k, chunkmode)
             else:
                plot_bedpick(port_fp, star_fp, (1/ft)*bed, dist_m, ft, shape_port, sonpath, 0, chunkmode)

          metadat['bed'] = bed[:nrec]

    else:
       metadat['bed'] = dep_m[:nrec]*ft

    metadat['heading'] = metadat['heading'][:nrec]
    metadat['lon'] = lon[:nrec]
    metadat['lat'] = lat[:nrec]
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
    metadat['es'] = metadat['es'][:nrec]
    metadat['ns'] = metadat['ns'][:nrec]
    metadat['caltime'] = metadat['caltime'][:nrec]

    savemat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')), metadat ,oned_as='row')

    f = open(os.path.normpath(os.path.join(sonpath,base+'rawdat.csv')), 'wt')
    writer = csv.writer(f)
    writer.writerow( ('longitude', 'latitude', 'easting', 'northing', 'depth (m)', 'distance (m)', 'instr. heading (deg)', 'heading (deg.)' ) )
    for i in range(0, nrec):
       writer.writerow(( float(lon[i]),float(lat[i]),float(es[i]),float(ns[i]),float(dep_m[i]),float(dist_m[i]), float(metadat['instr_heading'][i]), float(metadat['heading'][i]) ))
    f.close()

    del lat, lon, dep_m #, dist_m

    if doplot==1:

       plot_pos(sonpath, metadat, es, ns)

       if 'dwnlow_fp' in locals():

          plot_dwnlow(dwnlow_fp, chunkmode, sonpath)

       if 'dwnhi_fp' in locals():

          plot_dwnhi(dwnhi_fp, chunkmode, sonpath)

    if os.name=='posix': # true if linux/mac
       elapsed = (time.time() - start)
    else: # windows
       elapsed = (time.clock() - start)
    print "Processing took ", elapsed , "seconds to analyse"

    print "Done!"


# =========================================================
def getscans(sonfile, humfile, c, model, cs2cs_args):

   data = pyread_single.pyread(sonfile, humfile, c, model, cs2cs_args)

   a, b = data.getscan()

   if b == 'sidescan_port':
      dat = data.gethumdat()
      metadat = data.getmetadata()
   else:
      dat = None
      metadat = None

   return a, b, dat, metadat

# =========================================================
def plot_dwnhi(dwnhi_fp, chunkmode, sonpath):

    if chunkmode!=4:
       for k in xrange(len(dwnhi_fp)):
          fig = plt.figure()
          plt.imshow(dwnhi_fp[k],cmap='gray')
          plt.axis('normal'); plt.axis('tight')
          plt.xlabel('Ping Number (Time)')
          plt.ylabel('Range (Distance)')

          custom_save(sonpath,'raw_dwnhi'+str(k))
          plt.close(); del fig

    else:
       fig = plt.figure()
       plt.imshow(dwnhi_fp,cmap='gray')
       plt.axis('normal'); plt.axis('tight')
       plt.xlabel('Ping Number (Time)')
       plt.ylabel('Range (Distance)')

       custom_save(sonpath,'raw_dwnhi'+str(0))
       plt.close(); del fig

# =========================================================
def plot_dwnlow(dwnlow_fp, chunkmode, sonpath):

    if chunkmode!=4:
       for k in xrange(len(dwnlow_fp)):
          fig = plt.figure()
          plt.imshow(dwnlow_fp[k],cmap='gray')
          plt.axis('normal'); plt.axis('tight')
          plt.xlabel('Ping Number (Time)')
          plt.ylabel('Range (Distance)')

          custom_save(sonpath,'raw_dwnlow'+str(k))
          plt.close(); del fig
    else:
       fig = plt.figure()
       plt.imshow(dwnlow_fp,cmap='gray')
       plt.axis('normal'); plt.axis('tight')
       plt.xlabel('Ping Number (Time)')
       plt.ylabel('Range (Distance)')

       custom_save(sonpath,'raw_dwnlow'+str(0))
       plt.close(); del fig

# =========================================================
def plot_pos(sonpath, metadat, es, ns):

    fig = plt.figure()
    fig.subplots_adjust(wspace = 0.5, hspace=0.5)
    plt.subplot(221)
    plt.plot(metadat['e'],metadat['n'],'k')
    plt.plot(es,ns,'r.')
    plt.xlabel('Easting (m)')
    plt.ylabel('Northing (m)')
    plt.setp(plt.xticks()[1], rotation=30)
    plt.axis('normal'); plt.axis('tight')

    plt.subplot(222)
    plt.plot(metadat['lon'],metadat['lat'],'k')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.axis('normal'); plt.axis('tight')
    plt.setp(plt.xticks()[1], rotation=30)
    custom_save(sonpath,'raw_filt_pos')
    plt.close(); del fig

# =========================================================
def makechunks_scan(chunkmode, chunkval, metadat, data, flag):

    if chunkmode==1:
       nchunks = 0
       while nchunks<2:
          chunkval = chunkval-1
          tmp = metadat['dist_m']/chunkval #length_chunk
          nchunks = np.floor(tmp.max())
          del tmp
       if flag==0:
          print "port sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==1:
          print "starboard sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==2:
          print "low-freq. sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==3:
          print "high-freq. sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       chunkval = chunkval+1
       Zt, ind = makechunks_simple(data, nchunks)

    elif chunkmode==2:
       nchunks = 0
       while nchunks<2:
          chunkval = chunkval-1
          tmp = np.max(np.shape(data))/chunkval #length_chunk
          nchunks = np.floor(tmp)
          del tmp
       if flag==0:
          print "port sonar data will be parsed into %s, %s ping chunks" % (str(nchunks), str(chunkval))
       elif flag==1:
          print "starboard sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==2:
          print "low-freq. sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==3:
          print "high-freq. sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       chunkval = chunkval+1
       Zt, ind = makechunks_simple(data, nchunks)

    elif chunkmode==3:
       nchunks = 0
       while nchunks<2:
          chunkval = chunkval-1
          tmp = np.abs(metadat['heading']-metadat['heading'][0])/chunkval
          nchunks = np.floor(tmp.max())
          del tmp
       if flag==0:
          print "port sonar data will be parsed into %s, %s degree deviation chunks" % (str(nchunks), str(chunkval))
       elif flag==1:
          print "starboard sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==2:
          print "low-freq. sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       elif flag==3:
          print "high-freq. sonar data will be parsed into %s, %s m chunks" % (str(nchunks), str(chunkval))
       chunkval = chunkval+1
       Zt, ind = makechunks_simple(data, nchunks)

    elif chunkmode==4:
       Zt, ind = makechunks_simple(data, 1)

    return Zt, ind

# =========================================================
def custom_save(figdirec,root):
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=400)

# =========================================================
def makechunks_simple(dat, numchunks):
   Ny, Nx = np.shape(dat)
   # get windowed data
   try:
      return sliding_window(dat,(Ny,Nx/int(numchunks)))    #humutils.
   except:
      print "memory-mapping failed in sliding window - trying memory intensive version"
      return humutils.sliding_window_nomm(dat,(Ny,Nx/int(numchunks)))

# =========================================================
def plot_2bedpicks(dat_port, dat_star, Zbed, Zdist, Zx, ft, shape_port, sonpath, k, chunkmode):

   if chunkmode != 4:
      extent = shape_port[1] #np.shape(merge)[0]
   else:
      extent = shape_port[0]

   fig = plt.figure()
   fig.subplots_adjust(wspace = 0.1, hspace=0.1)
   plt.subplot(2,2,1)
   ax = plt.gca()
   ax.imshow(np.flipud(dat_port),cmap='gray',extent=[min(Zdist), max(Zdist), 0, extent*(1/ft)],origin='upper') #im =
   plt.ylabel('Range (m)'); #plt.xlabel('Distance along track (m)')
   plt.axis('normal'); plt.axis('tight')

   plt.subplot(2,2,3)
   ax = plt.gca()
   ax.imshow(dat_star,cmap='gray',extent=[min(Zdist), max(Zdist), extent*(1/ft), 0],origin='upper') #im =
   plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
   plt.axis('normal'); plt.axis('tight')

   axR=plt.subplot(1,2,2);
   axR.yaxis.tick_right()
   axR.yaxis.set_label_position("right")
   axR.imshow(dat_star,cmap='gray',extent=[min(Zdist), max(Zdist), extent*(1/ft), 0],origin='upper')
   plt.plot(Zdist,Zbed/ft,'k')
   plt.plot(Zdist,Zx[:len(Zdist)]/ft,'r')
   plt.axis('normal'); plt.axis('tight')
   plt.ylim(np.max(Zbed/ft)+5,0)
   plt.ylabel('Range (m)'); plt.xlabel('Distance along track (m)')
   custom_save(sonpath,'bed_2picks'+str(k))
   plt.close(); del fig

# =========================================================
def plot_bedpick(dat_port, dat_star, Zbed, Zdist, ft, shape_port, sonpath, k, chunkmode):

   if chunkmode != 4:
      extent = shape_port[1] #np.shape(merge)[0]
   else:
      extent = shape_port[0]

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
# Return a sliding window over a in any number of dimensions
def sliding_window(a,ws,ss = None,flatten = True):
   '''
   Return a sliding window over a in any number of dimensions
   '''
   if None is ss:
      # ss was not provided. the windows will not overlap in any direction.
      ss = ws
   ws = norm_shape(ws)
   ss = norm_shape(ss)
   # convert ws, ss, and a.shape to numpy arrays
   ws = np.array(ws)
   ss = np.array(ss)

   #import PyHum.io as io

   shape_tmp = io.set_mmap_data('', '', 'tmp.dat', 'float32', a)
   del a
   a = io.get_mmap_data('', '', 'tmp.dat', 'float32', shape_tmp)

   shap = np.array(a.shape)

   try:
      os.remove('tmp.dat')
   except:
      pass

   # ensure that ws, ss, and a.shape all have the same number of dimensions
   ls = [len(shap),len(ws),len(ss)]
   if 1 != len(set(ls)):
      raise ValueError(\
      'a.shape, ws and ss must all have the same length. They were %s' % str(ls))

   # ensure that ws is smaller than a in every dimension
   if np.any(ws > shap):
      raise ValueError(\
      'ws cannot be larger than a in any dimension.\
 a.shape was %s and ws was %s' % (str(a.shape),str(ws)))

   # how many slices will there be in each dimension?
   newshape = norm_shape(((shap - ws) // ss) + 1)
   # the shape of the strided array will be the number of slices in each dimension
   # plus the shape of the window (tuple addition)
   newshape += norm_shape(ws)
   # the strides tuple will be the array's strides multiplied by step size, plus

   try:
      # the array's strides (tuple addition)
      newstrides = norm_shape(np.array(a.strides) * ss) + a.strides
      a = ast(a,shape = newshape,strides = newstrides)
      if not flatten:
         return a
      # Collapse strided so that it has one more dimension than the window.  I.e.,
      # the new array is a flat list of slices.
      meat = len(ws) if ws.shape else 0
      firstdim = (int(np.product(newshape[:-meat])),) if ws.shape else ()
      dim = firstdim + (newshape[-meat:])
      # remove any dimensions with size 1
      dim = filter(lambda i : i != 1,dim)

      return a.reshape(dim), newshape

   except:

      from itertools import product
      print "memory error, windowing using slower method"
      # For each dimension, create a list of all valid slices
      slices = [[] for i in range(len(ws))]
      for i in xrange(len(ws)):
         nslices = ((shap[i] - ws[i]) // ss[i]) + 1
         for j in xrange(0,nslices):
            start = j * ss[i]
            stop = start + ws[i]
            slices[i].append(slice(start,stop))
      # Get an iterator over all valid n-dimensional slices of the input
      allslices = product(*slices)

      # Allocate memory to hold all valid n-dimensional slices
      nslices = np.product([len(s) for s in slices])
      #out = np.ndarray((nslices,) + tuple(ws),dtype = a.dtype)
      out=[]
      for i,s in enumerate(allslices):
         #out[i] = a[s]
         out.append(a[s])

      del a
      import dask.bag as db
      tmp = db.from_sequence(out, npartitions=1000)
      del out

      return tmp.compute(), newshape


# =========================================================
def norm_shape(shap):
   '''
   Normalize numpy array shapes so they're always expressed as a tuple,
   even for one-dimensional shapes.
   '''
   try:
      i = int(shap)
      return (i,)
   except TypeError:
      # shape was not a number
      pass

   try:
      t = tuple(shap)
      return t
   except TypeError:
      # shape was not iterable
      pass

   raise TypeError('shape must be an int, or a tuple of ints')

# =========================================================
# =========================================================
if __name__ == '__main__':

   read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, model, calc_bearing, filt_bearing, chunk) #cog

       #try:
       #except:
       #   print "memory error ... using dask array on high-freq. scan"
       #   dat = da.from_array(data_dwnhi, chunks=1000)
       #   del data_dwnhi
       #   Zt, ind_hi = makechunks_scan(chunkmode, chunkval, metadat, dat, 3)
       #   del dat


#    f : float, *optional* [Default=455]
#       frequency of sidescan transducer in kHz
#        f : frequency of sidescan sound (kHz)
#
#import matplotlib.colors as colors
#import simplekml

#import dask.array as da
   #import tkMessageBox

    #if f:
    #  f = np.asarray(f,int)
    #  print 'Frequency is %s kHz' % (str(f))
