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
#                             __            __                
#   ____ ___  ____ _____     / /____  _  __/ /___  __________ 
#  / __ `__ \/ __ `/ __ \   / __/ _ \| |/_/ __/ / / / ___/ _ \
# / / / / / / /_/ / /_/ /  / /_/  __/>  </ /_/ /_/ / /  /  __/
#/_/ /_/ /_/\__,_/ .___/   \__/\___/_/|_|\__/\__,_/_/   \___/ 
#               /_/                                           
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
from __future__ import division
from scipy.io import loadmat
import os, time, sys, getopt
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory
except:
   pass
from joblib import Parallel, delayed, cpu_count

# numerical
import numpy as np
import pyproj
import PyHum.utils as humutils
#from scipy.interpolate import griddata
#from scipy.spatial import cKDTree as KDTree
#from scipy.ndimage.filters import median_filter

import pyresample
import replace_nans
from scipy.ndimage import binary_dilation, binary_erosion, binary_fill_holes

# plotting
import matplotlib.pyplot as plt
try:
   from mpl_toolkits.basemap import Basemap
except:
   print "Error: Basemap could not be imported"
   pass
import simplekml

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

import warnings
warnings.filterwarnings("ignore")

__all__ = [
    'map_texture',
    'custom_save',
    'custom_save2',    
    'calc_beam_pos',
    ]

#################################################
def map_texture(humfile, sonpath, cs2cs_args = "epsg:26949", dogrid = 1, res = 0.5, dowrite = 0, mode=3, nn = 64, influence = 1, numstdevs=5):
         
    '''
    Create plots of the texture lengthscale maps made in PyHum.texture module 
    using the algorithm detailed by Buscombe et al. (forthcoming)
    This textural lengthscale is not a direct measure of grain size. Rather, it is a statistical 
    representation that integrates over many attributes of bed texture, of which grain size is the most important. 
    The technique is a physically based means to identify regions of texture within a sidescan echogram, 
    and could provide a basis for objective, automated riverbed sediment classification.

    Syntax
    ----------
    [] = PyHum.map_texture(humfile, sonpath, cs2cs_args, dogrid, res, dowrite, mode, nn, influence, numstdevs)

    Parameters
    ----------
    humfile : str
       path to the .DAT file
    sonpath : str
       path where the *.SON files are
    cs2cs_args : int, *optional* [Default="epsg:26949"]
       arguments to create coordinates in a projected coordinate system
       this argument gets given to pyproj to turn wgs84 (lat/lon) coordinates
       into any projection supported by the proj.4 libraries
    dogrid : float, *optional* [Default=1]
       if 1, textures will be gridded with resolution 'res'. 
       Otherwise, point cloud will be plotted
    res : float, *optional* [Default=0.5]
       grid resolution of output gridded texture map
    dowrite: int, *optional* [Default=0]
       if 1, point cloud data from each chunk is written to ascii file
       if 0, processing times are speeded up considerably but point clouds are not available for further analysis
    mode: int, *optional* [Default=3]
       gridding mode. 1 = nearest neighbour
                      2 = inverse weighted nearest neighbour
                      3 = Gaussian weighted nearest neighbour
    nn: int, *optional* [Default=64]
       number of nearest neighbours for gridding (used if mode > 1)
    influence: float, *optional* [Default=1]
       Radius of influence used in gridding. Cut off distance in meters   
    numstdevs: int, *optional* [Default = 4]
       Threshold number of standard deviations in texture lengthscale per grid cell up to which to accept 
           
    Returns
    -------
    sonpath+'x_y_class'+str(p)+'.asc' : text file
        contains the point cloud of easting, northing, and texture lengthscales
        of the pth chunk

    sonpath+'class_GroundOverlay'+str(p)+'.kml': kml file
        contains gridded (or point cloud) texture lengthscale map for importing into google earth
        of the pth chunk

    sonpath+'class_map'+str(p)+'.png' : 
        image overlay associated with the kml file

    sonpath+'class_map_imagery'+str(p)+'.png' : png image file
        gridded (or point cloud) texture lengthscale map
        overlain onto an image pulled from esri image server

    References
    ----------
      .. [1] Buscombe, D., Grams, P.E., and Smith, S.M.C., 2015, Automated riverbed sediment
       classification using low-cost sidescan sonar. Journal of Hydraulic Engineering, accepted
    '''

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

    if cs2cs_args:
       print 'cs2cs arguments are %s' % (cs2cs_args)

    if dogrid:
       dogrid = int(dogrid)
       if dogrid==1:
          print "Data will be gridded"         

    if res:
       res = np.asarray(res,float)
       print 'Gridding resolution: %s' % (str(res))      

    if dowrite:
       dowrite = int(dowrite)
       if dowrite==0:
          print "Point cloud data will be written to ascii file" 

    if mode:
       mode = int(mode)
       print 'Mode for gridding: %s' % (str(mode))      

    if nn:
       nn = int(nn)
       print 'Number of nearest neighbours for gridding: %s' % (str(nn))             

    if influence:
       influence = int(influence)
       print 'Radius of influence for gridding: %s (m)' % (str(influence))             

    if numstdevs:
       numstdevs = int(numstdevs)
       print 'Threshold number of standard deviations in texture lengthscale per grid cell up to which to accept: %s' % (str(numstdevs))             


    trans =  pyproj.Proj(init=cs2cs_args)

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

    meta = loadmat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')))

    esi = np.squeeze(meta['e'])
    nsi = np.squeeze(meta['n']) 

    pix_m = np.squeeze(meta['pix_m'])
    dep_m = np.squeeze(meta['dep_m'])
    c = np.squeeze(meta['c'])
    dist_m = np.squeeze(meta['dist_m'])

    theta = np.squeeze(meta['heading'])/(180/np.pi)

    # load memory mapped scans
    shape_port = np.squeeze(meta['shape_port'])
    if shape_port!='':
       #port_fp = np.memmap(sonpath+base+'_data_port_l.dat', dtype='float32', mode='r', shape=tuple(shape_port))
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat'))):
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat')), 'r') as ff:
             port_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_port))
       else:
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')), 'r') as ff:
             port_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_port))

    shape_star = np.squeeze(meta['shape_star'])
    if shape_star!='':
       #star_fp = np.memmap(sonpath+base+'_data_star_l.dat', dtype='float32', mode='r', shape=tuple(shape_star))
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat'))):
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat')), 'r') as ff:
             star_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_star))
       else:
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')), 'r') as ff:
             star_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_star))

    shape = shape_port.copy()
    shape[1] = shape_port[1] + shape_star[1]
    #class_fp = np.memmap(sonpath+base+'_data_class.dat', dtype='float32', mode='r', shape=tuple(shape))
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_class.dat')), 'r') as ff:
       class_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape))


    tvg = ((8.5*10**-5)+(3/76923)+((8.5*10**-5)/4))*c
    dist_tvg = ((np.tan(np.radians(25)))*dep_m)-(tvg)

    for p in xrange(len(class_fp)):

       e = esi[shape_port[-1]*p:shape_port[-1]*(p+1)]
       n = nsi[shape_port[-1]*p:shape_port[-1]*(p+1)]
       t = theta[shape_port[-1]*p:shape_port[-1]*(p+1)]
       d = dist_tvg[shape_port[-1]*p:shape_port[-1]*(p+1)]

       len_n = len(n)
   
       merge = class_fp[p].copy()

       merge[np.isnan(merge)] = 0
       merge[np.isnan(np.vstack((np.flipud(port_fp[p]),star_fp[p])))] = 0

       extent = shape_port[1]
       R1 = merge[extent:,:]
       R2 = np.flipud(merge[:extent,:])

       merge = np.vstack((R2,R1))
       del R1, R2

       # get number pixels in scan line
       extent = int(np.shape(merge)[0]/2)

       yvec = np.linspace(pix_m,extent*pix_m,extent)

       X, Y  = getXY(e,n,yvec,d,t,extent)

       merge[merge==0] = np.nan

       if len(merge.flatten()) != len(X):
          merge = merge[:,:len_n]

       merge = merge.T.flatten()

       index = np.where(np.logical_not(np.isnan(merge)))[0]

       X = X.flatten()[index]
       Y = Y.flatten()[index]
       merge = merge.flatten()[index]

       X = X[np.where(np.logical_not(np.isnan(Y)))]
       merge = merge.flatten()[np.where(np.logical_not(np.isnan(Y)))]
       Y = Y[np.where(np.logical_not(np.isnan(Y)))]

       Y = Y[np.where(np.logical_not(np.isnan(X)))]
       merge = merge.flatten()[np.where(np.logical_not(np.isnan(X)))]
       X = X[np.where(np.logical_not(np.isnan(X)))]

       X = X[np.where(np.logical_not(np.isnan(merge)))]
       Y = Y[np.where(np.logical_not(np.isnan(merge)))]
       merge = merge[np.where(np.logical_not(np.isnan(merge)))]

       if dowrite==1:
          # write raw bs to file
          #outfile = sonpath+'x_y_class'+str(p)+'.asc' 
          outfile = os.path.normpath(os.path.join(sonpath,'x_y_class'+str(p)+'.asc'))
          with open(outfile, 'w') as f:
             np.savetxt(f, np.hstack((humutils.ascol(X),humutils.ascol(Y), humutils.ascol(merge))), delimiter=' ', fmt="%8.6f %8.6f %8.6f")

       humlon, humlat = trans(X, Y, inverse=True)

       if dogrid==1:

          complete=0
          while complete==0:
             try:
                grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
                #del X, Y
                longrid, latgrid = trans(grid_x, grid_y, inverse=True)
                shape = np.shape(grid_x)
                #del grid_y, grid_x

                targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
                del longrid, latgrid

                orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())
                #del humlat, humlon
                if 'orig_def' in locals(): 
                   complete=1 
             except:
                print "memory error: trying grid resolution of %s" % (str(res*2))
                res = res*2

          if mode==1:

             complete=0
             while complete==0:
                try:
                   dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count()) 
                   if 'dat' in locals(): 
                      complete=1 
                except:
                   del grid_x, grid_y, targ_def, orig_def
                   dat, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)

          elif mode==2:
             # custom inverse distance 
             wf = lambda r: 1/r**2

             complete=0
             while complete==0:
                try:
                   dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())
                   if 'dat' in locals(): 
                      complete=1 
                except:
                   del grid_x, grid_y, targ_def, orig_def
                   dat, stdev, counts, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)


          elif mode==3:
             sigmas = 1 #m
             eps = 2

             complete=0
             while complete==0:
                try:
                   dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)
                   if 'dat' in locals(): 
                      complete=1 
                except:
                   del grid_x, grid_y, targ_def, orig_def
                   dat, stdev, counts, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)

          dat = dat.reshape(shape)

          if mode>1:
             stdev = stdev.reshape(shape)
             counts = counts.reshape(shape)

          mask = dat.mask.copy()

          dat[mask==1] = 0

          if mode>1:
             dat[(stdev>3) & (mask!=0)] = np.nan
             dat[(counts<nn) & (counts>0)] = np.nan

          dat2 = replace_nans.RN(dat.astype('float64'),1000,0.01,2,'localmean').getdata()
          dat2[dat==0] = np.nan

          # get a new mask
          mask = np.isnan(dat2)

          mask = ~binary_dilation(binary_erosion(~mask,structure=np.ones((15,15))), structure=np.ones((15,15)))
          #mask = binary_fill_holes(mask, structure=np.ones((15,15)))
          #mask = ~binary_fill_holes(~mask, structure=np.ones((15,15)))

          dat2[mask==1] = np.nan
          dat2[dat2<1] = np.nan

          del dat
          dat = dat2
          del dat2

       if dogrid==1:
          ## mask
          #dat[dist> 1 ] = np.nan 

          #del dist, tree

          dat[dat==0] = np.nan
          dat[np.isinf(dat)] = np.nan

          datm = np.ma.masked_invalid(dat)

          glon, glat = trans(grid_x, grid_y, inverse=True)
          del grid_x, grid_y

       try:
          print "drawing and printing map ..."
          fig = plt.figure(frameon=False)
          map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #26949,
           resolution = 'i', #h #f
           llcrnrlon=np.min(humlon)-0.0001, llcrnrlat=np.min(humlat)-0.0001,
           urcrnrlon=np.max(humlon)+0.0001, urcrnrlat=np.max(humlat)+0.0001)

          if dogrid==1:
             gx,gy = map.projtran(glon, glat)

          ax = plt.Axes(fig, [0., 0., 1., 1.], )
          ax.set_axis_off()
          fig.add_axes(ax)

          if dogrid==1:
             if datm.size > 25000000:
                print "matrix size > 25,000,000 - decimating by factor of 5 for display"
                map.pcolormesh(gx[::5,::5], gy[::5,::5], datm[::5,::5], cmap='YlOrRd', vmin=0.25, vmax=2)
             else:
                map.pcolormesh(gx, gy, datm, cmap='YlOrRd', vmin=0.25, vmax=2)

             #map.pcolormesh(gx, gy, datm, cmap='YlOrRd', vmin=0.25, vmax=2)
             del dat
          else: 
             ## draw point cloud
             x,y = map.projtran(humlon, humlat)
             map.scatter(x.flatten(), y.flatten(), 0.5, merge.flatten(), cmap='YlOrRd', linewidth = '0')

          custom_save(sonpath,'class_map'+str(p))
          del fig 

       except:
          print "error: map could not be created..."

       kml = simplekml.Kml()
       ground = kml.newgroundoverlay(name='GroundOverlay')
       ground.icon.href = 'class_map'+str(p)+'.png'
       ground.latlonbox.north = np.min(humlat)-0.00001
       ground.latlonbox.south = np.max(humlat)+0.00001
       ground.latlonbox.east =  np.max(humlon)+0.00001
       ground.latlonbox.west =  np.min(humlon)-0.00001
       ground.latlonbox.rotation = 0

       #kml.save(sonpath+'class_GroundOverlay'+str(p)+'.kml')
       kml.save(os.path.normpath(os.path.join(sonpath,'class_GroundOverlay'+str(p)+'.kml')))

    if dowrite==1:

       X = []; Y = []; S = [];
       for p in xrange(len(class_fp)):
          #dat = np.genfromtxt(sonpath+'x_y_class'+str(p)+'.asc', delimiter=' ')
          dat = np.genfromtxt(os.path.normpath(os.path.join(sonpath,'x_y_class'+str(p)+'.asc')), delimiter=' ')
          X.append(dat[:,0])
          Y.append(dat[:,1])
          S.append(dat[:,2])
          del dat

       # merge flatten and stack
       X = np.asarray(np.hstack(X),'float')
       X = X.flatten()

       # merge flatten and stack
       Y = np.asarray(np.hstack(Y),'float')
       Y = Y.flatten()

       # merge flatten and stack
       S = np.asarray(np.hstack(S),'float')
       S = S.flatten()

       humlon, humlat = trans(X, Y, inverse=True)

       if dogrid==1:

          complete=0
          while complete==0:
             try:
                grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
                #del X, Y
                longrid, latgrid = trans(grid_x, grid_y, inverse=True)
                shape = np.shape(grid_x)
                #del grid_y, grid_x

                targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
                del longrid, latgrid

                orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())
                #del humlat, humlon
                if 'orig_def' in locals(): 
                   complete=1 
             except:
                print "memory error: trying grid resolution of %s" % (str(res*2))
                res = res*2

          if mode==1:

             complete=0
             while complete==0:
                try:
                   dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count()) 
                   if 'dat' in locals(): 
                      complete=1 
                except:
                   del grid_x, grid_y, targ_def, orig_def
                   dat, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)

          elif mode==2:
             # custom inverse distance 
             wf = lambda r: 1/r**2

             complete=0
             while complete==0:
                try:
                   dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())
                   if 'dat' in locals(): 
                      complete=1 
                except:
                   del grid_x, grid_y, targ_def, orig_def
                   dat, stdev, counts, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)


          elif mode==3:
             sigmas = 1 #m
             eps = 2

             complete=0
             while complete==0:
                try:
                   dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)
                   if 'dat' in locals(): 
                      complete=1 
                except:
                   del grid_x, grid_y, targ_def, orig_def
                   dat, stdev, counts, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)

          dat = dat.reshape(shape)

          if mode>1:
             stdev = stdev.reshape(shape)
             counts = counts.reshape(shape)

          mask = dat.mask.copy()

          dat[mask==1] = 0

          if mode>1:
             dat[(stdev>5) & (mask!=0)] = np.nan
             dat[(counts<nn) & (counts>0)] = np.nan

          dat2 = replace_nans.RN(dat.astype('float64'),1000,0.01,2,'localmean').getdata()
          dat2[dat==0] = np.nan

          # get a new mask
          mask = np.isnan(dat2)

          mask = ~binary_dilation(binary_erosion(~mask,structure=np.ones((15,15))), structure=np.ones((15,15)))
          #mask = binary_fill_holes(mask, structure=np.ones((15,15)))
          #mask = ~binary_fill_holes(~mask, structure=np.ones((15,15)))

          dat2[mask==1] = np.nan
          dat2[dat2<1] = np.nan

          del dat
          dat = dat2
          del dat2


       if dogrid==1:
          ## mask
          #dat[dist> 1 ] = np.nan

          #el dist, tree

          dat[dat==0] = np.nan
          dat[np.isinf(dat)] = np.nan

          datm = np.ma.masked_invalid(dat)

          glon, glat = trans(grid_x, grid_y, inverse=True)
          del grid_x, grid_y

       levels = [0.5,0.75,1.25,1.5,1.75,2,3]

       try:
          print "drawing and printing map ..."
          fig = plt.figure()
          map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1],
           resolution = 'i',
           llcrnrlon=np.min(humlon)-0.00001, llcrnrlat=np.min(humlat)-0.00001,
           urcrnrlon=np.max(humlon)+0.00001, urcrnrlat=np.max(humlat)+0.00001)

          map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)
          if dogrid==1:
             gx,gy = map.projtran(glon, glat)

          if dogrid==1:

             if datm.size > 25000000:
                print "matrix size > 25,000,000 - decimating by factor of 5 for display"
                map.contourf(gx[::5,::5], gy[::5,::5], datm[::5,::5], levels, cmap='YlOrRd')
             else:
                map.contourf(gx, gy, datm, levels, cmap='YlOrRd')

             #map.contourf(gx, gy, datm, levels, cmap='YlOrRd')
          else: 
             ## draw point cloud
             x,y = map.projtran(humlon, humlat)
             map.scatter(x.flatten(), y.flatten(), 0.5, S.flatten(), cmap='YlOrRd', linewidth = '0')

          custom_save2(sonpath,'class_map_imagery'+str(p))
          del fig 
       except:
          print "error: map could not be created..."

# =========================================================
def getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res, mode):

   complete=0
   while complete==0:
      try:
         grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
         longrid, latgrid = trans(grid_x, grid_y, inverse=True)
         shape = np.shape(grid_x)
         targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
         del longrid, latgrid

         orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())

         if mode==1:
            dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count())
            stdev = None
            counts = None
         elif mode==2:
            dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())
         else:
            dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)
 
         if 'dat' in locals(): 
            complete=1 
      except:
         print "memory error: trying grid resolution of %s" % (str(res*2))
         res = res*2

   return dat, stdev, counts, res, complete

# =========================================================
def getmesh(minX, maxX, minY, maxY, res):

   complete=0
   while complete==0:
      try:
         grid_x, grid_y = np.meshgrid( np.arange(minX, maxX, res), np.arange(minY, maxY, res) )
         if 'grid_x' in locals(): 
            complete=1 
      except:
         print "memory error: trying grid resolution of %s" % (str(res*2))
         res = res*2
         
   return grid_x, grid_y, res

# =========================================================
def getxy(e, n, yvec, d, t,extent):
   x = np.concatenate((np.tile(e,extent) , np.tile(e,extent)))
   rangedist = np.sqrt(np.power(yvec, 2.0) - np.power(d, 2.0))
   y = np.concatenate((n+rangedist, n-rangedist))
   # Rotate line around center point
   xx = e - ((x - e) * np.cos(t)) - ((y - n) * np.sin(t))
   yy = n - ((x - e) * np.sin(t)) + ((y - n) * np.cos(t))
   xx, yy = calc_beam_pos(d, t, xx, yy)
   return xx, yy 


# =========================================================
def getXY(e,n,yvec,d,t,extent):
   print "getting point cloud ..." 

   o = Parallel(n_jobs = -1, verbose=0)(delayed(getxy)(e[k], n[k], yvec, d[k], t[k], extent) for k in xrange(len(n)))

   X, Y = zip(*o)

   # merge flatten and stack
   X = np.asarray(X,'float')
   X = X.flatten()

   # merge flatten and stack
   Y = np.asarray(Y,'float')
   Y = Y.flatten()

   return X, Y

# =========================================================
def custom_save(figdirec,root):
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=600,transparent=True)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=1000, transparent=True)

# =========================================================
def custom_save2(figdirec,root):
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=600)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=600)

# =========================================================
def calc_beam_pos(dist, bearing, x, y):

   dist_x, dist_y = (dist*np.sin(bearing), dist*np.cos(bearing))
   xfinal, yfinal = (x + dist_x, y + dist_y)
   return (xfinal, yfinal)

# =========================================================
# =========================================================
if __name__ == '__main__':

   map_texture(humfile, sonpath, cs2cs_args, dogrid, res, dowrite, mode, nn, influence, numstdevs)




#          complete=0
#          while complete==0:
#             try:
#                grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
#                #del X, Y
#                longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#                shape = np.shape(grid_x)
#                #del grid_y, grid_x

#                targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#                del longrid, latgrid

#                orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())
#                #del humlat, humlon
#                if 'dat' in locals(): 
#                   complete=1 
#             except:
#                print "memory error: trying grid resolution of %s" % (str(res*2))
#                res = res*2

#          if mode==1:

#             complete=0
#             while complete==0:
#                try:
#                   dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count())
#                   if 'dat' in locals(): 
#                      complete=1 
#                except:
#                   del grid_x, grid_y, longrid, latgrid, targ_def, orig_def
#                   grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
#                   longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#                   shape = np.shape(grid_x)

#                   targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#                   del longrid, latgrid

#                   orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())

#                   dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count()) 

#                   if 'dat' in locals(): 
#                      complete=1     


#          elif mode==2:
#             # custom inverse distance 
#             wf = lambda r: 1/r**2

#             complete=0
#             while complete==0:
#                try:
#                   dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())
#                   if 'dat' in locals(): 
#                      complete=1 
#                except:
#                   del grid_x, grid_y, longrid, latgrid, targ_def, orig_def
#                   grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
#                   longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#                   shape = np.shape(grid_x)

#                   targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#                   del longrid, latgrid

#                   orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())

#                   dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())

#                   if 'dat' in locals(): 
#                      complete=1     

#          elif mode==3:
#             sigmas = 1 #m
#             eps = 2

#             complete=0
#             while complete==0:
#                try:
#                   dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)
#                   if 'dat' in locals(): 
#                      complete=1 
#                except:
#                   del grid_x, grid_y, longrid, latgrid, targ_def, orig_def
#                   grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
#                   longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#                   shape = np.shape(grid_x)

#                   targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#                   del longrid, latgrid

#                   orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())

#                   dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)

#                   if 'dat' in locals(): 
#                      complete=1 



#             try:
#                # nearest neighbour
#                dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count())
#             except:
#                # nearest neighbour
#                dat, stdev, counts = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, with_uncert = True, nprocs = 1)


#          grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)

#          #grid_x, grid_y = np.meshgrid( np.arange(np.min(X), np.max(X), res), np.arange(np.min(Y), np.max(Y), res) )  

#          #del X, Y
#          longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#          shape = np.shape(grid_x)
#          #del grid_y, grid_x

#          targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#          del longrid, latgrid

#          orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())
#          #del humlat, humlon

#          #influence = 1 #m
#          #numneighbours = 64

#             try:
#                dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)
#             except:
#                dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = 1, epsilon = eps)

#             try:
#                dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())
#             except:
#                dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = 1)   

#             try:
#                # nearest neighbour
#                dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count())
#             except:
#                # nearest neighbour
#                dat, stdev, counts = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, with_uncert = True, nprocs = 1)


#          grid_x, grid_y = np.meshgrid( np.arange(np.min(X), np.max(X), res), np.arange(np.min(Y), np.max(Y), res) )  

#          #dat = griddata(np.c_[X.flatten(),Y.flatten()], S.flatten(), (grid_x, grid_y), method='nearest')
#          ### create mask for where the data is not
#          #tree = KDTree(np.c_[X.flatten(),Y.flatten()])
#          #dist, _ = tree.query(np.c_[grid_x.ravel(), grid_y.ravel()], k=1)
#          #dist = dist.reshape(grid_x.shape)

#          tree = KDTree(zip(X.flatten(), Y.flatten()))

#          #nearest neighbour
#          dist, inds = tree.query(zip(grid_x.flatten(), grid_y.flatten()), k = 1)
#          dat = S.flatten()[inds].reshape(grid_x.shape)

#          ## inverse distance weighting, using 10 nearest neighbours
#          #d, inds = tree.query(zip(grid_x.flatten(), grid_y.flatten()), k = 10)
#          #w = 1.0 / d**2
#          #dat = np.sum(w * merge.flatten()[inds], axis=1) / np.sum(w, axis=1)
#          #dat.shape = grid_x.shape

#          ## create mask for where the data is not
#          #dist, _ = tree.query(np.c_[grid_x.ravel(), grid_y.ravel()], k=1)
#          dist = dist.reshape(grid_x.shape)




#          #dat = griddata(np.c_[X.flatten(),Y.flatten()], S.flatten(), (grid_x, grid_y), method='nearest')
#          ### create mask for where the data is not
#          #tree = KDTree(np.c_[X.flatten(),Y.flatten()])
#          #dist, _ = tree.query(np.c_[grid_x.ravel(), grid_y.ravel()], k=1)
#          #dist = dist.reshape(grid_x.shape)

#          tree = KDTree(zip(X.flatten(), Y.flatten()))

#          #nearest neighbour
#          dist, inds = tree.query(zip(grid_x.flatten(), grid_y.flatten()), k = 1)
#          dat = merge.flatten()[inds].reshape(grid_x.shape)

#          ## inverse distance weighting, using 10 nearest neighbours
#          #d, inds = tree.query(zip(grid_x.flatten(), grid_y.flatten()), k = 10)
#          #w = 1.0 / d**2
#          #dat = np.sum(w * merge.flatten()[inds], axis=1) / np.sum(w, axis=1)
#          #dat.shape = grid_x.shape

#          ## create mask for where the data is not
#          #dist, _ = tree.query(np.c_[grid_x.ravel(), grid_y.ravel()], k=1)
#          dist = dist.reshape(grid_x.shape)

#       del X, Y

#    if not cs2cs_args:
#       # arguments to pass to cs2cs for coordinate transforms
#       cs2cs_args = "epsg:26949"
#       print '[Default] cs2cs arguments are %s' % (cs2cs_args)

#    if not dogrid:
#       if dogrid != 0:
#          dogrid = 1
#          print "[Default] Data will be gridded"

#    if not calc_bearing:
#       if calc_bearing != 1:
#          calc_bearing = 0
#          print "[Default] Heading recorded by instrument will be used"

#    if not filt_bearing:
#       if filt_bearing != 1:
#          filt_bearing = 0
#          print "[Default] Heading will not be filtered"

#    if not res:
#       res = 0.5
#       print '[Default] Grid resolution is %s m' % (str(res))

#    if not cog:
#       if cog != 0:
#          cog = 1
#          print "[Default] Heading based on course-over-ground"

#    if not dowrite:
#       if dowrite != 0:
#          dowrite = 1
#          print "[Default] Point cloud data will be written to ascii file"

#       print "getting point cloud ..."
#       # get the points by rotating the [x,y] vector so it lines up with boat heading
#       X=[]; Y=[]; 
#       for k in range(len(n)): 
#          x = np.concatenate((np.tile(e[k],extent) , np.tile(e[k],extent)))
#          #y = np.concatenate((n[k]+yvec, n[k]-yvec))
#          rangedist = np.sqrt(np.power(yvec, 2.0) - np.power(d[k], 2.0))
#          y = np.concatenate((n[k]+rangedist, n[k]-rangedist))
#          # Rotate line around center point
#          xx = e[k] - ((x - e[k]) * np.cos(t[k])) - ((y - n[k]) * np.sin(t[k]))
#          yy = n[k] - ((x - e[k]) * np.sin(t[k])) + ((y - n[k]) * np.cos(t[k]))
#          xx, yy = calc_beam_pos(d[k], t[k], xx, yy)
#          X.append(xx)
#          Y.append(yy) 

#       del e, n, t, x, y

#       # merge flatten and stack
#       X = np.asarray(X,'float')
#       X = X.flatten()

#       # merge flatten and stack
#       Y = np.asarray(Y,'float')
#       Y = Y.flatten()
