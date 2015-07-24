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
#                        
#   ____ ___  ____ _____ 
#  / __ `__ \/ __ `/ __ \
# / / / / / / /_/ / /_/ /
#/_/ /_/ /_/\__,_/ .___/ 
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
import pyproj

import write

# numerical
import numpy as np
import PyHum.utils as humutils
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
    'map',
    'custom_save',
    'custom_save2',    
    'calc_beam_pos',
    ]

#################################################
def map(humfile, sonpath, cs2cs_args = "epsg:26949", dogrid = 1, res = 0.1, dowrite = 0, mode=3, nn = 64, influence = 1, numstdevs=4):
         
    '''
    Create plots of the spatially referenced sidescan echograms

    Syntax
    ----------
    [] = PyHum.map(humfile, sonpath, cs2cs_args, dogrid, res, dowrite, mode, nn, influence, numstdevs)

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
    res : float, *optional* [Default=0.1]
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
       Threshold number of standard deviations in sidescan intensity per grid cell up to which to accept            

    
    Returns
    -------
    sonpath+'x_y_ss_raw'+str(p)+'.asc'  : text file
        contains the point cloud of easting, northing, and sidescan intensity
        of the pth chunk

    sonpath+'GroundOverlay'+str(p)+'.kml': kml file
        contains gridded (or point cloud) sidescan intensity map for importing into google earth
        of the pth chunk

    sonpath+'map'+str(p)+'.png' : 
        image overlay associated with the kml file

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
       print 'Threshold number of standard deviations in sidescan intensity per grid cell up to which to accept: %s' % (str(numstdevs))             


    # start timer
    if os.name=='posix': # true if linux/mac or cygwin on windows
       start = time.time()
    else: # windows
       start = time.clock()

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
    
    theta = np.squeeze(meta['heading'])/(180/np.pi)

    # load memory mapped scans
    shape_port = np.squeeze(meta['shape_port'])
    if shape_port!='':
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat'))):
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_lar.dat')), 'r') as ff:
             port_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_port))
       else:
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')), 'r') as ff:
             port_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_port))

    shape_star = np.squeeze(meta['shape_star'])
    if shape_star!='':
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat'))):
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_lar.dat')), 'r') as ff:
             star_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_star))
       else:
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')), 'r') as ff:
             star_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_star))

    # time varying gain
    tvg = ((8.5*10**-5)+(3/76923)+((8.5*10**-5)/4))*c
        
    # depth correction
    dist_tvg = ((np.tan(np.radians(25)))*dep_m)-(tvg)

    # read in range data
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_range.dat')), 'r') as ff:
       R_fp = np.memmap(ff, dtype='float32', mode='r', shape=tuple(shape_star))

    for p in xrange(len(star_fp)):
       res = make_map(esi[shape_port[-1]*p:shape_port[-1]*(p+1)], nsi[shape_port[-1]*p:shape_port[-1]*(p+1)], theta[shape_port[-1]*p:shape_port[-1]*(p+1)], dist_tvg[shape_port[-1]*p:shape_port[-1]*(p+1)], port_fp[p], star_fp[p], R_fp[p], pix_m, res, cs2cs_args, sonpath, p, dogrid, dowrite, mode, nn, influence, numstdevs)
       print "grid resolution is %s" % (str(res))

    if os.name=='posix': # true if linux/mac
       elapsed = (time.time() - start)
    else: # windows
       elapsed = (time.clock() - start)
    print "Processing took ", elapsed , "seconds to analyse"

    print "Done!"


# =========================================================
def custom_save(figdirec,root):
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=600,transparent=True)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=1000, transparent=True)

# =========================================================
def calc_beam_pos(dist, bearing, x, y):

   dist_x, dist_y = (dist*np.sin(bearing), dist*np.cos(bearing))
   xfinal, yfinal = (x + dist_x, y + dist_y)
   return (xfinal, yfinal)

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
def make_map(e, n, t, d, dat_port, dat_star, data_R, pix_m, res, cs2cs_args, sonpath, p, dogrid, dowrite, mode, nn, influence, numstdevs):
   
   trans =  pyproj.Proj(init=cs2cs_args)   
   
   merge = np.vstack((dat_port,dat_star))
   #merge = np.vstack((np.flipud(port_fp[p]),star_fp[p]))
   
   merge[np.isnan(merge)] = 0
   merge = merge[:,:len(n)]

   R = np.vstack((np.flipud(data_R),data_R))
   R = R[:np.shape(merge)[0],:np.shape(merge)[1]]
  
   # get number pixels in scan line
   extent = int(np.shape(merge)[0]/2)

   yvec = np.linspace(pix_m,extent*pix_m,extent)

   X, Y, D, h, t  = getXY(e,n,yvec,d,t,extent)
   
   D[np.isnan(D)] = 0
   h[np.isnan(h)] = 0
   t[np.isnan(t)] = 0
       
   X = X[np.where(np.logical_not(np.isnan(Y)))]
   merge = merge.flatten()[np.where(np.logical_not(np.isnan(Y)))]
   Y = Y[np.where(np.logical_not(np.isnan(Y)))]
   D = D[np.where(np.logical_not(np.isnan(Y)))]
   R = R.flatten()[np.where(np.logical_not(np.isnan(Y)))]
   h = h[np.where(np.logical_not(np.isnan(Y)))]
   t = t[np.where(np.logical_not(np.isnan(Y)))]   
         
   Y = Y[np.where(np.logical_not(np.isnan(X)))]
   merge = merge.flatten()[np.where(np.logical_not(np.isnan(X)))]
   X = X[np.where(np.logical_not(np.isnan(X)))]
   D = D[np.where(np.logical_not(np.isnan(X)))]
   R = R.flatten()[np.where(np.logical_not(np.isnan(X)))]
   h = h[np.where(np.logical_not(np.isnan(X)))]
   t = t[np.where(np.logical_not(np.isnan(X)))]   
         
   X = X[np.where(np.logical_not(np.isnan(merge)))]
   Y = Y[np.where(np.logical_not(np.isnan(merge)))]
   merge = merge[np.where(np.logical_not(np.isnan(merge)))]
   D = D[np.where(np.logical_not(np.isnan(merge)))]
   R = R[np.where(np.logical_not(np.isnan(merge)))]
   h = h[np.where(np.logical_not(np.isnan(merge)))]
   t = t[np.where(np.logical_not(np.isnan(merge)))]   
         
   if dowrite==1:
      ## write raw bs to file
      outfile = os.path.normpath(os.path.join(sonpath,'x_y_ss_raw'+str(p)+'.asc'))
      write.txtwrite( outfile, np.hstack((humutils.ascol(X.flatten()),humutils.ascol(Y.flatten()), humutils.ascol(merge.flatten()), humutils.ascol(D.flatten()), humutils.ascol(np.cos(R.flatten())), humutils.ascol(h.flatten()), humutils.ascol(t.flatten())  )) )
      
   del D, R, h, t
  
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


      #influence = 1 #m

      if mode==1:

         complete=0
         while complete==0:
            try:
               dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count()) 
               stdev = None
               counts = None
               if 'dat' in locals(): 
                  complete=1 
            except:
               del grid_x, grid_y, targ_def, orig_def
               dat, stdev, counts, res, complete = getgrid_lm(humlon, humlat, merge, influence, minX, maxX, minY, maxY, res*2, mode)         

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

      del X, Y

      dat = dat.reshape(shape)

      if mode>1:
         stdev = stdev.reshape(shape)
         counts = counts.reshape(shape)

      mask = dat.mask.copy()

      dat[mask==1] = 0

      if mode>1:
         dat[(stdev>numstdevs) & (mask!=0)] = np.nan
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
      ### mask
      #if np.floor(np.sqrt(1/res))-1 > 0.0:
      #   dat[dist> np.floor(np.sqrt(1/res))-1 ] = np.nan #np.floor(np.sqrt(1/res))-1 ] = np.nan
      #else:
      #   dat[dist> np.sqrt(1/res) ] = np.nan #np.floor(np.sqrt(1/res))-1 ] = np.nan

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
       llcrnrlon=np.min(humlon)-0.00001, llcrnrlat=np.min(humlat)-0.00001,
       urcrnrlon=np.max(humlon)+0.00001, urcrnrlat=np.max(humlat)+0.00001)

      if dogrid==1:
         gx,gy = map.projtran(glon, glat)

      ax = plt.Axes(fig, [0., 0., 1., 1.], )
      ax.set_axis_off()
      fig.add_axes(ax)

      if dogrid==1:
         if datm.size > 25000000:
            print "matrix size > 25,000,000 - decimating by factor of 5 for display"
            map.pcolormesh(gx[::5,::5], gy[::5,::5], datm[::5,::5], cmap='gray', vmin=np.nanmin(datm), vmax=np.nanmax(datm))
         else:
            map.pcolormesh(gx, gy, datm, cmap='gray', vmin=np.nanmin(datm), vmax=np.nanmax(datm))
         del datm, dat
      else: 
         ## draw point cloud
         x,y = map.projtran(humlon, humlat)
         map.scatter(x.flatten(), y.flatten(), 0.5, merge.flatten(), cmap='gray', linewidth = '0')

      custom_save(sonpath,'map'+str(p))
      del fig 

   except:
      print "error: map could not be created..."

   kml = simplekml.Kml()
   ground = kml.newgroundoverlay(name='GroundOverlay')
   ground.icon.href = 'map'+str(p)+'.png'
   ground.latlonbox.north = np.min(humlat)-0.00001
   ground.latlonbox.south = np.max(humlat)+0.00001
   ground.latlonbox.east =  np.max(humlon)+0.00001
   ground.latlonbox.west =  np.min(humlon)-0.00001
   ground.latlonbox.rotation = 0

   #kml.save(sonpath+'GroundOverlay'+str(p)+'.kml')
   kml.save(os.path.normpath(os.path.join(sonpath,'GroundOverlay'+str(p)+'.kml')))

   del humlat, humlon
   return res #return the new resolution

# =========================================================
def getxy(e, n, yvec, d, t,extent):
   x = np.concatenate((np.tile(e,extent) , np.tile(e,extent)))
   #y = np.concatenate((n[k]+yvec, n[k]-yvec))
   rangedist = np.sqrt(np.power(yvec, 2.0) - np.power(d, 2.0))
   y = np.concatenate((n+rangedist, n-rangedist))
   # Rotate line around center point
   xx = e - ((x - e) * np.cos(t)) - ((y - n) * np.sin(t))
   yy = n - ((x - e) * np.sin(t)) + ((y - n) * np.cos(t))
   xx, yy = calc_beam_pos(d, t, xx, yy)
   #x, y, eucl. dist, depth, theta 
   return xx, yy, np.sqrt((xx-e)**2 + (yy-n)**2), np.ones(len(xx))*d, np.ones(len(xx))*t

# =========================================================
def getXY(e,n,yvec,d,t,extent):
   print "getting point cloud ..." 

   o = Parallel(n_jobs = -1, verbose=0)(delayed(getxy)(e[k], n[k], yvec, d[k], t[k], extent) for k in xrange(len(n)))

   #eating, northing, distance to sonar, depth, heading
   X, Y, D, h, t = zip(*o)

   # merge flatten and stack
   X = np.asarray(X,'float').T
   X = X.flatten()

   # merge flatten and stack
   Y = np.asarray(Y,'float').T
   Y = Y.flatten()

   # merge flatten and stack
   D = np.asarray(D,'float').T
   D = D.flatten()

   # merge flatten and stack
   h = np.asarray(h,'float').T
   h = h.flatten()
   
   # merge flatten and stack
   t = np.asarray(t,'float').T
   t = t.flatten()
         
   return X, Y, D, h, t

# =========================================================
# =========================================================
if __name__ == '__main__':

   map(humfile, sonpath, cs2cs_args, dogrid, res, dowrite, mode, nn, influence, numstdevs)





#               grid_x, grid_y, res = getmesh(np.min(X), np.max(X), np.min(Y), np.max(Y), res)
#               longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#               shape = np.shape(grid_x)

#               targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#               del longrid, latgrid

#               orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())

#               dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count()) 

#               if 'dat' in locals(): 
#                  complete=1 

         #grid_x, grid_y = np.meshgrid( np.arange(minX, maxX, res), np.arange(minY, maxY, res*2) )
         #if 'grid_x' in locals(): 
         #   complete=1 

#         try:
#            dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)
#         except:

#            print "Memory error: trying a grid resolution twice as big"
#            res = res*2

#            grid_x, grid_y = np.meshgrid( np.arange(np.min(X), np.max(X), res), np.arange(np.min(Y), np.max(Y), res) )  

#            #del X, Y
#            longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#            shape = np.shape(grid_x)
#            #del grid_y, grid_x

#            targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#            del longrid, latgrid

#            orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten()) 

#            dat, stdev, counts = pyresample.kd_tree.resample_gauss(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, neighbours=nn, sigmas=sigmas, fill_value=None, with_uncert = np.nan, nprocs = cpu_count(), epsilon = eps)

#         # custom inverse distance 
#         wf = lambda r: 1/r**2

#         try:
#            dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count())

#         except:

#            print "Memory error: trying a grid resolution twice as big"
#            res = res*2

#            grid_x, grid_y = np.meshgrid( np.arange(np.min(X), np.max(X), res), np.arange(np.min(Y), np.max(Y), res) )  

#            #del X, Y
#            longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#            shape = np.shape(grid_x)
#            #del grid_y, grid_x

#            targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#            del longrid, latgrid

#            orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten()) 

#            dat, stdev, counts = pyresample.kd_tree.resample_custom(orig_def, merge.flatten(),targ_def, radius_of_influence=influence, neighbours=nn, weight_funcs=wf, fill_value=None, with_uncert = True, nprocs = cpu_count()) 


#         try:
#            # nearest neighbour
#            dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count())

#         except:

#            print "Memory error: trying a grid resolution twice as big"
#            res = res*2

#            grid_x, grid_y = np.meshgrid( np.arange(np.min(X), np.max(X), res), np.arange(np.min(Y), np.max(Y), res) )  

#            #del X, Y
#            longrid, latgrid = trans(grid_x, grid_y, inverse=True)
#            shape = np.shape(grid_x)
#            #del grid_y, grid_x

#            targ_def = pyresample.geometry.SwathDefinition(lons=longrid.flatten(), lats=latgrid.flatten())
#            del longrid, latgrid

#            orig_def = pyresample.geometry.SwathDefinition(lons=humlon.flatten(), lats=humlat.flatten())

#            dat = pyresample.kd_tree.resample_nearest(orig_def, merge.flatten(), targ_def, radius_of_influence=influence, fill_value=None, nprocs = cpu_count())


       #port_fp = np.memmap(sonpath+base+'_data_port_l.dat', dtype='float32', mode='r', shape=tuple(shape_port))
       #star_fp = np.memmap(sonpath+base+'_data_star_l.dat', dtype='float32', mode='r', shape=tuple(shape_star))
       
#    # over-ride measured bearing and calc from positions
#    if calc_bearing==1:
#       lat = np.squeeze(meta['lat'])
#       lon = np.squeeze(meta['lon']) 

#       #point-to-point bearing
#       bearing = np.zeros(len(lat))
#       for k in xrange(len(lat)-1):
#          bearing[k] = bearingBetweenPoints(lat[k], lat[k+1], lon[k], lon[k+1])
#       del lat, lon

#    else:
#       # reported bearing by instrument (Kalman filtered?)
#       bearing = np.squeeze(meta['heading'])

#    # if stdev in heading is large, there's probably noise that needs to be filtered out
#    if np.std(bearing)>180:
#       print "WARNING: large heading stdev - attempting filtering"
#       from sklearn.cluster import MiniBatchKMeans
#       # can have two modes
#       data = np.column_stack([bearing, bearing])
#       k_means = MiniBatchKMeans(2)
#       # fit the model
#       k_means.fit(data) 
#       values = k_means.cluster_centers_.squeeze()
#       labels = k_means.labels_

#       if np.sum(labels==0) > np.sum(labels==1):
#          bearing[labels==1] = np.nan
#       else:
#          bearing[labels==0] = np.nan

#       nans, y= humutils.nan_helper(bearing)
#       bearing[nans]= np.interp(y(nans), y(~nans), bearing[~nans])

#       # save this filtered version to file
#       #meta = loadmat(sonpath+base+'meta.mat')
#       meta['heading_filt'] = bearing
#       #savemat(sonpath+base+'meta.mat', meta ,oned_as='row')
#       savemat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')), meta ,oned_as='row')
#       #del meta   

#    if filt_bearing ==1:
#       bearing = humutils.runningMeanFast(bearing, len(bearing)/100)

#    theta = np.asarray(bearing, 'float')/(180/np.pi)

#    if cog==1:
#       #course over ground is given as a compass heading (ENU) from True north, or Magnetic north.
#       #To get this into NED (North-East-Down) coordinates, you need to rotate the ENU 
#       # (East-North-Up) coordinate frame. 
#       #Subtract pi/2 from your heading
#       theta = theta - np.pi/2
#       # (re-wrap to Pi to -Pi)
#       theta = np.unwrap(-theta)

