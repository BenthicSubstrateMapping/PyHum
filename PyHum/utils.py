'''
Part of PyHum software

INFO:


Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior.
For more information, see the official USGS copyright policy at
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
'''

from numpy.lib.stride_tricks import as_strided as ast
import os
import numpy as np
from numpy.matlib import repmat

import matplotlib.pyplot as plt

from simplekml import (Kml, OverlayXY, ScreenXY, Units, RotationXY,
                       AltitudeMode, Camera)

from sklearn.cluster import MiniBatchKMeans
from scipy.interpolate import RectBivariateSpline
import string, random
from scipy.ndimage.filters import median_filter

import dask.array as da

# suppress divide and invalid warnings
np.seterr(all='ignore')

__all__ = [
    'distBetweenPoints',
    'bearingBetweenPoints',
    'id_generator',
    'rm_spikes',
    'ascol',
    'rescale',
    'runningMeanFast',
    'nan_helper',
    'norm_shape',
    'sliding_window',
    'dpboundary',
    'cut_kmeans',
    'im_resize',
    'histeq',
    ]

#################################################

# =========================================================
def auto_bedpick(ft, dep_m, chunkmode, port_fp, c):
    #buff = 50#10

    # get bed from depth trace
    #bed = ft*dep_m
    tvg = ((8.5*10**-5)+(3/76923)+((8.5*10**-5)/4))*c
    bed = ft*(dep_m - tvg)

    bed = rm_spikes(bed,3)
    bed = runningMeanFast(bed, 100)

    buff = min(bed)
    imu = []

    if chunkmode!=4:
      for k in xrange(len(port_fp)):
         #imu.append(port_fp[k][int(np.min(bed)):int(np.max(bed)),:])
         imu.append(port_fp[k][np.max([0,int(np.min(bed))-buff]):int(np.max(bed))+buff,:])
      imu = np.hstack(imu)
    else:
      imu.append(port_fp[np.max([0,int(np.min(bed))-buff]):int(np.max(bed))+buff,:])

    imu = np.squeeze(np.asarray(imu, 'float64'))-buff

    try:
       imu = da.from_array(imu, chunks=1000)   #dask implementation
    except:
       pass

    #imu = median_filter(imu,(20,20))
    imu = median_filter(imu,(np.shape(imu)[0]/100,np.shape(imu)[1]/100))

    #autobed = dpboundary(-imu[buff:,:].T)+buff

    dx,dy = np.gradient(imu)
    lap = np.sqrt(dx**2 + dy**2)
    del dx, dy
    autobed = dpboundary(-lap[buff:,:].T)+buff
    del lap

    autobed = np.squeeze(autobed)

    if len(autobed) < len(bed):
       autobed2 = bed.copy()
       autobed2[:len(autobed)] = autobed
       del autobed
       autobed = autobed2

    ## narrow image to within range of estimated bed
    # use dynamic boundary tracing to get 2nd estimate of bed
    #x = np.squeeze(int(np.min(bed))+dpboundary(-imu.T))
    #x = np.squeeze(int(np.min(bed))+autobed)
    x = np.min(np.vstack((np.squeeze(bed),np.squeeze(autobed))), axis=0)
    del imu

    if len(x)<len(bed):
       x = np.append(x,x[-1]*np.ones(len(bed)-len(x)))
    elif len(x)>len(bed):
       bed = np.append(bed,bed[-1]*np.ones(len(x)-len(bed)))

    # if standard deviation of auto bed pick is too small, then use acoustic bed pick
    if np.std(x)<5:
       print "stdev of auto bed pick is low, using acoustic pick"
       x = bed.copy()

    return x, bed

# =========================================================
def make_trackline(lon,lat, sonpath, base):

    import simplekml
    # create kml for loading path into google earth
    kml = simplekml.Kml()
    ls = kml.newlinestring(name='trackline')
    ls.coords = zip(lon,lat)
    ls.extrude = 1
    ls.altitudemode = simplekml.AltitudeMode.relativetoground
    ls.style.linestyle.width = 5
    ls.style.linestyle.color = simplekml.Color.red
    kml.save(os.path.normpath(os.path.join(sonpath,base+'trackline.kml')))

# =========================================================
def get_depth(dep_m):

    dep_m = np.squeeze(dep_m)
    dep_m = rm_spikes(dep_m,2)
    return runningMeanFast(dep_m, 3)

# =========================================================
def get_dist(lat, lon):

    dist = np.zeros(len(lat))
    for k in xrange(len(lat)-1):
       dist[k] = distBetweenPoints(lat[k], lat[k+1], lon[k], lon[k+1])

    return np.cumsum(dist)

# =========================================================
def get_bearing(calc_bearing, filt_bearing, lat, lon, heading): #cog

    # over-ride measured bearing and calc from positions
    if calc_bearing==1:
       lat = np.squeeze(lat)
       lon = np.squeeze(lon)

       #point-to-point bearing
       bearing = np.zeros(len(lat))
       for k in xrange(len(lat)-1):
          bearing[k] = bearingBetweenPoints(lat[k], lat[k+1], lon[k], lon[k+1])
       #del lat, lon

    else:
       # reported bearing by instrument (Kalman filtered?)
       bearing = np.squeeze(heading)

    bearing = rm_spikes(bearing,3)

    # if stdev in heading is large, there's probably noise that needs to be filtered out
    if np.std(bearing)>180:
       print "WARNING: large heading stdev - attempting filtering"
       from sklearn.cluster import MiniBatchKMeans
       # can have two modes
       data = np.column_stack([bearing, bearing])

       try:
          data = da.from_array(data, chunks=1000)   #dask implementation
       except:
          pass

       k_means = MiniBatchKMeans(2)
       # fit the model
       k_means.fit(data)
       del data

       #values = k_means.cluster_centers_.squeeze()
       labels = k_means.labels_

       if np.sum(labels==0) > np.sum(labels==1):
          bearing[labels==1] = np.nan
       else:
          bearing[labels==0] = np.nan

       nans, y= nan_helper(bearing)
       bearing[nans]= np.interp(y(nans), y(~nans), bearing[~nans])

    if filt_bearing ==1:
       bearing = runningMeanFast(bearing, np.max((len(bearing)/10, 3)))

    ##if cog==1:
    theta = np.asarray(bearing, 'float')/(180/np.pi)
    ##course over ground is given as a compass heading (ENU) from True north, or Magnetic north.
    ##To get this into NED (North-East-Down) coordinates, you need to rotate the ENU
    ## (East-North-Up) coordinate frame.
    ##Subtract pi/2 from your heading
    theta = theta - np.pi/2
    # (re-wrap to Pi to -Pi)
    theta = np.unwrap(-theta)
    bearing = theta * (180/np.pi)

    #return (bearing + 360) % 360
    return bearing % 360

# =========================================================
def strip_base(base):
    # remove underscores, negatives and spaces from basename
    if base.find('_')>-1:
       base = base[:base.find('_')]

    if base.find('-')>-1:
       base = base[:base.find('-')]

    if base.find(' ')>-1:
       base = base[:base.find(' ')]

    if base.find('.')>-1:
       base = base[:base.find('.')]
    return base

# =========================================================
def distBetweenPoints(pos1_lat, pos2_lat, pos1_lon, pos2_lon):
   return 6378137.0 * 2.0 * np.arcsin(np.sqrt(np.power(np.sin((np.deg2rad(pos1_lat) - np.deg2rad(pos2_lat)) / 2.0), 2.0) + np.cos(np.deg2rad(pos1_lat)) * np.cos(np.deg2rad(pos2_lat)) * np.power(np.sin((np.deg2rad(pos1_lon) - np.deg2rad(pos2_lon)) / 2.0), 2.0)))


# # =========================================================
# def bearingBetweenPoints(pos1_lat, pos2_lat, pos1_lon, pos2_lon):
#    lat1 = np.deg2rad(pos1_lat)
#    lon1 = np.deg2rad(pos1_lon)
#    lat2 = np.deg2rad(pos2_lat)
#    lon2 = np.deg2rad(pos2_lon)
#
#    bearing = np.arctan2(np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(lon2 - lon1), np.sin(lon2 - lon1) * np.cos(lat2))
#
#    db = np.rad2deg(bearing)
#    return (90.0 - db + 360.0) % 360.0

# # =========================================================
def bearingBetweenPoints(pos1_lat, pos2_lat, pos1_lon, pos2_lon):
   lat1 = np.deg2rad(pos1_lat)
   lat2 = np.deg2rad(pos2_lat)

   diffLong = np.deg2rad(pos2_lon) - np.deg2rad(pos1_lon)
   #orig #bearing = np.arctan2(np.cos(lat1) * np.sin(lat2) - (np.sin(lat1) * np.cos(lat2) * np.cos(diffLong)), np.sin(diffLong) * np.cos(lat2))

   bearing = np.arctan2(np.sin(diffLong) * np.cos(lat2), np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(diffLong))

   db = np.rad2deg(bearing)
   #return (90.0 - db + 360.0) % 360.0
   return db % 360.0

# =========================================================
def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))

# =========================================================
def ascol( arr ):
   '''
   reshapes row matrix to be a column matrix (N,1).
   '''
   if len( arr.shape ) == 1: arr = arr.reshape( ( arr.shape[0], 1 ) )
   return arr

# =========================================================
def rm_spikes(dat,numstds):
   """
   remove spikes in dat
   """
   ht = np.mean(dat) + numstds*np.std(dat)
   lt = np.argmax(np.mean(dat) - numstds*np.std(dat),0)

   index = np.where(dat>ht);
   if index:
      dat[index] = np.nan

   index = np.where(dat<lt);
   if index:
      dat[index] = np.nan

   # fill nans using linear interpolation
   nans, y= nan_helper(dat)
   dat[nans]= np.interp(y(nans), y(~nans), dat[~nans])
   return dat

# =========================================================
def rescale(dat,mn,mx):
   """
   rescales an input dat between mn and mx
   """
   m = np.min(dat.flatten())
   M = np.max(dat.flatten())
   return (mx-mn)*(dat-m)/(M-m)+mn

# =========================================================
def runningMeanFast(x, N):
   '''
   flawed but fast running mean
   '''
   x = np.convolve(x, np.ones((N,))/N)[(N-1):]
   # the last N values will be crap, so they're set to the global mean
   x[-N:] = x[-N]
   return x

# =========================================================
def nan_helper(y):
   '''
   function to help manage indices of nans
   '''
   return np.isnan(y), lambda z: z.nonzero()[0]

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
# Return a sliding window over a in any number of dimensions
# version with no memory mapping
def sliding_window_nomm(a,ws,ss = None,flatten = True):
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
    shap = np.array(a.shape)
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
    # the array's strides (tuple addition)
    newstrides = norm_shape(np.array(a.strides) * ss) + a.strides
    a = ast(a,shape = newshape,strides = newstrides)
    if not flatten:
        return a
    # Collapse strided so that it has one more dimension than the window.  I.e.,
    # the new array is a flat list of slices.
    meat = len(ws) if ws.shape else 0
    firstdim = (np.product(newshape[:-meat]),) if ws.shape else ()
    dim = firstdim + (newshape[-meat:])
    # remove any dimensions with size 1
    dim = filter(lambda i : i != 1,dim)

    return a.reshape(dim), newshape


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

   import PyHum.io as io

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
# Return a sliding window over a in any number of dimensions
def sliding_window_sliced(a,density, ws,ss = None,flatten = True):
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

   r = np.arange(1,ws[0]-1,density, dtype=np.int)

   shap = np.array(a.shape)

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

   # the array's strides (tuple addition)
   newstrides = norm_shape(np.array(a.strides) * ss) + a.strides
   a = ast(a,shape = newshape,strides = newstrides)[:,:,:,r]
   if not flatten:
      return a
   # Collapse strided so that it has one more dimension than the window.  I.e.,
   # the new array is a flat list of slices.
   meat = len(ws) if ws.shape else 0
   firstdim = (int(np.product(newshape[:-meat])),) if ws.shape else ()
   dim = firstdim + (newshape[-meat:])

   dim = list(dim)
   dim[-1] = len(r)
   ## remove any dimensions with size 1
   dim = filter(lambda i : i != 1,dim)
   dim = tuple(dim)

   newshape = np.shape(a)

   return a.reshape(dim), newshape


# =========================================================
def dpboundary(imu):
   '''
   dynamic boundary tracing in an image
   (translated from matlab: CMP Vision Algorithms http://visionbook.felk.cvut.cz)
   '''
   m,n = np.shape(imu)
   c = np.zeros((m,n))
   p = np.zeros((m,n))
   c[0,:] = imu[0,:]

   for i in xrange(1,m):
      c0 = c[i-1,:]
      tmp1 = np.squeeze(ascol(np.hstack((c0[1:],c0[-1]))))
      tmp2 = np.squeeze(ascol(np.hstack((c0[0], c0[0:len(c0)-1]))))
      d = repmat( imu[i,:], 3, 1 ) + np.vstack( (c0,tmp1,tmp2) )
      del tmp1, tmp2
      p[i,:] =  np.argmin(d,axis=0)
      c[i,:] =  np.min(d,axis=0)

   p[p==0] = -1
   p = p+1

   x = np.zeros((m,1))
   #cost = np.min(c[-1,:])
   xpos = np.argmin( c[-1,:] )
   for i in reversed(range(1,m)):
      x[i] = xpos
      if p[i,xpos]==2 and xpos<n:
         xpos = xpos+1
      elif p[i,xpos]==3 and xpos>1:
         xpos = xpos-1
   x[0] = xpos
   return x

## =========================================================
def cut_kmeans(w,numclusters):
   '''
   perform a k-means segmentation of image
   '''
   wc = w.reshape((-1, 1)) # We need an (n_sample, n_feature) array

   try:
      wc = da.from_array(wc, chunks=1000)   #dask implementation
   except:
      pass

   k_means = MiniBatchKMeans(numclusters)
   # fit the model
   k_means.fit(wc)
   del wc
   values = k_means.cluster_centers_.squeeze()
   labels = k_means.labels_
   # make the cut and reshape
   wc = np.choose(labels, values)
   wc.shape = w.shape
   return wc, values

# =========================================================
def im_resize(im,Nx,Ny):
   '''
   resize array by bivariate spline interpolation
   '''
   ny, nx = np.shape(im)
   xx = np.linspace(0,nx,Nx)
   yy = np.linspace(0,ny,Ny)

   try:
      im = da.from_array(im, chunks=1000)   #dask implementation
   except:
      pass

   newKernel = RectBivariateSpline(np.r_[:ny],np.r_[:nx],im)
   return newKernel(yy,xx)

# =========================================================
def histeq(im,nbr_bins=256):

   im[np.isnan(im)] = 0
   #get image histogram
   imhist,bins = np.histogram(im.flatten(),nbr_bins,normed=True)
   cdf = imhist.cumsum() #cumulative distribution function
   cdf = 255 * cdf / cdf[-1] #normalize

   #use linear interpolation of cdf to find new pixel values
   im2 = np.interp(im.flatten(),bins[:-1],cdf)

   return im2.reshape(im.shape), cdf

## http://nbviewer.ipython.org/url/ocefpaf.github.com/python4oceanographers/downloads/notebooks/2014-03-10-google-earth.ipynb
# =========================================================
def gearth_fig(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat, pixels=1024):
    """Return a Matplotlib `fig` and `ax` handles for a Google-Earth Image."""
    aspect = np.cos(np.mean([llcrnrlat, urcrnrlat]) * np.pi/180.0)
    xsize = np.ptp([urcrnrlon, llcrnrlon]) * aspect
    ysize = np.ptp([urcrnrlat, llcrnrlat])
    aspect = ysize / xsize

    if aspect > 1.0:
        figsize = (10.0 / aspect, 10.0)
    else:
        figsize = (10.0, 10.0 * aspect)

    if False:
        plt.ioff()  # Make `True` to prevent the KML components from poping-up.
    fig = plt.figure(figsize=figsize,
                     frameon=False,
                     dpi=pixels//10)
    # KML friendly image.  If using basemap try: `fix_aspect=False`.
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(llcrnrlon, urcrnrlon)
    ax.set_ylim(llcrnrlat, urcrnrlat)
    return fig, ax


# =========================================================
def make_kml(llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat,
             figs, colorbar=None, **kw):
    """TODO: LatLon bbox, list of figs, optional colorbar figure,
    and several simplekml kw..."""

    kml = Kml()
    altitude = kw.pop('altitude', 2e7)
    roll = kw.pop('roll', 0)
    tilt = kw.pop('tilt', 0)
    altitudemode = kw.pop('altitudemode', AltitudeMode.relativetoground)
    camera = Camera(latitude=np.mean([urcrnrlat, llcrnrlat]),
                    longitude=np.mean([urcrnrlon, llcrnrlon]),
                    altitude=altitude, roll=roll, tilt=tilt,
                    altitudemode=altitudemode)

    kml.document.camera = camera
    draworder = 0
    for fig in figs:  # NOTE: Overlays are limited to the same bbox.
        draworder += 1
        ground = kml.newgroundoverlay(name='GroundOverlay')
        ground.draworder = draworder
        ground.visibility = kw.pop('visibility', 1)
        ground.name = kw.pop('name', 'overlay')
        ground.color = kw.pop('color', '9effffff')
        ground.atomauthor = kw.pop('author', 'ocefpaf')
        ground.latlonbox.rotation = kw.pop('rotation', 0)
        ground.description = kw.pop('description', 'Matplotlib figure')
        ground.gxaltitudemode = kw.pop('gxaltitudemode',
                                       'clampToSeaFloor')
        ground.icon.href = fig
        ground.latlonbox.east = llcrnrlon
        ground.latlonbox.south = llcrnrlat
        ground.latlonbox.north = urcrnrlat
        ground.latlonbox.west = urcrnrlon

    if colorbar:  # Options for colorbar are hard-coded (to avoid a big mess).
        screen = kml.newscreenoverlay(name='ScreenOverlay')
        screen.icon.href = colorbar
        screen.overlayxy = OverlayXY(x=0, y=0,
                                     xunits=Units.fraction,
                                     yunits=Units.fraction)
        screen.screenxy = ScreenXY(x=0.015, y=0.075,
                                   xunits=Units.fraction,
                                   yunits=Units.fraction)
        screen.rotationXY = RotationXY(x=0.5, y=0.5,
                                       xunits=Units.fraction,
                                       yunits=Units.fraction)
        screen.size.x = 0
        screen.size.y = 0
        screen.size.xunits = Units.fraction
        screen.size.yunits = Units.fraction
        screen.visibility = 1

    kmzfile = kw.pop('kmzfile', 'overlay.kmz')
    kml.savekmz(kmzfile)
