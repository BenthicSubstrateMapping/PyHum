'''
Part of PyHum software 

INFO:


Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.2.3      Revision: Apr, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
'''

from numpy.lib.stride_tricks import as_strided as ast
#from numpy import array, product, isnan, min, max, convolve, isnan, ones, mean, std, argmax, where, interp, shape, zeros, hstack, vstack, argmin, squeeze, choose, linspace, r_, cumsum, histogram, any, seterr, cos, sin, rad2deg, arctan2, deg2rad, arcsin, power, sqrt

import numpy as np

from numpy import nan as npnan
from numpy.matlib import repmat

from sklearn.cluster import MiniBatchKMeans
from scipy.interpolate import RectBivariateSpline
import string, random

# suppress divide and invalid warnings
#seterr(divide='ignore')
#seterr(invalid='ignore')
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
def get_bearing(calc_bearing, cog, filt_bearing, lat, lon, heading):

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

       nans, y= nan_helper(bearing)
       bearing[nans]= np.interp(y(nans), y(~nans), bearing[~nans]) 

    if filt_bearing ==1:
       bearing = runningMeanFast(bearing, len(bearing)/100)
       
    if cog==1:
       theta = np.asarray(bearing, 'float')/(180/np.pi)
       #course over ground is given as a compass heading (ENU) from True north, or Magnetic north.
       #To get this into NED (North-East-Down) coordinates, you need to rotate the ENU 
       # (East-North-Up) coordinate frame. 
       #Subtract pi/2 from your heading
       theta = theta - np.pi/2
       # (re-wrap to Pi to -Pi)
       theta = np.unwrap(-theta)
       bearing = theta * (180/np.pi)     
       
    return bearing
    

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
   dat[nans]= interp(y(nans), y(~nans), dat[~nans])
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
   firstdim = (int(np.product(newshape[:-meat])),) if ws.shape else ()
   dim = firstdim + (newshape[-meat:])
   # remove any dimensions with size 1
   dim = filter(lambda i : i != 1,dim) 
    
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
      d = np.repmat( imu[i,:], 3, 1 ) + np.vstack( (c0,tmp1,tmp2) )
      del tmp1, tmp2
      p[i,:] =  np.argmin(d,axis=0)
      c[i,:] =  np.min(d,axis=0)

   p[p==0] = -1
   p = p+1

   x = np.zeros((m,1))
   cost = np.min(c[-1,:])
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
   k_means = MiniBatchKMeans(numclusters)
   # fit the model
   k_means.fit(wc) 
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
    
