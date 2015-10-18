# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: nonecheck=False
"""
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
"""

from __future__ import division
import numpy as np
cimport numpy as np
cimport cython

#from pylab import fft2, ifft2, fftshift, real
from scipy.interpolate import RectBivariateSpline

# =========================================================
cdef class Noise:

   cdef object res

   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] im, float factor=1.25):
      cdef int cols, rows
      cols, rows = np.shape(im)

      cdef np.ndarray[np.float64_t, ndim=2] res = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] noise = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] filt = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] radius = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] mag = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] rr = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.complex128_t, ndim=2] ffrr = np.empty( (rows, cols), dtype=np.complex128)
      cdef np.ndarray[np.complex128_t, ndim=2] imfft = np.empty( (rows, cols), dtype=np.complex128)
      cdef np.ndarray[np.complex128_t, ndim=2] phase = np.empty( (rows, cols), dtype=np.complex128)
      cdef np.ndarray[np.int64_t, ndim=2] xi = np.empty((rows,cols),dtype=np.int64)
      cdef np.ndarray[np.int64_t, ndim=2] yi = np.empty((rows,cols),dtype=np.int64)

      rr = np.random.randn(cols,rows)
      ffrr = np.fft.fft2(rr)
      imfft = np.fft.fftshift(ffrr)
      mag = np.abs(imfft)  
      phase = imfft/mag  
      xi, yi = np.meshgrid(np.r_[:rows].astype('int64'),np.r_[:cols].astype('int64'))  
      radius = np.sqrt(xi**2 + yi**2)
      radius[int(cols/2 + 1), int(rows/2 + 1)] = 1
      radius[radius==0] = 1
      filt = np.divide(1,(radius**factor))
      noise = np.real(np.fft.ifft2(np.fft.fftshift(np.multiply(filt,phase)))) 
      noise = noise/noise.sum() 
      res = self._rescale(self._im_resize(noise[::2,::2],cols,rows),np.nanmin(im),np.nanmax(im))
      self.res = np.asarray(res,'float16').T
      return 

   # =========================================================
   cpdef np.ndarray _rescale(self, np.ndarray dat, float mn, float mx):
      """
      rescales an input dat between mn and mx
      """
      cdef float m = np.min(dat.flatten())
      cdef float M = np.max(dat.flatten())
      return (mx-mn)*(dat-m)/(M-m)+mn

   # =========================================================
   cpdef np.ndarray _im_resize(self, np.ndarray im, int Nx, int Ny):
      '''
      resize array by bivariate spline interpolation
      '''
      cdef int nx
      cdef int ny
      ny, nx = np.shape(im)
      cdef np.ndarray xx = np.linspace(0,nx,Nx)
      cdef np.ndarray yy = np.linspace(0,ny,Ny)
      cdef object newKernel = RectBivariateSpline(np.r_[:ny],np.r_[:nx],im) 
      return newKernel(yy,xx)

   # =========================================================
   cpdef np.ndarray getres(self):
      """
      return result
      """
      return self.res

