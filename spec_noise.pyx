'''
spec_noise.pyx
Part of PyHum software 

INFO:
Cython script to pad data edges with spectral colored noise

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.0      Revision: July, 2014

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government.
'''

from __future__ import division
import numpy as np
cimport numpy as np

from pylab import fft2, ifft2, fftshift, real
from scipy.interpolate import RectBivariateSpline

# =========================================================
cdef class Noise:

   cdef object res

   # =========================================================
   def __init__(self, np.ndarray im, float factor=1.25):
      cdef int cols
      cdef int rows
      cdef np.ndarray yi
      cdef np.ndarray xi

      cols, rows = np.shape(im)
      cdef np.ndarray rr = np.random.randn(cols,rows)
      cdef np.ndarray ffrr = fft2(rr)
      cdef np.ndarray imfft = fftshift(ffrr)
      cdef np.ndarray mag = abs(imfft)  
      cdef np.ndarray phase = imfft/mag  
      xi, yi = np.meshgrid(np.r_[:rows],np.r_[:cols])  
      cdef np.ndarray radius = np.sqrt(xi**2 + yi**2)
      radius[int(cols/2 + 1), int(rows/2 + 1)] = 1
      radius[radius==0] = 1
      cdef np.ndarray filter = np.divide(1,(radius**factor))
      cdef np.ndarray noise = real(ifft2(fftshift(np.multiply(filter,phase)))) 
      noise = noise/noise.sum() 
      cdef np.ndarray res = self._rescale(self._im_resize(noise[::2,::2],cols,rows),np.nanmin(im),np.nanmax(im))
      self.res = np.asarray(res,'float16').T
      return 

   # =========================================================
   def _rescale(self, np.ndarray dat, float mn, float mx):
      """
      rescales an input dat between mn and mx
      """
      cdef float m = np.min(dat.flatten())
      cdef float M = np.max(dat.flatten())
      return (mx-mn)*(dat-m)/(M-m)+mn

   # =========================================================
   def _im_resize(self, np.ndarray im, int Nx, int Ny):
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
   def getres(self):
      """
      return result
      """
      return self.res

