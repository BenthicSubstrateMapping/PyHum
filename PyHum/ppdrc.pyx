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

# =========================================================
cdef class ppdrc:

   cdef object res

   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] im, int wavelength=768, int n=2):


      # Reference:
      # Peter Kovesi, "Phase Preserving Tone Mapping of Non-Photographic High Dynamic 
      # Range Images".  Proceedings: Digital Image Computing: Techniques and
      # Applications 2012 (DICTA 2012). Available via IEEE Xplore

      # translated from matlab code posted on:
      # http://www.csse.uwa.edu.au/~pk/research/matlabfns/PhaseCongruency/
      cdef int cols, rows
      cdef float eps = 2.2204e-16

      rows,cols = np.shape(im)    

      cdef np.ndarray[np.float64_t, ndim=2] E = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] H = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] res = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] radius = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] u1 = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] u2 = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] ph = np.empty( [rows, cols], dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] h1f = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] h2f = np.empty( (rows, cols), dtype=np.float64)
      cdef np.ndarray[np.float64_t, ndim=2] f = np.empty( (rows, cols), dtype=np.float64)

      cdef np.ndarray[np.complex128_t, ndim=2] IM = np.empty( (rows, cols), dtype=np.complex128)
      cdef np.ndarray[np.complex128_t, ndim=2] H1 = np.empty( (rows, cols), dtype=np.complex128)
      cdef np.ndarray[np.complex128_t, ndim=2] H2 = np.empty( (rows, cols), dtype=np.complex128)

      IM = np.fft.fft2(im)

      # Generate horizontal and vertical frequency grids that vary from
      # -0.5 to 0.5 
      u1, u2 = np.meshgrid((np.r_[0:cols]-(np.fix(cols/2)+1))/(cols-np.mod(cols,2)),(np.r_[0:rows]-(np.fix(rows/2)+1))/(rows-np.mod(rows,2)))

      u1 = np.fft.ifftshift(u1)   # Quadrant shift to put 0 frequency at the corners
      u2 = np.fft.ifftshift(u2)
    
      radius = np.sqrt(u1*u1 + u2*u2)
      # Matrix values contain frequency values as a radius from centre (but quadrant shifted)
    
      # Get rid of the 0 radius value in the middle (at top left corner after
      # fftshifting) so that dividing by the radius, will not cause trouble.
      radius[1,1] = 1
    
      H1 = 1j*u1/radius   # The two monogenic filters in the frequency domain
      H2 = 1j*u2/radius
      H1[1,1] = 0
      H2[1,1] = 0
      radius[1,1] = 0  # undo fudge
 
      # High pass Butterworth filter
      H =  1.0 - 1.0 / (1.0 + (radius * wavelength)**(2*n))       
         
      f = np.real(np.fft.ifft2(H*IM))
      h1f = np.real(np.fft.ifft2(H*H1*IM))
      h2f = np.real(np.fft.ifft2(H*H2*IM))
    
      ph = np.arctan(f/np.sqrt(h1f*h1f + h2f*h2f + eps))
      E = np.sqrt(f*f + h1f*h1f + h2f*h2f)
      res = np.sin(ph)*np.log1p(E)
      self.res = res


   # =========================================================    
   cpdef np.ndarray getdata(self):
      return self.res
       
