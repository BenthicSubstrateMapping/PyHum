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

import RunningStats

# =========================================================
cdef class proc:
   """
   Returns an instance.
   """
   cdef object data

   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=1] points): 

      # pre-allocate arrays
      cdef list filled

      rs1 = RunningStats.RunningStats()
      # global stats, not detrended
      for k in points:
         rs1.Push(k)
 
      # compile all parameters into a list
      filled = [rs1.Mean(), rs1.StandardDeviation()]

      self.data = filled

      rs1.Clear()

      return

   # =========================================================
   cpdef list getdata(self):
      return self.data


