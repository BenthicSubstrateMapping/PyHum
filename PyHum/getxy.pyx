from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sin,cos

# =========================================================
cdef class GetXY:

    cdef object e
    cdef object n
    cdef object d
    cdef object t
    cdef object xx
    cdef object yy
                        
    # =========================================================
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    def __init__(self, float e, float n, np.ndarray yvec, float d, float t, int extent):

       cdef float sint
       cdef float cost
       cdef float dsq
       cdef float dist_y
       cdef float dist_x
                            
       self.e = e
       self.n = n
       self.d = d
       self.t = t                 

       with nogil:
          sint = sin(t)
          cost = cos(t)
          dsq = d*d
          dist_y = d*cost
          dist_x = d*sint 
                                
       rangedist = np.sqrt(yvec*yvec - dsq)               
       x = np.concatenate((np.tile(e,extent) , np.tile(e,extent)))
       y = np.concatenate((n+rangedist, n-rangedist))
       
       # Rotate line around center point
       xx = e - ((x - e) * cost) - ((y - n) * sint)
       yy = n - ((x - e) * sint) + ((y - n) * cost)
       
       self.xx = xx + dist_x
       self.yy = yy + dist_y
       return
  
       
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    def getdat(self):
        """
        returns data
        """  
        xxe = (self.xx-self.e)
        yyn = (self.yy-self.n)
        return self.xx, self.yy, np.sqrt((xxe*xxe) + (yyn*yyn)), np.ones(len(self.xx))*self.d, np.ones(len(self.xx))*self.t
        
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    def getdat2(self):
        """
        returns data
        """  
        return self.xx, self.yy
                
