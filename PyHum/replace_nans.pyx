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

#from cython.parallel cimport prange

# =========================================================
cdef class RN:
   """
   Returns an instance.
   """
   cdef object data

   # =========================================================
   def __init__(self, np.ndarray[np.float64_t, ndim=2] array, int max_iter, float tol, int kernel_size=1, str method='localmean'):
    """Replace NaN elements in an array using an iterative image inpainting algorithm.
    The algorithm is the following:
    1) For each element in the input array, replace it by a weighted average
    of the neighbouring elements which are not NaN themselves. The weights depends
    of the method type. If ``method=localmean`` weight are equal to 1/( (2*kernel_size+1)**2 -1 )
    2) Several iterations are needed if there are adjacent NaN elements.
    If this is the case, information is "spread" from the edges of the missing
    regions iteratively, until the variation is below a certain threshold.
    Parameters
    ----------
    array : 2d np.ndarray
    an array containing NaN elements that have to be replaced
    max_iter : int
    the number of iterations
    kernel_size : int
    the size of the kernel, default is 1
    method : str
    the method used to replace invalid values. Valid options are
    `localmean`, 'idw'.
    Returns
    -------
    filled : 2d np.ndarray
    a copy of the input array, where NaN elements have been replaced.
    """
     
    cdef int i, j, I, J, it, k, l
    cdef int n_invalids 
    cdef double n
    
    cdef np.ndarray[np.float64_t, ndim=2] filled = np.empty( [array.shape[0], array.shape[1]], dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=2] kernel = np.empty( (2*kernel_size+1, 2*kernel_size+1), dtype=np.float64)

    cdef np.ndarray[np.int64_t, ndim=1] inans = np.empty( [array.shape[0]], dtype=np.int64)
    cdef np.ndarray[np.int64_t, ndim=1] jnans = np.empty( [array.shape[1]], dtype=np.int64)

    # indices where array is NaN
    inans, jnans = np.nonzero( np.isnan(array) )
    
    # number of NaN elements
    cdef int n_nans 
    n_nans = len(inans)
    
    # arrays which contain replaced values to check for convergence
    cdef np.ndarray[np.float64_t, ndim=1] replaced_new = np.zeros( n_nans, dtype=np.float64)
    cdef np.ndarray[np.float64_t, ndim=1] replaced_old = np.zeros( n_nans, dtype=np.float64)
    
    # depending on kernel type, fill kernel array
    if method == 'localmean':
      
        #print 'kernel_size', kernel_size       
        for i in range(2*kernel_size+1):
            for j in range(2*kernel_size+1):
                kernel[i,j] = 1
        #print kernel, 'kernel'

    elif method == 'idw':
        kernel = np.array([[0, 0.5, 0.5, 0.5,0],
                  [0.5,0.75,0.75,0.75,0.5], 
                  [0.5,0.75,1,0.75,0.5],
                  [0.5,0.75,0.75,0.5,1],
                  [0, 0.5, 0.5 ,0.5 ,0]])
        #print kernel, 'kernel'      
    else:
        raise ValueError( 'method not valid. Should be one of `localmean`.')
    
    filled = array.copy()

    # make several passes
    # until we reach convergence
    #for it in xrange(max_iter):
    for it from 0 <= it < max_iter: 
        #print 'iteration', it
        # for each NaN element
        #for k in xrange(n_nans):
        #for k from 0 <= k < n_nans: 
        for k from 0 <= k < n_nans:                      
            i = inans[k]
            j = jnans[k]
            
            # initialize to zero
            filled[i,j] = 0.0
            n = 0.0
            
            # loop over the kernel
            #for I in xrange(2*kernel_size+1):
            for I from 0 <= I < 2*kernel_size+1:
                #for J in xrange(2*kernel_size+1):
                for J from 0 <= J < 2*kernel_size+1:                     
                    # if we are not out of the boundaries
                    if i+I-kernel_size < array.shape[0] and i+I-kernel_size >= 0:
                        if j+J-kernel_size < array.shape[1] and j+J-kernel_size >= 0:
                                                
                            # if the neighbour element is not NaN itself.
                            if filled[i+I-kernel_size, j+J-kernel_size] == filled[i+I-kernel_size, j+J-kernel_size] :
                                
                                # do not sum itself
                                if I-kernel_size != 0 and J-kernel_size != 0:
                                    
                                    # convolve kernel with original array
                                    filled[i,j] = filled[i,j] + filled[i+I-kernel_size, j+J-kernel_size]*kernel[I, J]
                                    n = n + 1*kernel[I,J]
                                    #print n

            # divide value by effective number of added elements
            if n != 0:
                filled[i,j] = filled[i,j] / n
                replaced_new[k] = filled[i,j]
            else:
                filled[i,j] = np.nan
                
        # check if mean square difference between values of replaced
        #elements is below a certain tolerance
        #print 'tolerance', np.mean( (replaced_new-replaced_old)**2 )
        if np.mean( (replaced_new-replaced_old)**2 ) < tol:
            break
        else:
            for l in xrange(n_nans):
                replaced_old[l] = replaced_new[l]
    
    self.data = filled
    return #filled

   # =========================================================    
   cpdef np.ndarray getdata(self):
      return self.data
       
       
    
