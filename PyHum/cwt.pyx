# cython module for wavelet computations
# Daniel Buscombe, June-Aug 2014
from __future__ import division
import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt,log,abs,floor

# =========================================================
cdef class Cwt:
    """
    continuous Morlet wavelet transform
    Implements via the Fourier transform
    Returns an instance.
    """
    cdef object data
    cdef object scale
    cdef object notes
    cdef object cwt
    cdef object fftdata
    cdef object nscale
    cdef object scales
    cdef object currentscale
    cdef object win
    #cdef object density
    cdef object r
    cdef object docalc
            
    # =========================================================
    @cython.boundscheck(False)
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    def __init__(self, np.ndarray[np.float32_t, ndim=2] matrix, int largestscale, int notes, int win):
    
        """
        Continuous Morlet wavelet transform of data

        data:    data in array to transform
        notes:   number of scale intervals per octave
        largestscale: largest scale as inverse fraction of length
                 of data array
                 scale = len(data)/largestscale
                 smallest scale should be >= 2 for meaningful data
        """
        self.win = win
        #self.density = density
        
        cdef float pi = 3.14159265
        
        #cdef np.ndarray r
        #r = np.arange(1,self.win-1,self.density, dtype=np.int)
        #self.r = r
        #cdef int lr = len(self.r)
        
        cdef int lr = np.shape(matrix)[1]
        self.r = lr
         
        cdef int i, scaleindex
        
        cdef np.ndarray[np.float64_t, ndim=0] currentscale 
        
        cdef double base2 
        #base2 = np.floor(log(self.win)/log(2) + 0.4999)
        with nogil:
           base2 = floor(log(win)/log(2) + 0.4999)
                           
        cdef int ndata = int(2**(base2+1)) #len(data)
        cdef int tmp = 0
        self.nscale = tmp
        self.scale = largestscale
        self._setscales(ndata,largestscale,notes)
        #cdef np.ndarray[np.complex64_t, ndim=3] cwt = np.zeros((self.nscale,ndata,len(self.r)),dtype=np.complex64)
        cdef np.ndarray[np.complex64_t, ndim=3] cwt = np.zeros((self.nscale,ndata,lr),dtype=np.complex64)        
        self.cwt = cwt
        cdef np.ndarray[np.float64_t, ndim=1] omega = np.empty(ndata, dtype=np.float64)
        omega = np.array(range(0,np.int(ndata/2))+range(-np.int(ndata/2),0))*(2.0*pi/ndata)
        
        cdef np.ndarray[np.float32_t,ndim=1] data = np.empty(self.win, dtype=np.float32)
        cdef np.ndarray[np.float64_t,ndim=1] data2 = np.empty(ndata, dtype=np.float64)
        cdef np.ndarray[np.complex128_t,ndim=1] datahat = np.empty(ndata, dtype=np.complex128)
        cdef np.ndarray[np.float64_t,ndim=1] s_omega = np.empty(ndata, dtype=np.float64)
        cdef np.ndarray[np.float64_t,ndim=1] psihat = np.empty(ndata, dtype=np.float64) 
        
        self.docalc = 0
        
        if np.sum(matrix)>0:
           self.docalc = 1

           for i from 0 <= i < lr:  
              #data = np.asarray( self._column(matrix, np.int(self.r[i]) ) )
              data = np.asarray( self._column(matrix, i ) )           
              data2 = self._pad2nxtpow2(data - np.mean(data), base2) 
                      
              datahat = np.fft.fft(data2)
              self.fftdata = datahat
      
              for scaleindex from 0 <= scaleindex < self.nscale:
                 currentscale = np.asarray(self.scales[scaleindex])
                 self.currentscale = currentscale  # for internal use
                 s_omega = omega*currentscale
                 psihat = self._wf(s_omega) * sqrt(2.0*pi*currentscale)
                 self.cwt[scaleindex,0:ndata,i] = np.fft.ifft(psihat * datahat)        
           return
        
        else:
           return        
        

    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef int _log2(self, double x):
        '''
        utility function to return (integer) log2
        '''
        with nogil:
           return int(log(x+0.0001)/ log(2.0)+0.0001)
        
    
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef int _setscales(self, int ndata, int largestscale, int notes):
        """
        returns a log scale based on notes per ocave
        """
        cdef int noctave = self._log2( ndata/largestscale/2 )
        self.nscale = notes*noctave
        cdef np.ndarray[np.float64_t, ndim=1] scales = np.empty(self.nscale,np.float64)        
        
        self.scales = scales
        for j from 0 <= j < self.nscale:
             self.scales[j] = ndata/(self.scale*(2.0**(self.nscale-1-j)/notes))
        return 0            
 
        
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef np.ndarray _wf(self, np.ndarray s_omega):
       """
       Morlet mother wavelet
       """    
       cdef np.ndarray[np.int64_t, ndim=1] H = np.ones(np.shape(s_omega), np.int64)
               
       return 0.75112554*( np.exp(-(s_omega-6.0)**2/2.0))*H
      
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef np.ndarray _pad2nxtpow2(self, np.ndarray data, double base2):
       """
       zero pad numpy array up to next power 2
       """
       cdef np.ndarray[np.float64_t, ndim=2] Y = np.zeros((1, np.int(2**(base2+1)) ), np.float64)
       
       Y.flat[np.arange(self.win)] = data
       return np.squeeze(Y)
      
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef list _column(self, np.ndarray matrix, int i):
       """
       return a column from a matrix
       """
       return [row[i] for row in matrix]          

    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef np.ndarray _getwave(self):
        """
        return power spectra
        """
        #cdef np.ndarray[np.float64_t, ndim=3] wave = np.empty((self.nscale,self.win,len(self.r)), np.float64)
        cdef np.ndarray[np.float64_t, ndim=3] wave = np.empty((self.nscale,self.win,self.r), np.float64)        
        #for i from 0 <= i < len(self.r):
        for i from 0 <= i < self.r:            
           wave[:,:,i] = np.tile(self.scales**-1, (self.win,1)).T*(abs(self.cwt[:,0:self.win,i])**2)
        return wave
        
    # =========================================================
    @cython.boundscheck(False)   
    @cython.cdivision(True)
    @cython.wraparound(False)
    @cython.nonecheck(False)
    cpdef double getvar(self):
        """
        return variance of wave
        """
        cdef float pi = 3.14159265
        cdef np.ndarray n
        cdef np.ndarray[np.float64_t, ndim=1] dat= np.empty(len(self.scales), np.float64)  
                                 
        if self.docalc == 1:

           n = np.r_[0:len(self.scales)]-(len(self.scales)-1)/2
           wave = self._getwave()
        
           dat = np.var(np.var(wave.T,axis=1),axis=0)

           dat = dat/np.sum(dat) * np.exp(-(0.5)*((pi/2)*n/((len(self.scales)-1)/2))**2)           
           dat = dat/np.sum(dat)

           dat = dat/(self.scales**2)
           dat = dat/np.sum(dat)
           
           #return sqrt(np.sum(dat*((self.scales- np.sum(dat*self.scales) )**2)))
           return np.sum(dat*self.scales)
        else:
           return np.nan
        
