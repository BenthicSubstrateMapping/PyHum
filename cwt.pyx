'''
cwt.pyx
Part of PyHum software 

INFO:
Cython script to carry out Morlet wavelet computations

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

# =========================================================
cdef class Cwt:
    """
    continuous Morelet wavelet transform
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
    cdef object density
    cdef object r
                   
    # =========================================================
    def __init__(self, np.ndarray matrix, np.int largestscale, np.int notes, np.int win, np.int density):
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
        self.density = density
        
        cdef np.ndarray r = np.arange(1,self.win-1,self.density, dtype=np.int)
        self.r = r
        
        cdef np.ndarray currentscale 
        cdef np.ndarray s_omega         
        cdef np.ndarray psihat 
        cdef np.ndarray datahat          
        
        cdef np.float base2 = np.float(np.floor(np.log(self.win)/np.log(2) + 0.4999))
                   
        cdef np.int ndata = np.int(2**(base2+1)) #len(data)
        cdef np.int tmp = 0
        self.nscale = tmp
        self.scale = largestscale
        self._setscales(ndata,largestscale,notes)
        cdef np.ndarray cwt = np.zeros((self.nscale,ndata,len(self.r)), np.complex64)
        self.cwt = cwt
        cdef np.ndarray omega = np.array(range(0,np.int(ndata/2))+range(-np.int(ndata/2),0))*(2.0*np.pi/ndata)
        
        for i from 0 <= i < len(self.r):  
           data = np.asarray( self._column(matrix, np.int(self.r[i]) ) )
           data = self._pad2nxtpow2(data - np.mean(data), base2) 
                      
           datahat = np.fft.fft(data)
           self.fftdata = datahat
      
           for scaleindex from 0 <= scaleindex < self.nscale:
              currentscale = np.asarray(self.scales[scaleindex])
              self.currentscale = currentscale  # for internal use
              s_omega = omega*currentscale
              psihat = self._wf(s_omega) * np.sqrt(2.0*np.pi*currentscale)
              self.cwt[scaleindex,0:ndata,i] = np.fft.ifft(psihat * datahat)        
        return

    # =========================================================
    def _log2(self, np.float x):
        '''
        utility function to return (integer) log2
        '''
        return np.int(np.log(x+0.0001)/ np.log(2.0)+0.0001)
        
    # =========================================================        
    def _setscales(self, np.int ndata, np.int largestscale, np.int notes):
        """
        returns a log scale based on notes per ocave
        """
        cdef np.int noctave = self._log2( ndata/largestscale/2 )
        self.nscale = notes*noctave
        cdef np.ndarray scales = np.zeros(self.nscale,np.float)
        self.scales = scales
        for j from 0 <= j < self.nscale:
             self.scales[j] = ndata/(self.scale*(2.0**(np.float(self.nscale-1-j)/notes)))
        return
        
    # =========================================================         
    def _wf(self, np.ndarray s_omega):
       """
       Morlet mother wavelet
       """    
       cdef np.ndarray H = np.ones(np.shape(s_omega), np.int)
       return 0.75112554*( np.exp(-(s_omega-6.0)**2/2.0))*H
      
    # =========================================================      
    def _pad2nxtpow2(self, np.ndarray data, np.float base2):
       """
       zero pad numpy array up to next power 2
       """
       #cdef np.float base2 = np.float(np.floor(np.log(self.win)/np.log(2) + 0.4999))
       cdef np.ndarray Y = np.zeros((1, np.int(2**(base2+1)) ), np.float)
       Y.flat[np.arange(self.win)] = data
       return np.squeeze(Y)
      
    # =========================================================
    def _column(self, np.ndarray matrix, np.int i):
       """
       return a column from a matrix
       """
       return [row[i] for row in matrix]          

    # =========================================================
    def _getwave(self):
        """
        return power spectra
        """
        cdef np.ndarray wave = np.zeros((self.nscale,self.win,len(self.r)), np.float)
        for i from 0 <= i < len(self.r):  
           wave[:,:,i] = np.tile(self.scales**-1, (self.win,1)).T*(np.absolute(self.cwt[:,0:self.win,i])**2)
        return wave
        
    # =========================================================
    def getvar(self):
        """
        return variance of wave
        """
        cdef np.ndarray n = np.r_[0:len(self.scales)]-(len(self.scales)-1)/2
        wave = self._getwave()
        
        cdef np.ndarray dat = np.var(np.var(wave.T,axis=1),axis=0)

        dat = dat/np.sum(dat) * np.exp(-(0.5)*((np.pi/2)*n/((len(self.scales)-1)/2))**2)
        dat = dat/np.sum(dat)
           
        return np.sqrt(np.sum(dat*((self.scales- np.sum(dat*self.scales) )**2)))
        
        
