
import os
import numpy as np
#dtype = 'int16'
#string = '_data_port.dat'
    

# =========================================================
def get_mmap_data(sonpath, base, string, dtype, shape):
    #we are only going to access the portion of memory required
    with open(os.path.normpath(os.path.join(sonpath,base+string)), 'r') as ff:
       fp = np.memmap(ff, dtype=dtype, mode='r', shape=shape)
    return fp      
  
# =========================================================  
def set_mmap_data(sonpath, base, string, dtype, Zt):
    # create memory mapped file for Z
    #with open(os.path.normpath(os.path.join(sonpath,base+string)), 'w+') as ff:
    #   fp = np.memmap(ff, dtype=dtype, mode='w+', shape=np.shape(Zt))
    try:
       os.remove(os.path.normpath(os.path.join(sonpath,base+string)))
    except:
       pass

    try:
       with open(os.path.normpath(os.path.join(sonpath,base+string)), 'w+') as ff:
          fp = np.memmap(ff, dtype=dtype, mode='readwrite', shape=np.shape(Zt))
       fp[:] = Zt[:]

    except:
       with open(os.path.normpath(os.path.join(sonpath,base+string)), 'w+') as ff:
          fp = np.memmap(ff, dtype=dtype, mode='copyonwrite', shape=np.shape(Zt))
       fp[:] = Zt[:]

    del fp
    shape = np.shape(Zt)
    del Zt
    return shape      
