
import os

#dtype = 'int16'
#string = '_data_port.dat'

def get_mmap_data(sonpath, base, string, dtype, shape):
    #we are only going to access the portion of memory required
    with open(os.path.normpath(os.path.join(sonpath,base+string)), 'r') as ff:
       fp = np.memmap(ff, dtype=dtype, mode='r', shape=shape)
    return fp      
