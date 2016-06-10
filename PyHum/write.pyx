# cython: boundscheck=False
# cython: cdivision=True
# cython: wraparound=False
# cython: nonecheck=False
import numpy as np
cimport numpy as np
cimport cython

# =========================================================
def txtwrite(str outfile, np.ndarray[np.float32_t, ndim=2] towrite):
   '''
   Custom fast numpy array to comma-delimited ASCII txt file

   Syntax
   ----------
   () = write.txtwrite(infile, towrite)

   Parameters
   ------------
   outfile : str
   	name of file to write to
   towrite : ndarray
   	ndarray containing Nx4 point cloud

   Returns
   ----------
   None

   '''
   with open(outfile, 'wb') as f:

      np.savetxt(f, towrite, fmt="%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f") 

