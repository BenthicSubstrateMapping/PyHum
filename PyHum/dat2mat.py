## PyHum (Python program for Humminbird(R) data processing) 
## has been developed at the Grand Canyon Monitoring & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
## Project homepage: <https://github.com/dbuscombe-usgs/PyHum>
##
##This software is in the public domain because it contains materials that originally came from 
##the United States Geological Survey, an agency of the United States Department of Interior. 
##For more information, see the official USGS copyright policy at 
##http://www.usgs.gov/visual-id/credit_usgs.html#copyright
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.


from Tkinter import Tk
from tkFileDialog import askopenfilename#, askdirectory
#import tkMessageBox
import os
from scipy.io import loadmat, savemat
import numpy as np

def get_mmap_data(filename, dtype, shape):
    with open(filename, 'r') as ff:
       fp = np.memmap(ff, dtype=dtype, mode='r', shape=shape)
    return fp  

def do_savemat(fp, datfile):
   for k in xrange(len(fp)):
      savemat(datfile.split('.dat')[0]+'_'+str(k)+'.mat', mdict = {os.path.basename(datfile).split('.dat')[0]:np.squeeze(fp[k])})


Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
humfile = askopenfilename(title="Select the .DAT file?", filetypes=[("DAT file:","*.DAT")]) 

datfiles = askopenfilename(title="Which file do you wish to convert?", filetypes=[("dat file:","*.dat")], multiple=True) 

sonpath = os.path.dirname(datfiles[0])

# if son path name supplied has no separator at end, put one on
if sonpath[-1]!=os.sep:
   sonpath = sonpath + os.sep

base = humfile.split('.DAT') # get base of file name for output
base = base[0].split(os.sep)[-1]

# remove underscores, negatives and spaces from basename
if base.find('_')>-1:
   base = base[:base.find('_')]
if base.find('-')>-1:
   base = base[:base.find('-')]
if base.find(' ')>-1:
   base = base[:base.find(' ')]
if base.find('.')>-1:
   base = base[:base.find('.')]

# load metadata array for our data array shapes
meta = loadmat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')))

for datfile in datfiles:
   if datfile.find("star")>0:
      shape_star = np.squeeze(meta['shape_star'])
      if datfile.find("_l")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_star))
      elif datfile.find("_la")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_star))
      elif datfile.find("_lw")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_star))
      elif datfile.find("_lar")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_star))
      else:
         fp = get_mmap_data(datfile, 'int16', tuple(shape_star))
      do_savemat(fp, datfile)

   elif datfile.find("port")>0:
      shape_port = np.squeeze(meta['shape_port']) 
      if datfile.find("_l")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_port))
      elif datfile.find("_la")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_port))
      elif datfile.find("_lw")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_port))
      elif datfile.find("_lar")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_port))
      else:
         fp = get_mmap_data(datfile, 'int16', tuple(shape_port))
      do_savemat(fp, datfile)

   elif datfile.find("dwnhi")>0:
      shape_hi = np.squeeze(meta['shape_hi'])
      if datfile.find("_l")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_hi))
      elif datfile.find("_la")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_hi))
      else:
         fp = get_mmap_data(datfile, 'int16', tuple(shape_hi))

   elif datfile.find("dwnlow")>0: 
      shape_low = np.squeeze(meta['shape_low'])
      if datfile.find("_l")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_low))
      elif datfile.find("_la")>0:
         fp = get_mmap_data(datfile, 'float32', tuple(shape_low))
      else:
         fp = get_mmap_data(datfile, 'int16', tuple(shape_low))
      do_savemat(fp, datfile)

   elif datfile.find("incidentangle")>0: 
      shape_star = np.squeeze(meta['shape_star'])
      fp = get_mmap_data(datfile, 'float32', tuple(shape_star))
      do_savemat(fp, datfile)

   elif datfile.find("range")>0: 
      shape_star = np.squeeze(meta['shape_star'])
      fp = get_mmap_data(datfile, 'float32', tuple(shape_star))
      do_savemat(fp, datfile)

   elif datfile.find("class")>0: 
      shape_star = np.squeeze(meta['shape_star'])
      shape_port = np.squeeze(meta['shape_port'])

      if len(shape_star)>2:
         shape = shape_port.copy()
         shape[1] = shape_port[1] + shape_star[1]
      else:
         shape = []
         shape.append(1)
         shape.append(shape_port[0])
         shape.append(shape_port[1])
         shape[1] = shape_port[0] + shape_star[0]

      fp = get_mmap_data(datfile, 'float32', tuple(shape))
      do_savemat(fp, datfile)



