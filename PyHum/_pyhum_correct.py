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

#"""
# ____        _   _                         
#|  _ \ _   _| | | |_   _ _ __ ___    _   _ 
#| |_) | | | | |_| | | | | '_ ` _ \  (_) (_)
#|  __/| |_| |  _  | |_| | | | | | |  _   _ 
#|_|    \__, |_| |_|\__,_|_| |_| |_| (_) (_)
#       |___/                               
#
#                                  __ 
#  _________  _____________  _____/ /_
# / ___/ __ \/ ___/ ___/ _ \/ ___/ __/
#/ /__/ /_/ / /  / /  /  __/ /__/ /_  
#\___/\____/_/  /_/   \___/\___/\__/  
#                                     
#
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
#+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#|d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
#|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+

#"""

# =========================================================
# ====================== libraries ======================
# =========================================================

#operational
from __future__ import division
from scipy.io import savemat, loadmat
import os, time, sys, getopt
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory
except:
   pass
from joblib import Parallel, delayed, cpu_count

#numerical
import numpy as np
import PyHum.utils as humutils
#from scipy.stats import nanmean, nanmedian
import ppdrc
#from scipy.ndimage.filters import median_filter

#plotting
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

import warnings
warnings.filterwarnings("ignore")


# =========================================================
# =============== begin program ======================
# ========================================================

__all__ = [
    'correct',
    'custom_save',
    'remove_water',
    'correct_scans',
    ]

#################################################
def correct(humfile, sonpath, maxW=1000, doplot=1, dofilt=0, correct_withwater=0, ph = 7, temp = 10, salinity = 0, dconcfile = None):

    '''
    Remove water column and carry out some rudimentary radiometric corrections, 
    accounting for directivity and attenuation with range

    Syntax
    ----------
    [] = PyHum.correct(humfile, sonpath, maxW, doplot, correct_withwater, ph, temp, salinity, dconcfile)

    Parameters
    ----------
    humfile : str
       path to the .DAT file

    sonpath : str
       path where the *.SON files are

    maxW : int, *optional* [Default=1000]
       maximum transducer power

    doplot : int, *optional* [Default=1]
       1 = make plots, otherwise do not

    dofilt : int, *optional* [Default=0]
       1 = apply a phase preserving filter to the scans

    correct_withwater : int, *optional* [Default=0]
       1 = apply radiometric correction but don't remove water column from scans

    ph : float, *optional* [Default=7.0]
       water acidity in pH

    temp : float, *optional* [Default=10.0]
       water temperature in degrees Celsius

    salinity : float, *optional* [Default=0.0]
       salinity of water in parts per thousand

    dconcfile : str, *optional* [Default=None]
       file path of a text file containing sediment concentration data
       this file must contain the following fields separated by spaces:
       size (microns) conc (mg/L) dens (kg/m3)
       with one row per grain size, for example:
       30 1700 2200
       100 15 2650

    Returns
    -------
    sonpath+base+'_data_star_l.dat': memory-mapped file
        contains the starboard scan with water column removed

    sonpath+base+'_data_port_l.dat': memory-mapped file
        contains the portside scan with water column removed

    sonpath+base+'_data_star_la.dat': memory-mapped file
        contains the starboard scan with water column removed and 
        radiometrically corrected

    sonpath+base+'_data_port_la.dat': memory-mapped file
        contains the portside scan with water column removed and
        radiometrically corrected

    sonpath+base+'_data_range.dat': memory-mapped file
        contains the cosine of the range which is used to correct
        for attenuation with range

    sonpath+base+'_data_dwnlow_l.dat': memory-mapped file
        contains the low freq. downward scan with water column removed

    sonpath+base+'_data_dwnhi_l.dat': memory-mapped file
        contains the high freq. downward  scan with water column removed

    sonpath+base+'_data_dwnlow_la.dat': memory-mapped file
        contains the low freq. downward  scan with water column removed and 
        radiometrically corrected

    sonpath+base+'_data_dwnhi_la.dat': memory-mapped file
        contains the high freq. downward  scan with water column removed and
        radiometrically corrected
    
    if correct_withwater == 1:
    
       sonpath+base+'_data_star_lw.dat': memory-mapped file
           contains the starboard scan with water column retained and 
           radiometrically corrected

       sonpath+base+'_data_port_lw.dat': memory-mapped file
           contains the portside scan with water column retained and
           radiometrically corrected

    '''

    # prompt user to supply file if no input file given
    if not humfile:
      print 'An input file is required!!!!!!'
      Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
      inputfile = askopenfilename(filetypes=[("DAT files","*.DAT")]) 

    # prompt user to supply directory if no input sonpath is given
    if not sonpath:
      print 'A *.SON directory is required!!!!!!'
      Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
      sonpath = askdirectory() 

    # print given arguments to screen and convert data type where necessary
    if humfile:
      print 'Input file is %s' % (humfile)

    if sonpath:
      print 'Sonar file path is %s' % (sonpath)

    if maxW:
      maxW = np.asarray(maxW,float)
      print 'Max. transducer power is %s W' % (str(maxW))

    if doplot:
      doplot = int(doplot)
      if doplot==0:
         print "Plots will not be made"

    if dofilt:
      dofilt = int(dofilt)
      if dofilt==0:
         print "Phase preserving filter will not be applied"
      else:
         print "Phase preserving filter will be applied"

    if correct_withwater:
      correct_withwater = int(correct_withwater)
      if correct_withwater==1:
         print "Correction will be applied without removing water column"

    if salinity:
       salinity = np.asarray(salinity,float)
       print 'Salinity is %s ppt' % (str(salinity))

    if ph:
       ph = np.asarray(ph,float)
       print 'pH is %s' % (str(ph))

    if temp:
       temp = np.asarray(temp,float)
       print 'Temperature is %s' % (str(temp))

    if dconcfile is not None:
       print 'Suspended sediment size/conc. file is %s' % (dconcfile)
       dconc = np.genfromtxt(dconcfile).T
       conc = dconc[1]
       dens = dconc[2]
       d = dconc[0]

    #================================
    # start timer
    if os.name=='posix': # true if linux/mac or cygwin on windows
       start = time.time()
    else: # windows
       start = time.clock()

    # if son path name supplied has no separator at end, put one on
    if sonpath[-1]!=os.sep:
       sonpath = sonpath + os.sep

    base = humfile.split('.DAT') # get base of file name for output
    base = base[0].split(os.sep)[-1]

    # remove underscores, negatives and spaces from basename
    base = humutils.strip_base(base)   

    # add wattage to metadata dict 
    meta = loadmat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')))

    dep_m = meta['dep_m'][0]
    pix_m = meta['pix_m'][0]

    meta['maxW'] = maxW
    savemat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')), meta ,oned_as='row')

    bed = np.squeeze(meta['bed'])
    ft = 1/(meta['pix_m'])
    dist_m = np.squeeze(meta['dist_m'])

    if dconcfile is not None:
       # sediment attenuation
       alpha = sed_atten(meta['f'],conc,dens,d,meta['c'])
    else:
       alpha = 0

    # load memory mapped scans
    shape_port = np.squeeze(meta['shape_port'])
    if shape_port!='':
       
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_port2.dat'))):
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), 'r') as ff:
             port_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_port))

       else:
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port.dat')), 'r') as ff:
             port_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_port))

    shape_star = np.squeeze(meta['shape_star'])
    if shape_star!='':
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_star2.dat'))):
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star2.dat')), 'r') as ff:
             star_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_star))
       else:
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star.dat')), 'r') as ff:
             star_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_star))

    if len(shape_star)==2:
       extent = shape_star[0] 
    else:
       extent = shape_star[1] #np.shape(data_port)[0]

    bed = np.asarray(bed,'int')+int(0.25*ft)

    # calculate in dB
    ######### star
    Zt, R, A = remove_water(star_fp, bed, shape_star, dep_m, pix_m, 1,  maxW)

    Zt = np.squeeze(Zt)
    R = np.squeeze(R)
    A = np.squeeze(A)

    R[np.isnan(R)] = 0

    alpha_w = water_atten(R,meta['f'],meta['c'], ph, temp, salinity)

    # compute transmission losses
    TL = (40 * np.log10(R) + alpha_w + (2*alpha)*R/1000)/255

    TL[np.isnan(TL)] = 0
    TL[TL<0] = 0

    # create memory mapped file for Z)
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_l.dat')), 'w+') as ff:
       fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
    fp[:] = Zt[:]
    del fp
    shape_star = np.shape(Zt)
    del Zt

    # create memory mapped file for R
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_range.dat')), 'w+') as ff:
       fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(R))
    fp[:] = R[:]
    del fp
    del R


    # create memory mapped file for A
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_incidentangle.dat')), 'w+') as ff:
       fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(A))
    fp[:] = A[:]
    del fp
    del A

    with open(os.path.normpath(os.path.join(sonpath,base+'_data_incidentangle.dat')), 'r') as ff:
       A_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_star)

    if correct_withwater == 1:
       Zt = correct_scans(star_fp, A_fp, TL, dofilt)

       # create memory mapped file for Z)
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_lw.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_star = np.shape(Zt)
       del Zt

    #we are only going to access the portion of memory required
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_l.dat')), 'r') as ff:
       star_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_star)

    Zt = correct_scans(star_fp, A_fp, TL, dofilt)

    Zt = np.squeeze(Zt)

    # create memory mapped file for Z
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')), 'w+') as ff:
       fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
    fp[:] = Zt[:]
    del fp
    shape_star = np.shape(Zt)
    del Zt
    #we are only going to access the portion of memory required
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_la.dat')), 'r') as ff:
       star_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_star)


    ######### port
    if correct_withwater == 1:
       Zt = correct_scans(port_fp, A_fp, TL, dofilt)

       # create memory mapped file for Z)
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_lw.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_port = np.shape(Zt)
       del Zt

    Zt = remove_water(port_fp, bed, shape_port, dep_m, pix_m, 0,  maxW)

    Zt = np.squeeze(Zt)

    # create memory mapped file for Z
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_l.dat')), 'w+') as ff:
       fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
    fp[:] = Zt[:]
    del fp
    shape_port = np.shape(Zt)
    del Zt
    #we are only going to access the portion of memory required
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_l.dat')), 'r') as ff:
       port_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_port)

    Zt = correct_scans(port_fp, A_fp, TL, dofilt)

    # create memory mapped file for Z
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')), 'w+') as ff:
       fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
    fp[:] = Zt[:]
    del fp
    shape_port = np.shape(Zt)
    del Zt
    #we are only going to access the portion of memory required
    with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_la.dat')), 'r') as ff:
       port_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_port)

    ## do plots of merged scans
    if doplot==1:
       if correct_withwater == 1:

          with open(os.path.normpath(os.path.join(sonpath,base+'_data_port_lw.dat')), 'r') as ff:
             port_fpw = np.memmap(ff, dtype='float32', mode='r', shape=shape_port)


          with open(os.path.normpath(os.path.join(sonpath,base+'_data_star_lw.dat')), 'r') as ff:
              star_fpw = np.memmap(ff, dtype='float32', mode='r', shape=shape_star)

          if len(np.shape(star_fpw))>2:
             for p in xrange(len(star_fpw)):
                plot_merged_scans(port_fpw[p], star_fpw[p], dist_m, shape_port, ft, sonpath, p)
          else:
             plot_merged_scans(port_fpw, star_fpw, dist_m, shape_port, ft, sonpath, 0)

       else:

          if len(np.shape(star_fp))>2:
             for p in xrange(len(star_fp)):
                plot_merged_scans(port_fp[p], star_fp[p], dist_m, shape_port, ft, sonpath, p)
          else:
             plot_merged_scans(port_fp, star_fp, dist_m, shape_port, ft, sonpath, 0)


    # load memory mapped scans
    shape_low = np.squeeze(meta['shape_low'])

    if shape_low!='':
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow2.dat'))):
          try:
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow2.dat')), 'r') as ff:
                low_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_low))

          except:
             if 'shape_hi' in locals():
                with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow2.dat')), 'r') as ff:
                   low_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_hi))

       else:

          try:
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')), 'r') as ff:
                low_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_low))

          except:
             if 'shape_hi' in locals():
                with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow.dat')), 'r') as ff:
                   low_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_hi))

    shape_hi = np.squeeze(meta['shape_hi'])

    if shape_hi!='':
       if os.path.isfile(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi2.dat'))):
          try:
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
                hi_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_hi))

          except:
             if 'shape_low' in locals():
                with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
                   hi_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_low))

       else:
          try:
             with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
                hi_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_hi))

          except:
             if 'shape_low' in locals():
                with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
                   hi_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_low))


    if 'low_fp' in locals():
       ######### low
       Zt = remove_water(low_fp, bed, shape_low, dep_m, pix_m, 0,  maxW)
       Zt = np.squeeze(Zt)

       # create memory mapped file for Z
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow_l.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_low = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow_l.dat')), 'r') as ff:
          low_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_low)

       Zt = correct_scans2(low_fp, TL)

       # create memory mapped file for Z
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow_la.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_low = np.shape(Zt)
       del Zt
       #we are only going to access the lowion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnlow_la.dat')), 'r') as ff:
          low_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_low)


       if doplot==1:
          if len(np.shape(low_fp))>2:
             for p in xrange(len(low_fp)):
                plot_dwnlow_scans(low_fp[p], dist_m, shape_low, ft, sonpath, p)
          else:
             plot_dwnlow_scans(low_fp, dist_m, shape_low, ft, sonpath, 0)

          #for p in xrange(len(low_fp)):
          #   plot_dwnlow_scans(low_fp[p], dist_m, shape_low, ft, sonpath, p)

    if 'hi_fp' in locals():
       ######### hi
       Zt = remove_water(hi_fp, bed, shape_hi, dep_m, pix_m, 0,  maxW)
       Zt = np.squeeze(Zt)

       # create memory mapped file for Z
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi_l.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_hi = np.shape(Zt)
       del Zt
       #we are only going to access the portion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi_l.dat')), 'r') as ff:
          hi_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_hi)

       Zt = correct_scans2(hi_fp, TL)

       # create memory mapped file for Z
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi_la.dat')), 'w+') as ff:
          fp = np.memmap(ff, dtype='float32', mode='w+', shape=np.shape(Zt))
       fp[:] = Zt[:]
       del fp
       shape_hi = np.shape(Zt)
       del Zt
       #we are only going to access the hiion of memory required
       with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi_la.dat')), 'r') as ff:
          hi_fp = np.memmap(ff, dtype='float32', mode='r', shape=shape_hi)

       if doplot==1:
          if len(np.shape(hi_fp))>2:
             for p in xrange(len(hi_fp)):
                plot_dwnhi_scans(hi_fp[p], dist_m, shape_hi, ft, sonpath, p)
          else:
             plot_dwnhi_scans(hi_fp, dist_m, shape_hi, ft, sonpath, 0)

    if os.name=='posix': # true if linux/mac
       elapsed = (time.time() - start)
    else: # windows
       elapsed = (time.clock() - start)
    print "Processing took ", elapsed , "seconds to analyse"

    print "Done!"


# =========================================================
def water_atten(H,f,c,pH,T,S):
   '''
   calculate absorption of sound in water.
   '''
   H = np.abs(H)
   P1 = 1 # cosntant
   A1 = (8.86/c)*(10**(0.78*pH - 5))
   f1 = 2.8*(S/35)**0.5 * 10**(4 - 1245/(T + 273))
   A2 = 21.44*(S/c)*(1 + 0.025*T)
   A3 = (4.937 *10**-4) - (2.59 * 10**-5)*T + (9.11* 10**-7)*T**2- (1.5 * 10**-8)*T**3
   f2 = (8.17 * 10**(8 - 1990/(T + 273))) / (1 + 0.0018*(S - 35))
   P2= 1 - (1.37*10**-4) * H + (6.2 * 10**-9)* H**2
   P3 = 1 - (3.83 * 10**-5)*H + (4.9 *10**(-10) )* H**2
   # absorption sound water dB/km
   alphaw = ( (A1*P1*f1*f**2)/(f**2 + f1**2) ) + ( (A2*P2*f2*f**2)/(f**2 + f2**2) ) + (A3*P3*f**2)
   return 2*(alphaw/1000)*H # depth(m) * dB/m = dB


# =========================================================
def sed_atten(f,conc,dens,d,c):
   '''
   calculates the attenuation due to sediment, in dB/m
   given sediment concentration, density, grain size, frequency and speed of sound
   according to Urick (1948), JASA
   http://www.rdinstruments.com/pdfs/Use-ADCP-suspended-sediment%20discharge-NY.pdf
   http://rspa.royalsocietypublishing.org/content/459/2037/2153.full.pdf
   example values
   f = 400  freq, kHz
   c = 1490 speed sound in water, m/s
   d = [40, 100] microns
   dens = [2000, 2650] sediment density, kg/m^3
   conc = [1000, 100] mg/L
   '''

   if np.any(conc)>0:
      f = f * 1000 # frequency, Hz
      sigma = dens/1000 # ratio sediment to fluid density
      d = d/1e6 # particle diameter, m
      nu = 1.004e-6 # viscosity fresh water, m^2/s
      lam = c/f # acoustic wavelength, m
      k = (2*np.pi)/lam # acoustic wavenumber 
      w = (2*np.pi)*f # radian frequency
      delta_v = np.sqrt(2*nu/w)
      phi = (conc/1e6)/dens #sediment volume fraction
      a = d/2 # particle radius, m
      tau = (1/2) + (9/4)*(delta_v/a)
      s = (9/4)*(delta_v/a)*(1+(delta_v/a))
      alpha = phi*( (1/6) *k**4 *a**3 + k*(sigma-1)**2 *( s/( s**2+(sigma+tau)**2 ) ) )*1e4
      return np.sum(alpha) # times 2 because 2 way travel
   else:
      return np.nan

# =========================================================
def custom_save(figdirec,root):
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=400)

# =========================================================
def remove_water(fp,bed,shape, dep_m, pix_m, calcR,  maxW):
    Zt = []
    if calcR==1:
       R = []
       A = []

    if  len(np.shape(fp))>2:
       for p in xrange(len(fp)):
          data_dB = fp[p]*(10*np.log10(maxW)/255)

          Zbed = np.squeeze(bed[shape[-1]*p:shape[-1]*(p+1)])

          # shift proportionally depending on where the bed is
          for k in xrange(np.shape(data_dB)[1]):
             try:
                data_dB[:,k] = np.r_[data_dB[Zbed[k]:,k], np.zeros( (np.shape(data_dB)[0] -  np.shape(data_dB[Zbed[k]:,k])[0] ,) )]
             except:
                data_dB[:,k] = np.ones(np.shape(data_dB)[0])

          Zt.append(data_dB)    

          if calcR ==1:
             extent = shape[1]
             yvec = np.linspace(pix_m,extent*pix_m,extent)
             d = dep_m[shape[-1]*p:shape[-1]*(p+1)]

             a = np.ones(np.shape(fp[p]))
             for k in range(len(d)): 
                a[:,k] = d[k]/yvec

             r = np.ones(np.shape(fp[p]))
             for k in range(len(d)): 
                r[:,k] = np.sqrt(yvec**2 - d[k]**2)

             # shift proportionally depending on where the bed is
             for k in xrange(np.shape(r)[1]):
                try:
                   r[:,k] = np.r_[r[Zbed[k]:,k], np.zeros( (np.shape(r)[0] -  np.shape(r[Zbed[k]:,k])[0] ,) )]
                   a[:,k] = np.r_[a[Zbed[k]:,k], np.zeros( (np.shape(a)[0] -  np.shape(a[Zbed[k]:,k])[0] ,) )]
                except:
                   r[:,k] = np.ones(np.shape(r)[0])
                   a[:,k] = np.ones(np.shape(a)[0])

             R.append(r)
             A.append(a)

    else:

       data_dB = fp*(10*np.log10(maxW)/255)

       Zbed = np.squeeze(bed)

       # shift proportionally depending on where the bed is
       for k in xrange(np.shape(data_dB)[1]):
          try:
            data_dB[:,k] = np.r_[data_dB[Zbed[k]:,k], np.zeros( (np.shape(data_dB)[0] -  np.shape(data_dB[Zbed[k]:,k])[0] ,) )]
          except:
             data_dB[:,k] = np.ones(np.shape(data_dB)[0])

       Zt.append(data_dB)    

       if calcR ==1:
          extent = shape[0]
          yvec = np.linspace(pix_m,extent*pix_m,extent)
          d = dep_m

          a = np.ones(np.shape(fp))
          for k in range(len(d)): 
             a[:,k] = d[k]/yvec

          r = np.ones(np.shape(fp))
          for k in range(len(d)): 
             r[:,k] = np.sqrt(yvec**2 - d[k]**2)

          # shift proportionally depending on where the bed is
          for k in xrange(np.shape(r)[1]):
             try:
                r[:,k] = np.r_[r[Zbed[k]:,k], np.zeros( (np.shape(r)[0] -  np.shape(r[Zbed[k]:,k])[0] ,) )]
                a[:,k] = np.r_[a[Zbed[k]:,k], np.zeros( (np.shape(a)[0] -  np.shape(a[Zbed[k]:,k])[0] ,) )]
             except:
                r[:,k] = np.ones(np.shape(r)[0])
                a[:,k] = np.ones(np.shape(a)[0])

          R.append(r)
          A.append(a)

    if calcR ==1:
       return Zt, R, np.pi/2 - np.arctan(A)
    else:
       return Zt
 
# =========================================================
def correct_scans(fp, a_fp, TL, dofilt):
    return Parallel(n_jobs = cpu_count(), verbose=0)(delayed(c_scans)(fp[p], a_fp[p], TL[p], dofilt) for p in xrange(len(fp)))

# =========================================================
def c_scans(fp, a_fp, TL, dofilt):
   if dofilt==1:
      fp = do_ppdrc(fp, np.shape(fp)[-1]/2)
   #mg = 10**np.log10(np.asarray(fp*np.cos(a_fp),'float32')+0.001)
   mg = 10**np.log10(np.asarray(fp * 1-np.cos(a_fp),'float32')+0.001 + TL)
   mg[fp==0] = np.nan
   mg[mg<0] = np.nan
   return mg

# =========================================================
def correct_scans2(fp, TL):
    return Parallel(n_jobs = cpu_count(), verbose=0)(delayed(c_scans2)(fp[p], TL) for p in xrange(len(fp)))

# =========================================================
def c_scans2(fp, TL):
   try:
      mg = 10**np.log10(np.asarray(fp,'float32')+0.001 + TL[:,::2] )
   except:
      mg = 10**np.log10(np.asarray(fp,'float32')+0.001 )

   mg[fp==0] = np.nan
   mg[mg<0] = np.nan
   return mg

# =========================================================
def do_ppdrc(fp, filtsize):
   dat = fp.astype('float64')
   dat[np.isnan(dat)] = 0
   dat1 = ppdrc.ppdrc(dat, filtsize)
   dat1 = humutils.rescale(dat1.getdata(),np.min(dat),np.max(dat))
   dat1[np.isnan(fp)] = np.nan
   return dat1

# =========================================================
def plot_merged_scans(dat_port, dat_star, dist_m, shape_port, ft, sonpath, p):

   if len(shape_port)>2:
      Zdist = dist_m[shape_port[-1]*p:shape_port[-1]*(p+1)]
      extent = shape_port[1] #np.shape(merge)[0]
   else:
      Zdist = dist_m
      extent = shape_port[0] #np.shape(merge)[0]

   fig = plt.figure()
   plt.imshow(np.vstack((np.flipud(dat_port), dat_star)), cmap='gray', extent=[min(Zdist), max(Zdist), -extent*(1/ft), extent*(1/ft)])
   plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

   plt.axis('normal'); plt.axis('tight')
   custom_save(sonpath,'merge_corrected_scan'+str(p))
   del fig

# =========================================================
def plot_dwnlow_scans(dat_dwnlow, dist_m, shape_low, ft, sonpath, p):

    if len(shape_low)>2:
       Zdist = dist_m[shape_low[-1]*p:shape_low[-1]*(p+1)]
       extent = shape_low[1] #np.shape(merge)[0]
    else:
      Zdist = dist_m
      extent = shape_low[0] #np.shape(merge)[0]  
 
    fig = plt.figure()
    plt.imshow(dat_dwnlow, cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0])
    plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

    plt.axis('normal'); plt.axis('tight')
    custom_save(sonpath,'dwnlow_corrected_scan'+str(p))
    del fig

# =========================================================
def plot_dwnhi_scans(dat_dwnhi, dist_m, shape_hi, ft, sonpath, p):

    if len(shape_hi)>2:
       Zdist = dist_m[shape_hi[-1]*p:shape_hi[-1]*(p+1)]
       extent = shape_hi[1] #np.shape(merge)[0]
    else:
      Zdist = dist_m
      extent = shape_hi[0] #np.shape(merge)[0]  
    
    fig = plt.figure()
    plt.imshow(dat_dwnhi, cmap='gray', extent=[min(Zdist), max(Zdist), extent*(1/ft), 0])
    plt.ylabel('Range (m)'), plt.xlabel('Distance along track (m)')

    plt.axis('normal'); plt.axis('tight')
    custom_save(sonpath,'dwnhi_corrected_scan'+str(p))
    del fig


# =========================================================
# =========================================================
if __name__ == '__main__':

   correct(humfile, sonpath, maxW, doplot)

