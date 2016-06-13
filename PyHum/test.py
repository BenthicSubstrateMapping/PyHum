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

#python -c "import PyHum; PyHum.test.dotest()"

import PyHum
import os
import shutil
import errno
 
def dircopy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dest)
        else:
            print('Directory not copied. Error: %s' % e)

__all__ = [
    'dotest',
    ]

def dotest():

   # copy files over to somewhere read/writeable
   dircopy(PyHum.__path__[0], os.path.expanduser("~")+os.sep+'pyhum_test')
   shutil.copy(PyHum.__path__[0]+os.sep+'test.DAT', os.path.expanduser("~")+os.sep+'pyhum_test'+os.sep+'test.DAT')

   # general settings   
   humfile = os.path.normpath(os.path.join(os.path.expanduser("~"),'pyhum_test','test.DAT'))
   sonpath = os.path.normpath(os.path.join(os.path.expanduser("~"),'pyhum_test'))

   doplot = 1 #yes

   # reading specific settings
   cs2cs_args = "epsg:26949" #arizona central state plane
   bedpick = 1 # auto bed pick
   c = 1450 # speed of sound fresh water
   t = 0.108 # length of transducer
   draft = 0.3 # draft in metres
   flip_lr = 1 # flip port and starboard
   model = 998 # humminbird model
   calc_bearing = 1 #1=yes
   filt_bearing = 1 #1=yes
   chunk = 'd100' # distance, 100m
   #chunk = 'p1000' # pings, 1000
   #chunk = 'h10' # heading deviation, 10 deg
          
   # correction specific settings
   maxW = 1000 # rms output wattage
   dofilt = 0 # 1 = apply a phase preserving filter (WARNING!! takes a very long time for large scans)
   correct_withwater = 0 # don't retain water column in radiometric correction (1 = retains water column for radiomatric corrections)
   ph = 7.0 # acidity on the pH scale
   temp = 10.0 # water temperature in degrees Celsius
   salinity = 0.0

   # for shadow removal
   shadowmask = 0 #automatic shadow removal
   win = 31

   # for texture calcs
   shift = 10 # pixel shift
   density =win/2 # win/2 
   numclasses = 4 # number of discrete classes for contouring and k-means
   maxscale = 20 # Max scale as inverse fraction of data length (for wavelet analysis)
   notes = 4 # Notes per octave (for wavelet analysis)

   # for mapping
   res = 0.25 #99 # grid resolution in metres
   # if res==99, the program will automatically calc res from the spatial res of the scans
   mode = 1 # gridding mode (simple nearest neighbour)
   #mode = 2 # gridding mode (inverse distance weighted nearest neighbour)
   #mode = 3 # gridding mode (gaussian weighted nearest neighbour)
   use_uncorrected = 0

   nn = 64 #number of nearest neighbours for gridding (used if mode > 1)
   ##influence = 1 #Radius of influence used in gridding. Cut off distance in meters 
   numstdevs = 5 #Threshold number of standard deviations in sidescan intensity per grid cell up to which to accept 

   # for downward-looking echosounder echogram (e1-e2) analysis
   beam = 20.0
   transfreq = 200.0 # frequency (kHz) of downward looking echosounder
   integ = 5
   numclusters = 3 # number of acoustic classes to group observations

   ## read data in SON files into PyHum memory mapped format (.dat)
   PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, bedpick, flip_lr, model, calc_bearing, filt_bearing, chunk) #cog

   ## correct scans and remove water column
   PyHum.correct(humfile, sonpath, maxW, doplot, dofilt, correct_withwater, ph, temp, salinity)

   ## remove acoustic shadows (caused by distal acoustic attenuation or sound hitting shallows or shoreline)
   PyHum.rmshadows(humfile, sonpath, win, shadowmask, doplot)

   win = 100 # pixel window
   
   ## Calculate texture lengthscale maps using the method of Buscombe et al. (2015)
   PyHum.texture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

   ## grid and map the scans
   PyHum.map(humfile, sonpath, cs2cs_args, res, mode, nn, numstdevs, use_uncorrected) #dowrite, 

   res = 1 # grid resolution in metres
   numstdevs = 5
   
   ## grid and map the texture lengthscale maps
   PyHum.map_texture(humfile, sonpath, cs2cs_args, res, mode, nn, numstdevs)

   ## calculate and map the e1 and e2 acoustic coefficients from the downward-looking sonar
   PyHum.e1e2(humfile, sonpath, cs2cs_args, ph, temp, salinity, beam, transfreq, integ, numclusters, doplot)
   
   #res = 0
   #nn = 5
   #noisefloor = 10
   
   ## create mosaic out of all chunks with weighting according to distance from nadir, grazing angle, or both
   #PyHum.mosaic(humfile, sonpath, cs2cs_args, res, nn, noisefloor)

if __name__ == '__main__':
   dotest()

