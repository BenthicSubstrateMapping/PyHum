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
            shutil.copy(src, dst)
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
   humfile = os.path.expanduser("~")+os.sep+'pyhum_test'+os.sep+'test.DAT'+os.sep #PyHum.__path__[0]+os.sep+'test.DAT'
   sonpath = os.path.expanduser("~")+os.sep+'pyhum_test'+os.sep #PyHum.__path__[0]
   doplot = 1 #yes

   # reading specific settings
   cs2cs_args = "epsg:26949" #arizona central state plane
   bedpick = 1 # auto bed pick
   c = 1450 # speed of sound fresh water
   t = 0.108 # length of transducer
   f = 455 # frequency kHz
   draft = 0.3 # draft in metres
   flip_lr = 1 # flip port and starboard

   # correction specific settings
   maxW = 1000 # rms output wattage

   # for texture calcs
   win = 50 # pixel window
   shift = 10 # pixel shift
   density = win/2 
   numclasses = 4 # number of discrete classes for contouring and k-means
   maxscale = 20 # Max scale as inverse fraction of data length (for wavelet analysis)
   notes = 4 # Notes per octave (for wavelet analysis)

   # for mapping
   #imagery = 1 # server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery'
   dogrid = 1 # yes
   calc_bearing = 0 #no
   filt_bearing = 1 #yes
   res = 0.05 # grid resolution in metres
   chunk_size = 0 # auto chunk size

   PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunk_size)

   PyHum.correct(humfile, sonpath, maxW, doplot)

   PyHum.texture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

   PyHum.map(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res)

   res = 0.5 # grid resolution in metres
   
   PyHum.map_texture(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res)

if __name__ == '__main__':
   dotest()

