"""
Part of PyHum software 

INFO:


Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.1.4      Revision: Mar, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    
"""
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
   humfile = os.path.expanduser("~")+os.sep+'test.DAT' #PyHum.__path__[0]+os.sep+'test.DAT'
   sonpath = os.path.expanduser("~")+os.sep+'pyhum_test' #PyHum.__path__[0]
   doplot = 1

   # reading specific settings
   cs2cs_args = "epsg:26949"
   draft = 0
   bedpick = 1
   c = 1450
   t = 0.108
   f = 455
   draft = 0.3
   flip_lr = 1

   # correction specific settings
   maxW = 1000

   # for texture calcs
   win = 100
   shift = 10
   density = win/2
   numclasses = 4
   maxscale = 20
   notes = 4
   shorepick = 0
   do_two = 0

   # for mapping
   imagery = 1 # server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery'

   PyHum.humread(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr)

   PyHum.humcorrect(humfile, sonpath, maxW, doplot)

   PyHum.humtexture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes, shorepick, do_two)

   PyHum.domap(humfile, sonpath, cs2cs_args, imagery)

if __name__ == '__main__':
   dotest()

