'''
pyhum_map.py
Part of PyHum software 

INFO:

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.1.2      Revision: Mar, 2015

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 and 1198 series instruments. 
'''

# =========================================================
# ====================== libraries ======================
# =========================================================

# operational
from __future__ import division
from scipy.io import loadmat
import os, time, sys, getopt
from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory

# numerical
import numpy as np
import pyproj
import PyHum.utils as humutils
#from pyhum_utils import ascol, runningMeanFast, rescale

# plotting
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import simplekml

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')

__all__ = [
    'domap',
    'custom_save',
    ]

#################################################
def domap(humfile, sonpath, cs2cs_args, imagery):
         
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
      if cs2cs_args:
         print 'cs2cs arguments are %s' % (cs2cs_args)
      if imagery:
         imagery = int(imagery)
         if imagery==0:
            print "ESRI_Imagery_World_2D will be used"

      if not cs2cs_args:
         # arguments to pass to cs2cs for coordinate transforms
         cs2cs_args = "epsg:26949"
         print '[Default] cs2cs arguments are %s' % (cs2cs_args)

      if not imagery:
         if imagery != 0:
            imagery = 1
            print "[Default] World imagery will be used"

      trans =  pyproj.Proj(init=cs2cs_args)

      # if son path name supplied has no separator at end, put one on
      if sonpath[-1]!=os.sep:
         sonpath = sonpath + os.sep

      base = humfile.split('.DAT') # get base of file name for output
      base = base[0].split('/')[-1]

      try:
         esi = np.squeeze(loadmat(sonpath+base+'meta.mat')['e']) #+ 395
         nsi = np.squeeze(loadmat(sonpath+base+'meta.mat')['n']) #- 58
         pix_m = np.squeeze(loadmat(sonpath+base+'meta.mat')['pix_m'])
         bearing = np.squeeze(loadmat(sonpath+base+'meta.mat')['heading'])
         port_la = np.asarray(np.squeeze(loadmat(sonpath+base+'port_la.mat')['port_mg_la']),'float16')
         star_la = np.asarray(np.squeeze(loadmat(sonpath+base+'star_la.mat')['star_mg_la']),'float16')
      except:
         esi = np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['e']) #+ 395
         nsi = np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['n']) #- 58
         pix_m = np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['pix_m'])
         bearing = np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'meta.mat')['heading'])
         port_la = np.asarray(np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'port_la.mat')['port_mg_la']),'float16')
         star_la = np.asarray(np.squeeze(loadmat(os.path.expanduser("~")+os.sep+base+'star_la.mat')['star_mg_la']),'float16')

      # reported bearing by instrument (Kalman filtered?)

      bearing = humutils.runningMeanFast(bearing, len(bearing)/100)
      theta = np.asarray(bearing, 'float')/(180/np.pi)


      #merge = np.vstack((np.flipud(port_la),star_la))
      merge = np.vstack((port_la,star_la))
      merge[np.isnan(merge)] = 0

      del port_la, star_la

      # get number pixels in scan line
      extent = int(np.shape(merge)[0]/2)

      yvec = np.linspace(pix_m,extent*pix_m,extent)

      print "getting point cloud ..."
      # get the points by rotating the [x,y] vector so it lines up with boat heading, assumed to be the same as the curvature of the [e,n] trace
      X=[]; Y=[];
      for k in range(len(nsi)): 
         x = np.concatenate((np.tile(esi[k],extent) , np.tile(esi[k],extent)))
         y = np.concatenate((nsi[k]+yvec, nsi[k]-yvec))
         # Rotate line around center point
         X.append(esi[k] - ((x - esi[k]) * np.cos(theta[k])) - ((y - nsi[k]) * np.sin(theta[k])))
         Y.append(nsi[k] - ((x - esi[k]) * np.sin(theta[k])) + ((y - nsi[k]) * np.cos(theta[k])))   

      del esi, nsi, theta, x, y

      # merge flatten and stack
      X = np.asarray(X,'float').T
      X = X.flatten()

      # merge flatten and stack
      Y = np.asarray(Y,'float').T
      Y = Y.flatten()

      # write raw bs to file
      try:
         outfile = sonpath+'x_y_ss_raw.asc' 
         with open(outfile, 'w') as f:
            np.savetxt(f, np.hstack((humutils.ascol(X.flatten()), humutils.ascol(Y.flatten()), humutils.ascol(merge.flatten()))), delimiter=' ', fmt="%8.6f %8.6f %8.6f")
      except:
         outfile = os.path.expanduser("~")+os.sep+'x_y_ss_raw.asc' 
         with open(outfile, 'w') as f:
            np.savetxt(f, np.hstack((humutils.ascol(X.flatten()), humutils.ascol(Y.flatten()), humutils.ascol(merge.flatten()))), delimiter=' ', fmt="%8.6f %8.6f %8.6f")

      humlon, humlat = trans(X, Y, inverse=True)

      del X, Y, bearing, pix_m, yvec

      print "drawing and printing map ..."
      fig = plt.figure()
      map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #26949,
          resolution = 'i', #h #f
          llcrnrlon=np.min(humlon)-0.001, llcrnrlat=np.min(humlat)-0.001,
          urcrnrlon=np.max(humlon)+0.001, urcrnrlat=np.max(humlat)+0.001)

      # draw point cloud
      x,y = map.projtran(humlon, humlat)

      map.scatter(x.flatten(), y.flatten(), 0.01, merge.flatten(), cmap='gray', linewidth = '0')

      try:
         if imagery==1:
            map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)
         else:
            map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=1000, ypixels=None, dpi=300)
      except:
         print "servor error: trying again"
         if imagery==1:
            map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)
         else:
            map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=1000, ypixels=None, dpi=300)

#plt.show()

      custom_save(sonpath,base+'map')
      del fig 

      #============================
      kml = simplekml.Kml()
      ground = kml.newgroundoverlay(name='GroundOverlay')
      ground.icon.href = sonpath+base+'map.png'

      ground.latlonbox.north = np.min(humlat)-0.001
      ground.latlonbox.south = np.max(humlat)+0.001
      ground.latlonbox.east =  np.max(humlon)+0.001
      ground.latlonbox.west =  np.min(humlon)-0.001
      ground.latlonbox.rotation = 0

      try:
         kml.save(sonpath+base+"GroundOverlay.kml")
      except:
         kml.save(os.path.expanduser("~")+os.sep+base+"GroundOverlay.kml")

# =========================================================
def custom_save(figdirec,root):
   try:
      plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)
   except:
      plt.savefig(os.path.expanduser("~")+os.sep+root,bbox_inches='tight',dpi=400)      




