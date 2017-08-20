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
#       ___       ___ 
#  ___ <  /  ___ |__ \
# / _ \/ /  / _ \__/ /
#/  __/ /  /  __/ __/ 
#\___/_/   \___/____/ 
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

# operational
from __future__ import division
from scipy.io import loadmat
import os #, time, sys, getopt
try:
   from Tkinter import Tk
   from tkFileDialog import askopenfilename, askdirectory
except:
   pass
from joblib import Parallel, delayed #, cpu_count
#import pyproj

# numerical
import numpy as np
#import pyproj
import csv
#from math import sqrt as sqrt
from numpy import power as pow
#from math import sin as sin
from math import pi as pi
from math import log10 as log10
from sklearn.cluster import MiniBatchKMeans

import PyHum.utils as humutils #runningMeanFast, nan_helper

# plotting
import matplotlib.pyplot as plt
try:
   from mpl_toolkits.basemap import Basemap
except:
   print "Error: Basemap could not be imported"
   pass
import simplekml

# suppress divide and invalid warnings
np.seterr(divide='ignore')
np.seterr(invalid='ignore')


#################################################
def e1e2(humfile, sonpath, cs2cs_args="epsg:26949", ph=7.0, temp=10.0, salinity=0.0, beam=20.0, transfreq=200.0, integ=5, numclusters=3, doplot=1):
         
    '''
    Analysis of first (e1, 'roughness') and second (e2, 'hardness') echo returns from the high-frequency downward looking echosounder
    Generates generalised acoustic parameters
    for the purposes of point classification of submerged substrates/vegetation
    Accounts for the absorption of sound in water
    Does a basic k-means cluster of e1 and e2 coefficients into specified number of 'acoustic classes'
    based on code by Barb Fagetter (blueseas@oceanecology.ca)

    Syntax
    ----------
    [] = PyHum.e1e2(humfile, sonpath, cs2cs_args, ph, temp, salinity, beam, transfreq, integ, numclusters, doplot)

    Parameters
    ----------
    humfile : str
       path to the .DAT file

    sonpath : str
       path where the *.SON files are

    cs2cs_args : int, *optional* [Default="epsg:26949"]
       arguments to create coordinates in a projected coordinate system
       this argument gets given to pyproj to turn wgs84 (lat/lon) coordinates
       into any projection supported by the proj.4 libraries

    ph : float, *optional* [Default=7.0]
       water acidity in pH

    temp : float, *optional* [Default=10.0]
       water temperature in degrees Celsius

    salinity : float, *optional* [Default=0.0]
       salinity of water in parts per thousand

    beam : float, *optional* [Default=20.0]
       beam width in degrees

    transfreq : float, *optional* [Default=200.0]
       transducer frequency in kHz

    integ : int, *optional* [Default=5]
       number of pings over which to integrate

    numclusters : int, *optional* [Default=3]
       number of acoustic classes to classify all the data into

    doplot : int, *optional* [Default=1]
       1 = make plots, otherwise do not

    Returns
    -------
    sonpath+base+'rough_and_hard'+str(p)+'.csv'  : csv file
        contains the following fields: 'longitude', 'latitude', 'easting', 'northing', 'depth', 
        'roughness', 'hardness', 'average roughness', 'average hardness','k-mean label'
        of the pth chunk
        'average' implies average over 'integ' successive pings

    The following are returned if doplot==1:

    sonpath+'e1e2_scan'+str(p).png : png image file
       png image file showing the downward echosounder echogram overlain with the locations of the start and 
       end of the first and second echo region envelope 

    sonpath+'e1e2_kmeans'+str(p).png: png image file
        png image file showing 1) (left) volume scattering coefficient 1 versus volume scattering coefficient 2, colour-coded
        by k-means acoustic class, and
        2) (right) e1 versus e2, colour-coded
        by k-means acoustic class

    sonpath+'rgh_hard_kmeans'+str(p).png : png image file
        png image file showing scatter plot of easting versus northing colour-coded by k-means acoustic class 

    sonpath+'map_rgh'+str(p).png : png image file
        png image file showing scatter plot of 'roughness' (e1) overlying an aerial image pulled from an ESRI image server 

    sonpath+'map_hard'+str(p).png : png image file
        png image file showing scatter plot of 'hardness' (e2) overlying an aerial image pulled from an ESRI image server 

    sonpath,'Rough'+str(p).png : png image file 
        png image overlay associated with the kml file, sonpath,'Hard'+str(p).kml

    sonpath,'Rough'+str(p).kml : kml file
        kml overlay for showing roughness scatter plot (sonpath,'Rough'+str(p).png)

    sonpath,'Hard'+str(p).png : png image file
        png image overlay associated with the kml file, sonpath,'Hard'+str(p).kml
    
    sonpath,'Hard'+str(p).kml : kml file
        kml overlay for showing harness scatter plot (sonpath,'Hard'+str(p).png)

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

    if cs2cs_args:
       print 'cs2cs arguments are %s' % (cs2cs_args)

    if beam:
       beam = np.asarray(beam,float)
       print 'Beam is %s deg' % (str(beam))

    if salinity:
       salinity = np.asarray(salinity,float)
       print 'Salinity is %s ppt' % (str(salinity))

    if ph:
       ph = np.asarray(ph,float)
       print 'pH is %s' % (str(ph))

    if temp:
       temp = np.asarray(temp,float)
       print 'Temperature is %s' % (str(temp))

    if transfreq:
       transfreq = np.asarray(transfreq,float)
       print 'Dwnward sonar freq. is %s' % (str(transfreq))

    if integ:
       integ = np.asarray(integ,int)
       print 'number of records for integration is %s' % (str(integ))

    if numclusters:
       numclusters = np.asarray(numclusters,int)
       print 'number of returned acoustic clusters is %s' % (str(numclusters))

    if doplot:
      doplot = int(doplot)
      if doplot==0:
         print "Plots will not be made"


    # if son path name supplied has no separator at end, put one on
    if sonpath[-1]!=os.sep:
       sonpath = sonpath + os.sep

    base = humfile.split('.DAT') # get base of file name for output
    base = base[0].split(os.sep)[-1]

    # remove underscores, negatives and spaces from basename
    base = humutils.strip_base(base)

    meta = loadmat(os.path.normpath(os.path.join(sonpath,base+'meta.mat')))

    beamwidth = beam*(np.sqrt(0.5))
    equivbeam = (5.78/(np.power(1.6,2)))*(np.power((np.sin((beamwidth*np.pi)/(2*180))),2))

    meta = loadmat(sonpath+base+'meta.mat')

    c = np.squeeze(meta['c'])
    t = np.squeeze(meta['t'])
    f = np.squeeze(meta['f'])
    maxW = np.squeeze(meta['maxW'])

    lat = np.squeeze(meta['lat'])
    lon = np.squeeze(meta['lon'])
    es = np.squeeze(meta['e'])
    ns = np.squeeze(meta['n'])
    dep = np.squeeze(meta['dep_m'])
    #del meta

    # load memory mapped scans
    shape_hi= np.squeeze(meta['shape_hi'])
    if shape_hi!='':
       try:
          #dwnhi_fp = np.memmap(sonpath+base+'_data_dwnhi.dat', dtype='int16', mode='r', shape=tuple(shape_hi))
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
             dwnhi_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_hi))

       except:
          shape_lo= np.squeeze(meta['shape_low'])
          #dwnhi_fp = np.memmap(sonpath+base+'_data_dwnhi.dat', dtype='int16', mode='r', shape=tuple(shape_lo))
          with open(os.path.normpath(os.path.join(sonpath,base+'_data_dwnhi.dat')), 'r') as ff:
             dwnhi_fp = np.memmap(ff, dtype='int16', mode='r', shape=tuple(shape_lo))
    

    if 'dwnhi_fp' in locals():

       theta3dB = np.arcsin(c/(t*(f*1000))) # *(180/pi) # to see in degs
       ft = (np.pi/2)*(1/theta3dB)
       bed = ft*dep

       if len(shape_hi)>2:    
          i = np.linspace(1,shape_hi[0]*shape_hi[2], len(bed)) 
          #np.shape(beam_data)[1],len(bed))
          #bedi = np.interp(np.linspace(1,np.shape(beam_data)[1],np.shape(beam_data)[1]), i, bed)
          bedi = np.interp(np.linspace(1,shape_hi[0]*shape_hi[2],shape_hi[0]*shape_hi[2]), i, bed)  
          ei = np.interp(np.linspace(1,shape_hi[0]*shape_hi[2],shape_hi[0]*shape_hi[2]), i, es)    
          ni = np.interp(np.linspace(1,shape_hi[0]*shape_hi[2],shape_hi[0]*shape_hi[2]), i, ns)    
          lati = np.interp(np.linspace(1,shape_hi[0]*shape_hi[2],shape_hi[0]*shape_hi[2]), i, lat)   
          loni = np.interp(np.linspace(1,shape_hi[0]*shape_hi[2],shape_hi[0]*shape_hi[2]), i, lon)    
          del i
       else:
          i = np.linspace(1,shape_hi[1], len(bed)) 
          #np.shape(beam_data)[1],len(bed))
          #bedi = np.interp(np.linspace(1,np.shape(beam_data)[1],np.shape(beam_data)[1]), i, bed)
          bedi = np.interp(np.linspace(1,shape_hi[1],shape_hi[1]), i, bed)  
          ei = np.interp(np.linspace(1,shape_hi[1],shape_hi[1]), i, es)    
          ni = np.interp(np.linspace(1,shape_hi[1],shape_hi[1]), i, ns)    
          lati = np.interp(np.linspace(1,shape_hi[1],shape_hi[1]), i, lat)   
          loni = np.interp(np.linspace(1,shape_hi[1],shape_hi[1]), i, lon)    
          del i 
          
       bedi = np.asarray(bedi,'int')

       depi = ((1/ft)*bedi) 

       # near-field region
       nf = int(ft*(1000*(0.105**2)*f/(4*1500)))

       #absorption = calcAb(c, ph, salinity, temp, np.asarray(depi), transfreq)
       absorption = water_atten(np.asarray(depi), transfreq, c, ph, temp, salinity)

       if len(shape_hi)>2:    
          for p in xrange(len(dwnhi_fp)):
             #make an index of every other record
             ind = range(0,np.shape(dwnhi_fp[p])[1])

             Zdepi = depi[shape_hi[2]*p:shape_hi[2]*(p+1)]
             Zabsorp = absorption[shape_hi[2]*p:shape_hi[2]*(p+1)]
             Zlat = lati[shape_hi[2]*p:shape_hi[2]*(p+1)]
             Zlon = loni[shape_hi[2]*p:shape_hi[2]*(p+1)]       
             Zes = ei[shape_hi[2]*p:shape_hi[2]*(p+1)]
             Zns = ni[shape_hi[2]*p:shape_hi[2]*(p+1)]       
                
             try: #parallel processing with all available cores
               w = Parallel(n_jobs=-1, verbose=0)(delayed(get_rgh_hrd)(dwnhi_fp[p][:,i],Zdepi[i],Zabsorp[i],c,nf,transfreq,equivbeam,maxW,pi,ft) for i in ind)
             except: #fall back to serial
               w = Parallel(n_jobs=1, verbose=0)(delayed(get_rgh_hrd)(dwnhi_fp[p][:,i],Zdepi[i],Zabsorp[i],c,nf,transfreq,equivbeam,maxW,pi,ft) for i in ind)

             rough, hard, sv_e1, sv_e2, e1a, e1b, e2a, e2b = zip(*w) 

             rough = np.array(rough,'float')
             rough[rough==0.0] = np.nan

             hard = np.array(hard,'float')
             hard[hard==0.0] = np.nan

             sv_e1 = np.array(sv_e1,'float')
             sv_e1[sv_e1==0.0] = np.nan

             sv_e2 = np.array(sv_e2,'float')
             sv_e2[sv_e2==0.0] = np.nan

             try:
                nans, y= humutils.nan_helper(rough)
                rough[nans]= np.interp(y(nans), y(~nans), rough[~nans])
             except:
                pass

             try:
                nans, y= humutils.nan_helper(hard)
                hard[nans]= np.interp(y(nans), y(~nans), hard[~nans])
             except:
                pass

             try:
                nans, y= humutils.nan_helper(sv_e1)
                sv_e1[nans]= np.interp(y(nans), y(~nans), sv_e1[~nans])
             except:
                pass

             try:
                nans, y= humutils.nan_helper(sv_e2)
                sv_e2[nans]= np.interp(y(nans), y(~nans), sv_e2[~nans])
             except:
                pass

             data = np.column_stack([sv_e1, sv_e2])
             k_means = MiniBatchKMeans(numclusters)
             # fit the model
             k_means.fit(data) 
             values = k_means.cluster_centers_.squeeze()
             labels = k_means.labels_
    
             hardav = humutils.runningMeanFast(hard,integ)
             roughav = humutils.runningMeanFast(rough,integ)

             #f = open(sonpath+base+'rough_and_hard'+str(p)+'.csv', 'wt')
             f = open(os.path.normpath(os.path.join(sonpath,base+'rough_and_hard'+str(p)+'.csv')), 'wt')
             writer = csv.writer(f)
             writer.writerow( ('longitude', 'latitude', 'easting', 'northing', 'depth', 'roughness', 'hardness', 'average roughness', 'average hardness','k-mean label') )
             for i in range(0, len(rough)):
                writer.writerow(( float(Zlon[i]),float(Zlat[i]),float(Zes[i]),float(Zns[i]),float(Zdepi[i]),float(rough[i]),float(hard[i]),float(roughav[i]),float(hardav[i]), labels[i].astype(int) ))
             f.close()

             if doplot==1:
                try:

                   fig = plt.figure()
                   plt.imshow(dwnhi_fp[p], cmap='gray')
                   plt.plot(e1a,'r');
                   plt.plot(e1b,'y');
                   plt.plot(e2a,'c');
                   plt.plot(e2b,'m');
                   plt.axis('tight')
                   #plt.show()
                   custom_save(sonpath,'e1e2_scan'+str(p))
                   del fig

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:

                   fig = plt.figure()
                   fig.subplots_adjust(wspace = 0.4, hspace=0.4)
                   plt.subplot(221)   
                   plt.plot(sv_e1[labels==0],sv_e2[labels==0],'ko');
                   plt.plot(sv_e1[labels==1],sv_e2[labels==1],'ro');
                   plt.plot(sv_e1[labels==2],sv_e2[labels==2],'bo');
                   plt.xlabel('SV1'); plt.ylabel('SV2')
                   plt.xlim(0,1); plt.ylim(0,1)

                   plt.subplot(222)   
                   plt.plot(rough[labels==0],hard[labels==0],'ko');
                   plt.plot(rough[labels==1],hard[labels==1],'ro');
                   plt.plot(rough[labels==2],hard[labels==2],'bo');
                   plt.xlabel('E1'); plt.ylabel('E2')
                   plt.xlim(1,8); plt.ylim(1,8)
                   #plt.show()
                   custom_save(sonpath,'e1e2_kmeans'+str(p))
                   del fig

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:

                   fig = plt.figure()
                   s=plt.scatter(Zes[labels==0],Zns[labels==0],marker='o',c='k', s=10, linewidth=0, vmin=0, vmax=8);
                   s=plt.scatter(Zes[labels==1],Zns[labels==1],marker='o',c='r', s=10, linewidth=0, vmin=0, vmax=8);
                   s=plt.scatter(Zes[labels==2],Zns[labels==2],marker='o',c='b', s=10, linewidth=0, vmin=0, vmax=8);
                   custom_save(sonpath,'rgh_hard_kmeans'+str(p))
                   del fig

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:

                   print "drawing and printing map ..."
                   fig = plt.figure(frameon=False)
                   #fig.subplots_adjust(wspace = 0.4, hspace=0.4)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #epsg=26949,
                      resolution = 'i', #h #f
                      llcrnrlon=np.min(Zlon)-0.0001, llcrnrlat=np.min(Zlat)-0.0001,
                      urcrnrlon=np.max(Zlon)+0.0001, urcrnrlat=np.max(Zlat)+0.0001)

                   # draw point cloud
                   x,y = map.projtran(Zlon, Zlat)

                   cs = map.scatter(x.flatten(), y.flatten(), 1, rough.flatten(), linewidth=0, vmin=0, vmax=8)

                   try:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=1000, ypixels=None, dpi=300)
                   except:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)                   

                   cbar = map.colorbar(cs,location='bottom',pad="5%")
                   cbar.set_label('E1')
                   cbar.set_ticks([0,2,4,6,8])

                   custom_save(sonpath,'map_rgh'+str(p))
                   del fig 

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:
                   fig = plt.figure()
                   #fig.subplots_adjust(wspace = 0.4, hspace=0.4)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1],
                      resolution = 'i', #h #f
                      llcrnrlon=np.min(Zlon)-0.0001, llcrnrlat=np.min(Zlat)-0.0001,
                      urcrnrlon=np.max(Zlon)+0.0001, urcrnrlat=np.max(Zlat)+0.0001)

                   # draw point cloud
                   x,y = map.projtran(Zlon, Zlat)

                   cs = map.scatter(x.flatten(), y.flatten(), 1, hard.flatten(), linewidth=0, vmin=0, vmax=8)

                   try:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=1000, ypixels=None, dpi=300)
                   except:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)
                      
                   cbar = map.colorbar(cs,location='bottom',pad="5%")
                   cbar.set_label('E2')
                   cbar.set_ticks([0,2,4,6,8])

                   custom_save(sonpath,'map_hard'+str(p))
                   del fig 

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:
                   print "drawing and printing map ..."
                   fig = plt.figure(frameon=False)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #26949,
                    resolution = 'i', #h #f
                    llcrnrlon=np.min(Zlon)-0.001, llcrnrlat=np.min(Zlat)-0.001,
                    urcrnrlon=np.max(Zlon)+0.001, urcrnrlat=np.max(Zlat)+0.001)

                   ax = plt.Axes(fig, [0., 0., 1., 1.], )
                   ax.set_axis_off()
                   fig.add_axes(ax)

                   ## draw point cloud
                   x,y = map.projtran(Zlon, Zlat)
                   map.scatter(x.flatten(), y.flatten(), 1, rough.flatten(), linewidth = '0', vmin=0, vmax=8)

                   custom_save(sonpath,'Rough'+str(p))
                   del fig 

                   kml = simplekml.Kml()
                   ground = kml.newgroundoverlay(name='GroundOverlay')
                   ground.icon.href = 'Rough'+str(p)+'.png'
                   ground.latlonbox.north = np.min(Zlat)-0.001
                   ground.latlonbox.south = np.max(Zlat)+0.001
                   ground.latlonbox.east =  np.max(Zlon)+0.001
                   ground.latlonbox.west =  np.min(Zlon)-0.001
                   ground.latlonbox.rotation = 0

                   #kml.save(sonpath+'Rough'+str(p)+'.kml')
                   kml.save(os.path.normpath(os.path.join(sonpath,'Rough'+str(p)+'.kml')))

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:
                   print "drawing and printing map ..."
                   fig = plt.figure(frameon=False)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #26949,
                    resolution = 'i', #h #f
                    llcrnrlon=np.min(Zlon)-0.001, llcrnrlat=np.min(Zlat)-0.001,
                    urcrnrlon=np.max(Zlon)+0.001, urcrnrlat=np.max(Zlat)+0.001)

                   ax = plt.Axes(fig, [0., 0., 1., 1.], )
                   ax.set_axis_off()
                   fig.add_axes(ax)

                   ## draw point cloud
                   x,y = map.projtran(Zlon, Zlat)
                   map.scatter(x.flatten(), y.flatten(), 1, hard.flatten(), linewidth = '0', vmin=0, vmax=8)

                   custom_save(sonpath,'Hard'+str(p))
                   del fig 

                   kml = simplekml.Kml()
                   ground = kml.newgroundoverlay(name='GroundOverlay')
                   ground.icon.href = 'Hard'+str(p)+'.png'
                   ground.latlonbox.north = np.min(Zlat)-0.001
                   ground.latlonbox.south = np.max(Zlat)+0.001
                   ground.latlonbox.east =  np.max(Zlon)+0.001
                   ground.latlonbox.west =  np.min(Zlon)-0.001
                   ground.latlonbox.rotation = 0

                   #kml.save(sonpath+'Hard'+str(p)+'.kml')
                   kml.save(os.path.normpath(os.path.join(sonpath,'Hard'+str(p)+'.kml')))

                except:
                   print "plot could not be produced"

       else:
          if 2 > 1: # need to tiday all this up later!!
             #make an index of every other record
             ind = range(0,np.shape(dwnhi_fp)[1])

             Zdepi = depi
             Zabsorp = absorption
             Zlat = lati
             Zlon = loni       
             Zes = ei
             Zns = ni      
                
             try: #parallel processing with all available cores
               w = Parallel(n_jobs=-1, verbose=0)(delayed(get_rgh_hrd)(dwnhi_fp[:,i],Zdepi[i],Zabsorp[i],c,nf,transfreq,equivbeam,maxW,pi,ft) for i in ind)
             except: #fall back to serial
               w = Parallel(n_jobs=1, verbose=0)(delayed(get_rgh_hrd)(dwnhi_fp[:,i],Zdepi[i],Zabsorp[i],c,nf,transfreq,equivbeam,maxW,pi,ft) for i in ind)

             rough, hard, sv_e1, sv_e2, e1a, e1b, e2a, e2b = zip(*w) 

             rough = np.array(rough,'float')
             rough[rough==0.0] = np.nan

             hard = np.array(hard,'float')
             hard[hard==0.0] = np.nan

             sv_e1 = np.array(sv_e1,'float')
             sv_e1[sv_e1==0.0] = np.nan

             sv_e2 = np.array(sv_e2,'float')
             sv_e2[sv_e2==0.0] = np.nan

             try:
                nans, y= humutils.nan_helper(rough)
                rough[nans]= np.interp(y(nans), y(~nans), rough[~nans])
             except:
                pass

             try:
                nans, y= humutils.nan_helper(hard)
                hard[nans]= np.interp(y(nans), y(~nans), hard[~nans])
             except:
                pass

             try:
                nans, y= humutils.nan_helper(sv_e1)
                sv_e1[nans]= np.interp(y(nans), y(~nans), sv_e1[~nans])
             except:
                pass

             try:
                nans, y= humutils.nan_helper(sv_e2)
                sv_e2[nans]= np.interp(y(nans), y(~nans), sv_e2[~nans])
             except:
                pass

             data = np.column_stack([sv_e1, sv_e2])
             k_means = MiniBatchKMeans(numclusters)
             # fit the model
             k_means.fit(data) 
             values = k_means.cluster_centers_.squeeze()
             labels = k_means.labels_
    
             hardav = humutils.runningMeanFast(hard,integ)
             roughav = humutils.runningMeanFast(rough,integ)

             #f = open(sonpath+base+'rough_and_hard'+str(p)+'.csv', 'wt')
             f = open(os.path.normpath(os.path.join(sonpath,base+'rough_and_hard'+str(0)+'.csv')), 'wt')
             writer = csv.writer(f)
             writer.writerow( ('longitude', 'latitude', 'easting', 'northing', 'depth', 'roughness', 'hardness', 'average roughness', 'average hardness','k-mean label') )
             for i in range(0, len(rough)):
                writer.writerow(( float(Zlon[i]),float(Zlat[i]),float(Zes[i]),float(Zns[i]),float(Zdepi[i]),float(rough[i]),float(hard[i]),float(roughav[i]),float(hardav[i]), labels[i].astype(int) ))
             f.close()

             if doplot==1:
                try:

                   fig = plt.figure()
                   plt.imshow(dwnhi_fp, cmap='gray')
                   plt.plot(e1a,'r');
                   plt.plot(e1b,'y');
                   plt.plot(e2a,'c');
                   plt.plot(e2b,'m');
                   plt.axis('tight')
                   #plt.show()
                   custom_save(sonpath,'e1e2_scan'+str(0))
                   del fig

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:

                   fig = plt.figure()
                   fig.subplots_adjust(wspace = 0.4, hspace=0.4)
                   plt.subplot(221)   
                   plt.plot(sv_e1[labels==0],sv_e2[labels==0],'ko');
                   plt.plot(sv_e1[labels==1],sv_e2[labels==1],'ro');
                   plt.plot(sv_e1[labels==2],sv_e2[labels==2],'bo');
                   plt.xlabel('SV1'); plt.ylabel('SV2')
                   plt.xlim(0,1); plt.ylim(0,1)

                   plt.subplot(222)   
                   plt.plot(rough[labels==0],hard[labels==0],'ko');
                   plt.plot(rough[labels==1],hard[labels==1],'ro');
                   plt.plot(rough[labels==2],hard[labels==2],'bo');
                   plt.xlabel('E1'); plt.ylabel('E2')
                   plt.xlim(1,8); plt.ylim(1,8)
                   #plt.show()
                   custom_save(sonpath,'e1e2_kmeans'+str(0))
                   del fig

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:

                   fig = plt.figure()
                   s=plt.scatter(Zes[labels==0],Zns[labels==0],marker='o',c='k', s=10, linewidth=0, vmin=0, vmax=8);
                   s=plt.scatter(Zes[labels==1],Zns[labels==1],marker='o',c='r', s=10, linewidth=0, vmin=0, vmax=8);
                   s=plt.scatter(Zes[labels==2],Zns[labels==2],marker='o',c='b', s=10, linewidth=0, vmin=0, vmax=8);
                   custom_save(sonpath,'rgh_hard_kmeans'+str(0))
                   del fig

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:

                   print "drawing and printing map ..."
                   fig = plt.figure(frameon=False)
                   #fig.subplots_adjust(wspace = 0.4, hspace=0.4)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #epsg=26949,
                      resolution = 'i', #h #f
                      llcrnrlon=np.min(Zlon)-0.0001, llcrnrlat=np.min(Zlat)-0.0001,
                      urcrnrlon=np.max(Zlon)+0.0001, urcrnrlat=np.max(Zlat)+0.0001)

                   # draw point cloud
                   x,y = map.projtran(Zlon, Zlat)

                   cs = map.scatter(x.flatten(), y.flatten(), 1, rough.flatten(), linewidth=0, vmin=0, vmax=8)

                   try:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=1000, ypixels=None, dpi=300)
                   except:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)
  
                   cbar = map.colorbar(cs,location='bottom',pad="5%")
                   cbar.set_label('E1')
                   cbar.set_ticks([0,2,4,6,8])

                   custom_save(sonpath,'map_rgh'+str(0))
                   del fig 

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:
                   fig = plt.figure()
                   #fig.subplots_adjust(wspace = 0.4, hspace=0.4)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1],
                      resolution = 'i', #h #f
                      llcrnrlon=np.min(Zlon)-0.0001, llcrnrlat=np.min(Zlat)-0.0001,
                      urcrnrlon=np.max(Zlon)+0.0001, urcrnrlat=np.max(Zlat)+0.0001)

                   # draw point cloud
                   x,y = map.projtran(Zlon, Zlat)

                   cs = map.scatter(x.flatten(), y.flatten(), 1, hard.flatten(), linewidth=0, vmin=0, vmax=8)

                   try:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='ESRI_Imagery_World_2D', xpixels=1000, ypixels=None, dpi=300)
                   except:
                      map.arcgisimage(server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery', xpixels=1000, ypixels=None, dpi=300)

                   cbar = map.colorbar(cs,location='bottom',pad="5%")
                   cbar.set_label('E2')
                   cbar.set_ticks([0,2,4,6,8])

                   custom_save(sonpath,'map_hard'+str(0))
                   del fig 

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:
                   print "drawing and printing map ..."
                   fig = plt.figure(frameon=False)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #26949,
                    resolution = 'i', #h #f
                    llcrnrlon=np.min(Zlon)-0.001, llcrnrlat=np.min(Zlat)-0.001,
                    urcrnrlon=np.max(Zlon)+0.001, urcrnrlat=np.max(Zlat)+0.001)

                   ax = plt.Axes(fig, [0., 0., 1., 1.], )
                   ax.set_axis_off()
                   fig.add_axes(ax)

                   ## draw point cloud
                   x,y = map.projtran(Zlon, Zlat)
                   map.scatter(x.flatten(), y.flatten(), 1, rough.flatten(), linewidth = '0', vmin=0, vmax=8)

                   custom_save(sonpath,'Rough'+str(0))
                   del fig 

                   kml = simplekml.Kml()
                   ground = kml.newgroundoverlay(name='GroundOverlay')
                   ground.icon.href = 'Rough'+str(0)+'.png'
                   ground.latlonbox.north = np.min(Zlat)-0.001
                   ground.latlonbox.south = np.max(Zlat)+0.001
                   ground.latlonbox.east =  np.max(Zlon)+0.001
                   ground.latlonbox.west =  np.min(Zlon)-0.001
                   ground.latlonbox.rotation = 0

                   #kml.save(sonpath+'Rough'+str(p)+'.kml')
                   kml.save(os.path.normpath(os.path.join(sonpath,'Rough'+str(0)+'.kml')))

                except:
                   print "plot could not be produced"

             if doplot==1:
                try:
                   print "drawing and printing map ..."
                   fig = plt.figure(frameon=False)
                   map = Basemap(projection='merc', epsg=cs2cs_args.split(':')[1], #26949,
                    resolution = 'i', #h #f
                    llcrnrlon=np.min(Zlon)-0.001, llcrnrlat=np.min(Zlat)-0.001,
                    urcrnrlon=np.max(Zlon)+0.001, urcrnrlat=np.max(Zlat)+0.001)

                   ax = plt.Axes(fig, [0., 0., 1., 1.], )
                   ax.set_axis_off()
                   fig.add_axes(ax)

                   ## draw point cloud
                   x,y = map.projtran(Zlon, Zlat)
                   map.scatter(x.flatten(), y.flatten(), 1, hard.flatten(), linewidth = '0', vmin=0, vmax=8)

                   custom_save(sonpath,'Hard'+str(0))
                   del fig 

                   kml = simplekml.Kml()
                   ground = kml.newgroundoverlay(name='GroundOverlay')
                   ground.icon.href = 'Hard'+str(0)+'.png'
                   ground.latlonbox.north = np.min(Zlat)-0.001
                   ground.latlonbox.south = np.max(Zlat)+0.001
                   ground.latlonbox.east =  np.max(Zlon)+0.001
                   ground.latlonbox.west =  np.min(Zlon)-0.001
                   ground.latlonbox.rotation = 0

                   #kml.save(sonpath+'Hard'+str(p)+'.kml')
                   kml.save(os.path.normpath(os.path.join(sonpath,'Hard'+str(0)+'.kml')))

                except:
                   print "plot could not be produced"         

    else:
       print "high-frequency downward echosounder data not available"


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
   return H*(alphaw/1000) # depth(m) * dB/m = dB

# =========================================================
def custom_save(figdirec,root):
    #plt.savefig(figdirec+root,bbox_inches='tight',dpi=400, transparent=True)
    plt.savefig(os.path.normpath(os.path.join(figdirec,root)),bbox_inches='tight',dpi=400, transparent=True)

# =========================================================
def get_rgh_hrd(beamdat,dep,absorp,c,nf,transfreq,equivbeam,maxW,pi,ft):
    peakstart = (int)((float(dep)/c)*76923)
    noiseend = (int)(round(0.9*peakstart))
    E1start = peakstart+int(ft/2) #27
    E1end = peakstart+(nf*3) #131
    E2start = (int)(2*peakstart)
    E2end = (int)(2*peakstart)+(nf*3) #131
    sum = 0

    try:
       for k in range(nf,noiseend): #80
          backstrength = ((30.0/255.0)*(float((np.squeeze(beamdat))[k])))+(20*log10(float(dep)))+(2*(absorp/1000)*(float(dep)))-(10*(log10(((maxW*(pow((c/(transfreq/1000)),2))*c*0.0007*equivbeam)/(32*(pow(pi,2)))))))
          backcoeff = pow(10,(backstrength/10))
          sum = sum + backcoeff
       n = noiseend - nf + 1 #80 + 1
       noise = (4*pi*(pow(1852.0,2))*(2*sum))/max(n,1)
       sum = 0
       for k in range(E1start,E1end):
          backstrength = ((30.0/255.0)*(float((np.squeeze(beamdat))[k])))+(20*log10(float(dep)))+(2*(absorp/1000)*(float(dep)))-(10*(log10(((maxW*(pow((c/(transfreq/1000)),2))*c*0.0007*equivbeam)/(32*(pow(pi,2)))))))
          backcoeff = pow(10,(backstrength/10))
          sum = sum + backcoeff
       sv_e1 = sum
       n = E1end - E1start + 1
       energy = (4*pi*(pow(1852.0,2))*(2*sum))-(max(n,1)*noise)
       if energy < 0:
          energy = 1.0
       rough = log10(energy)
    except:
       rough = np.nan
       sv_e1 = np.nan

    try:    
       sum = 0
       for k in range(E2start,E2end):
          backstrength = ((30.0/255.0)*(float((np.squeeze(beamdat))[k])))+(20*log10(float(dep)))+(2*(absorp/1000)*(float(dep)))-(10*(log10(((maxW*(pow((c/(transfreq/1000)),2))*c*0.0007*equivbeam)/(32*(pow(pi,2)))))))
          backcoeff = pow(10,(backstrength/10))
          sum = sum + backcoeff
       sv_e2 = sum
       n = E2end - E2start + 1
       energy = (4*pi*(pow(1852.0,2))*(2*sum))-(max(n,1)*noise)
       if energy < 0:
          energy = 1.0
       hard = log10(energy)
       sum = 0
    except:
       hard = np.nan
       sv_e2 = np.nan

    return rough, hard, sv_e1, sv_e2, E1start, E1end, E2start, E2end


# =========================================================
# =========================================================
if __name__ == '__main__':
   
   e1e2(humfile, sonpath, cs2cs_args, ph, temp, salinity, beam, transfreq, integ, numclusters, doplot)


#    if not beam:
#       beam = 20.0
#       print '[Default] Beam is %s deg' % (str(beam))

#    if not salinity:
#       if salinity != 0.0:
#          salinity = 0.0
#          print '[Default] Salinity is %s ppt' % (str(salinity))

#    if not ph:
#       ph = 7.0
#       print '[Default] pH is %s' % (str(ph))

#    if not integ:
#       integ = 5
#       print '[Default] Number of records for integration is %s' % (str(ph))

#    if not numclusters:
#       numclusters = 3
#       print '[Default] Number of acoustic clusters is %s' % (str(ph))

#    if not temp:
#       temp = 10.0
#       print '[Default] Temperature is %s degC' % (str(temp))

#    if not transfreq:
#       transfreq = 200.0
#       print '[Default] Dwnward freq. is %s kHz' % (str(transfreq))

#    if not cs2cs_args:
#       # arguments to pass to cs2cs for coordinate transforms
#       cs2cs_args = "epsg:26949"
#       print '[Default] cs2cs arguments are %s' % (cs2cs_args)

#    if not doplot:
#      if doplot != 0:
#         doplot = 1
#         print "[Default] Plots will be made"

