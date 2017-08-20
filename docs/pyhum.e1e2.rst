.. pyhum.e1e2:

pyhum.e1e2 module
=========================

Analysis of first (e1, 'roughness') and second (e2, 'hardness') echo returns from the high-frequency downward looking echosounder

Generates generalised acoustic parameters for the purposes of point classification of submerged substrates/vegetation

Accounts for the absorption of sound in water

Does a basic k-means cluster of e1 and e2 coefficients into specified number of 'acoustic classes'

Based on code by Barb Fagetter (blueseas@oceanecology.ca)

Syntax
----------
You call the function like this::

  [] = PyHum.e1e2(humfile, sonpath, cs2cs_args, ph, temp, salinity, beam, transfreq, integ, numclusters, doplot)


Parameters
------------
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
       transducer frequency in kHz

    doplot : int, *optional* [Default=1]
       1 = make plots, otherwise do not

Returns
--------
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


  .. image:: _static/pyhum_logo_colour_sm.png


