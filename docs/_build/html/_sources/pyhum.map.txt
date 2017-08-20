.. pyhum.map:

pyhum.map module
=========================

    Create plots of the spatially referenced sidescan echograms

Syntax
----------
You call the function like this::

    [] = PyHum.map(humfile, sonpath, cs2cs_args, res, dowrite, mode, nn, numstdevs)


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
    res : float, *optional* [Default=0]
       grid resolution of output gridded texture map
       if res=0, res will be determined automatically from the spatial resolution of 1 pixel
    dowrite: int, *optional* [Default=0]
       if 1, point cloud data from each chunk is written to ascii file
       if 0, processing times are speeded up considerably but point clouds are not available for further analysis
    mode: int, *optional* [Default=3]
       gridding mode. 1 = nearest neighbour
                      2 = inverse weighted nearest neighbour
                      3 = Gaussian weighted nearest neighbour
    nn: int, *optional* [Default=64]
       number of nearest neighbours for gridding (used if mode > 1)   
    numstdevs: int, *optional* [Default = 4]
       Threshold number of standard deviations in sidescan intensity per grid cell up to which to accept            


Returns
----------
    sonpath+'x_y_ss_raw'+str(p)+'.asc'  : text file
        contains the point cloud of easting, northing, and sidescan intensity
        of the pth chunk

    sonpath+'GroundOverlay'+str(p)+'.kml': kml file
        contains gridded (or point cloud) sidescan intensity map for importing into google earth
        of the pth chunk

    sonpath+'map'+str(p)+'.png' : 
        image overlay associated with the kml file

  .. image:: _static/pyhum_logo_colour_sm.png


