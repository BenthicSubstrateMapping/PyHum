.. pyhum.map:

pyhum.map module
=========================

    Create plots of the spatially referenced sidescan echograms

Syntax
----------
You call the function like this::

  [] = PyHum.map(humfile, sonpath, cs2cs_args, calc_bearing, filt_bearing, res, cog, dowrite)


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
    calc_bearing : float, *optional* [Default=0]
       if 1, bearing will be calculated from coordinates
    filt_bearing : float, *optional* [Default=0]
       if 1, bearing will be filtered
    res : float, *optional* [Default=0.1]
       grid resolution of output gridded texture map
    cog : int, *optional* [Default=1]
       if 1, heading calculated assuming GPS course-over-ground rather than
       using a compass
    dowrite: int, *optional* [Default=1]
       if 1, point cloud data from each chunk is written to ascii file
       if 0, processing times are speeded up considerably but point clouds are not available for further analysis


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


