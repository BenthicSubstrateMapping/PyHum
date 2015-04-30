.. pyhum.map_texture:

pyhum.map_texture module
=========================

    Create plots of the texture lengthscale maps made in PyHum.texture module 
    using the algorithm detailed by Buscombe et al. (forthcoming)
    This textural lengthscale is not a direct measure of grain size. Rather, it is a statistical 
    representation that integrates over many attributes of bed texture, of which grain size is the most important. 
    The technique is a physically based means to identify regions of texture within a sidescan echogram, 
    and could provide a basis for objective, automated riverbed sediment classification.

Syntax
----------
You call the function like this::

  [] = PyHum.map_texture(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res)

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
    dogrid : float, *optional* [Default=1]
       if 1, textures will be gridded with resolution 'res'. 
       Otherwise, point cloud will be plotted
    calc_bearing : float, *optional* [Default=1]
       if 1, bearing will be calculated from coordinates
    filt_bearing : float, *optional* [Default=1]
       if 1, bearing will be filtered
    res : float, *optional* [Default=1]
       grid resolution of output gridded texture map

Returns
---------
    sonpath+'x_y_class'+str(p)+'.asc' : text file
        contains the point cloud of easting, northing, and texture lengthscales
        of the pth chunk

    sonpath+'class_GroundOverlay'+str(p)+'.kml': kml file
        contains gridded (or point cloud) texture lengthscale map for importing into google earth
        of the pth chunk

    sonpath+'class_map'+str(p)+'.png' : 
        image overlay associated with the kml file

    sonpath+'class_map_imagery'+str(p)+'.png' : png image file
        gridded (or point cloud) texture lengthscale map
        overlain onto an image pulled from esri image server

References
-------------
     .. [1] Buscombe, D., Grams, P.E., and Smith, S.M.C., Automated riverbed sediment
       classification using low-cost sidescan sonar. submitted to
       Journal of Hydraulic Engineering

  .. image:: _static/pyhum_logo_colour_sm.png


