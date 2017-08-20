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

    [] = PyHum.map_texture(humfile, sonpath, cs2cs_args, res, mode, nn, numstdevs)

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
    res : float, *optional* [Default=0.5]
       grid resolution of output gridded texture map
    mode: int, *optional* [Default=3]
       gridding mode. 1 = nearest neighbour
                      2 = inverse weighted nearest neighbour
                      3 = Gaussian weighted nearest neighbour
    nn: int, *optional* [Default=64]
       number of nearest neighbours for gridding (used if mode > 1) 
    numstdevs: int, *optional* [Default = 4]
       Threshold number of standard deviations in texture lengthscale per grid cell up to which to accept 
       
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
-----------

     [1] Buscombe, D., Grams, P.E., and Smith, S.M.C., Automated riverbed sediment
     classification using low-cost sidescan sonar. 
     Journal of Hydraulic Engineering,  10.1061/(ASCE)HY.1943-7900.0001079, 06015019.

     [2] Buscombe, D., 2017, Shallow water benthic imaging and substrate characterization using recreational-grade sidescan-sonar. 
         ENVIRONMENTAL MODELLING & SOFTWARE 89, 1-18.


  .. image:: _static/pyhum_logo_colour_sm.png


