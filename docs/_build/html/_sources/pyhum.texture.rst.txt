.. pyhum.texture:

pyhum.texture module
======================

     Create a texture lengthscale map using the algorithm detailed by Buscombe et al. (2015) Journal of Hydraulic Engineering

     This textural lengthscale is not a direct measure of grain size. Rather, it is a statistical 

     representation that integrates over many attributes of bed texture, of which grain size is the most important. 

     The technique is a physically based means to identify regions of texture within a sidescan echogram, 

     and could provide a basis for objective, automated riverbed sediment classification.

Syntax
----------

You call the function like this::

  [] = PyHum.texture2(humfile, sonpath, win, doplot, numclasses)

Parameters
------------

     humfile : str
       path to the .DAT file
     sonpath : str
       path where the *.SON files are
     win : int, *optional* [Default=100]
       pixel in pixels of the moving window
     shift : int, *optional* [Default=10]
       shift in pixels for moving window operation
     doplot : int, *optional* [Default=1]
       if 1, make plots, otherwise do not make plots
     numclasses : int, *optional* [Default=4]
       number of 'k means' that the texture lengthscale will be segmented into

Returns
----------

     sonpath+base+'_data_class.dat': memory-mapped file
        contains the texture lengthscale map

     sonpath+base+'_data_kclass.dat': memory-mapped file
        contains the k-means segmented texture lengthscale map

References
-----------

     [1] Buscombe, D., Grams, P.E., and Smith, S.M.C., Automated riverbed sediment
     classification using low-cost sidescan sonar. 
     Journal of Hydraulic Engineering,  10.1061/(ASCE)HY.1943-7900.0001079, 06015019.

  .. image:: _static/pyhum_logo_colour_sm.png


