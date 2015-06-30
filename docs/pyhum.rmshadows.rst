.. pyhum.rmshadows:

pyhum.rmshadows module
======================

    Remove dark shadows in scans caused by shallows, shorelines, and attenuation of acoustics with distance
    Manual or automated processing options available
    Works on the radiometrically corrected outputs of the correct module

Syntax
----------

You call the function like this::

    [] = PyHum.rmshadows(humfile, sonpath, win, shadowmask, kvals, doplot)

Parameters
-------------
    humfile : str
       path to the .DAT file
    sonpath : str
       path where the *.SON files are
    win : int, *optional* [Default=100]
       window size (pixels) for the automated shadow removal algorithm
    shadowmask : int, *optional* [Default=0]
       1 = do manual shadow masking, otherwise do automatic shadow masking
    kvals : int, *optional* [Default=8]
       if automatic shadowmask, this parameter sets the number of k-means to calculate 
       (the one with the lowest value will be the shadow which is removed)
    doplot : int, *optional* [Default=1]
       1 = make plots, otherwise do not

Returns
---------
    sonpath+base+'_data_star_la.dat': memory-mapped file
        contains the starboard scan with water column removed and 
        radiometrically corrected, and shadows removed

    sonpath+base+'_data_port_la.dat': memory-mapped file
        contains the portside scan with water column removed and
        radiometrically corrected, and shadows removed


  .. image:: _static/pyhum_logo_colour_sm.png


