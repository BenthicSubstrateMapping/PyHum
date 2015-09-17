.. pyhum.read:

pyhum.read module
======================

    Read a .DAT and associated set of .SON files recorded by a Humminbird(R)
    instrument. 
    
    Parse the data into a set of memory mapped files that will
    subsequently be used by the other functions of the PyHum module. 
    
    Export time-series data and metadata in other formats. 
    
    Create a kml file for visualising boat track

Syntax
----------

You call the function like this::

   [] = PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunksize, model, calc_bearing, filt_bearing, cog, chunk)

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
    c : float, *optional* [Default=1450.0]
       speed of sound in water (m/s). Defaults to a value of freshwater
    draft : float, *optional* [Default=0.3]
       draft from water surface to transducer face (m)
    doplot : float, *optional* [Default=1]
       if 1, plots will be made
    t : float, *optional* [Default=0.108]
       length of transducer array (m).
       Default value is that of the 998 series Humminbird(R)
    f : float, *optional* [Default=455]
       frequency of sidescan transducer in kHz
    bedpick : int, *optional* [Default=1]
       if 1, bedpicking with be carried out automatically
       if 0, user will be prompted to pick the bed location on screen
    flip_lr : int, *optional* [Default=0]
       if 1, port and starboard scans will be flipped
       (for situations where the transducer is flipped 180 degrees)
    model: int, *optional* [Default=998]
       A 3 or 4 number code indicating the model number 
       Examples: 998, 997, 1198, 1199
    cog : int, *optional* [Default=1]
       if 1, heading calculated assuming GPS course-over-ground rather than
       using a compass
    calc_bearing : float, *optional* [Default=0]
       if 1, bearing will be calculated from coordinates
    filt_bearing : float, *optional* [Default=0]
       if 1, bearing will be filtered
    chunk : str, *optional* [Default='d100' (distance, 100 m)]
       letter, followed by a number.

       There are the following letter options:

       'd' - parse chunks based on distance, then number which is distance in m

       'p' - parse chunks based on number of pings, then number which is number of pings 

       'h' - parse chunks based on change in heading, then number which is the change in heading in degrees
       '1' - process just 1 chunk

Returns
----------

    sonpath+base+'_data_port.dat': memory-mapped file
        contains the raw echogram from the port side
        sidescan sonar (where present)

    sonpath+base+'_data_port.dat': memory-mapped file
        contains the raw echogram from the starboard side
        sidescan sonar (where present)

    sonpath+base+'_data_dwnhi.dat': memory-mapped file
        contains the raw echogram from the high-frequency
        echosounder (where present)

    sonpath+base+'_data_dwnlow.dat': memory-mapped file
        contains the raw echogram from the low-frequency
        echosounder (where present)
        
    sonpath+base+"trackline.kml": google-earth kml file
        contains the trackline of the vessel during data
        acquisition
     
    sonpath+base+'rawdat.csv': comma separated value file
        contains time-series data. columns corresponding to
        longitude
        latitude
        easting (m)
        northing (m)
        depth to bed (m)
        alongtrack cumulative distance (m)
        vessel heading (deg.)
     
    sonpath+base+'meta.mat': .mat file
        matlab format file containing a dictionary object
        holding metadata information. Fields are:
        e : ndarray, easting (m)
        n : ndarray, northing (m)
        es : ndarray, low-pass filtered easting (m)
        ns : ndarray, low-pass filtered northing (m)
        lat : ndarray, latitude
        lon : ndarray, longitude
        shape_port : tuple, shape of port scans in memory mapped file
        shape_star : tuple, shape of starboard scans in memory mapped file
        shape_hi : tuple, shape of high-freq. scans in memory mapped file
        shape_low : tuple, shape of low-freq. scans in memory mapped file
        dep_m : ndarray, depth to bed (m)
        dist_m : ndarray, distance along track (m)
        heading : ndarray, heading of vessel (deg. N)
        pix_m: float, size of 1 pixel in across-track dimension (m)
        bed : ndarray, depth to bed (m)
        c : float, speed of sound in water (m/s)
        t : length of sidescan transducer array (m)
        f : frequency of sidescan sound (kHz)
        spd : ndarray, vessel speed (m/s)
        time_s : ndarray, time elapsed (s)
        caltime : ndarray, unix epoch time (s)

  .. image:: _static/pyhum_logo_colour_sm.png


