.. pyhum.correct:

pyhum.correct module
======================

    Remove water column and carry out some rudimentary radiometric corrections, 
    accounting for directivity and attenuation with range

Syntax
----------

You call the function like this::

  [] = PyHum.correct(humfile, sonpath, maxW, doplot)

Parameters
-------------
    humfile : str
       path to the .DAT file
    sonpath : str
       path where the *.SON files are
    maxW : int, *optional* [Default=1000]
       maximum transducer power
    doplot : int, *optional* [Default=1]
       1 = make plots, otherwise do not

Returns
---------
    sonpath+base+'_data_star_l.dat': memory-mapped file
        contains the starboard scan with water column removed

    sonpath+base+'_data_port_l.dat': memory-mapped file
        contains the portside scan with water column removed

    sonpath+base+'_data_star_la.dat': memory-mapped file
        contains the starboard scan with water column removed and 
        radiometrically corrected

    sonpath+base+'_data_port_la.dat': memory-mapped file
        contains the portside scan with water column removed and
        radiometrically corrected

    sonpath+base+'_data_range.dat': memory-mapped file
        contains the cosine of the range which is used to correct
        for attenuation with range

    sonpath+base+'_data_dwnlow_l.dat': memory-mapped file
        contains the low freq. downward scan with water column removed

    sonpath+base+'_data_dwnhi_l.dat': memory-mapped file
        contains the high freq. downward  scan with water column removed

    sonpath+base+'_data_dwnlow_la.dat': memory-mapped file
        contains the low freq. downward  scan with water column removed and 
        radiometrically corrected

    sonpath+base+'_data_dwnhi_la.dat': memory-mapped file
        contains the high freq. downward  scan with water column removed and
        radiometrically corrected


  .. image:: _static/pyhum_logo_colour_sm.png


