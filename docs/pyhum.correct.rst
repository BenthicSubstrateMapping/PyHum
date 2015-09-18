.. pyhum.correct:

pyhum.correct module
======================

    Remove water column and carry out some rudimentary radiometric corrections, 
    accounting for directivity and attenuation with range

Syntax
----------

You call the function like this::

  [] = PyHum.correct(humfile, sonpath, maxW, doplot, correct_withwater, ph, temp, salinity, dconcfile)

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

    dofilt : int, *optional* [Default=0]
       1 = apply a phase preserving filter to the scans

    correct_withwater : int, *optional* [Default=0]
       1 = apply radiometric correction but don't remove water column from scans

    ph : float, *optional* [Default=7.0]
       water acidity in pH

    temp : float, *optional* [Default=10.0]
       water temperature in degrees Celsius

    salinity : float, *optional* [Default=0.0]
       salinity of water in parts per thousand

    dconcfile : str, *optional* [Default=None]
       file path of a text file containing sediment concentration data
       this file must contain the following fields separated by spaces:
       size (microns) conc (mg/L) dens (kg/m3)
       with one row per grain size, for example:
       30 1700 2200
       100 15 2650

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
    
    if correct_withwater == 1:
    
       sonpath+base+'_data_star_lw.dat': memory-mapped file
           contains the starboard scan with water column retained and 
           radiometrically corrected

       sonpath+base+'_data_port_lw.dat': memory-mapped file
           contains the portside scan with water column retained and
           radiometrically corrected

  .. image:: _static/pyhum_logo_colour_sm.png


