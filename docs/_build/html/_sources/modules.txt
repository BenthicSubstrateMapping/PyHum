.. _modules:

***************
Modules
***************

.. _overview:

Overview
=========

The programs in this package are as follows:

**read**: read Humminbird DAT and associated SON files, and export data in various formats

**correct**: read output **read**, perform some radiometric corrections and produce some rudimentary plots

**texture**: read radiometrically corrected Humminbird data (output from **correct**), perform a textural analysis using the spectral method of Buscombe et al (forthcoming) and produce some rudimentary plots

**map**: script to generate a point cloud (X,Y,sidescan intensity), save it to ascii format file, grid it and make a raster overlay for a kml file for google-earth

**map_texture**: script to generate a point cloud (X,Y,texture lengthscale), save it to ascii format file, grid it and make a raster overlay for a kml file for google-earth

These are all command-line/modular programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options


  .. image:: _static/pyhum_logo_colour_sm.png

