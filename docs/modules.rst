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

**e1e2**: script to analyse the first (e1, 'roughness') and second (e2, 'hardness') echo returns from the high-frequency downward looking echosounder, and generate generalised acoustic parameters for the purposes of point classification of submerged substrates/vegetation. The processing accounts for the absorption of sound in water, and does a basic k-means cluster of e1 and e2 coefficients into specified number of 'acoustic classes'. This code is based on code by Barb Fagetter (blueseas@oceanecology.ca). Georeferenced parameters are saved in csv form, and optionally plots and kml files are generated

These are all command-line/modular programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options


  .. image:: _static/pyhum_logo_colour_sm.png

