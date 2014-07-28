PyHum

INFO:
Python/Cython scripts to read Humminbird DAT and associated SON files, export data, carry out rudimentary radiometric corrections to data, and classify bed texture using the algorithm detailed in Buscombe et al (in prep), "Automated riverbed sediment classification using low-cost sidescan sonar", to be submitted to Journal of Hydraulic Engineering

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.0      Revision: July, 2014

For latest code version please visit:
https://github.com/dbuscombe-usgs

====================================
   This function is part of PyHum software
   This software is in the public domain because it contains materials that originally came 
   from the United States Geological Survey, an agency of the United States Department of Interior. 
   For more information, see the official USGS copyright policy at 
   http://www.usgs.gov/visual-id/credit_usgs.html#copyright

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government. 
====================================

thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 series instruments. 

The programs in this package are as follows:
1) pyhum_read.py
Python script to read Humminbird DAT and associated SON files, and export data in MAT format

2) pyhum_correct.py
Python script to read Humminbird data in MAT format (output from pyhum_read.py) and perform some radiometric corrections and produce some rudimentary plots

3) pyhum_texture.py
Python script to read radiometrically corrected Humminbird data in MAT format (output from pyhum_correct.py) and perform a textural analysis using the spectral method of Buscombe et al (in prep) and produce some rudimentary plots

These are all command-line programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options

Installation:

PYTHON LIBRARIES YOU MAY NEED TO INSTALL TO USE PyHum:
1) Joblib: http://pythonhosted.org/joblib/
2) Pyproj: http://code.google.com/p/pyproj/
3) SciPy: http://www.scipy.org/scipylib/download.html
4) Numpy: http://www.scipy.org/scipylib/download.html
5) Matplotlib: http://matplotlib.org/downloads.html
6) Scikit-learn: http://scikit-learn.org/stable/
7) Python Image LIbrary (PIL) http://www.pythonware.com/products/pil/

All of the above are available through pip (https://pypi.python.org/pypi/pip) and easy_install (https://pythonhosted.org/setuptools/easy_install.html)

OTHER LIBRARIES (CYTHON) NEED TO BE COMPILED FOR SPEED:
1) pyread.pyx
2) ppdrc.pyx
3) cwt.pyx
4) replace_nans.pyx
5) spec_noise.pyx

This compilation can be carried out by running the supplied bash script, 'compile_pyhum.sh'

To run the example using the data within 'test.DAT' and files within the folder 'test_data', run the bash script 'run_example.sh' 

This is a new project written and maintained by Daniel Buscombe. Thus far extensive testing has not been possible so bugs are expected. 

Please download, try, report bugs, fork, modify, evaluate, discuss. Please address all suggestions, comments and queries to: dbuscombe@usgs.gov. Thanks for stopping by! 



