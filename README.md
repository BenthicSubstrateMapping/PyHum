pyhum.py

INFO:
Python script to read Humminbird DAT and associated SON files, and export data

Author:    Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.0      Revision: June, 2014

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'PyHum' software
This software is in the public domain because it contains materials that originally came from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4, and Windows 7.
This software has (so far) been used only with Humminbird 998 series instruments. 

SYNTAX:
python pyhum_read.py -i datfile -s sonpath
where datfile is the .DAT file associated with the survey, and sonpath is the (absolute or relative) path to where the associated .SON files are

Optional arguments:
-d measured draft of instrument in metres [Default = 0]
-c coordinate transformation. By default coordinates are output in WGS84. If you would like an additional coordinate transformation, specify the EPSG ID number for use in cs2cs (pyproj). See: http://cs2cs.mygeodata.eu/; http://www.eye4software.com/resources/stateplane/ [Default is Arizona State Plane: -c "epsg:26949"]
-p make simple plots of data [Default = 1 (yes); 0 = no]

EXAMPLES:
1) show help
python pyhum_read.py -h

2) run the provided test case with all defaults
python pyhum_read.py -i ./test.DAT -s ./test_data/ (linux)
python pyhum_read.py -i test.DAT -s \test_data\ (windows)

3) run a file and output eastings/northings in NAD Colorado Central (26954) with a draft of 0.4m
python pyhum_read.py -i ./test.DAT -s ./test_data/ -c "epsg:26954" -d 0.4

4) run a file and output eastings/northings in OSGB36 British National Grid (27700) with a draft of 1m 
python pyhum_read.py -i ./test.DAT -s ./test_data/ -c "epsg:27700" -d 1 

5) run the provided test case with all defaults except no plots
python pyhum_read.py -i ./test.DAT -s ./test_data/ -p 0

OUTPUTS:
Files are created which contain the raw and parsed meta data. They are prefixed by the root of the input file (*) followed by:
1) *.mat = all raw data (port side scan, starboard side scan, downward low freq, downward high freq)
2) *meta.mat = longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)

These are python/matlab/octave .mat data format. To read use, for example:
data = loadmat('test.mat')

If doplot =1 (see above) the program will also create some rudimentary plots of the data (mainly to check everything is ok). These are stored in the same directory as the .son files and are hopefully self explanatory

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


Please download, try, report bugs, fork, modify, evaluate, discuss. Thanks for stopping by!
