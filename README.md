pyhum.py

INFO:
Python script to read Humminbird DAT and associated SON files, and export data

Version: 1.0
Author:  Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.0      Revision: January, 2014

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of 'pyhum' software
This software is in the public domain because it contains materials that originally came 
from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 
http://www.usgs.gov/visual-id/credit_usgs.html#copyright

thanks to Barb Fagetter (blueseas@oceanecology.ca) for format info

This software has been tested with Python 2.7 on Linux Fedora and Windows 7.
This software has been used only with Humminbird 998 series instruments. 

SYNTAX:
python pyhum.py -i datfile -s sonpath
where datfile is the .DAT file associated with the survey, and sonpath is the (absolute or relative) path to where the associated .SON files are

Optional arguments:
-n number of CPU processors to use [Default = 8]
-d measured draft of instrument in metres [Default = 0]
-c coordinate transformation. By default coordinates are output in WGS84. If you would like an additional coordinate transformation, specify the EPSG ID number for use in cs2cs (pyproj). See: http://cs2cs.mygeodata.eu/; http://www.eye4software.com/resources/stateplane/ [Default is Arizona State Plane: -c "epsg:26949"]
-p make simple plots of data [Default = 1 (yes); 0 = no]

EXAMPLES:
1) run the provided test case with all defaults
python pyhum.py -i ./test.DAT -s ./test_data/

2) run a file and output eastings/northings in NAD Colorado Central (26954) with a draft of 0.4m
python pyhum.py -i ./test.DAT -s ./test_data/ -c "epsg:26954" -d 0.4

3) run a file and output eastings/northings in OSGB36 British National Grid (27700) with a draft of 1m using 4 processors
python pyhum.py -i ./test.DAT -s ./test_data/ -c "epsg:27700" -d 1 -n 4

4) run the provided test case with all defaults except no plots
python pyhum.py -i ./test.DAT -s ./test_data/ -p 0


OUTPUTS:
Several files are created which contain the raw and parsed data. They are prefixed by the root of the input file (*) followed by:
1) *.pkl = all raw data
2) *raw_port.pkl = port side scan, longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)
3) *raw_star.pkl = starboard side scan, longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)
4) *raw_low.pkl = low freq. downward sonar, longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)
5) *raw_hi.pkl = high freq. downward sonar, longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)
6) *meta.pkl = longitude, latitude, depth (m), speed (mph), easting (m), northing (m), time (unix epoch)

These are python pickled data format. To read use, for example:

with open('testraw_hi.pkl') as f:
      c_hi,lon,lat,spd,e,n,time = cPickle.load(f)

If doplot =1 (see above) the program will also create some rudimentary plots of the data (mainly to check everything is ok). These are stored in the same directory as the .son files and are hopefully self explanatory

PYTHON LIBRARIES YOU MAY NEED TO INSTALL:
1) Joblib: http://pythonhosted.org/joblib/
2) Pyproj: http://code.google.com/p/pyproj/
3) cPickle: http://docs.python.org/release/2.5/lib/module-cPickle.html
4) Numpy: http://www.scipy.org/scipylib/download.html
5) Matplotlib: http://matplotlib.org/downloads.html

OTHER LIBRARIES USED WHICH ARE USUALLY INSTALLED WITH PYTHON:
1) sys
2) getopt
3) glob
4) struct
5) time
6) array
7) tkinter
8) os

Please download, try, report bugs, fork, modify, evaluate, discuss. Thanks for stopping by!
