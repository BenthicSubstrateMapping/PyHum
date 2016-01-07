### PyHum


PyHum - a Python framework for reading and processing data from a Humminbird low-cost sidescan sonar

PyHum is an open-source project dedicated to provide a generic Python framework 
for reading and exporting data from Humminbird(R) instruments, carrying out rudimentary radiometric corrections to the data,
classify bed texture, and produce some maps on aerial photos and kml files for google-earth


1. read Humminbird DAT and associated SON files
2. export data
3. carry out rudimentary radiometric corrections to data, and 
4. classify bed texture using the algorithm detailed in Buscombe, Grams, Smith, (2015) "Automated riverbed sediment classification using low-cost sidescan sonar", Journal of Hydraulic Engineering, 10.1061/(ASCE)HY.1943-7900.0001079, 06015019.
5. produce some maps on aerial photos and kml files for google-earth

The software is designed to read Humminbird data (.SON, .IDX, and .DAT files) and works on both sidescan and downward-looking echosounder data, where available.

Some aspects of the program are detailed in:
Buscombe, D., Grams, P.E., and Smith, S. (2015) "Automated riverbed sediment classification using low-cost sidescan sonar", Journal of Hydraulic Engineering, 10.1061/(ASCE)HY.1943-7900.0001079, 06015019.

Full documentation of the program is forthcoming


![alt tag](http://dbuscombe-usgs.github.io/figs/class_R01560.png)
*Sand dunes on the bed of the Colorado River in Grand Canyon, Arizona*

![alt tag](http://dbuscombe-usgs.github.io/figs/Texas_reef_merged_cropped.png)
*Submerged reef and a school of fish, in the Gulf of Mexico off the coast of Texas. Credit: Richard Kline and Michael Bollinger from the University of Texas Rio Grande Valley.  Funded by the Texas Parks and Wildlife Department, Artificial Reef Program.*

![alt tag](http://dbuscombe-usgs.github.io/figs/DogRiver_AllenAven.png)
*Boat marina in Dog River, Alabama. Credit: Allen Aven, Dauphin Island Sea Lab*

![alt tag](http://dbuscombe-usgs.github.io/figs/PyHum_glencanyon.png)
*Fine/coarse transition on the bed of the Colorado River in Glen Canyon, Arizona*

![alt tag](http://dbuscombe-usgs.github.io/figs/SuwaneeRiver_Florida.png)
*Lower Suwannee River, Florida. Credit: Florida Fish and Wildlife Commission: Information Science and Management. Funded by U.S. Fish and Wildlife Service State Wildlife Grant (SWG)*

### Contributing & Credits

 Primary Developer |    Daniel Buscombe 
 ------ | ---------------
         |  Grand Canyon Monitoring and Research Center
          | United States Geological Survey
          | Flagstaff, AZ 86001
          | dbuscombe@usgs.gov

 Co-Developer |    Daniel Hamill
 ------ | ---------------
         |  Department of Watershed Sciences
          | Utah State University
          | Logan, UT 84322
          | dhamill@usgs.gov

Version: 1.3.4    |  Revision: Jan, 2016

For latest code version please visit:
https://github.com/dbuscombe-usgs

This function is part of PyHum software
This software is in the public domain because it contains materials that originally came 
from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 

```
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
```

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government. 

Thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info, Dan Hamill (Utah State University) and Paul Anderson (Quest Geophysical Asia) for debugging and suggestions for improvements

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4 & 14.4, Windows 7.
This software has (so far) been used only with Humminbird 798, 998, 1198 and 1199 series instruments. 

### Contents

The programs in this package are as follows:

1. read
script to read Humminbird DAT and associated SON files, export data, and produce some rudimentary plots

2. correct
script to read Humminbird data (output from 'read') and perform some radiometric corrections and produce some rudimentary plots

3. rmshadows
read output 'correct', and remove dark shadows in scans caused by shallows, shorelines, and attenuation of acoustics with distance

4. texture
script to read radiometrically corrected Humminbird data in MAT format (output from pyhum_correct.py) and perform a textural analysis using the spectral method of Buscombe et al (forthcoming) and produce some rudimentary plots

5. map
script to generate a point cloud (X,Y,sidescan intensity), save it to ascii format file, grid it and make a raster overlay on an aerial image (pulled automatically from the ESRI GIS image server), and a kml file for showing the same thing in google-earth

6. map_texture
script to generate a point cloud (X,Y,texture lengthscale - calculated using pyhum_texture), save it to ascii format file, grid it and make a raster overlay on an aerial image (pulled automatically from the ESRI GIS image server), and a kml file for showing the same thing in google-earth

7. e1e2
script to analyse the first (e1, 'roughness') and second (e2, 'hardness') echo returns from the high-frequency downward looking echosounder, and generate generalised acoustic parameters for the purposes of point classification of submerged substrates/vegetation. The processing accounts for the absorption of sound in water, and does a basic k-means cluster of e1 and e2 coefficients into specified number of 'acoustic classes'. This code is based on code by Barb Fagetter (blueseas@oceanecology.ca). Georeferenced parameters are saved in csv form, and optionally plots and kml files are generated

These are all command-line programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options

![alt tag](http://dbuscombe-usgs.github.io/figs/PyHum_glencanyon_class.png)
*Automated bed sediment measurement, Colorado River in Glen Canyon*

## Setup

### Notes for Windows/Anaconda users

PyHum currently has only been tested with Python 2.7, so you'll need that version of Anaconda

Step 1. Before installing PyHum, install Basemap using:

```
conda install basemap
```

Step 2. Install pyproj. pip seems to have a bug with pyproj depending on what c-compiler your python distribution uses. Therefore, you may have to install pyproj (and other dependencies) from `here <http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyproj>`_

a) download the .whl file. Choose the file associated with python 2.7 ("cp27") and the architecture you are using, i.e. 32-bit (win32) or 64-bit (amd64)
b) then move that file to your root Anaconda directory. This is probably C:\Users\yourusername\AppData\Local\Continuum\Anaconda (when you open an Anaconda command prompt it's the directory that's listed before the prompt '>')
c) then use pip to install it, e.g.:


```
pip install pyproj-1.9.4-cp27-none-win_amd64.whl
```

Step 4. Install simplekml, using:

```
pip install simplekml
```

Step 3. Assuming a Anaconda distribution which comes with almost all required program dependencies:

```
pip uninstall PyHum (removes any previous installation)
pip install PyHum
```

If you get C++ compiler errors (such as "Unable to find vcvarsall.bat"), you will need to install the Microsoft Visual C++ compiler from `here <http://aka.ms/vcpython27>`_


(Advanced) If you have git installed (from `here <https://git-scm.com/download/win>`_), you can install the latest 'bleeding edge' (pre-release) version directly from github:

```
pip install git+https://github.com/dbuscombe-usgs/PyHum.git
```


### Notes for Linux users

```
pip uninstall PyHum (removes previous installation)
pip install PyHum
```

or:

```
pip install git+https://github.com/dbuscombe-usgs/PyHum.git
```

Automatic Installation from github:

```
git clone git@github.com:dbuscombe-usgs/PyHum.git
cd PyHum
python setup.py install
```

or a local installation:

```
python setup.py install --user
```

or with admin privileges, e.g.:

```
sudo python setup.py install
```


You could try before you install, using a virtual environment:

```
virtualenv venv
source venv/bin/activate
pip install numpy
pip install Cython
pip install scipy
pip install simplekml
pip install pyproj
pip install scikit-learn
pip install Pillow
pip install matplotlib
pip install basemap --allow-external basemap --allow-unverified basemap
pip install pyresample
pip install PyHum
python -c "import PyHum; PyHum.test.dotest()"
deactivate (or source venv/bin/deactivate)
```

The results will live in "venv/lib/python2.7/site-packages/PyHum"

Note for Fedora linux users: you need the geos-devel package for basemap, and the blas and libpack libraries for scipy


###Manual Installation

PYTHON LIBRARIES YOU MAY NEED TO INSTALL TO USE PyHum:

1. Pyproj (http://code.google.com/p/pyproj/)
2. SciPy (http://www.scipy.org/scipylib/download.html)
3. Numpy (http://www.scipy.org/scipylib/download.html)
4. Matplotlib (http://matplotlib.org/downloads.html)
5. Scikit-learn (http://scikit-learn.org/stable/)
6. Python Image LIbrary (PIL) (http://www.pythonware.com/products/pil/)
7. simplekml (http://simplekml.readthedocs.org/en/latest/index.html)
8. pyproj (https://pypi.python.org/pypi/pyproj)
9. basemap (http://matplotlib.org/basemap/)
10. pyresample (http://pyresample.readthedocs.org/en/latest/index.html#)

All of the above are available through pip (https://pypi.python.org/pypi/pip) and easy_install (https://pythonhosted.org/setuptools/easy_install.html)

### Test

A test can be carried out by running the supplied script:

```
python -c "import PyHum; PyHum.test.dotest()"
```

which carries out the following operations:

```
   # copy files over to somewhere read/writeable
   dircopy(PyHum.__path__[0], os.path.expanduser("~")+os.sep+'pyhum_test')
   shutil.copy(PyHum.__path__[0]+os.sep+'test.DAT', os.path.expanduser("~")+os.sep+'pyhum_test'+os.sep+'test.DAT')

   # general settings   
   humfile = os.path.normpath(os.path.join(os.path.expanduser("~"),'pyhum_test','test.DAT'))
   sonpath = os.path.normpath(os.path.join(os.path.expanduser("~"),'pyhum_test'))

   doplot = 1 #yes

   # reading specific settings
   cs2cs_args = "epsg:26949" #arizona central state plane
   bedpick = 1 # auto bed pick
   c = 1450 # speed of sound fresh water
   t = 0.108 # length of transducer
   draft = 0.3 # draft in metres
   flip_lr = 1 # flip port and starboard
   model = 998 # humminbird model
   cog = 1 # GPS course-over-ground used for heading
   calc_bearing = 0 #no
   filt_bearing = 0 #no
   chunk = 'd100' # distance, 100m
   #chunk = 'p1000' # pings, 1000
   #chunk = 'h10' # heading deviation, 10 deg
          
   # correction specific settings
   maxW = 1000 # rms output wattage
   dofilt = 0 # 1 = apply a phase preserving filter (WARNING!! takes a very long time for large scans)
   correct_withwater = 0 # don't retain water column in radiometric correction (1 = retains water column for radiomatric corrections)
   ph = 7.0 # acidity on the pH scale
   temp = 10.0 # water temperature in degrees Celsius
   salinity = 0.0

   # for shadow removal
   shadowmask = 0 #automatic shadow removal

   # for texture calcs
   win = 100 # pixel window
   shift = 10 # pixel shift
   density = win/2 
   numclasses = 4 # number of discrete classes for contouring and k-means
   maxscale = 20 # Max scale as inverse fraction of data length (for wavelet analysis)
   notes = 4 # Notes per octave (for wavelet analysis)

   # for mapping
   res = 0.1 # grid resolution in metres
   mode = 1 # gridding mode (simple nearest neighbour)
   #mode = 2 # gridding mode (inverse distance weighted nearest neighbour)
   #mode = 3 # gridding mode (gaussian weighted nearest neighbour)
   dowrite = 0 #disable writing of point cloud data to file

   nn = 64 #number of nearest neighbours for gridding (used if mode > 1)
   influence = 1 #Radius of influence used in gridding. Cut off distance in meters 
   numstdevs = 4 #Threshold number of standard deviations in sidescan intensity per grid cell up to which to accept 

   # for downward-looking echosounder echogram (e1-e2) analysis
   beam = 20.0
   transfreq = 200.0 # frequency (kHz) of downward looking echosounder
   integ = 5
   numclusters = 3 # number of acoustic classes to group observations

   # read data in SON files into PyHum memory mapped format (.dat)
   PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, bedpick, flip_lr, model, calc_bearing, filt_bearing, cog, chunk)

   # correct scans and remove water column
   PyHum.correct(humfile, sonpath, maxW, doplot, dofilt, correct_withwater, ph, temp, salinity)

   # remove acoustic shadows (caused by distal acoustic attenuation or sound hitting shallows or shoreline)
   PyHum.rmshadows(humfile, sonpath, win, shadowmask, doplot)

   # Calculate texture lengthscale maps using the method of Buscombe et al. (2015)
   PyHum.texture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

   # grid and map the scans
   PyHum.map(humfile, sonpath, cs2cs_args, res, dowrite, mode, nn, influence, numstdevs)

   res = 0.5 # grid resolution in metres
   numstdevs = 5
   
   # grid and map the texture lengthscale maps
   PyHum.map_texture(humfile, sonpath, cs2cs_args, res, dowrite, mode, nn, influence, numstdevs)

   # calculate and map the e1 and e2 acoustic coefficients from the downward-looking sonar
   PyHum.e1e2(humfile, sonpath, cs2cs_args, ph, temp, salinity, beam, transfreq, integ, numclusters, doplot)

```

on the following files:

1. test.DAT
2. B003.SON
3. B002.SON
4. B001.SON
5. B000.SON

and results in a set of outputs such as csv, mat and kml files, and including some rudimentary figures such as:

![alt tag](http://dbuscombe-usgs.github.io/figs/bed_pick.png)
*port and starboard scans showing automated bed picks*

![alt tag](http://dbuscombe-usgs.github.io/figs/merge_corrected_scan_ppdrc.png)
*a merged port/starboard scan*

![alt tag](http://dbuscombe-usgs.github.io/figs/raw_dwnhi.png)
*a raw 200 kHz downward sonar scan*

![alt tag](http://dbuscombe-usgs.github.io/figs/raw_dwnlow.png)
*a raw 83 kHz downward sonar scan*

![alt tag](http://dbuscombe-usgs.github.io/figs/testclass1.png)
*radiometrically corrected scan (top) and wavelet lengthscale classification (bottom)*

![alt tag](http://dbuscombe-usgs.github.io/figs/testclass_kmeans1.png)
*k=4 means lengthscale classification*


### Getting Started

Inputs to the program are a .DAT file (e.g. R0089.DAT) and a folder of .SON and .IDX files (e.g. /my/folder/R0089). The program will read the .SON files with or without the accompanying .IDX files, but will be faster if the .IDX files are present. 

PyHum is modular so can be called from within a python or ipython console, from an IDE (such as IDLE or Spyder), or by running a script.

The following example script could be saved as, for example "proc_mysidescandata.py" and run from the command line using

```
python proc_mysidescandata.py -i C:\MyData\R0087.DAT -s C:\MyData\R0087
```


```
import sys, getopt

from Tkinter import Tk
from tkFileDialog import askopenfilename, askdirectory

import PyHum
import os

if __name__ == '__main__': 

    argv = sys.argv[1:]
    humfile = ''; sonpath = ''
    
    # parse inputs to variables
    try:
       opts, args = getopt.getopt(argv,"hi:s:")
    except getopt.GetoptError:
         print 'error'
         sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
         print 'help'
         sys.exit()
       elif opt in ("-i"):
          humfile = arg
       elif opt in ("-s"):
          sonpath = arg

    # prompt user to supply file if no input file given
    if not humfile:
       print 'An input file is required!!!!!!'
       Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
       humfile = askopenfilename(filetypes=[("DAT files","*.DAT")]) 

    # prompt user to supply directory if no input sonpath is given
    if not sonpath:
       print 'A *.SON directory is required!!!!!!'
       Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
       sonpath = askdirectory() 

    # print given arguments to screen and convert data type where necessary
    if humfile:
       print 'Input file is %s' % (humfile)

    if sonpath:
       print 'Son files are in %s' % (sonpath)
                 
    doplot = 1 #yes

    # reading specific settings
    cs2cs_args = "epsg:26949" #arizona central state plane
    bedpick = 1 # auto bed pick
    c = 1450 # speed of sound fresh water
    t = 0.108 # length of transducer
    draft = 0.3 # draft in metres
    flip_lr = 1 # flip port and starboard
    model = 998 # humminbird model
    cog = 1 # GPS course-over-ground used for heading
    calc_bearing = 0 #no
    filt_bearing = 0 #no
    chunk = 'd100' # distance, 100m
     
    # correction specific settings
    maxW = 1000 # rms output wattage
    dofilt = 0 # 1 = apply a phase preserving filter (WARNING!! takes a very long time for large scans)
    correct_withwater = 0 # don't retain water column in radiometric correction (1 = retains water column for radiomatric corrections)
    ph = 7.0 # acidity on the pH scale
    temp = 10.0 # water temperature in degrees Celsius
    #salinity = 0.0

    # for shadow removal
    shadowmask = 1 #manual shadow removal

    # for mapping
    calc_bearing = 0 #no
    filt_bearing = 0 #no
    res = 0.05 # grid resolution in metres
    cog = 1 # GPS course-over-ground used for heading

    PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, bedpick, flip_lr, chunk_size, model)

    PyHum.correct(humfile, sonpath, maxW, doplot)

    PyHum.rmshadows(humfile, sonpath, win, shadowmask, kvals, doplot)

    PyHum.map(humfile, sonpath, cs2cs_args, calc_bearing, filt_bearing, res, cog, dowrite)
```

or from within ipython (with a GUI prompt to navigate to the files):

```
   run proc_mysidescandata.py
```


### Trouble Shooting

1. Problem: pyhum read hangs for a long time (several minutes) on the test script. 
Try this: uninstall joblib and install an older version::

```
   pip uninstall joblib
   pip install joblib==0.7.1
```

2. Problem: you get an "invalid mode or file name" error.
Try this: construct file paths using raw strings e.g.:: 

```
   r'C:\Users\me\mydata\R0089' 
```

or using os, e.g.::

```
   import os
   os.path.abspath(os.path.join('C:\Users','me','mydata','R0089'))
```


### Support

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please use the 'Issues' tab in github

https://github.com/dbuscombe-usgs/PyHum

for all bug reports, suggestions, comments and queries. Thanks for stopping by! 



