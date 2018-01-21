## PyHum

![pyhum_logo_colour_sm](https://user-images.githubusercontent.com/3596509/35187745-ed45a05a-fde5-11e7-8f7e-5142b59dc772.png)

A Python framework for reading and processing data from a Humminbird low-cost sidescan sonar

Project website [here](http://dbuscombe-usgs.github.io/PyHum/) for more details

PyHum is an open-source project dedicated to provide a generic Python framework 
for reading and exporting data from Humminbird(R) instruments, carrying out rudimentary radiometric corrections to the data,
classify bed texture, and produce some maps on aerial photos and kml files for google-earth

1. read Humminbird DAT and associated SON files
2. export data
3. carry out rudimentary radiometric corrections to data, and 
4. classify bed texture using the algorithm detailed in Buscombe, Grams, Smith, (2015) "Automated riverbed sediment classification using low-cost sidescan sonar", Journal of Hydraulic Engineering, 10.1061/(ASCE)HY.1943-7900.0001079, 06015019.
5. produce some maps on aerial photos and kml files for google-earth

The software is designed to read Humminbird data (.SON, .IDX, and .DAT files) and works on both sidescan and downward-looking echosounder data, where available.

### Please cite! 

If you use PyHum in your published work, please cite the following papers:

1. Buscombe, D., Grams, P.E., and Smith, S. (2015) "Automated riverbed sediment classification using low-cost sidescan sonar", Journal of Hydraulic Engineering, 10.1061/(ASCE)HY.1943-7900.0001079, 06015019.

2. Buscombe, D., 2017, Shallow water benthic imaging and substrate characterization using recreational-grade sidescan-sonar. ENVIRONMENTAL MODELLING & SOFTWARE 89, 1-18.


### Contributing & Credits

 Primary Developer |    Daniel Buscombe 
 ------ | ---------------
         |  Northern Arizona University
          | Flagstaff, AZ 86011
          | daniel.buscombe@nau.edu

 Co-Developer |    Daniel Hamill
 ------ | ---------------
         |  Department of Watershed Sciences
          | Utah State University
          | Logan, UT 84322
          | hamill.daniel@gmail.com

Version: 1.4.2    |  Revision: Jan, 2018


### Please Read

1. This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12 -- 17, Windows 7. Python 3 is not yet supported

2. This software has (so far) been used only with Humminbird 700, 800, 900, 1100, HELIX, MEGA and ONIX series instruments. 

3. PyHum is not yet tested with SOLIX, ICE, ION, PMAX systems. Please make example data available and we'll see what we can do!

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

8. gui
A graphical user interface which essentially serves as a 'wrapper' to the above functions, allowing graphical input of processing options and sequential analysis of the data using PyHum modules

These are all command-line programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options

<!--![alt tag](http://dbuscombe-usgs.github.io/figs/PyHum_glencanyon_class.png)-->
<!--*Automated bed sediment measurement, Colorado River in Glen Canyon*-->

## Setup

### PyHum only works in python 2.X. Python 3 is not yet supported. 

### Installing in a conda virtual env (recommended)

In a conda (miniconda/anaconda) python 2 environment: 

Linux:

```
conda create --name pyhum python=2
source activate pyhum
conda install -c conda-forge basemap-data-hires -y
conda install scipy numpy scikit-image
pip install simplekml sklearn pandas dask
pip install joblib toolz cython
pip install pyresample
pip install PyHum --no-deps #(or pip install git+https://github.com/dbuscombe-usgs/PyHum.git --no-deps)
python -c"import PyHum;PyHum.dotest()" 
```

Windows:

```
conda create --name pyhum python=2
activate pyhum
conda install -c conda-forge basemap-data-hires -y
conda install scipy numpy scikit-image
pip install simplekml sklearn pandas dask
pip install joblib toolz cython
pip install pyresample==1.1.4
pip install PyHum --no-deps #(or pip install git+https://github.com/dbuscombe-usgs/PyHum.git --no-deps)
python -c"import PyHum;PyHum.dotest()" 
```


### Installing as a library accessible outside of virtual env

1. From PyPI::

```
pip install PyHum
```

2. the latest 'bleeding edge' (pre-release) version directly from github::

```
pip install git+https://github.com/dbuscombe-usgs/PyHum.git
```

(Windows users) install git from here: https://git-scm.com/download/win


3. from github repo clone::

```
git clone git@github.com:dbuscombe-usgs/PyHum.git
cd PyHum
python setup.py install
```

or a local installation:

```
python setup.py install --user
```


4. linux users, using a virtual environment:

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
pip install toolz
pip install dask
pip install pandas
pip install PyHum
python -c "import PyHum; PyHum.test.dotest()"
deactivate (or source venv/bin/deactivate)
```

The results will live in "venv/lib/python2.7/site-packages/PyHum". Note for Fedora linux users: you need the geos-devel package for basemap, and the blas and libpack libraries for scipy


### Running the test
A test can be carried out by running the supplied script. From the command line (terminal)::

```
python -c"import PyHum;PyHum.dotest()" 
```

or (if python3 is your default python)::

```
python2 -c"import PyHum;PyHum.dotest()" 
```

### Using the GUI
From the command line (terminal)::

```
python -c "import PyHum; PyHum.gui()"
```

## Using PyHum

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
    #salinity = 0.0

    # for shadow removal
    shadowmask = 1 #manual shadow removal

    # for mapping
    res = 99 # grid resolution in metres
    # if res==99, the program will automatically calc res from the spatial res of the scans
    mode = 1 # gridding mode (simple nearest neighbour)
    #mode = 2 # gridding mode (inverse distance weighted nearest neighbour)
    #mode = 3 # gridding mode (gaussian weighted nearest neighbour)
    dowrite = 0 #disable writing of point cloud data to file
    scalemax = 60 # max color scale value (60 is a good place to start)

    ## read data in SON files into PyHum memory mapped format (.dat)
    PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, bedpick, flip_lr, model, calc_bearing, filt_bearing, chunk) #cog

    ## correct scans and remove water column
    PyHum.correct(humfile, sonpath, maxW, doplot, dofilt, correct_withwater, ph, temp, salinity)

    ## remove acoustic shadows (caused by distal acoustic attenuation or sound hitting shallows or shoreline)
    PyHum.rmshadows(humfile, sonpath, win, shadowmask, doplot)
   
    ## Calculate texture lengthscale maps using the method of Buscombe et al. (2015)
    win = 10
    PyHum.texture2(humfile, sonpath, win, doplot, numclasses)

    ## grid and map the scans
    PyHum.map(humfile, sonpath, cs2cs_args, res, mode, nn, numstdevs, use_uncorrected, scalemax) #dowrite, 


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

If you get C++ compiler errors (such as "Unable to find vcvarsall.bat"), you will need to install the Microsoft Visual C++ compiler from here: http://aka.ms/vcpython27


### Support

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them. Please use the 'Issues' tab in github

https://github.com/dbuscombe-usgs/PyHum

Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. 

Project website [here](http://dbuscombe-usgs.github.io/PyHum/) for more details

Thanks for stopping by! 

### Acknowledgements

This function is part of PyHum software
This software is in the public domain because it contains materials that originally came 
from the United States Geological Survey, an agency of the United States Department of Interior. 
For more information, see the official USGS copyright policy at 

```
http://www.usgs.gov/visual-id/credit_usgs.html#copyright
```

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government. 

Thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info, Dan Hamill (Utah State University), Paul Anderson (Quest Geophysical Asia) and various others for debugging and suggestions for improvements


