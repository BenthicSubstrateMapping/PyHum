### PyHum

Python/Cython scripts to: 

1. read Humminbird DAT and associated SON files
2. export data
3. carry out rudimentary radiometric corrections to data, and 
4. classify bed texture using the algorithm detailed in Buscombe, Grams, Smith, "Automated riverbed sediment classification using low-cost sidescan sonar", forthcoming.
5. produce some maps on aerial photos and kml files for google-earth


![alt tag](http://dbuscombe-usgs.github.io/figs/class_R01560.png)
*Sand dunes on the bed of the Colorado River in Grand Canyon*


![alt tag](http://dbuscombe-usgs.github.io/figs/R00426_map.jpg)
*Dog River, Alabama. Credit: Allen Aven, Dauphin Island Sea Lab*

![alt tag](http://dbuscombe-usgs.github.io/figs/PyHum_glencanyon.png)
*Fine/coarse transition on the bed of the Colorado River in Glen Canyon*

### Contributing & Credits

 Author |    Daniel Buscombe 
 ------ | ---------------
         |  Grand Canyon Monitoring and Research Center
          | United States Geological Survey
          | Flagstaff, AZ 86001
          | dbuscombe@usgs.gov
Version: 1.1.9    |  Revision: May, 2015

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

thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info

This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4 & 14.4, Windows 7.
This software has (so far) been used only with Humminbird 798, 998, 1198 and 1199 series instruments. 

### Contents

The programs in this package are as follows:

1. read
script to read Humminbird DAT and associated SON files, export data, and produce some rudimentary plots

2. correct
script to read Humminbird data (output from 'read') and perform some radiometric corrections and produce some rudimentary plots

3. texture
script to read radiometrically corrected Humminbird data in MAT format (output from pyhum_correct.py) and perform a textural analysis using the spectral method of Buscombe et al (forthcoming) and produce some rudimentary plots

4. map
script to generate a point cloud (X,Y,sidescan intensity), save it to ascii format file, grid it and make a raster overlay on an aerial image (pulled automatically from the ESRI GIS image server), and a kml file for showing the same thing in google-earth

5. map_texture
script to generate a point cloud (X,Y,texture lengthscale - calculated using pyhum_texture), save it to ascii format file, grid it and make a raster overlay on an aerial image (pulled automatically from the ESRI GIS image server), and a kml file for showing the same thing in google-earth

6. e1e2
script to analyse the first (e1, 'roughness') and second (e2, 'hardness') echo returns from the high-frequency downward looking echosounder, and generate generalised acoustic parameters for the purposes of point classification of submerged substrates/vegetation. The processing accounts for the absorption of sound in water, and does a basic k-means cluster of e1 and e2 coefficients into specified number of 'acoustic classes'. This code is based on code by Barb Fagetter (blueseas@oceanecology.ca). Georeferenced parameters are saved in csv form, and optionally plots and kml files are generated

These are all command-line programs which take a number of input (some required, some optional). Please see the individual files for a comprehensive list of input options

![alt tag](http://dbuscombe-usgs.github.io/figs/PyHum_glencanyon_class.png)
*Automated bed sediment measurement, Colorado River in Glen Canyon*

## Setup

### Automatic Installation from PyPI 

```
pip uninstall PyHum (removes previous installation)
pip install PyHum
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

### Notes for Windows/Anaconda users

Before installing PyHum, install Basemap using::

```
conda install basemap
```

Assuming a Anaconda distribution which comes with almost all required program dependencies:

```
pip install simplekml
pip uninstall PyHum (removes previous installation)
pip install PyHum
```

pip seems to have a bug with pyproj depending on what c-compiler your python distribution uses. Therefore, you may have to install pyproj (and other dependencies) from here:

http://www.lfd.uci.edu/~gohlke/pythonlibs/

Download the .whl file, then use pip to install it, e.g.

```
pip install pyproj-1.9.4-cp27-none-win_amd64.whl
```


### Notes for Linux users

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
   f = 455 # frequency kHz
   draft = 0.3 # draft in metres
   flip_lr = 1 # flip port and starboard
   model = 998 # humminbird model
   chunk_size = 1000 # chunk size = 1000 pings
   dowrite = 0 #disable writing of point cloud data to file

   # correction specific settings
   maxW = 1000 # rms output wattage

   # for texture calcs
   win = 50 # pixel window
   shift = 10 # pixel shift
   density = win/2 
   numclasses = 4 # number of discrete classes for contouring and k-means
   maxscale = 20 # Max scale as inverse fraction of data length (for wavelet analysis)
   notes = 4 # Notes per octave (for wavelet analysis)

   # for mapping
   dogrid = 1 # yes
   calc_bearing = 0 #no
   filt_bearing = 1 #yes
   res = 0.2 # grid resolution in metres
   cog = 1 # GPS course-over-ground used for heading

   # for downward-looking echosounder echogram (e1-e2) analysis
   ph = 7.0 # acidity on the pH scale
   temp = 10.0 # water temperature in degrees Celsius
   salinity = 0.0
   beam = 20.0
   transfreq = 200.0
   integ = 5
   numclusters = 3

   PyHum.read(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunk_size, model)

   PyHum.correct(humfile, sonpath, maxW, doplot)

   PyHum.texture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

   PyHum.map(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res, cog, dowrite)

   res = 0.5 # grid resolution in metres
   
   PyHum.map_texture(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res, cog, dowrite)

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

### Support

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please address all suggestions, comments and queries to: dbuscombe@usgs.gov. Thanks for stopping by! 



