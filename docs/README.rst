.. _getting_started:


***************
Getting started
***************

.. _about:

About
======

PyHum - a Python framework for reading and processing data from low-cost sidescan sonar

PyHum is an open-source project dedicated to provide a generic Python framework 
for reading and exporting data from Humminbird(R) instruments, carrying out rudimentary radiometric corrections to the data,
classify bed texture, and produce some maps on aerial photos and kml files for google-earth

Some aspects of the program are detailed in:
Buscombe, D., Grams, P.E., and Smith, S. "Automated riverbed sediment classification using low-cost sidescan sonar", forthcoming.

For the source code visit `the project github site <https://github.com/dbuscombe-usgs/PyHum/>`_


.. _license:

License
========

This software is in the public domain because it contains materials that
originally came from the United States Geological Survey, an agency of the
United States Department of Interior. For more information, 
see `the official USGS copyright policy <http://www.usgs.gov/visual-id/credit_usgs.html#copyright>`_

Any use of trade, product, or firm names is for descriptive purposes only 
and does not imply endorsement by the U.S. government.

This software is issued under the `GNU Lesser General Public License, Version 3 <http://www.gnu.org/copyleft/lesser.html>`_

Any use of trade, product, or firm names is for descriptive purposes only and does not imply endorsement by the U.S. government. 

Thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info

.. _setup:

Setup
========

Automatic Installation from PyPI::


  pip uninstall PyHum (removes previous installation)
  pip install PyHum


Automatic Installation from github::


  git clone git@github.com:dbuscombe-usgs/PyHum.git
  cd PyHum
  python setup.py install


or a local installation::


  python setup.py install --user


or with admin privileges, e.g.::


  sudo python setup.py install


Setup on Anaconda for Windows
===============================
Before installing PyHum, install Basemap using::

  conda install basemap


Assuming a Anaconda distribution which comes with almost all required program dependencies::


  pip install simplekml
  pip uninstall PyHum (removes previous installation)
  pip install PyHum


pip seems to have a bug with pyproj depending on what c-compiler your python distribution uses. Therefore, you may have to install pyproj (and other dependencies) from `here <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`_

Download the .whl file, then use pip to install it, e.g.::


  pip install pyproj-1.9.4-cp27-none-win_amd64.whl


This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12.4 & 13.4 & 14.4, Windows 7. This software has (so far) been used only with Humminbird 998 and 1198 series instruments. 


.. _virtualenv:

Virtual environment
====================

You could try before you install, using a virtual environment::

  virtualenv venv
  source venv/bin/activate
  pip install numpy
  pip install cython
  pip install scipy
  pip install joblib
  pip install simplekml
  pip install pyproj
  pip install scikit-learn
  pip install Pillow
  pip install matplotlib
  pip install basemap --allow-external basemap --allow-unverified basemap
  pip install PyHum
  python -c "import PyHum; PyHum.test()"
  deactivate #(or source venv/bin/deactivate)

The results will live in "venv/lib/python2.7/site-packages/PyHum"


.. _manualinstall:

Manual installation
====================

Python libraries you need to have installed to use PyHum:

1. `Nifty <http://www.mpa-garching.mpg.de/ift/nifty/index.html>`_
2. `SciPy <http://www.scipy.org/scipylib/download.html>`_
3. `Numpy <http://www.scipy.org/scipylib/download.html>`_
4. `Matplotlib <http://matplotlib.org/downloads.html>`_
5. `cython <http://cython.org/>`_
6. `joblib <https://pythonhosted.org/joblib/>`_
7. `Scikit-learn <http://scikit-learn.org/stable/>`_
8. `Python Image LIbrary (PIL) <http://www.pythonware.com/products/pil/>`_
9. `simplekml <http://simplekml.readthedocs.org/en/latest/index.html>`_
10. `pyproj <https://pypi.python.org/pypi/pyproj>`_
11. `basemap <http://matplotlib.org/basemap/>`_

All of the above are available through `pip <https://pypi.python.org/pypi/pip>`_ and `easy_install <https://pythonhosted.org/setuptools/easy_install.html>`_


Installation on Amazon Linux EC-2 instance
============================================

It's best to install numpy, scipy, cython and matplotlib through the OS package manager::

  sudo yum install gcc gcc-c++
  sudo yum install python27-numpy python27-Cython python27-scipy python27-matplotlib

Then install geos libraries using yum and Basemap using pip::
   
  sudo yum install geos geos-devel geos-python27
  sudo pip install basemap --allow-external basemap --allow-unverified basemap

Then PyHum using pip (which will install Pillow, pyproj, simplekml, joblib and scikit-learn)::

  sudo pip install PyHum


.. _test:

Test
======

A test can be carried out by running the supplied script::

  python -c "import PyHum; PyHum.test()"

which carries out the following operations::

   # general settings   
   humfile = os.path.expanduser("~")+os.sep+'pyhum_test'+os.sep+'test.DAT' #PyHum.__path__[0]+os.sep+'test.DAT'
   sonpath = os.path.expanduser("~")+os.sep+'pyhum_test' #PyHum.__path__[0]
   doplot = 1 #yes

   # reading specific settings
   cs2cs_args = "epsg:26949" #arizona central state plane
   bedpick = 1 # auto bed pick
   c = 1450 # speed of sound fresh water
   t = 0.108 # length of transducer
   f = 455 # frequency kHz
   draft = 0.3 # draft in metres
   flip_lr = 1 # flip port and starboard

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
   #imagery = 1 # server='http://server.arcgisonline.com/ArcGIS', service='World_Imagery'
   dogrid = 1 # yes
   calc_bearing = 0 #no
   filt_bearing = 1 #yes
   res = 0.05 # grid resolution in metres
   chunk_size = 0 # auto chunk size

   PyHum.humread(humfile, sonpath, cs2cs_args, c, draft, doplot, t, f, bedpick, flip_lr, chunk_size)

   PyHum.humcorrect(humfile, sonpath, maxW, doplot)

   PyHum.humtexture(humfile, sonpath, win, shift, doplot, density, numclasses, maxscale, notes)

   PyHum.domap(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res)

   res = 0.5 # grid resolution in metres
   
   PyHum.domap_texture(humfile, sonpath, cs2cs_args, dogrid, calc_bearing, filt_bearing, res)

.. _support:

Support
=========

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please address all suggestions, comments and queries to: dbuscombe@usgs.gov. Thanks for stopping by! 

  .. image:: _static/pyhum_logo_colour_sm.png

