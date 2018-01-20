.. _getting_started:


***************
Getting started
***************

.. _about:

About
======

PyHum - a Python framework for reading and processing data from a Humminbird low-cost sidescan sonar

PyHum is an open-source project dedicated to provide a generic Python framework 
for reading and exporting data from Humminbird(R) instruments, carrying out rudimentary radiometric corrections to the data,
classify bed texture, and produce some maps on aerial photos and kml files for google-earth

The software is designed to read Humminbird data (.SON, .IDX, and .DAT files) and works on both sidescan and downward-looking echosounder data, where available.

The following Humminbird units are currently supported:
1. 700, 800, 900 and 1100 series
2. ONIX series
    
.. _info:

Technical Info
==============
Full documentation of the procedures behind the program is in the following publication:

Buscombe, D., 2017, Shallow water benthic imaging and substrate characterization using recreational-grade sidescan-sonar. ENVIRONMENTAL MODELLING & SOFTWARE 89, 1-18.

.. _cite:

Please Cite!
============

If you use PyHum in your published work, please cite the following papers:

1. Buscombe, D., Grams, P.E., and Smith, S. (2015) "Automated riverbed sediment classification using low-cost sidescan sonar", Journal of Hydraulic Engineering, 10.1061/(ASCE)HY.1943-7900.0001079, 06015019.

2. Buscombe, D., 2017, Shallow water benthic imaging and substrate characterization using recreational-grade sidescan-sonar. ENVIRONMENTAL MODELLING & SOFTWARE 89, 1-18.


For the source code visit `the project github site <https://github.com/dbuscombe-usgs/PyHum/>`_

.. _credits:

Credits
========

 Primary Developer |    Daniel Buscombe 
 ------ | ---------------
         |  Northern Arizona University
          | Flagstaff, AZ 86001
          | daniel.buscombe@nau.edu

 Co-Developer |    Daniel Hamill
 ------ | ---------------
         |  Department of Watershed Sciences
          | Utah State University
          | Logan, UT 84322
          | dhamill@usgs.gov

Version: 1.4.3    |  Revision: Jan, 2018


.. _notes:

Notes
========

1. This software has been tested with Python 2.7 on Linux Fedora 16 & 20, Ubuntu 12 -- 17, Windows 7. Python 3 is not yet supported

2. This software has (so far) been used only with Humminbird 798, 898, 998, 1198, 1199 and ONIX series instruments. 

3. PyHum is not yet compatable with newer HELIX/MEGA systems


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

Thanks to Barb Fagetter (blueseas@oceanecology.ca) for some format info and Paul Anderson (Quest Geophysical Asia) for debugging and suggestions for improvements

.. _setup:


Installation
=============

Installing in a conda virtual env (recommended)

In a conda (miniconda/anaconda) python 2 environment: 

Linux::

  conda create --name pyhum python=2
  source activate pyhum
  conda install -c conda-forge basemap-data-hires -y
  conda install scipy numpy scikit-image
  pip install simplekml sklearn pandas dask
  pip install joblib toolz cython
  pip install pyresample
  pip install PyHum --no-deps #(or pip install git+https://github.com/dbuscombe-usgs/PyHum.git --no-deps)
  python -c"import PyHum;PyHum.dotest()" 

Windows::

  conda create --name pyhum python=2
  activate pyhum
  conda install -c conda-forge basemap-data-hires -y
  conda install scipy numpy scikit-image
  pip install simplekml sklearn pandas dask
  pip install joblib toolz cython
  pip install pyresample==1.1.4
  pip install PyHum --no-deps #(or pip install git+https://github.com/dbuscombe-usgs/PyHum.git --no-deps)
  python -c"import PyHum;PyHum.dotest()" 


Installing as a library accessible outside of virtual env

From PyPI::

  pip install PyHum

The latest 'bleeding edge' (pre-release) version directly from github::

  pip install git+https://github.com/dbuscombe-usgs/PyHum.git

(Windows users) install git from here: https://git-scm.com/download/win

From github repo clone::

  git clone git@github.com:dbuscombe-usgs/PyHum.git
  cd PyHum
  python setup.py install

A local installation::

  python setup.py install --user


If you get C++ compiler errors (such as "Unable to find vcvarsall.bat"), you will need to install the Microsoft Visual C++ compiler from `here <http://aka.ms/vcpython27>`_


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

  python -c "import PyHum; PyHum.test.dotest()"


.. _gui:

Using the GUI
==============

From the command line (terminal)::


   python -c "import PyHum; PyHum.gui()"

.. _gettingstarted:

Getting Started
================

Inputs to the program are a .DAT file (e.g. R0089.DAT) and a folder of .SON and .IDX files (e.g. /my/folder/R0089). The program will read the .SON files with or without the accompanying .IDX files, but will be faster if the .IDX files are present. 

PyHum is modular so can be called from within a python or ipython console, from an IDE (such as IDLE or Spyder), or by running a script.

The following example script::
 
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
    PyHum.map(humfile, sonpath, cs2cs_args, res, mode, nn, numstdevs, use_uncorrected) #dowrite, 


could be saved as, for example "proc_mysidescandata.py" and run from the command line using::


   python proc_mysidescandata.py -i C:\MyData\R0087.DAT -s C:\MyData\R0087


or from within ipython (with a GUI prompt to navigate to the files)::

   %run proc_mysidescandata.py
   
If you are in bash (or git bash) you might want to automate through a folder of subfolders like this::

   for k in $(find $PWD -type d -maxdepth 1 -mindepth 1); do python proc_mysidescandata.py -i "$k/${k##*/}.DAT" -s $k; done

which assumes the .DAT file is in the folder with the same root (such as a folder called R00123 which contains SON and IDX files as well as a file called R00123.DAT)


.. _support:

Support
=========

This is a new project written and maintained by Daniel Buscombe. Bugs are expected - please report them, I will fix them quickly. Feedback and suggestions for improvements are *very* welcome

Please download, try, report bugs, fork, modify, evaluate, discuss, collaborate. Please use the 'Issues' tab in github `here <https://github.com/dbuscombe-usgs/PyHum>`_

Thanks for stopping by! 


.. _troubleshooting:

Trouble Shooting
================

1. Problem: pyhum read hangs for a long time (several minutes) on the test script. 
Try this: uninstall joblib and install an older version::

   pip uninstall joblib
   pip install joblib==0.7.1

2. Problem: you get an "invalid mode or file name" error.
Try this: construct file paths using raw strings e.g.:: 

   r'C:\Users\me\mydata\R0089' 

or using os, e.g.::

   import os
   os.path.abspath(os.path.join('C:\Users','me','mydata','R0089'))

3. Problem: on Linux, PyHum is using an older version of scipy than 0.16, as revealed by::

   python -c 'import scipy;print(scipy.__version__)'

Try this: remove a system installed file e.g.::

   sudo apt-get remove python-scipy ##(Debian based)
   yum remove scipy ##(Fedora based)

4. Problem: do I have the latest version of PyHum installed? Check your version using this::

   python -c 'import PyHum;print(PyHum.__version__)'

Check this against the latest `bleeding-edge' version `here <https://github.com/dbuscombe-usgs/PyHum/blob/master/PyHum/__init__.py>`_ (line 47)


.. image:: _static/pyhum_logo_colour_sm.png

