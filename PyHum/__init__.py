# encoding: utf-8
"""
PyHum - a Python framework for sidescan data texture classification.

PyHum is an open-source project dedicated to provide a Python framework for
processing low-cost sidescan data. It provides parsers for Humminbird file formats,
and signal processing routines which allow the manipulation of sidescan data and automated texture classification (see Buscombe et al., forthcoming).

For more information visit http://dbuscombe-usgs.github.io/PyHum/

:license:
    GNU Lesser General Public License, Version 3
    (http://www.gnu.org/copyleft/lesser.html)
    
    This software is in the public domain because it contains materials that
    originally came from the United States Geological Survey, an agency of the
    United States Department of Interior. For more information, 
    see the official USGS copyright policy at
    http://www.usgs.gov/visual-id/credit_usgs.html#copyright
    Any use of trade, product, or firm names is for descriptive purposes only 
    and does not imply endorsement by the U.S. government.
    
"""

__version__ = '1.2.3'

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

from PyHum._pyhum_read_class import humread
from PyHum._pyhum_correct_class import humcorrect
from PyHum._pyhum_texture_class import humtexture
from PyHum._pyhum_map import domap
from PyHum._pyhum_map_texture import domap_texture
from PyHum.utils import *
from PyHum.test import *


