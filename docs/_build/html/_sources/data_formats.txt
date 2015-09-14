.. _data_formats:

***************
Data Formats
***************

.. _overview:

Overview
=========

PyHum reads data from binary files created by a Humminbird unit (using the 'Record' function). Humminbird units write out the following file formats:

1. *.DAT files which contain basic information about the sonar, time, position and sonar settings

2. *.SON files which contain the sonar data (echograms)

3. *.IDX files (1 per SON file)

One set of data in PyHum consists of up to 4 *.SON files and 1 *.DAT file. If present, the software will use the *.IDX files to more efficiently read in the echogram data from the *.SON files.
