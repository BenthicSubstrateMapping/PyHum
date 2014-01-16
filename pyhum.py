
"""
pyhum.py

INFO:
Python script to read Humminbird DAT and associated SON files, and export data

Author:  Daniel Buscombe
           Grand Canyon Monitoring and Research Center
           United States Geological Survey
           Flagstaff, AZ 86001
           dbuscombe@usgs.gov
Version: 1.1      Revision: January, 2014

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
"""

# =========================================================
# =============== import libraries ======================
# =========================================================
from __future__ import generators
import sys, getopt, pyproj, cPickle
import glob, struct, os, time
from array import array as arr
from Tkinter import Tk
from tkFileDialog import askopenfilename
from joblib import Parallel, delayed
import random as rand
import numpy as np
import matplotlib.pyplot as plt

# =========================================================
# =============== begin subfunctions ======================
# ========================================================
####################################################################
def fread(file, num, typ):
    dat = arr(typ)
    dat.fromfile(file, num)
    if typ == 'c': #character
        return(''.join(dat.tolist()))
    elif num == 1: # only 1 byte
        return(dat[0])
    else: 
        return(dat)

####################################################################
def decode_humdat(humfile, trans, transWGS84): 

    dat=[] #pre-allocate list

    fid = open(humfile,'r')
    head = fread(fid, 1, 'B')
    dat.append(fread(fid,1,'B')) # water
    dummy = fread(fid,2,'B')
    #print fid.tell()
    dat.append(struct.unpack('>i', fread(fid,4,'c'))[0]) #sonar_name
    dummy = struct.unpack('>iii', fread(fid,3*4,'c'))
    dat.append(struct.unpack('>i',fread(fid,4,'c'))[0]) # unix time
    dat.append(struct.unpack('>i',fread(fid,4,'c'))[0]) # utm x 
    dat.append(struct.unpack('>i',fread(fid,4,'c'))[0]) # utm y
    dat.append(fread(fid,10,'c')) #filename
    dummy = fread(fid,2,'B')
    dat.append(struct.unpack('>i',fread(fid,4,'c'))[0]) #numrecords
    #Total records in all SON file
    ##%If there are 4 SON files, Records per SON file
    dat.append(struct.unpack('>i',fread(fid,4,'c'))[0]) #recordlen_ms
    dat.append(struct.unpack('>i',fread(fid,4,'c'))[0]) #linesize
    dummy = fread(fid,1,'i')
    tails = fread(fid,4,'B')

    fid.close()

    # water type
    if dat[0]==0:
        dat.append('fresh')
    elif dat[0]==1:
        dat.append('deep salt')
    elif dat[0]==2:
        dat.append('shallow salt')
    else:
        dat.append('unknown')

    print "\t head: %s" % head
    print "\t water type: %s" % dat[0]
    print "\t sonar name: %s" % dat[1]
    print "\t unix time: %s" % dat[2]
    print "\t x utm: %s" % dat[3]
    print "\t y utm: %s" % dat[4]
    print "\t file name: %s" % dat[5]
    print "\t number of records: %s" % dat[6]
    print "\t record length (ms): %s" % dat[7]
    print "\t line size: %s" % dat[8]
    print "\t year: %s" % time.strftime('%Y', time.localtime(dat[2]) )
    print "\t month: %s" % time.strftime('%m', time.localtime(dat[2]) )
    print "\t day: %s" % time.strftime('%d', time.localtime(dat[2]) )
    print "\t hour: %s" % time.strftime('%H', time.localtime(dat[2]) )
    print "\t minute: %s" % time.strftime('%M', time.localtime(dat[2]) )
    print "\t second: %s" % time.strftime('%S', time.localtime(dat[2]) )
    print "\t time string: %s" % time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(dat[2]) )

    lon, lat = transWGS84(dat[3],dat[4], inverse=True)
    print "\t Longitude is: %s" % lon
    print "\t Latitude is: %s" % lat

    lon, lat = trans(lon, lat)
    dat.append(lat)
    dat.append(lon)
    print "\t Transformed coordinate x: %s" % lon
    print "\t Transformed coordinate y: %s" % lat

    return dat

####################################################################
def gethead(fid,trans, transWGS84): #cs2cs_args):
    hd = fread(fid, 3, 'B')

    head=[] #pre-allocate list

    # error catches - if headers not in right place, roll back
    if hd[0]==222:
        fid.seek(fid,-4, 1) #'cof')
        hd = fread(fid, 3, 'B')
    elif hd[0]==171:
        fid.seek(fid,-5, 1) #'cof')
        hd = fread(fid, 3, 'B')
    elif hd[0]==33:
        fid.seek(fid,-6, 1) #'cof')
        hd = fread(fid, 3, 'B')

    if hd[0]!=192 & hd[1]!=222 & hd[2]!=171:
        flag=1
    else:
        flag=0

    unit = fread(fid, 1, 'B')
    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>i', fread(fid,4,'c'))[0]) #recnum
    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>i', fread(fid,4,'c'))[0]) #time_ms
    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>i', fread(fid,4,'c'))[0]) # x_utm
    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>i', fread(fid,4,'c'))[0]) # y_utm
    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>h', fread(fid,2,'c'))[0]) # gps1
    head.append(float(struct.unpack('>h', fread(fid,2,'c'))[0])/10) # heading_deg

    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>h', fread(fid,2,'c'))[0]) # gps2
    head.append(float(struct.unpack('>h', fread(fid,2,'c'))[0])/10) # speed_ms
    spacer = fread(fid, 1, 'B')
    head.append(float(struct.unpack('>i', fread(fid,4,'c'))[0])/10) # depth_m
    spacer = fread(fid, 1, 'B')
    #%0 (50 or 83 kHz), 1 (200 kHz), 2 (SI Poort), 3 (SI Starboard)
    head.append(fread(fid, 1, 'B')) #beam
    spacer = fread(fid, 1, 'B')
    head.append(fread(fid, 1, 'B')) #voltscale
    spacer = fread(fid, 1, 'B')
    head.append(struct.unpack('>i', fread(fid,4,'c'))[0]/1000) # freq_khz
    spacer = fread(fid, 5, 'B')
    const = struct.unpack('>i', fread(fid,4,'c'))
    dummy = fread(fid, 5, 'B')
    head.append(struct.unpack('>i', fread(fid,4,'c'))[0]) #sentlen
    endhead = fread(fid, 1, 'B')

    # channel name
    if head[9]==0:
        head.append('down_lowfreq')
    elif head[9]==1:
        head.append('down_highfreq')
    elif head[9]==2:
        head.append('sidescan_port')
    elif head[9]==3:
        head.append('sidescan_starboard')
    else:
        head.append('unknown')

    lon, lat = transWGS84(head[2],head[3], inverse=True)

    head.append(lat)
    head.append(lon)

    lon, lat = trans(lon, lat)

    head.append(lat)
    head.append(lon)

    return head

####################################################################
def KnuthMorrisPratt(text, pattern):
# Knuth-Morris-Pratt string matching
# David Eppstein, UC Irvine, 1 Mar 2002
# http://code.activestate.com/recipes/117214/

    '''Yields all starting positions of copies of the pattern in the text.
Calling conventions are similar to string.find, but its arguments can be
lists or iterators, not just strings, it returns all matches, not just
the first one, and it does not need the whole text in memory at once.
Whenever it yields, it will have read the text exactly up to and including
the match that caused the yield.'''

    # allow indexing into pattern and protect against change during yield
    pattern = list(pattern)

    # build table of shift amounts
    shifts = [1] * (len(pattern) + 1)
    shift = 1
    for pos in range(len(pattern)):
        while shift <= pos and pattern[pos] != pattern[pos-shift]:
            shift += shifts[pos-shift]
        shifts[pos+1] = shift

    # do the actual search
    startPos = 0
    matchLen = 0
    for c in text:
        while matchLen == len(pattern) or \
              matchLen >= 0 and pattern[matchLen] != c:
            startPos += shifts[matchLen]
            matchLen -= shifts[matchLen]
        matchLen += 1
        if matchLen == len(pattern):
            yield startPos

####################################################################
def decode_son(sonfile, headbytes, trans, transWGS84):
    fid = open(sonfile,'r')
    dump = fread(fid,os.path.getsize(sonfile), 'c');
    fid.close()

    # unpack into integers
    int_list = [] # pre-allocate 
    for i in xrange(0, len(dump)):
        int_list.append(struct.unpack('>B', dump[i])[0])

    # find the start sequences in the list of integers
    fbreak=[]
    for s in KnuthMorrisPratt(int_list, [192,222,171,33,128]): 
        fbreak.append(s)

    # find the difference in indices of start sequences
    dfbreak = [ x-y for (x,y) in zip(fbreak[1:],fbreak[:-1]) ]

    # open the son file again and get each packet and data 
    fid = open(sonfile,'r')
    data=[]
    for i in xrange(0, len(dfbreak)):
        data.append(gethead(fid,trans,transWGS84)) # get header for packet
        int_list = []
        for j in xrange(0, dfbreak[i]-headbytes):
            int_list.append(struct.unpack('>B', fread(fid,1,'c'))[0])
        data.append(int_list) # grab the sonar data
    return data


####################################################################
def get_scans_hi(hi,packet): 
   # a function to obtain port scans 

   d = np.zeros( (packet,) )
   tmp = np.squeeze(hi) 
   d[:len(tmp)] = tmp[:packet]
   c_hi = d[:packet]
   return c_hi

####################################################################
def get_scans_low(low,packet): 
   # a function to obtain port scans 

   d = np.zeros( (packet,) )
   tmp = np.squeeze(low) 
   d[:len(tmp)] = tmp[:packet]
   c_low = d[:packet] 
   return c_low 

####################################################################
def get_scans_port(port,packet):
   # a function to obtain port scans 

   d = np.zeros( (packet,) )
   tmp = np.squeeze(port) 
   d[:len(tmp)] = tmp[:packet]
   c_port = d[:packet] 
   return c_port 

####################################################################
def get_scans_star(star,packet): 
   # a function to obtain straboard scans 

   d = np.zeros( (packet,) )
   tmp = np.squeeze(star) 
   d[:len(tmp)] = tmp[:packet]
   c_star = d[:packet] 
   return c_star

####################################################################
def get_scans_meta(k): 
   # a function to obtain other data
   # lat, lon, spd, time_s, dep_m, e, n 

   return get_meta(k,15), get_meta(k,14), get_meta(k,7), get_meta(k,1)/1000, get_meta(k,8), get_meta(k,17), get_meta(k,16) 

####################################################################
def get_meta(k,index):
   # function to return the indexth value from the kth data packet

   tmp=data_port[k-1]
   return float(tmp[0][index])

####################################################################
def custom_save(figdirec,root):
    plt.savefig(figdirec+root,bbox_inches='tight',dpi=400)

# =========================================================
# =============== begin program ======================
# ========================================================
print "####################################"
print "############### PYHUM ##############"
print "####################################"
print "########## H-B PROCESSING ##########"
print "############# READ RAW #############"
print "################ by ################"
print "########## DANIEL BUSCOMBE #########"
print "######## dbuscombe@usgs.gov ########"
print "############ G.C.M.R.C #############"
print "###### U.S. GEOLOGICAL SURVEY ######"
print "####################################"

# get list of input arguments and pre-allocate arrays
argv = sys.argv[1:]
humfile = ''; sonpath = ''
cs2cs_args = ''; numproc = ''
draft = ''; doplot = ''

# parse inputs to variables
try:
   opts, args = getopt.getopt(argv,"hi:s:c:n:d:p:")
except getopt.GetoptError:
     print 'pyhum.py -i <.DAT file> -s <sonpath (/where/*.SON/are/)> -c <cs2cs coordinate trandform arguments> -n <number of available processors> -d <draft (m)> -p <do plots, 0=no or 1=yes>'
     sys.exit(2)
for opt, arg in opts:
   if opt == '-h':
     print 'pyhum.py -i <.DAT file> -s <sonpath (/where/*.SON/are/)> -c <cs2cs coordinate trandform arguments> -n <number of available processors> -d <draft (m)> -p <do plots, 0=no or 1=yes>'
     sys.exit()
   elif opt in ("-i"):
      humfile = arg
   elif opt in ("-s"):
      sonpath = arg
   elif opt in ("-c"):
      cs2cs_args = arg
   elif opt in ("-n"):
      numproc = arg
   elif opt in ("-d"):
      draft = arg
   elif opt in ("-p"):
      doplot = arg

# prompt user to supply file if no input file given
if not humfile:
   print 'An input file is required!!!!!!'
   Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
   inputfile = askopenfilename(filetypes=[("DAT files","*.DAT")]) 

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
if cs2cs_args:
   print 'cs2cs arguments are %s' % (cs2cs_args)
if numproc:
   numproc = int(numproc)
   print 'Number of processors: %s' % (str(numproc))
if draft:
   draft = float(draft)
   print 'Draft: %s' % (str(draft))
if doplot:
   doplot = int(doplot)

if not numproc:
   numproc = 8
   print '[Default] Number of processors: %s' % (str(numproc))
if not draft:
   draft = 0
   print '[Default] Draft = %s metres' % (str(draft))
if not cs2cs_args:
   # arguments to pass to cs2cs for coordinate transforms
   cs2cs_args = "epsg:26949"
   print '[Default] cs2cs arguments are %s' % (cs2cs_args)
if not doplot:
   doplot = 1

##############################################################

# get the transformation matrix of desired output coordinates
try:
   trans =  pyproj.Proj(init=cs2cs_args)
except:
   trans =  pyproj.Proj(cs2cs_args)

cs2cs_args = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
transWGS84 = pyproj.Proj(cs2cs_args)

# number of bytes in a header packet in SON file
headbytes = 67

# if son path name supplied has no separator at end, put one on
if sonpath[-1]!=os.sep:
   sonpath = sonpath + os.sep

# get dat header from DAT file
dat = decode_humdat(humfile, trans, transWGS84) 

base = humfile.split('.DAT') # get base of file name for output
base = base[0].split('/')[-1]

# get the SON files from this directory
sonfiles = glob.glob(sonpath+'*.SON')

print "WARNING: Because files have to be read in byte by byte,"
print "this could take a very long time ..."

data = Parallel(n_jobs=4, verbose=100)(delayed(decode_son)(sonfiles[k], headbytes, trans, transWGS84) for k in range(len(sonfiles)))

# sometimes not all 4 sonars are recorded, either by accident or design, so here we have to make checks to see what data channel is what
if len(data)==4:
   # unpack into 4 data streams
   data1 = zip(data[0])
   data2 = zip(data[1])
   data3 = zip(data[2])
   data4 = zip(data[3])
   del data

   # we need to ascertain which data stream is which
   if data1[0][0][13] == 'sidescan_port':
      data_port = data1
   elif data1[0][0][13] == 'sidescan_starboard':
      data_star = data1
   elif data1[0][0][13] == 'down_lowfreq':
      data_dwnlow = data1
   elif data1[0][0][13] == 'down_highfreq':
      data_dwnhi = data1
   del data1

   if data2[0][0][13] == 'sidescan_port':
      data_port = data2
   elif data2[0][0][13] == 'sidescan_starboard':
      data_star = data2
   elif data2[0][0][13] == 'down_lowfreq':
      data_dwnlow = data2
   elif data2[0][0][13] == 'down_highfreq':
      data_dwnhi = data2
   del data2

   if data3[0][0][13] == 'sidescan_port':
      data_port = data3
   elif data3[0][0][13] == 'sidescan_starboard':
      data_star = data3
   elif data3[0][0][13] == 'down_lowfreq':
      data_dwnlow = data3
   elif data3[0][0][13] == 'down_highfreq':
      data_dwnhi = data3
   del data3

   if data4[0][0][13] == 'sidescan_port':
      data_port = data4
   elif data4[0][0][13] == 'sidescan_starboard':
      data_star = data4
   elif data4[0][0][13] == 'down_lowfreq':
      data_dwnlow = data4
   elif data4[0][0][13] == 'down_highfreq':
      data_dwnhi = data4
   del data4

elif len(data)==3:
   # unpack into 4 data streams
   data1 = zip(data[0])
   data2 = zip(data[1])
   data3 = zip(data[2])
   del data

   # we need to ascertain which data stream is which
   if data1[0][0][13] == 'sidescan_port':
      data_port = data1
   elif data1[0][0][13] == 'sidescan_starboard':
      data_star = data1
   elif data1[0][0][13] == 'down_lowfreq':
      data_dwnlow = data1
   elif data1[0][0][13] == 'down_highfreq':
      data_dwnhi = data1
   del data1

   if data2[0][0][13] == 'sidescan_port':
      data_port = data2
   elif data2[0][0][13] == 'sidescan_starboard':
      data_star = data2
   elif data2[0][0][13] == 'down_lowfreq':
      data_dwnlow = data2
   elif data2[0][0][13] == 'down_highfreq':
      data_dwnhi = data2
   del data2

   if data3[0][0][13] == 'sidescan_port':
      data_port = data3
   elif data3[0][0][13] == 'sidescan_starboard':
      data_star = data3
   elif data3[0][0][13] == 'down_lowfreq':
      data_dwnlow = data3
   elif data3[0][0][13] == 'down_highfreq':
      data_dwnhi = data3
   del data3

if 'data_dwnlow' in locals():
   if 'data_dwnhi' not in locals():
      # pickle raw data to file
      flag = 1
      with open(sonpath+base+'.pkl', 'wb') as fid:
         cPickle.dump([dat, data_port, data_star, data_dwnlow], fid) 
   else:
      flag = 2
      with open(sonpath+base+'.pkl', 'wb') as fid:
         cPickle.dump([dat, data_port, data_star, data_dwnlow, data_dwnhi], fid) 
else:
   if 'data_dwnhi' in locals():
      flag = 3
      # pickle raw data to file
      with open(sonpath+base+'.pkl', 'wb') as fid:
         cPickle.dump([dat, data_port, data_star, data_dwnhi], fid) 
   else:
      flag = 4
      with open(sonpath+base+'.pkl', 'wb') as fid:
         cPickle.dump([dat, data_port, data_star], fid) 

# here we just want the number of data values in the first packet. 
# The rest of the record will be cropped to this
packet = data_port[0]
packet = packet[0][12]

# make an index of every other record
ind = range(0,len(data_port))
ind = ind[1::2]

start_time = np.asarray(dat[2],'float')

verbosity=0
print "Dealing out data into variables ..."

# delete variables to save memory
del data_star
if 'data_dwnlow' in locals():
   del data_dwnlow
if 'data_dwnhi' in locals():
   del data_dwnhi

# dole out the raw data into port and starboard scans
# first port side
c_port = Parallel(n_jobs = numproc, verbose=verbosity)(delayed(get_scans_port)(data_port[i][0],packet) for i in ind)
c_port = np.squeeze(c_port).T
# pickle raw data
with open(sonpath+base+'raw_port.pkl', 'wb') as fid:
   cPickle.dump([c_port], fid) 

del data_port

# make a plot if requested 
if doplot==1:
   fig = plt.figure()
   plt.imshow(c_port,cmap='gray', origin = 'upper'); plt.colorbar(); 
   plt.title('Portside Raw Scan'); plt.xlabel('Ping Number (Time)'); plt.ylabel('Range (Distance)')
   custom_save(sonpath,'raw_port')
   del fig
   plt.close()

del c_port

# we're only reading in the data we need for the next segment to save memory
# get only the starboard data
if flag==1:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow = cPickle.load(f)
      del dat, data_port, data_dwnlow
elif flag==2:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow, data_dwnhi = cPickle.load(f)
      del dat, data_port, data_dwnlow, data_dwnhi
elif flag==3:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnhi = cPickle.load(f)
      del dat, data_port, data_dwnhi
elif flag==4:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star = cPickle.load(f)
      del dat, data_port

# now starboard side
c_star = Parallel(n_jobs = numproc, verbose=verbosity)(delayed(get_scans_star)(data_star[i][0],packet) for i in ind)
c_star = np.squeeze(c_star).T
with open(sonpath+base+'raw_star.pkl', 'wb') as fid:
   cPickle.dump([c_star], fid)

del data_star

# make a plot if requested 
if doplot==1:
   fig = plt.figure()
   plt.imshow(c_star,cmap='gray', origin = 'upper'); plt.colorbar(); 
   plt.title('Starboardside Raw Scan'); plt.xlabel('Ping Number (Time)'); plt.ylabel('Range (Distance)')
   custom_save(sonpath,'raw_star')
   del fig
   plt.close()

del c_star

# we're only reading in the data we need for the next segment to save memory
# get only the port data
if flag==1:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow = cPickle.load(f)
      del dat, data_star, data_dwnlow
elif flag==2:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow, data_dwnhi = cPickle.load(f)
      del dat, data_star, data_dwnlow, data_dwnhi
elif flag==3:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnhi = cPickle.load(f)
      del dat, data_star, data_dwnhi
elif flag==4:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star = cPickle.load(f)
      del dat, data_star

# now we deal with the meta data (coordinates, speeds, and depths)
d = Parallel(n_jobs = numproc, verbose=verbosity)(delayed(get_scans_meta)(i) for i in ind)
del ind
lon, lat, spd, time_s, dep_m, e, n = zip(*d)

del data_port

# make a plot if requested 
if doplot==1:
   # plot longitude vs latitude
   fig = plt.figure()
   plt.plot(lon,lat,'ko')
   plt.title('Boat Course'); plt.xlabel('Longitude'); plt.ylabel('Latitude')
   custom_save(sonpath,'raw_trace')
   del fig
   plt.close()

   # plot depth
   fig = plt.figure()
   plt.plot(dep_m,'k')
   plt.title('Depth'); plt.xlabel('Ping Number (Time)'); plt.ylabel('Depth (m)')
   custom_save(sonpath,'raw_depth')
   del fig
   plt.close()

del d

caltime = np.asarray(start_time + time_s,'float')
dep_m = np.asarray(dep_m,'float')+draft

with open(sonpath+base+'meta.pkl', 'wb') as fid:
   cPickle.dump([np.asarray(lat,'float'),np.asarray(lon,'float'),np.asarray(spd,'float'),np.asarray(time_s,'float'),np.asarray(e,'float'),np.asarray(n,'float'),dep_m,caltime], fid) 

# add meta data to port scan raw and rewrite file
try:
   with open(sonpath+base+'raw_port.pkl') as f:
      c_port = cPickle.load(f)
   with open(sonpath+base+'raw_port.pkl', 'wb') as fid:
      cPickle.dump([c_port,np.asarray(lat,'float'),np.asarray(lon,'float'),np.asarray(spd,'float'),np.asarray(time_s,'float'),np.asarray(e,'float'),np.asarray(n,'float'),dep_m,caltime], fid) 
   del c_port
except:
   print "meta data not written to %s" % (sonpath+base+'raw_port.pkl')

try:
   # add meta data to star scan raw and rewrite file
   with open(sonpath+base+'raw_star.pkl') as f:
      c_star = cPickle.load(f)
   with open(sonpath+base+'raw_star.pkl', 'wb') as fid:
      cPickle.dump([c_star,np.asarray(lat,'float'),np.asarray(lon,'float'),np.asarray(spd,'float'),np.asarray(time_s,'float'),np.asarray(e,'float'),np.asarray(n,'float'),dep_m,caltime], fid) 
   del c_star
except:
   print "meta data not written to %s" % (sonpath+base+'raw_star.pkl')


# we're only reading in the data we need for the next segment to save memory
# get only the low freq. data
if flag==1:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow = cPickle.load(f)
      del dat, data_port, data_star
elif flag==2:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow, data_dwnhi = cPickle.load(f)
      del dat, data_port, data_star, data_dwnhi

# if there is a low frequency downward looking sonar recorded, we parse that data too
if 'data_dwnlow' in locals():
   # make an index of every other record
   ind = range(0,len(data_dwnlow))
   ind = ind[1::2]
   c_low = Parallel(n_jobs = numproc, verbose=verbosity)(delayed(get_scans_low)(data_dwnlow[i][0],packet) for i in ind)
   del data_dwnlow
   c_low = np.squeeze(c_low).T

   with open(sonpath+base+'raw_low.pkl', 'wb') as fid:
      cPickle.dump([c_low], fid)

   # make a plot if requested 
   if doplot==1:
      fig = plt.figure()
      plt.imshow(c_low,cmap='gray', origin = 'upper'); plt.colorbar(); 
      plt.title('Low Freq. Downward Raw Sonar'); plt.xlabel('Ping Number (Time)'); plt.ylabel('Range (Distance)')
      custom_save(sonpath,'raw_dwnlow')
      del fig
      plt.close()

   del c_low

   with open(sonpath+base+'raw_low.pkl') as f:
      c_low = cPickle.load(f)
   with open(sonpath+base+'raw_low.pkl', 'wb') as fid:
      cPickle.dump([c_low,np.asarray(lat,'float'),np.asarray(lon,'float'),np.asarray(spd,'float'),np.asarray(time_s,'float'),np.asarray(e,'float'),np.asarray(n,'float'),dep_m,caltime], fid) 
   del c_low

# we're only reading in the data we need for the next segment to save memory
# get only the high freq. data
if flag==2:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnlow, data_dwnhi = cPickle.load(f)
      del dat, data_star, data_port, data_dwnlow
elif flag==3:
   with open(sonpath+base+'.pkl') as f:
      dat, data_port, data_star, data_dwnhi = cPickle.load(f)
      del dat, data_star, data_port

# if there is a high frequency downward looking sonar recorded, we parse that data too
if 'data_dwnhi' in locals():
   # make an index of every other record
   ind = range(0,len(data_dwnhi))
   ind = ind[1::2]
   c_hi = Parallel(n_jobs = numproc, verbose=verbosity)(delayed(get_scans_hi)(data_dwnhi[i][0],packet) for i in ind)
   del data_dwnhi
   c_hi = np.squeeze(c_hi).T

   with open(sonpath+base+'raw_hi.pkl', 'wb') as fid:
      cPickle.dump([c_hi], fid)

   # make a plot if requested 
   if doplot==1:
      fig = plt.figure()
      plt.imshow(c_hi,cmap='gray', origin = 'upper'); plt.colorbar(); 
      plt.title('High Freq. Downward Raw Sonar'); plt.xlabel('Ping Number (Time)'); plt.ylabel('Range (Distance)')
      custom_save(sonpath,'raw_dwnhi')
      del fig
      plt.close()

   del c_hi

   # add meta data to high sonar raw
   with open(sonpath+base+'raw_hi.pkl') as f:
      c_hi = cPickle.load(f)
   with open(sonpath+base+'raw_hi.pkl', 'wb') as fid:
      cPickle.dump([c_hi,np.asarray(lat,'float'),np.asarray(lon,'float'),np.asarray(spd,'float'),np.asarray(time_s,'float'),np.asarray(e,'float'),np.asarray(n,'float'),dep_m,caltime], fid) 
   del c_hi

print "Done!"


#with open(base[0]+'meta.pkl') as f:
#   lat, lon, spd, time_s, e, n, dep_m, caltime = cPickle.load(f)
# add meta data to low sonar raw

## save data to mat file
#savemat(base[0]+'.mat', mdict={'head':dat, 'data_port': data_port, 'data_star': data_star, 'data_dwnlow': data_dwnlow, 'data_dwnhi': data_dwnhi})
#savemat(outbase+'raw_low.mat', mdict={'c_low':c_low})
#savemat(outbase+'raw_hi.mat', mdict={'c_hi':c_hi})
#savemat(outbase+'raw.mat', mdict={'lat':lat, 'lon':lon, 'spd':spd, 'time_s':time_s, 'e':e, 'n':n, 'dep_m':dep_m, 'caltime':caltime}) 
#savemat(outbase+'raw_star.mat', mdict={'c_star':c_star})
#savemat(outbase+'raw_port.mat', mdict={'c_port':c_port})

