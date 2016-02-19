## PyHum (Python program for Humminbird(R) data processing) 
## has been developed at the Grand Canyon Monitoring & Research Center,
## U.S. Geological Survey
##
## Author: Daniel Buscombe
## Project homepage: <https://github.com/dbuscombe-usgs/PyHum>
##
##This software is in the public domain because it contains materials that originally came from 
##the United States Geological Survey, an agency of the United States Department of Interior. 
##For more information, see the official USGS copyright policy at 
##http://www.usgs.gov/visual-id/credit_usgs.html#copyright
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see <http://www.gnu.org/licenses/>.

#"""
# ____        _   _                         
#|  _ \ _   _| | | |_   _ _ __ ___    _   _ 
#| |_) | | | | |_| | | | | '_ ` _ \  (_) (_)
#|  __/| |_| |  _  | |_| | | | | | |  _   _ 
#|_|    \__, |_| |_|\__,_|_| |_| |_| (_) (_)
#       |___/                               
#                   _ 
#   ____ _ __  __  (_)
#  / __ `// / / / / / 
# / /_/ // /_/ / / /  
# \__, (_)__,_(_)_/   
#/____/                         
#
##+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#|b|y| |D|a|n|i|e|l| |B|u|s|c|o|m|b|e|
#+-+-+ +-+-+-+-+-+-+ +-+-+-+-+-+-+-+-+
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#|d|b|u|s|c|o|m|b|e|@|u|s|g|s|.|g|o|v|
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+
#|U|.|S|.| |G|e|o|l|o|g|i|c|a|l| |S|u|r|v|e|y|
#+-+-+-+-+ +-+-+-+-+-+-+-+-+-+-+ +-+-+-+-+-+-+

#"""
#=======================
# python -c "import PyHum; PyHum.gui()"

import Tkinter
from Tix import *

import ttk
from tkFileDialog import askopenfilename
import os
import PyHum
from PIL import Image, ImageTk

#import webbrowser
import tkMessageBox
from ScrolledText import ScrolledText

import os

#################################################
def gui():

	#=======================
	# NOTE: Frame will make a top-level window if one doesn't already exist which
	# can then be accessed via the frame's master attribute
	# make a Frame whose parent is root, named "pyhum"
	master = Tkinter.Frame(name='pyhum')

	self = master.master  # short-cut to top-level window
	master.pack()  # pack the Frame into root, defaults to side=TOP
	self.title('PyHum GUI')  # name the window

	# field for DAT filename
	self.DATfilename = Tkinter.StringVar()
	self.DATfilename.set(u" ")

	# field for sed filename
	self.sedfilename = Tkinter.StringVar()
	self.sedfilename.set(None)

	# some defaults 
	self.model = 998

	self.calc_heading = 0
	self.filt_heading = 0
	self.bedpick = 1

	self.doplot = 1
	self.cog = 1
	self.chunk = 'd100'

	self.mode = 1
	self.dogrid = 1
	self.integ = 5       
	       
	# create notebook
	demoPanel = Tkinter.Frame(master, name='demo')  # create a new frame slaved to master
	demoPanel.pack()  # pack the Frame into root

	# create (notebook) demo panel
	nb = ttk.Notebook(demoPanel, name='notebook')  # create the ttk.Notebook widget

	# extend bindings to top level window allowing
	#   CTRL+TAB - cycles thru tabs
	#   SHIFT+CTRL+TAB - previous tab
	#   ALT+K - select tab using mnemonic (K = underlined letter)
	nb.enable_traversal()

	nb.pack(fill=Tkinter.BOTH, expand=Tkinter.Y, padx=2, pady=3)  # add margin

	#==============================================================
	#==============================================================
	#========START about tab

	# create description tab
	# frame to hold (tab) content
	frame = Tkinter.Frame(nb, name='descrip')

	frame.configure(background='black')

	# widgets to be displayed on 'Description' tab
	About_msg = [
	    "PyHum is a program to read and process data collected ",
	    "by a Humminbird fishfinder / sidescan sonar\n",
	    "\n",
	    "The program is written and maintained by Daniel Buscombe\n",
	    "U.S. Geological Survey. Email dbuscombe@usgs.gov\n",
	    "\n",    
	    "Aspects of this program are documented in the following papers:\n",
	    "\n",
	    "1. Buscombe, D., Grams, P.E., Smith, S. (2015) ",
	    "Automated riverbed sediment classification using low-cost sidescan sonar. ",
	    "JOURNAL OF HYDRAULIC ENGINEERING, 06015019\n",
	    "\n",
	    "2. Buscombe, D., (2016, submitted) ",
	    "Processing and georeferencing recreational-grade sidescan-sonar data to support the democratization of acoustic imaging in shallow water",
	    "Intended for LIMNOLOGY AND OCEANOGRAPHY: METHODS\n",
	    "\n",
	    "The tabs are to be navigated in order (read, filter, etc)\n",
	    "\n",
	    "Please visit the website for more info: http://dbuscombe-usgs.github.io/PyHum/"]

	lbl = Tkinter.Label(frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(About_msg), fg="white")
	lbl.configure(background='black')
		        
	#neatVar = Tkinter.StringVar()
	#btn = Tkinter.Button(frame, text='Neat!', underline=0,
	#                 command=lambda v=neatVar: _say_neat(master, v))
	#neat = Tkinter.Label(frame, textvariable=neatVar, name='neat')

	# position and set resize behavior
	lbl.grid(row=0, column=0, columnspan=2, sticky='new', pady=5)
	#btn.grid(row=1, column=0, pady=(2,4))
	#neat.grid(row=1, column=1,  pady=(2,4))
	frame.rowconfigure(1, weight=1)
	frame.columnconfigure((0,1), weight=1, uniform=1)

	# make panel for logo
	image_panel = Tkinter.Canvas(frame, width = 250, height =232)#, cursor = "cross")
	image_panel.grid(column = 0, row = 2)
	image_panel.configure(background='black')
	 
	show_image = ImageTk.PhotoImage(Image.open(PyHum.__path__[0]+os.sep+"pyhum_logo_white_sm.png"))
	# show the panel
	image_panel.create_image(0, 0, anchor=Tkinter.NW, image=show_image) 
		
	# add to notebook (underline = index for short-cut character)
	nb.add(frame, text='About', underline=0, padding=2)

	#==============================================================
	#==============================================================
	#========END about tab

	#==============================================================
	#==============================================================
	#========START read tab

	# Populate the second pane. Note that the content doesn't really matter
	read_frame = Tkinter.Frame(nb)
	nb.add(read_frame, text='Read')#, state='disabled')

	read_frame.configure(background='MediumPurple1')

	Read_msg = [
	    "Read a .DAT and associated set of .SON files recorded by a Humminbird(R) instrument.\n\n",
	    "Parse the data into a set of memory mapped files that will",
	    "subsequently be used by the other functions of the PyHum module.\n\n" ,  
	    "Export time-series data and metadata in other formats.\n\n",    
	    "Create a kml file for visualising boat track be selected.\n\n"
	    "Create rudimentary plots of the data"]

	lbl2 = Tkinter.Label(read_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(Read_msg))

	lbl2.configure(background='thistle3', fg="black")
		        
	# position and set resize behavior
	lbl2.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello1():
	   tkMessageBox.showinfo("Read Data Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n Model: 3 or 4 number code indicating model of sidescan unit\n\n Bed Pick: 1 = auto, 0 = manual, 3=auto with manual override\n\n Sound speed: typically, 1450 m/s in freshwater, 1500 in saltwater\n\n Transducer length: m\n\n Draft: m\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n Chunk: partition data to chunks. 'd' - distance, m. 'p' - number of pings. 'h' - change in heading, degs. '1' - just 1 chunk\n\n Flip Port/Star: flip port and starboard sidescans\n\n Filter bearing: spike removing filter to bearing\n\n Calculate bearing: recaluclate bearing from positions")

	def hello1_alt():
	   try:
	      root = Tk()
	      root.wm_title("Read module")
	      S = Scrollbar(root)
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Read Data Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n Model: 3 or 4 number code indicating model of sidescan unit\n\n Bed Pick: 1 = auto, 0 = manual, 3=auto with manual override\n\n Sound speed: typically, 1450 m/s in freshwater, 1500 in saltwater\n\n Transducer length: m\n\n Draft: m\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n Chunk: partition data to chunks. 'd' - distance, m. 'p' - number of pings. 'h' - change in heading, degs. '1' - just 1 chunk\n\n Flip Port/Star: flip port and starboard sidescans\n\n Filter bearing: spike removing filter to bearing\n\n Calculate bearing: recalculate bearing from positions"""
	      T.insert(END, quote)
	   except:
	      hello1()   

	MSG1_btn = Tkinter.Button(read_frame, text = "Instructions", command = hello1_alt)
	MSG1_btn.grid(row=0, column=1, pady=(2,4))
	MSG1_btn.configure(background='thistle3', fg="black")

	read_frame.rowconfigure(1, weight=1)
	read_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.read_dat_btn = Tkinter.Button(read_frame, text='Get DAT file', underline=0,
		         command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(read_frame, textvariable=datVar, name='dat')
	self.read_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.read_dat_btn.configure(background='thistle3', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.read_son_btn = Tkinter.Button(read_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(read_frame, textvariable=sonVar, name='dat')
	self.read_son_btn.grid(row=1, column=1, pady=(2,4))
	self.read_son_btn.configure(background='thistle3', fg="black")

	#=======================
	#menu for model
	self.mb=  Tkinter.Menubutton ( read_frame, text="Humminbird model", relief=Tkinter.RAISED)
	self.mb.grid(column = 0, row = 2, pady=(2,4))
	self.mb.menu  =  Tkinter.Menu ( self.mb, tearoff = 0 , background='PaleVioletRed1', fg="black" )
	      
	self.mb.menu.add_command(label="997", command = lambda v=1: _SetModel(master, v))
	self.mb.menu.add_command(label="998", command = lambda v=2: _SetModel(master, v))
	self.mb.menu.add_command(label="1198", command = lambda v=3: _SetModel(master, v))
	self.mb.menu.add_command(label="1199", command = lambda v=4: _SetModel(master, v))                
	self.mb["menu"]  =  self.mb.menu
	self.mb.configure(background='thistle3', fg="black")

	#=======================
	#menu for bedpick
	self.bb=  Tkinter.Menubutton ( read_frame, text="Bed pick type", relief=Tkinter.RAISED)
	self.bb.grid(column = 1, row = 2, pady=(2,4))
	self.bb.menu  =  Tkinter.Menu ( self.bb, tearoff = 0 , background='PaleVioletRed1', fg="black" )
		
	self.bb.menu.add_command(label="Auto", command = lambda v=1: _SetBedPick(master, v))
	self.bb.menu.add_command(label="Manual", command = lambda v=2: _SetBedPick(master, v))
	self.bb.menu.add_command(label="Hybrid", command = lambda v=3: _SetBedPick(master, v))
	self.bb["menu"]  =  self.bb.menu
	self.bb.configure(background='thistle3', fg="black")

	#=======================
	# sound speed
	self.cvar = Tkinter.DoubleVar()
	cscale = Tkinter.Scale( read_frame, variable = self.cvar, from_=1440, to=1510, resolution=5, tickinterval=10, label = 'Sound Velocity [m/s]' )
	cscale.set(1450)
	cscale.grid(row=3, column=0,  pady=(2,4))
	cscale.configure(background='thistle3', fg="black")

	##=======================
	## frequency
	#self.fvar = Tkinter.DoubleVar()
	#fscale = Tkinter.Scale( read_frame, variable = self.fvar, from_=400, to=800, resolution=5, tickinterval=100, label = 'SS Frequency [kHz]' )
	#fscale.set(455)
	#fscale.grid(row=3, column=1,  pady=(2,4))
	#fscale.configure(background='thistle3', fg="black")

	#=======================
	# transducer length
	self.tvar = Tkinter.DoubleVar()
	tscale = Tkinter.Scale( read_frame, variable = self.tvar, from_=.1, to=.5, resolution=.001, tickinterval=.1, label = 'Transducer Length [m]' )
	tscale.set(0.108)
	tscale.grid(row=4, column=0,  pady=(2,4))
	tscale.configure(background='thistle3', fg="black")

	#=======================
	# draft
	self.dvar = Tkinter.DoubleVar()
	dscale = Tkinter.Scale( read_frame, variable = self.dvar, from_=0, to=1, resolution=.01, tickinterval=.2, label = 'Draft [m]' )
	dscale.set(0.3)
	dscale.grid(row=4, column=1,  pady=(2,4))
	dscale.configure(background='thistle3', fg="black")

	#=======================        
	# epsg
	self.epsg = Tkinter.StringVar()
	self.epsg_entry = Tkinter.Entry(read_frame, width = 30, textvariable = self.epsg)
	self.epsg_entry.grid(column = 0, row = 5, pady=(2,4)) #, columnspan = 2,sticky = 'EW')
	self.epsg_entry.bind("<Return>", lambda epsg=self.epsg.get(): _OnPressEnter1(self))
	self.epsg.set(u"epsg:26949")       
	self.epsg_entry.configure(background='thistle3', fg="black")

	#=======================
	#menu for chunk
	self.cb=  Tkinter.Menubutton ( read_frame, text="chunk argument", relief=Tkinter.RAISED)
	self.cb.grid(column = 1, row = 5, pady=(2,4))
	self.cb.menu  =  Tkinter.Menu ( self.cb, tearoff = 0, background='PaleVioletRed1', fg="black" )
	      
	self.cb.menu.add_command(label="d50", command = lambda v=1: _SetChunk(master, v))
	self.cb.menu.add_command(label="d100", command = lambda v=2: _SetChunk(master, v))
	self.cb.menu.add_command(label="d200", command = lambda v=3: _SetChunk(master, v))
	self.cb.menu.add_command(label="d400", command = lambda v=4: _SetChunk(master, v)) 

	self.cb.menu.add_command(label="p1000", command = lambda v=5: _SetChunk(master, v))
	self.cb.menu.add_command(label="p5000", command = lambda v=6: _SetChunk(master, v))
	self.cb.menu.add_command(label="p10000", command = lambda v=7: _SetChunk(master, v)) 
	self.cb.menu.add_command(label="p20000", command = lambda v=8: _SetChunk(master, v)) 
	   
	self.cb.menu.add_command(label="h5", command = lambda v=9: _SetChunk(master, v))
	self.cb.menu.add_command(label="h10", command = lambda v=10: _SetChunk(master, v))
	self.cb.menu.add_command(label="h20", command = lambda v=11: _SetChunk(master, v))
	self.cb.menu.add_command(label="h40", command = lambda v=12: _SetChunk(master, v)) 

	self.cb["menu"]  =  self.cb.menu
	self.cb.configure(background='thistle3', fg="black")

	#=======================        
	# check button for flip_lr        
	self.flipvar = Tkinter.IntVar()      
	flip_entry = Tkinter.Checkbutton(read_frame, text='Flip Port/Star', variable=self.flipvar)
	flip_entry.config(indicatoron=1, bd=4, width=12) 
	flip_entry.grid(column = 0, row = 6, pady=(2,4))
	flip_entry.configure(background='thistle3', fg="black")

	#=======================
	# check button for calc_heading        
	self.calcheadvar = Tkinter.IntVar()      
	calchead_entry = Tkinter.Checkbutton(read_frame, text='Calculate Heading', variable=self.calcheadvar)
	calchead_entry.config(indicatoron=1, bd=4, width=12) 
	calchead_entry.grid(column = 1, row = 6, pady=(2,4))
	calchead_entry.configure(background='thistle3', fg="black")

	#=======================
	# check button for filt_heading        
	self.filtheadvar = Tkinter.IntVar()      
	filthead_entry = Tkinter.Checkbutton(read_frame, text='Filter Heading', variable=self.filtheadvar)
	filthead_entry.config(indicatoron=1, bd=4, width=12) 
	filthead_entry.grid(column = 0, row = 7, pady=(2,4))
	filthead_entry.configure(background='thistle3', fg="black")

	#=======================
	# process button
	proc_btn = Tkinter.Button(read_frame, text='Process!', underline=0,
		         command=lambda filt_heading=self.filtheadvar.get(): _proc(self))
	proc_btn.grid(row=7, column=1, pady=(2,4))
	proc_btn.configure(background='thistle3', fg="black")

	##=======================
	## bind for button short-cut key # (must be bound to toplevel window)
	#master.winfo_toplevel().bind('<Alt-n>', lambda e, v=datVar: self._get_DAT(v))

	#==============================================================
	#========START functions for read tab

	#=======================
	# must press enter to set density
	def _OnPressEnter1(self):
	    """
	    sets epsg code on Enter press
	    """
	    self.epsg.set( self.epsg.get() )
	    self.epsg_entry.focus_set()
	    self.epsg_entry.selection_range(0, Tkinter.END)
	    print 'epsg code set to %s ' % (str(self.epsg.get()))
	#    if int(self.c.get())>1500:
	#       tkMessageBox.showinfo("High value", "Consider 1450 for freshwater and 1500 for salt water")

	#=======================
	def _proc(self):
	    # function to invoke PyHum.read
	    
	    # build error checking into here
	    print 'Processing ...'
	#    print "filt_head: " + str(self.filtheadvar.get())
	#    print "calc_head: " + str(self.calcheadvar.get())    
	#    print "humfile: " + self.DATfilename.get()
	#    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	#    print "model: " + str(self.model)
	#    print "bedpick: " + str(self.bedpick)
	#    print "flip_lr: " + str(self.flipvar.get())
	#    print "c: " + str(self.cvar.get())    
	#    print "f: " + str(self.fvar.get())    
	#    print "t: " + str(self.tvar.get())    
	#    print "draft: " + str(self.dvar.get())                
	#    print "cs2cs arguments: " + str(self.epsg.get())                
	#    print "chunk argument: " + str(self.chunk)                
	    # do stuff here
	    PyHum.read(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), str(self.epsg.get()), self.cvar.get(), self.dvar.get(), self.doplot, self.tvar.get(), self.bedpick, self.flipvar.get(), self.model, self.calcheadvar.get(), self.filtheadvar.get(), self.cog, self.chunk)
	    self.update()
	    tkMessageBox.showinfo("Done!", "Read module finished") 

	#=======================        
	def _get_DAT(master, v):
	    self.DATfile = askopenfilename(filetypes=[("DAT files","*.DAT")], multiple=False)

	    self.DATfilename.set(self.DATfile)

	    print 'You chose: %s' % (self.DATfilename.get())
	    
	    self.read_dat_btn.configure(fg='thistle3', background="black")
	    self.correct_dat_btn.configure(fg='lightgoldenrod', background="black")     
	    self.map_dat_btn.configure(fg='IndianRed1', background="black")
	    self.rmshadows_dat_btn.configure(fg='LightBlue2', background="black")
	    self.texture_dat_btn.configure(fg='SeaGreen1', background="black")   
	    self.maptexture_dat_btn.configure(fg='gold2', background="black")
	    self.e1e2_dat_btn.configure(fg='SlateGray1', background="black")
		    
	    self.update()    

	#=======================
	def _get_SON(master, v):
	    self.SONfiles = askopenfilename(filetypes=[("SON files","*.SON")], multiple=True)

	    for k in xrange(len(self.SONfiles)):
	       print 'You chose: %s' % (self.SONfiles[k])
	    self.folder = os.path.dirname(self.SONfiles[0])
	    
	    #self.son_btn.configure(fg='thistle3', background="black")

	    self.read_son_btn.configure(fg='thistle3', background="black")
	    self.correct_son_btn.configure(fg='lightgoldenrod', background="black")     
	    self.map_son_btn.configure(fg='IndianRed1', background="black")
	    self.rmshadows_son_btn.configure(fg='LightBlue2', background="black")
	    self.texture_son_btn.configure(fg='SeaGreen1', background="black")   
	    self.maptexture_son_btn.configure(fg='gold2', background="black")
	    self.e1e2_son_btn.configure(fg='SlateGray1', background="black")
		    
	    self.update()    

	#=======================
	def _SetModel(master, v):
	   if v==1:
	      self.model=997
	      print "model is 997"
	   elif v==2:
	      self.model=998
	      print "model is 998"
	   elif v==3:
	      self.model=1198
	      print "model is 1198"
	   elif v==4:
	      self.model=1199
	      print "model is 1199"
	      
	   self.mb.configure(fg='thistle3', background="black")
	      
	   self.update() 
	   
	#=======================   
	def _SetBedPick(master, v):
	   if v==1:
	      self.bedpick=1
	      print "bed picking is auto"
	   elif v==2:
	      self.bedpick=2
	      print "bed picking is manual"
	   elif v==3:
	      self.bedpick=3
	      print "bed picking is hybrid"
	      
	   self.bb.configure(fg='thistle3', background="black")
		 
	   self.update()          

	#=======================   
	def _SetChunk(master, v):
	   if v==1:
	      self.chunk='d50'; print "chunk is 50 m"
	   elif v==2:
	      self.chunk='d100'; print "chunk is 100 m"
	   elif v==3:
	      self.chunk='d200'; print "chunk is 200 m"
	   elif v==4:
	      self.chunk='d400'; print "chunk is 400 m"
	   elif v==5:
	      self.chunk='p1000'; print "chunk is 1000 pings"
	   elif v==6:
	      self.chunk='p5000'; print "chunk is 5000 pings"
	   elif v==7:
	      self.chunk='p10000'; print "chunk is 10,000 pings"
	   elif v==8:
	      self.chunk='p20000'; print "chunk is 20,000 pings"            
	   elif v==9:
	      self.chunk='h5'; print "chunk is a 5 degree heading change" 
	   elif v==10:
	      self.chunk='h10'; print "chunk is a 10 degree heading change" 
	   elif v==11:
	      self.chunk='h20'; print "chunk is a 20 degree heading change" 
	   elif v==12:
	      self.chunk='h40'; print "chunk is a 40 degree heading change"
	      
	   self.cb.configure(fg='thistle3', background="black")
	       
	   self.update()  
	   
	#==============================================================
	#========END functions for read tab

	#==============================================================
	#==============================================================
	#========END read tab

	#==============================================================
	#==============================================================
	#========START filter tab

	filter_frame = Tkinter.Frame(nb)
	nb.add(filter_frame, text='Filter', state='disabled')

	#==============================================================
	#==============================================================
	#========END filter tab

	#==============================================================
	#==============================================================
	#========START correct tab

	correct_frame = Tkinter.Frame(nb)
	nb.add(correct_frame, text='Correct')#, state='disabled')

	correct_frame.configure(background='GoldenRod1')

	Correct_msg = [
	    "Remove water column and carry out some rudimentary radiometric corrections,\n",
	    "accounting for directivity and attenuation of sound by water.\n\n"
	    "Radiometric correction carried out using a simple Cosine model (Lambert's Law).\n\n"
	    "Data converted from arbitrary to acoustic units (dB W).\n\n"
	    "Create some rudimentary plots of the data.\n\n"
	    "Optionally, also correct for attenuation of sound by sediment.\n\n"
	    "Optionally, carry out a phase-preserving denoising filter on the data"]
	    
	lbl3 = Tkinter.Label(correct_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(Correct_msg))

	lbl3.configure(background='lightgoldenrod', fg="black")
		        
	# position and set resize behavior
	lbl3.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello2():
	   tkMessageBox.showinfo("Scan Correction Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n RMS Output: maximum transducer power in W\n\n pH: water acidity\n\n Temp: water temperature in deg C\n\n Salinity: water salinity in ppt\n\n Do Filter: check to apply phase-preserving filter to scans (warning: can take a very long time)\n\n Retain water column: check if you don't want to remove water column from scans\n\n Get sediment file: select text file containing sediment concentration data (must contain the following fields separated by spaces: size (microns) conc (mg/L) dens (kg/m), with one row per grain size)")

	def hello2_alt():
	   try:
	      root = Tk()
	      root.wm_title("Correct module")      
	      S = Scrollbar(root)
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Scan Correction Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n RMS Output: maximum transducer power in W\n\n pH: water acidity\n\n Temp: water temperature in deg C\n\n Salinity: water salinity in ppt\n\n Do Filter: check to apply phase-preserving filter to scans (warning: can take a very long time)\n\n Retain water column: check if you don't want to remove water column from scans\n\n Get sediment file: select text file containing sediment concentration data (must contain the following fields separated by spaces: size (microns) conc (mg/L) dens (kg/m), with one row per grain size)"""
	      T.insert(END, quote)
	   except:
	      hello2()
	      
	MSG2_btn = Tkinter.Button(correct_frame, text = "Instructions", command = hello2_alt)
	MSG2_btn.grid(row=0, column=1, pady=(2,4))
	MSG2_btn.configure(background='lightgoldenrod', fg="black")

	correct_frame.rowconfigure(1, weight=1)
	correct_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.correct_dat_btn = Tkinter.Button(correct_frame, text='Get DAT file', underline=0,
		         command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(correct_frame, textvariable=datVar, name='dat')
	self.correct_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.correct_dat_btn.configure(background='lightgoldenrod', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.correct_son_btn = Tkinter.Button(correct_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(correct_frame, textvariable=sonVar, name='dat')
	self.correct_son_btn.grid(row=1, column=1, pady=(2,4))
	self.correct_son_btn.configure(background='lightgoldenrod', fg="black")

	#=======================
	# rms output wattage
	self.Wvar = Tkinter.DoubleVar()
	Wscale = Tkinter.Scale( correct_frame, variable = self.Wvar, from_=500, to=1500, resolution=100, tickinterval=200, label = 'RMS Output [W]' )
	Wscale.set(1000)
	Wscale.grid(row=2, column=0,  pady=(2,4))
	Wscale.configure(background='lightgoldenrod', fg="black")

	#=======================
	# pH
	self.Pvar = Tkinter.DoubleVar()
	Pscale = Tkinter.Scale( correct_frame, variable = self.Pvar, from_=4, to=8, resolution=.1, tickinterval=1, label = 'pH' )
	Pscale.set(7)
	Pscale.grid(row=2, column=1,  pady=(2,4))
	Pscale.configure(background='lightgoldenrod', fg="black")

	#=======================
	# temp
	self.Tvar = Tkinter.DoubleVar()
	Tscale = Tkinter.Scale( correct_frame, variable = self.Tvar, from_=3, to=25, resolution=1, tickinterval=5, label = 'Temp [degC]' )
	Tscale.set(10)
	Tscale.grid(row=3, column=0,  pady=(2,4))
	Tscale.configure(background='lightgoldenrod', fg="black")

	#=======================
	# salinity
	self.Svar = Tkinter.DoubleVar()
	Sscale = Tkinter.Scale( correct_frame, variable = self.Svar, from_=0, to=35, resolution=1, tickinterval=5, label = 'Salinity [ppt]' )
	Sscale.set(0)
	Sscale.grid(row=3, column=1,  pady=(2,4))
	Sscale.configure(background='lightgoldenrod', fg="black")

	#=======================
	# check button for dofilt        
	self.dofiltvar = Tkinter.IntVar()      
	dofilt_entry = Tkinter.Checkbutton(correct_frame, text='Do Filter', variable=self.dofiltvar)
	dofilt_entry.config(indicatoron=1, bd=4, width=12) 
	dofilt_entry.grid(column = 0, row = 4, pady=(2,4))
	dofilt_entry.configure(background='lightgoldenrod', fg="black")

	#=======================
	# check button for correct_withwater       
	self.watervar = Tkinter.IntVar()      
	water_entry = Tkinter.Checkbutton(correct_frame, text='Retain Water Column', variable=self.watervar)
	water_entry.config(indicatoron=1, bd=4, width=20) 
	water_entry.grid(column = 1, row = 4, pady=(2,4))
	water_entry.configure(background='lightgoldenrod', fg="black")

	#=======================
	# get sediment file           
	sedVar = Tkinter.StringVar()
	self.correct_sed_btn = Tkinter.Button(correct_frame, text='Get sediment file', underline=0,
		         command=lambda v=sedVar: _get_sedfile(master, v))
	sed = Tkinter.Label(correct_frame, textvariable=sedVar, name='sed')
	self.correct_sed_btn.grid(row=5, column=0, pady=(2,4))
	self.correct_sed_btn.configure(background='lightgoldenrod', fg="black")

	#=======================
	# process button
	proc2_btn = Tkinter.Button(correct_frame, text='Process!', underline=0,
		         command=lambda watervar=self.watervar.get(): _proc2(self))
	proc2_btn.grid(row=5, column=1, pady=(2,4))
	proc2_btn.configure(background='lightgoldenrod', fg="black")

	#==============================================================
	#========START functions for correct tab

	#=======================
	def _proc2(self):
	    # function to invoke PyHum.correct

	    #dofilt = 0 # 1 = apply a phase preserving filter (WARNING!! takes a very long time for large scans)    
	    # build error checking into here
	    print 'Processing ...'   
	#    print "humfile: " + self.DATfilename.get()
	#    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	#    print "max wattage: " + str(self.Wvar.get())        
	#    print "pH: " + str(self.Pvar.get())    
	#    print "temperature: " + str(self.Tvar.get())    
	#    print "salinity: " + str(self.Svar.get())    
	#    print "dofilt: " + str(self.dofiltvar.get())
	#    print "correct_withwater: " + str(self.watervar.get())                 
	    # do stuff here
	    PyHum.correct(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), self.Wvar.get(), self.doplot, self.dofiltvar.get(), self.watervar.get(), self.Pvar.get(), self.Tvar.get(), self.Svar.get(), self.sedfilename.get())    
	    self.update() 
	    tkMessageBox.showinfo("Done!", "Correct module finished") 
	    
	#=======================        
	def _get_sedfile(master, v):
	    self.sedfile = askopenfilename(filetypes=[("Sediment file","*.txt")], multiple=False)

	    self.sedfilename.set(self.sedfile)

	    print 'You chose: %s' % (self.sedfilename.get())
	    
	    self.correct_sed_btn.configure(fg='lightgoldenrod', background="black")     
		    
	    self.update()  
	#==============================================================
	#========END functions for correct tab

	#==============================================================
	#==============================================================
	#========END correct tab

	#==============================================================
	#==============================================================
	#========START rmshadows tab

	rmshadows_frame = Tkinter.Frame(nb)
	nb.add(rmshadows_frame, text='Remove Shadows') #, state='disabled')

	rmshadows_frame.configure(background='SteelBlue1')

	rmshadows_msg = [
	    "Remove dark shadows in scans caused by shallows, shorelines, and attenuation of acoustics with distance\n\n",
	    "Manual or automated processing options available\n\n",
	    "Works on the radiometrically corrected outputs of the correct module\n\n",
	    "If manual, the user will be prompted to delineate a path with the cursor,\n", 
	    "on each starbaord and port scan chunk in turn, beyond which the data will be removed\n\n",
	    "If automatic, portions of scans will be asessed for shadows (and subsequently removed)\n",
	    "using an approach based on grey-level co-occurrence matrices"]
	    

	lbl4 = Tkinter.Label(rmshadows_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(rmshadows_msg))

	lbl4.configure(background='LightBlue2', fg="black")
		        
	# position and set resize behavior
	lbl4.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello3():
	   tkMessageBox.showinfo("Shadow Removal Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n win: window size (pixels) for the automated shadow removal algorithm. Larger windows means less well resolved blocking out of shadows\n\n Manual mask: check for manual masking of echograms, otherwise the shaow removal will be carried out automatically")

	def hello3_alt():
	   try:
	      root = Tk()
	      root.wm_title("Shadow removal module")      
	      S = Scrollbar(root)
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Shadow Removal Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n win: window size (pixels) for the automated shadow removal algorithm. Larger windows means less well resolved blocking out of shadows\n\n Manual mask: check for manual masking of echograms, otherwise the shaow removal will be carried out automatically"""
	      T.insert(END, quote)
	   except:
	      hello3()
		  
	MSG3_btn = Tkinter.Button(rmshadows_frame, text = "Instructions", command = hello3_alt)
	MSG3_btn.grid(row=0, column=1, pady=(2,4))
	MSG3_btn.configure(background='LightBlue2', fg="black")

	rmshadows_frame.rowconfigure(1, weight=1)
	rmshadows_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.rmshadows_dat_btn = Tkinter.Button(rmshadows_frame, text='Get DAT file', underline=0,command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(rmshadows_frame, textvariable=datVar, name='dat')
	self.rmshadows_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.rmshadows_dat_btn.configure(background='LightBlue2', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.rmshadows_son_btn = Tkinter.Button(rmshadows_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(rmshadows_frame, textvariable=sonVar, name='dat')
	self.rmshadows_son_btn.grid(row=1, column=1, pady=(2,4))
	self.rmshadows_son_btn.configure(background='LightBlue2', fg="black")

	#=======================
	# window size
	self.Winvar = Tkinter.DoubleVar()
	Winscale = Tkinter.Scale( rmshadows_frame, variable = self.Winvar, from_=5, to=500, resolution=1, tickinterval=50, label = 'Window Size [pixels]' )
	Winscale.set(31)
	Winscale.grid(row=2, column=0,  pady=(2,4))
	Winscale.configure(background='LightBlue2', fg="black")

	#=======================
	# check button for manual mask        
	self.manmaskvar = Tkinter.IntVar()      
	manmask_entry = Tkinter.Checkbutton(rmshadows_frame, text='Manual Mask', variable=self.manmaskvar)
	manmask_entry.config(indicatoron=1, bd=4, width=12) 
	manmask_entry.grid(column = 1, row = 2, pady=(2,4))
	manmask_entry.configure(background='LightBlue2', fg="black")

	#=======================
	# process button
	proc3_btn = Tkinter.Button(rmshadows_frame, text='Process!', underline=0,
		         command=lambda watervar=self.watervar.get(): _proc3(self))
	proc3_btn.grid(row=5, column=1, pady=(2,4))
	proc3_btn.configure(background='LightBlue2', fg="black")


	#==============================================================
	#========START functions for rmshadows tab

	#=======================
	def _proc3(self):
	    # function to invoke PyHum.rmshadows

	    # build error checking into here
	    print 'Processing ...'   
	#    print "humfile: " + self.DATfilename.get()
	#    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	#    print "window size: " + str(self.Winvar.get())        
	#    print "manually mask: " + str(self.manmaskvar.get())                    
	    # do stuff here
	    PyHum.rmshadows(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), self.Winvar.get(), self.manmaskvar.get(), self.doplot)  
	    self.update() 
	    tkMessageBox.showinfo("Done!", "Shadow removal module finished")     

	#==============================================================
	#========END functions for rmshadows tab

	#==============================================================
	#==============================================================
	#========END rmshadows tab

	#==============================================================
	#==============================================================
	#========START map tab

	map_frame = Tkinter.Frame(nb)
	nb.add(map_frame, text='Map Sidescan')#, state='disabled')

	map_frame.configure(background='firebrick3')

	map_msg = [
	    "Load the sidescan chunks and create point clouds of spatially referenced sidescan intensity\n",
	    "By positioning and georectifying each sidescan pixel\n\n",
	    "Create a regular grid (raster) the point cloud data according to user-specified parameters\n\n"
	    "Create simple kml images and plots of the spatially referenced sidescan echograms"]
	    
	lbl5 = Tkinter.Label(map_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(map_msg))

	lbl5.configure(background='IndianRed1', fg="black")
		        
	# position and set resize behavior
	lbl5.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello4():
	   tkMessageBox.showinfo("Map Sidescan Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n res: grid resolution of the output gridded sidescan map. If 0, an appropriate grid resolution will be automatically found\n\n gridding mode: 1 = simple nearest neighbour; 2 = nearest neighbour, weighted using inverse distance metric; 3 = Gaussian-weighted nearest neighbour\n\n Nearest neighbours: the number of nearst neighbours to use in the gridding (larger number means more spatial smoothing of the data)\n\n Number of standard deviationss: threshold number of standard deviations in sidescan intensity per grid cell up to which to accept\n\n")    

	def hello4_alt():
	   try:
	      root = Tk()
	      root.wm_title("Map module")      
	      S = Scrollbar(root)
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Map Sidescan Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n res: grid resolution of the output gridded sidescan map. If 0, an appropriate grid resolution will be automatically found\n\n gridding mode: 1 = simple nearest neighbour; 2 = nearest neighbour, weighted using inverse distance metric; 3 = Gaussian-weighted nearest neighbour\n\n Nearest neighbours: the number of nearst neighbours to use in the gridding (larger number means more spatial smoothing of the data)\n\n Number of standard deviationss: threshold number of standard deviations in sidescan intensity per grid cell up to which to accept\n\n"""
	      T.insert(END, quote)
	   except:
	      hello4()    

	MSG4_btn = Tkinter.Button(map_frame, text = "Instructions", command = hello4_alt)
	MSG4_btn.grid(row=0, column=1, pady=(2,4))
	MSG4_btn.configure(background='IndianRed1', fg="black")

	map_frame.rowconfigure(1, weight=1)
	map_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.map_dat_btn = Tkinter.Button(map_frame, text='Get DAT file', underline=0,
		         command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(map_frame, textvariable=datVar, name='dat')
	self.map_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.map_dat_btn.configure(background='IndianRed1', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.map_son_btn = Tkinter.Button(map_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(map_frame, textvariable=sonVar, name='dat')
	self.map_son_btn.grid(row=1, column=1, pady=(2,4))
	self.map_son_btn.configure(background='IndianRed1', fg="black")

	#=======================
	#menu for mode
	self.Mb=  Tkinter.Menubutton ( map_frame, text="Gridding mode", relief=Tkinter.RAISED)
	self.Mb.grid(column = 0, row = 2, pady=(2,4))
	self.Mb.menu  =  Tkinter.Menu ( self.Mb, tearoff = 0 , background='orchid2', fg="black" )
		
	self.Mb.menu.add_command(label="NN", command = lambda v=1: _SetMode(master, v))
	self.Mb.menu.add_command(label="IDW", command = lambda v=2: _SetMode(master, v))
	self.Mb.menu.add_command(label="Gaussian", command = lambda v=3: _SetMode(master, v))
	self.Mb["menu"]  =  self.Mb.menu
	self.Mb.configure(background='IndianRed1', fg="black")

	#=======================        
	# epsg
	self.epsg = Tkinter.StringVar()
	self.epsg_entry = Tkinter.Entry(map_frame, width = 30, textvariable = self.epsg)
	self.epsg_entry.grid(column = 1, row = 2, pady=(2,4)) #, columnspan = 2,sticky = 'EW')
	self.epsg_entry.bind("<Return>", lambda epsg=self.epsg.get(): _OnPressEnter1(self))
	self.epsg.set(u"epsg:26949")       
	self.epsg_entry.configure(background='IndianRed1', fg="black")

	#=======================
	# res
	self.resvar = Tkinter.DoubleVar()
	resscale = Tkinter.Scale( map_frame, variable = self.resvar, from_=0, to=2, resolution=0.01, tickinterval=0.25, label = 'Resolution [m] 0 = auto' )
	resscale.set(0)
	resscale.grid(row=3, column=0,  pady=(2,4))
	resscale.configure(background='IndianRed1', fg="black")

	#=======================
	# nn
	self.nnvar = Tkinter.DoubleVar()
	nnscale = Tkinter.Scale( map_frame, variable = self.nnvar, from_=1, to=512, resolution=1, tickinterval=50, label = 'Nearest neighbours' )
	nnscale.set(64)
	nnscale.grid(row=3, column=1,  pady=(2,4))
	nnscale.configure(background='IndianRed1', fg="black")

	#=======================
	## influence
	#self.infvar = Tkinter.DoubleVar()
	#infscale = Tkinter.Scale( map_frame, variable = self.infvar, from_=0.1, to=10, resolution=0.1, tickinterval=1, label = 'Influence [m]' )
	#infscale.set(1)
	#infscale.grid(row=4, column=0,  pady=(2,4))
	#infscale.configure(background='IndianRed1', fg="black")

	#=======================
	# numstdevs
	self.nstdvar = Tkinter.DoubleVar()
	nstdscale = Tkinter.Scale( map_frame, variable = self.nstdvar, from_=1, to=10, resolution=1, tickinterval=2, label = 'Numb. stan. dev.' )
	nstdscale.set(5)
	nstdscale.grid(row=4, column=1,  pady=(2,4))
	nstdscale.configure(background='IndianRed1', fg="black")

	#=======================
	# check button for dowrite       
	self.dowritevar = Tkinter.IntVar()      
	dowrite_entry = Tkinter.Checkbutton(map_frame, text='Write Point Cloud', variable=self.dowritevar)
	dowrite_entry.config(indicatoron=1, bd=4, width=12) 
	dowrite_entry.grid(column = 0, row = 5, pady=(2,4))
	dowrite_entry.configure(background='IndianRed1', fg="black")

	#=======================
	# process button
	proc4_btn = Tkinter.Button(map_frame, text='Process!', underline=0,
		         command=lambda watervar=self.watervar.get(): _proc4(self))
	proc4_btn.grid(row=5, column=1, pady=(2,4))
	proc4_btn.configure(background='IndianRed1', fg="black")


	#==============================================================
	#========START functions for map tab

	#=======================
	# must press enter to set density
	def _OnPressEnter1(self):
	    """
	    sets epsg code on Enter press
	    """
	    self.epsg.set( self.epsg.get() )
	    self.epsg_entry.focus_set()
	    self.epsg_entry.selection_range(0, Tkinter.END)
	    print 'epsg code set to %s ' % (str(self.epsg.get()))
	#    if int(self.c.get())>1500:
	#       tkMessageBox.showinfo("High value", "Consider 1450 for freshwater and 1500 for salt water")


	#=======================
	def _proc4(self):
	    # function to invoke PyHum.map
	    # build error checking into here
	    print 'Processing ...'   
	#    print "humfile: " + self.DATfilename.get()
	#    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	#    print "mode: " + str(self.mode)        
	#    print "cs2cs arguments: " + str(self.epsg.get())   
	#    print "resolution: " + str(self.resvar.get())
	#    print "max. number of nearest neighbours: " + str(self.nnvar.get()) 
	#    print "gridding influence [m]: " + str(self.infvar.get()) 
	#    print "number std. dev. to accept: " + str(self.nstdvar.get())             
	#    print "Write point cloud to file: " + str(self.dowritevar.get())             
		   
	    if self.resvar.get()==0:
	       self.resvar=99                 
	    # do stuff here
	    PyHum.map(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), self.epsg.get(), self.resvar.get(), self.dowritevar.get(), self.mode, self.nnvar.get(), self.nstdvar.get()) #self.infvar.get(),
	    self.update() 
	    tkMessageBox.showinfo("Done!", "Map module finished") 
	    
	#=======================   
	def _SetMode(master, v):
	   if v==1:
	      self.mode=1
	      print "grid mode is nearest neighbour"
	   elif v==2:
	      self.mode=2
	      print "grid mode is inverse weight distance nearest neighbour"
	   elif v==3:
	      self.mode=3
	      print "grid mode is gaussian"
	      
	   self.Mb.configure(fg='IndianRed1', background="black")
	   self.Mbt.configure(fg='gold2', background="black")
	      
	   self.update() 
	   
	#==============================================================
	#========END functions for map tab         

	#==============================================================
	#==============================================================
	#========END map tab

	#==============================================================
	#==============================================================
	#========START texture tab

	texture_frame = Tkinter.Frame(nb)
	nb.add(texture_frame, text='Texture')#, state='disabled')

	texture_frame.configure(background='SpringGreen3')
	 
	texture_msg = [
	    "Read radiometrically corrected Humminbird data (output from the correct module),\n",
	    "perform a textural analysis using the spectral method of Buscombe et al (2015).\n\n",
	    "Data are split into windows of specified size and with specified overlap.\n"
	    "The data in each window are analysed spectrally to compute a texture lengthscale per window.\n"
	    "These per-window texture lengthscales are re-mapped according to the centroids,\n"
	    "with interpolation of values between centroids where appropriate, to the same extent as the original scans.\n\n"
	    "Plots of the texture lengthscale maps, as well as k-mean segmentations of the texture lengthscale maps, are created"
	]

	#    "Note that this textural lengthscale is is unit length (m) but is not a direct measure of grain size.\n",
	#    "Rather, it is a statistical representation that integrates over many attributes of bed texture,\n",
	#    "of which grain size is the most important. The technique is a physically based means to identify\n", 
	#    "regions of texture within a sidescan echogram, and could provide a basis for objective, automated riverbed sediment classification.\n\n",

	lbl6 = Tkinter.Label(texture_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,text=''.join(texture_msg))

	lbl6.configure(background='SeaGreen1', fg="black")
		        
	# position and set resize behavior
	lbl6.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello5():
	   tkMessageBox.showinfo("Texture Calculation Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n Window size: size of square window in pixels\n\n shift size: length of pixel shift (in 2 dimensions) between windows, which dictates the degree of overlap (no overlap would be shift size equalling window size)\n\n density: Tells the program how much of the image to use. If sample density is 10, then every 10 lines of the echogram will be used. If set to 2, then every other row will be used. If win=100 and density=50, 2 slices would be used to calculate texture lengthscale. If density=25, 4 slices would be used.\n\n Number of classes: number of classes to classify the texture lengthscale map into, using k-means segmentation \n\n Max. scale: maximum scale as an inverse fraction of the window size (e.g. if max scale is 10, and window size is 100, the maximum texture lengthscale would be 10 pixels)\n\n notes: notes (scale increments) per octave to consider in continuous wavelet transform. The number of texture lengthscales considered is notes * octaves, where octaves = log2( win/max scale/2 )")    

	def hello5_alt():
	   try:
	      root = Tk()
	      S = Scrollbar(root)
	      root.wm_title("Texture module")      
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Texture Calculation Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n Window size: size of square window in pixels\n\n shift size: length of pixel shift (in 2 dimensions) between windows, which dictates the degree of overlap (no overlap would be shift size equalling window size)\n\n density: Tells the program how much of the image to use. If sample density is 10, then every 10 lines of the echogram will be used. If set to 2, then every other row will be used. If win=100 and density=50, 2 slices would be used to calculate texture lengthscale. If density=25, 4 slices would be used.\n\n Number of classes: number of classes to classify the texture lengthscale map into, using k-means segmentation \n\n Max. scale: maximum scale as an inverse fraction of the window size (e.g. if max scale is 10, and window size is 100, the maximum texture lengthscale would be 10 pixels)\n\n notes: notes (scale increments) per octave to consider in continuous wavelet transform. The number of texture lengthscales considered is notes * octaves, where octaves = log2( win/max scale/2 )"""
	      T.insert(END, quote)
	   except:
	      hello5()
	      
	MSG5_btn = Tkinter.Button(texture_frame, text = "Instructions", command = hello5_alt)
	MSG5_btn.grid(row=0, column=1, pady=(2,4))
	MSG5_btn.configure(background='SeaGreen1', fg="black")

	texture_frame.rowconfigure(1, weight=1)
	texture_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.texture_dat_btn = Tkinter.Button(texture_frame, text='Get DAT file', underline=0,
		         command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(texture_frame, textvariable=datVar, name='dat')
	self.texture_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.texture_dat_btn.configure(background='SeaGreen1', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.texture_son_btn = Tkinter.Button(texture_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(texture_frame, textvariable=sonVar, name='dat')
	self.texture_son_btn.grid(row=1, column=1, pady=(2,4))
	self.texture_son_btn.configure(background='SeaGreen1', fg="black")

	#=======================
	# window size
	self.Winvar = Tkinter.DoubleVar()
	Winscale = Tkinter.Scale( texture_frame, variable = self.Winvar, from_=5, to=500, resolution=1, tickinterval=50, label = 'Window Size [pixels]' )
	Winscale.set(100)
	Winscale.grid(row=2, column=0,  pady=(2,4))
	Winscale.configure(background='SeaGreen1', fg="black")

	#=======================
	# shift size
	self.shiftvar = Tkinter.DoubleVar()
	shiftscale = Tkinter.Scale( texture_frame, variable = self.shiftvar, from_=2, to=50, resolution=1, tickinterval=10, label = 'Shift Size [pixels]' )
	shiftscale.set(10)
	shiftscale.grid(row=2, column=1,  pady=(2,4))
	shiftscale.configure(background='SeaGreen1', fg="black")

	#=======================
	# density
	self.densvar = Tkinter.DoubleVar()
	densscale = Tkinter.Scale( texture_frame, variable = self.densvar, from_=1, to=100, resolution=1, tickinterval=10, label = 'Density' )
	densscale.set(10)
	densscale.grid(row=3, column=0,  pady=(2,4))
	densscale.configure(background='SeaGreen1', fg="black")

	#=======================
	# maxscale
	self.maxscalevar = Tkinter.DoubleVar()
	maxscalescale = Tkinter.Scale( texture_frame, variable = self.maxscalevar, from_=2, to=50, resolution=1, tickinterval=10, label = 'Max scale' )
	maxscalescale.set(20)
	maxscalescale.grid(row=3, column=1,  pady=(2,4))
	maxscalescale.configure(background='SeaGreen1', fg="black")

	#=======================
	# notes
	self.notesvar = Tkinter.DoubleVar()
	notesscale = Tkinter.Scale( texture_frame, variable = self.notesvar, from_=2, to=20, resolution=1, tickinterval=5, label = 'Notes' )
	notesscale.set(4)
	notesscale.grid(row=4, column=0,  pady=(2,4))
	notesscale.configure(background='SeaGreen1', fg="black")

	#=======================
	# numclasses
	self.ncvar = Tkinter.DoubleVar()
	ncscale = Tkinter.Scale( texture_frame, variable = self.ncvar, from_=2, to=10, resolution=1, tickinterval=2, label = 'Number classes' )
	ncscale.set(4)
	ncscale.grid(row=4, column=1,  pady=(2,4))
	ncscale.configure(background='SeaGreen1', fg="black")

	#=======================
	# process button
	proc5_btn = Tkinter.Button(texture_frame, text='Process!', underline=0,
		         command=lambda watervar=self.watervar.get(): _proc5(self))
	proc5_btn.grid(row=5, column=1, pady=(2,4))
	proc5_btn.configure(background='SeaGreen1', fg="black")

	#==============================================================
	#========START functions for texture tab         

	#=======================
	def _proc5(self):
	    # function to invoke PyHum.texture
	    # build error checking into here
	    print 'Processing ...'   
	#    print "humfile: " + self.DATfilename.get()
	#    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	#    print "window size [pixels]: " + str(self.Winvar.get())        
	#    print "shift size [pixels]: " + str(self.shiftvar.get())        
	#    print "density [pixels]: " + str(self.densvar.get())
	#    print "maxscale: " + str(self.maxscalevar.get())
	#    print "notes: " + str(self.notesvar.get())
	#    print "number of texture classes: " + str(self.ncvar.get())
		             
	    # do stuff here
	    PyHum.texture(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), self.Winvar.get(), self.shiftvar.get(), self.doplot, self.densvar.get(), self.ncvar.get(), self.maxscalevar.get(), self.notesvar.get())
	    self.update() 
	    tkMessageBox.showinfo("Done!", "Texture module finished") 
	    
	#==============================================================
	#========END functions for texture tab         

	#==============================================================
	#==============================================================
	#========END texture tab

	#==============================================================
	#==============================================================
	#========START map texture tab

	map_texture_frame = Tkinter.Frame(nb)
	nb.add(map_texture_frame, text='Map Texture')#, state='disabled')

	map_texture_frame.configure(background='dark orange')
	 
	map_texture_msg = [
	    "Load the sidescan chunks and create point clouds of spatially referenced texture lengthscale values\n",
	    "By positioning and georectifying each texture lengthscale value at a specified resolution\n\n",
	    "Save the point cloud to ascci format for use by this and other programs\n\n"
	    "Create a regular grid (raster) the point cloud data according to user-specified parameters\n\n"
	    "Create simple kml images and plots of the spatially referenced texture lengthscales computed from sidescan echograms"]
	    
	lbl7 = Tkinter.Label(map_texture_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(map_texture_msg))

	lbl7.configure(background='gold2', fg="black")
		        
	# position and set resize behavior
	lbl7.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello6():
	   tkMessageBox.showinfo("Map Texture Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n res: grid resolution of the output gridded sidescan map. If 0, an appropriate grid resolution will be automatically found\n\n gridding mode: 1 = simple nearest neighbour; 2 = nearest neighbour, weighted using inverse distance metric; 3 = Gaussian-weighted nearest neighbour\n\n Nearest neighbours: the number of nearst neighbours to use in the gridding (larger number means more spatial smoothing of the data)\n\n Number of standard deviationss: threshold number of standard deviations in texture lengthscale per grid cell up to which to accept\n\n")    

	def hello6_alt():
	   try:
	      root = Tk()
	      S = Scrollbar(root)
	      root.wm_title("Map texture module")      
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Map Texture Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n res: grid resolution of the output gridded sidescan map. If 0, an appropriate grid resolution will be automatically found\n\n gridding mode: 1 = simple nearest neighbour; 2 = nearest neighbour, weighted using inverse distance metric; 3 = Gaussian-weighted nearest neighbour\n\n Nearest neighbours: the number of nearst neighbours to use in the gridding (larger number means more spatial smoothing of the data)\n\n Number of standard deviationss: threshold number of standard deviations in texture lengthscale per grid cell up to which to accept"""
	      T.insert(END, quote)
	   except:
	      hello6()
	    
	MSG6_btn = Tkinter.Button(map_texture_frame, text = "Instructions", command = hello6_alt)
	MSG6_btn.grid(row=0, column=1, pady=(2,4))
	MSG6_btn.configure(background='gold2', fg="black")

	map_texture_frame.rowconfigure(1, weight=1)
	map_texture_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.maptexture_dat_btn = Tkinter.Button(map_texture_frame, text='Get DAT file', underline=0,
		         command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(map_texture_frame, textvariable=datVar, name='dat')
	self.maptexture_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.maptexture_dat_btn.configure(background='gold2', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.maptexture_son_btn = Tkinter.Button(map_texture_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(map_texture_frame, textvariable=sonVar, name='dat')
	self.maptexture_son_btn.grid(row=1, column=1, pady=(2,4))
	self.maptexture_son_btn.configure(background='gold2', fg="black")


	#=======================
	#menu for mode
	self.Mbt=  Tkinter.Menubutton ( map_texture_frame, text="Gridding mode", relief=Tkinter.RAISED)
	self.Mbt.grid(column = 0, row = 2, pady=(2,4))
	self.Mbt.menu  =  Tkinter.Menu ( self.Mbt, tearoff = 0 , background='coral1', fg="black" )
		
	self.Mbt.menu.add_command(label="NN", command = lambda v=1: _SetMode(master, v))
	self.Mbt.menu.add_command(label="IDW", command = lambda v=2: _SetMode(master, v))
	self.Mbt.menu.add_command(label="Gaussian", command = lambda v=3: _SetMode(master, v))
	self.Mbt["menu"]  =  self.Mbt.menu
	self.Mbt.configure(background='gold2', fg="black")

	#=======================        
	# epsg
	self.epsg = Tkinter.StringVar()
	self.epsg_entry = Tkinter.Entry(map_texture_frame, width = 30, textvariable = self.epsg)
	self.epsg_entry.grid(column = 1, row = 2, pady=(2,4)) #, columnspan = 2,sticky = 'EW')
	self.epsg_entry.bind("<Return>", lambda epsg=self.epsg.get(): _OnPressEnter1(self))
	self.epsg.set(u"epsg:26949")       
	self.epsg_entry.configure(background='gold2', fg="black")

	#=======================
	# res
	self.resvar = Tkinter.DoubleVar()
	#resscale = Tkinter.Scale( map_texture_frame, variable = self.resvar, from_=0.2, to=10, resolution=0.01, tickinterval=1, label = 'Resolution [m]' )
	resscale = Tkinter.Scale( map_texture_frame, variable = self.resvar, from_=0, to=10, resolution=0.1, tickinterval=0.25, label = 'Resolution [m] 0 = auto' )
	resscale.set(0.5)
	resscale.grid(row=3, column=0,  pady=(2,4))
	resscale.configure(background='gold2', fg="black")

	#=======================
	# nn
	self.nnvar = Tkinter.DoubleVar()
	nnscale = Tkinter.Scale( map_texture_frame, variable = self.nnvar, from_=1, to=512, resolution=1, tickinterval=50, label = 'Nearest neighbours' )
	nnscale.set(64)
	nnscale.grid(row=3, column=1,  pady=(2,4))
	nnscale.configure(background='gold2', fg="black")

	#=======================
	# influence
	#self.infvar = Tkinter.DoubleVar()
	#infscale = Tkinter.Scale( map_texture_frame, variable = self.infvar, from_=0.1, to=10, resolution=0.1, tickinterval=1, label = 'Influence [m]' )
	#infscale.set(1)
	#infscale.grid(row=4, column=0,  pady=(2,4))
	#infscale.configure(background='gold2', fg="black")

	#=======================
	# numstdevs
	self.nstdvar = Tkinter.DoubleVar()
	nstdscale = Tkinter.Scale( map_texture_frame, variable = self.nstdvar, from_=1, to=10, resolution=1, tickinterval=2, label = 'Numb. stan. dev.' )
	nstdscale.set(5)
	nstdscale.grid(row=4, column=1,  pady=(2,4))
	nstdscale.configure(background='gold2', fg="black")

	#=======================
	# process button
	proc6_btn = Tkinter.Button(map_texture_frame, text='Process!', underline=0,
		         command=lambda watervar=self.watervar.get(): _proc6(self))
	proc6_btn.grid(row=5, column=1, pady=(2,4))
	proc6_btn.configure(background='gold2', fg="black")

	#==============================================================
	#========START functions for map_texture tab         

	#=======================
	def _proc6(self):
	    # function to invoke PyHum.map_texture
	    # build error checking into here
	    print 'Processing ...'   
	#    print "humfile: " + self.DATfilename.get()
	#    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	#    print "mode: " + str(self.mode)        
	#    print "cs2cs arguments: " + str(self.epsg.get())   
	#    print "resolution: " + str(self.resvar.get())
	#    print "max. number of nearest neighbours: " + str(self.nnvar.get()) 
	#    print "gridding influence [m]: " + str(self.infvar.get()) 
	#    print "number std. dev. to accept: " + str(self.nstdvar.get())             
		             
	    # do stuff here
	    PyHum.map_texture(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), self.epsg.get(),self.resvar.get(), self.mode, self.nnvar.get(), self.nstdvar.get())  #self.infvar.get(),  
	    self.update() 
	    tkMessageBox.showinfo("Done!", "Map texture module finished") 
		
	#==============================================================
	#==============================================================
	#========END map_texture tab


	#==============================================================
	#==============================================================
	#========START e1e2 tab

	e1e2_frame = Tkinter.Frame(nb)
	nb.add(e1e2_frame, text='Bed Class')#, state='disabled')

	e1e2_frame.configure(background='PeachPuff2')
	 
	e1e2_msg = [
	    "Analysis of first (e1, 'roughness') and second (e2, 'hardness') echo returns\n",
	    "from the high-frequency downward looking echosounder. Generates\n",
	    "generalised acoustic parameters for the purposes of point classification\n",
	    "of submerged substrates/vegetation\n\n",
	    "Accounts for the absorption of sound in water\n\n",
	    "Does a basic k-means cluster of e1 and e2 coefficients into specified number of 'acoustic classes'\n\n",
	    "and creates some rudimentary plots, kml files and text outputs\n\n",
	    "Based on code by Barb Fagetter (blueseas@oceanecology.ca)"]
	    
	lbl8 = Tkinter.Label(e1e2_frame, wraplength='4i', justify=Tkinter.LEFT, anchor=Tkinter.N,
		        text=''.join(e1e2_msg))

	lbl8.configure(background='SlateGray1', fg="black")
		        
	# position and set resize behavior
	lbl8.grid(row=0, column=0, columnspan=1, sticky='new', pady=5)

	def hello7():
	   tkMessageBox.showinfo("Bed class instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n H: water acidity\n\n Temp: water temperature in deg C\n\n salinity: in ppt\n\n beam: downward echosounder beam width in degrees\n\n Transducer frequency: downward echosounder frequency in kHz\n\n Integ: number of pings over which to integrate\n\n Number of clusters: number of acoustic classes that will be generated")
	    
	def hello7_alt():
	   try:
	      root = Tk()
	      root.wm_title("Bed class module")      
	      S = Scrollbar(root)
	      T = Text(root, height=40, width=60, wrap=WORD)
	      S.pack(side=RIGHT, fill=Y)
	      T.pack(side=LEFT, fill=Y)
	      S.config(command=T.yview)
	      T.config(yscrollcommand=S.set)

	      T.tag_configure('bold_italics', font=('Arial', 12, 'bold', 'italic'))
	      T.tag_configure('big', font=('Verdana', 20, 'bold'))
	      T.tag_configure('color', foreground='#476042', font=('Tempus Sans ITC', 12, 'bold'))

	      quote = """Bed class Instructions", "DAT file: path to the .DAT file\n\n SON files: path to *.SON files\n\n cs2cs_args: argument given to pyproj to turn wgs84 coords. to projection supported by proj.4. Default='epsg:26949'\n\n pH: water acidity\n\n Temp: water temperature in deg C\n\n salinity: in ppt\n\n beam: downward echosounder beam width in degrees\n\n Beam frequency: downward echosounder frequency in kHz\n\n Number of classes: number of acoustic classes that will be generated"""
	      T.insert(END, quote)  
	   except:
	      hello7()
	      
	MSG7_btn = Tkinter.Button(e1e2_frame, text = "Instructions", command = hello7_alt)
	MSG7_btn.grid(row=0, column=1, pady=(2,4))
	MSG7_btn.configure(background='SlateGray1', fg="black")

	e1e2_frame.rowconfigure(1, weight=1)
	e1e2_frame.columnconfigure((0,1), weight=1, uniform=1)

	#=======================
	# get dat file           
	datVar = Tkinter.StringVar()
	self.e1e2_dat_btn = Tkinter.Button(e1e2_frame, text='Get DAT file', underline=0,
		         command=lambda v=datVar: _get_DAT(master, v))
	dat = Tkinter.Label(e1e2_frame, textvariable=datVar, name='dat')
	self.e1e2_dat_btn.grid(row=1, column=0, pady=(2,4))
	self.e1e2_dat_btn.configure(background='SlateGray1', fg="black")

	#=======================
	# get son files
	sonVar = Tkinter.StringVar()
	self.e1e2_son_btn = Tkinter.Button(e1e2_frame, text='Get SON files', underline=0,
		         command=lambda v=datVar: _get_SON(master, v))
	son = Tkinter.Label(e1e2_frame, textvariable=sonVar, name='dat')
	self.e1e2_son_btn.grid(row=1, column=1, pady=(2,4))
	self.e1e2_son_btn.configure(background='SlateGray1', fg="black")

	#=======================
	# num acoustic classes
	self.Nvar = Tkinter.DoubleVar()
	Nscale = Tkinter.Scale( e1e2_frame, variable = self.Nvar, from_=2, to=10, resolution=1, tickinterval=2, label = 'Num. classes' )
	Nscale.set(3)
	Nscale.grid(row=2, column=0,  pady=(2,4))
	Nscale.configure(background='SlateGray1', fg="black")

	#=======================
	# pH
	self.Pvar = Tkinter.DoubleVar()
	Pscale = Tkinter.Scale( e1e2_frame, variable = self.Pvar, from_=4, to=8, resolution=.1, tickinterval=1, label = 'pH' )
	Pscale.set(7)
	Pscale.grid(row=2, column=1,  pady=(2,4))
	Pscale.configure(background='SlateGray1', fg="black")

	#=======================
	# temp
	self.Tvar = Tkinter.DoubleVar()
	Tscale = Tkinter.Scale( e1e2_frame, variable = self.Tvar, from_=3, to=25, resolution=1, tickinterval=5, label = 'Temp [degC]' )
	Tscale.set(10)
	Tscale.grid(row=3, column=0,  pady=(2,4))
	Tscale.configure(background='SlateGray1', fg="black")

	#=======================
	# salinity
	self.Svar = Tkinter.DoubleVar()
	Sscale = Tkinter.Scale( e1e2_frame, variable = self.Svar, from_=0, to=35, resolution=1, tickinterval=5, label = 'Salinity [ppt]' )
	Sscale.set(0)
	Sscale.grid(row=3, column=1,  pady=(2,4))
	Sscale.configure(background='SlateGray1', fg="black")

	#=======================
	# beam freq
	self.BFvar = Tkinter.DoubleVar()
	BFscale = Tkinter.Scale( e1e2_frame, variable = self.BFvar, from_=100, to=500, resolution=5, tickinterval=100, label = 'Beam freq. [kHz]' )
	BFscale.set(200)
	BFscale.grid(row=4, column=0,  pady=(2,4))
	BFscale.configure(background='SlateGray1', fg="black")

	#=======================
	# beam angle
	self.BAvar = Tkinter.DoubleVar()
	BAscale = Tkinter.Scale( e1e2_frame, variable = self.BAvar, from_=5, to=35, resolution=1, tickinterval=5, label = 'Beam angle [deg]' )
	BAscale.set(20)
	BAscale.grid(row=4, column=1,  pady=(2,4))
	BAscale.configure(background='SlateGray1', fg="black")

	#=======================        
	# epsg
	self.epsg = Tkinter.StringVar()
	self.epsg_entry = Tkinter.Entry(e1e2_frame, width = 30, textvariable = self.epsg)
	self.epsg_entry.grid(column = 0, row = 5, pady=(2,4)) #, columnspan = 2,sticky = 'EW')
	self.epsg_entry.bind("<Return>", lambda epsg=self.epsg.get(): _OnPressEnter1(self))
	self.epsg.set(u"epsg:26949")       
	self.epsg_entry.configure(background='SlateGray1', fg="black")

	#=======================
	# process button
	proc7_btn = Tkinter.Button(e1e2_frame, text='Process!', underline=0,
		         command=lambda watervar=self.watervar.get(): _proc7(self))
	proc7_btn.grid(row=5, column=1, pady=(2,4))
	proc7_btn.configure(background='SlateGray1', fg="black")

	#==============================================================
	#========START functions for e1e2 tab         

	#=======================
	def _proc7(self):
	    # function to invoke PyHum.e1e2
	    # build error checking into here
	    print 'Processing ...'   
	    print "humfile: " + self.DATfilename.get()
	    print "sonpath: " + os.path.dirname(self.SONfiles[0])
	    print "pH: " + str(self.Pvar.get())    
	    print "temperature: " + str(self.Tvar.get())    
	    print "salinity: " + str(self.Svar.get())  
	    print "beam frequency: " + str(self.BFvar.get())  
	    print "beam angle: " + str(self.BAvar.get())            
	    print "number of acoustic classes: " + str(self.Nvar.get())            
	    # do stuff here
	    PyHum.e1e2(self.DATfilename.get(), os.path.dirname(self.SONfiles[0]), self.epsg.get(), self.Pvar.get(), self.Tvar.get(), self.Svar.get(), self.BAvar.get(), self.BFvar.get(), self.integ, self.Nvar.get(), self.doplot)       
	    self.update() 
	    tkMessageBox.showinfo("Done!", "Bed class finished") 
		
	#    # for downward-looking echosounder echogram (e1-e2) analysis
	#   beam = 20.0
	#   transfreq = 200.0 # frequency (kHz) of downward looking echosounder
	#   numclusters = 3 # number of acoustic classes to group observations
	   
	    
	#==============================================================
	#==============================================================
	#========END e1e2 tab

	# start app
	master.mainloop()


# =========================================================
# =========================================================
if __name__ == '__main__':
   
   gui()

