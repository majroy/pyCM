#!/usr/bin/env python
'''
VTK and QT to allow for preprocessing FEAs associated with the contour method. 
See documentation for the following main methods in this package:
-point_cloud
-align_average
-fit_surface
-preprocess
-------------------------------------------------------------------------------
ver 1.0 17-01-04 - contains all elements except for a post-processor
'''
__author__ = "M.J. Roy"
__version__ = "1.0"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2018"

import sys
import os.path
from pkg_resources import Requirement, resource_filename
import numpy as np
import scipy.io as sio
import vtk
import vtk.util.numpy_support as vtk_to_numpy
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from .pyCMcommon import *
import point_cloud as pc
import align_average as aa
import fit_surface as sf
import preprocess as pre

class Ui_MainWindow(object):
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.setWindowTitle("pyCM v%s" %__version__)
		MainWindow.setEnabled(True)
		MainWindow.resize(1280, 720)

		self.centralwidget = QtWidgets.QWidget(MainWindow)
		self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
		self.mainLayout = QtWidgets.QHBoxLayout(self.centralwidget)
		self.mainLayout.addWidget(self.tabWidget)
		MainWindow.setCentralWidget(self.centralwidget)
		self.tabLayout = QtWidgets.QHBoxLayout(self.tabWidget)

		#add tabs
		self.pctab = QtWidgets.QWidget(self.tabWidget)
		self.tabWidget.addTab(self.pctab, "Point cloud editor")

		self.aatab = QtWidgets.QWidget()
		self.tabWidget.addTab(self.aatab, "Alignment/averaging")
		self.sftab = QtWidgets.QWidget()
		self.tabWidget.addTab(self.sftab, "Surface fitting")
		self.pretab = QtWidgets.QWidget()
		self.tabWidget.addTab(self.pretab, "Preprocessing")
		self.posttab = QtWidgets.QWidget()
		
		MainWindow.setCentralWidget(self.centralwidget)
		#detect changes between tabs
		self.tabWidget.currentChanged.connect(self.on_change)
		#add menubar
		self.menubar = QtWidgets.QMenuBar(MainWindow)
		#do file menu first
		fileMenu = self.menubar.addMenu('&File')
		
		loadButton = QtWidgets.QAction('Load', MainWindow)
		loadButton.setShortcut('Ctrl+L')
		loadButton.setStatusTip('Load pyCM data file')
		loadButton.triggered.connect(self.populate)
		
		exitButton = QtWidgets.QAction('Exit', MainWindow)
		exitButton.setShortcut('Ctrl+Q')
		exitButton.setStatusTip('Exit application')
		exitButton.triggered.connect(MainWindow.close)
		#add buttons to menus
		fileMenu.addAction(loadButton)
		fileMenu.addAction(exitButton)

		
		#options menu
		optMenu = self.menubar.addMenu('&Options')
		setFEA = QtWidgets.QAction('Set FEA configuration', MainWindow)
		setFEA.setShortcut('Ctrl+E')
		setFEA.setStatusTip('Identify/set FEA executable locations')
		setFEA.triggered.connect(self.getFEAconfig)
		
		setWorkDir = QtWidgets.QAction('Set FEA working directory', MainWindow)
		setWorkDir.setShortcut('Ctrl+D')
		setWorkDir.setStatusTip('Select working directory to perform FEA')
		setWorkDir.triggered.connect(self.setFEAworkdir)
		
		optMenu.addAction(setWorkDir)
		optMenu.addAction(setFEA)
		
		##To be implemented
		#reporting menu
		optMenu = self.menubar.addMenu('&Reporting')

		MainWindow.setMenuBar(self.menubar)
		self.statusbar = QtWidgets.QStatusBar(MainWindow)
		self.statusbar.setObjectName("statusbar")
		MainWindow.setStatusBar(self.statusbar)

		self.tabWidget.setCurrentIndex(0)
		QtCore.QMetaObject.connectSlotsByName(MainWindow)

		#populate tabs
		self.setup_pc()
		self.setup_aa()
		self.setup_sf()
		self.setup_pre()
		
	def setup_pc(self):
		lhLayout=QtWidgets.QHBoxLayout(self.pctab)
		self.pctab.setLayout(lhLayout)
		self.pcui=pc.pnt_interactor(self.centralwidget)
		lhLayout.addWidget(self.pcui)
		self.pcui.iren.Initialize()
		self.pcui.unsaved_changes=False

	def setup_aa(self):
		lhLayout=QtWidgets.QHBoxLayout(self.aatab)
		self.aatab.setLayout(lhLayout)
		self.aaui=aa.aa_interactor(self.centralwidget)
		lhLayout.addWidget(self.aaui)
		self.aaui.iren.Initialize()
		self.aaui.averaged=False
		self.aaui.unsaved_changes=False
		
	def setup_sf(self):
		lhLayout=QtWidgets.QHBoxLayout(self.sftab)
		self.sftab.setLayout(lhLayout)
		self.sfui=sf.surf_int(self.centralwidget)
		lhLayout.addWidget(self.sfui)
		self.sfui.iren.Initialize()
		self.sfui.unsaved_changes=False

	def setup_pre(self):
		lhLayout=QtWidgets.QHBoxLayout(self.pretab)
		self.pretab.setLayout(lhLayout)
		self.preui=pre.msh_interactor(self.centralwidget)
		lhLayout.addWidget(self.preui)
		self.preui.iren.Initialize()
		self.preui.unsaved_changes=False

	def populate(self):
		"""
		Populates all ui's with contents of mat file
		"""
		self.activeFile, _, = get_file('*.mat')
		MainWindow.setWindowTitle("%s  -  pyCM v%s" %(self.activeFile,__version__))
		self.pcui.fileo=self.activeFile
		self.pcui.load_mat()
		self.aaui.get_input_data(self.activeFile)
		self.sfui.get_input_data(self.activeFile)
		self.preui.get_input_data(self.activeFile)
		
	def on_change(self):
		#check if there are unsaved changes from the point editor
		if self.pcui.unsaved_changes:
			ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
				"There are unsaved changes pending to a point cloud. Ignore?", \
				QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
			if ret == QtWidgets.QMessageBox.Yes: #don't incorporate
				#reload relevant parts of results file & clear unsaved changes.
				self.pcui.load_mat()
				self.pcui.unsaved_changes=False
				return
			else: 
				#change unsaved flag temporarily to suppress on_change dialog
				self.pcui.unsaved_changes=False
				self.tabWidget.setCurrentIndex(0)
				self.pcui.unsaved_changes=True

		#check if there's an activeFile
		if not hasattr(self,'activeFile'):
			#then this is the first tab change since running the point editor
			try:
				self.activeFile=self.pcui.fileo
				MainWindow.setWindowTitle("%s  -  pyCM v%s" %(self.activeFile,__version__))
				self.aaui.get_input_data(self.activeFile)
			except:
				ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
				"No saved data to progress analysis.", \
				QtWidgets.QMessageBox.Ok)

		#check if there are unsaved changes from alignment/averaging
		if self.aaui.unsaved_changes:
			ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
				"There are unsaved changes pending to a point cloud. Ignore?", \
				QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
			if ret == QtWidgets.QMessageBox.Yes: #don't incorporate
				#reload relevant parts of results file & clear unsaved changes.
				self.aaui.get_input_data(self.activeFile)
				return
			else: 
				#change unsaved flag temporarily to suppress on_change dialog
				self.aaui.unsaved_changes=False
				self.tabWidget.setCurrentIndex(1)
				self.aaui.unsaved_changes=True
				
		#check the status of alignment/averaging against the surface input file
		if self.aaui.averaged and not hasattr(self.sfui,"fileo"):
			try:
				self.sfui.get_input_data(self.activeFile)

			except:
				ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
				"No saved data to progress analysis.", \
				QtWidgets.QMessageBox.Ok)
				
	def getFEAconfig(self):
		#check config file
		try:
			self.filec = resource_filename("pyCM","pyCMconfig.yml")#needs to be pyCM
		except: #resource_filename will inform if the directory doesn't exist
			print("Did not find config file in the pyCM installation directory.")
		try:
			with open(self.preui.filec,'r') as ymlfile:
				readcfg = yaml.load(ymlfile)	
				l=[readcfg['FEA']['abaqusExec'],readcfg['FEA']['gmshExec'],readcfg['FEA']['ccxExec']]
				self.preui.cfg=pre.GetFEAconfig(l,self.filec)
		except:
			try:
				pre.GetFEAconfig(['','',''],self.filec)
			except:
				print("Failed to set config file.")
				return
	
	def setFEAworkdir(self):
		#sets directory to run FEA from - writes outputd variable to preprocessor interactor
		self.preui.outputd = str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select directory", "",QtWidgets.QFileDialog.ShowDirsOnly))

if __name__ == "__main__":
	import sys
	app = QtWidgets.QApplication(sys.argv)

	
	
	spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
	splash_pix = QtGui.QPixmap(spl_fname,'PNG')
	splash = QtWidgets.QSplashScreen(splash_pix)
	splash.setMask(splash_pix.mask())

	splash.show()
	app.processEvents()
	
	MainWindow = QtWidgets.QMainWindow()
	ui = Ui_MainWindow()
	ui.setupUi(MainWindow)
	
	MainWindow.show()
	splash.finish(MainWindow)
	
	sys.exit(app.exec_())

