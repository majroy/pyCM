#!/usr/bin/env python
'''
Uses VTK python to allow for editing point clouds associated with the contour 
method. Full interaction requires a 3-button mouse and keyboard.
-------------------------------------------------------------------------------
Current mapping is as follows:
LMB - rotate about point cloud centroid.
MMB - pan
RMB - zoom/refresh window extents
1 - view 1, default, looks down z axis onto xy plane
2 - view 2, looks down x axis onto zy plane
3 - view 3, looks down y axis onto zx plane
z - increase z-aspect ratio
x - decrease z-aspect ratio
c - return to default z-aspect
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
r - remove compass/axes
h - flip on x plane
k - flip on y
d - flip visibility of reference and floating point clouds/outlines
a - align
A - align using perimeters
v - average
e - write output and exit
-------------------------------------------------------------------------------
ver 1.1 16-11-06
1.1 - Initial release
1.2 - Refactored to use pyQt interface and eliminated global variables
'''
__author__ = "M.J. Roy"
__version__ = "1.2"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2017"

import sys
import os.path
import vtk
import vtk.util.numpy_support as VN
import numpy as np
import numpy.matlib
import scipy.io as sio
from scipy.interpolate import griddata
from scipy.spatial.distance import pdist, squareform
from matplotlib import path
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt4 import QtCore, QtGui
from pkg_resources import Requirement, resource_filename
from pyCMcommon import *


def ali_avg_interactor(*args,**kwargs):
	"""
	Main function, builds qt interaction
	"""	
	app = QtGui.QApplication.instance()
	if app is None:
		app = QApplication(sys.argv)
	
	spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
	splash_pix = QtGui.QPixmap(spl_fname,'PNG')
	splash = QtGui.QSplashScreen(splash_pix)
	splash.setMask(splash_pix.mask())

	splash.show()
	app.processEvents()
	
	window = aa_interactor()

	if len(args)==1: 
		aa_interactor.get_input_data(window,args[0])
	else: 
		aa_interactor.get_input_data(window,None)

	window.show()
	splash.finish(window)
	window.iren.Initialize() # Need this line to actually show the render inside Qt

	ret = app.exec_()
	
	if sys.stdin.isatty() and not hasattr(sys,'ps1'):
		sys.exit(ret)
	else:
		return window

class ali_avg(object):
	"""
	Class to build qt interaction, including VTK widget
	setupUi builds, initialize starts VTK widget
	"""
	
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.setWindowTitle("pyCM - Alignment and averaging tool v%s" %__version__)
		MainWindow.resize(1280, 720)
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.Boxlayout = QtGui.QHBoxLayout(self.centralWidget)
		self.Subtendlayout=QtGui.QVBoxLayout()
		mainUiBox = QtGui.QFormLayout()

		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
		self.vtkWidget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		self.vtkWidget.setMinimumSize(1150, 700); #leave 100 px on the size for i/o

		self.Subtendlayout.addWidget(self.vtkWidget)
		self.activeFileLabel=QtGui.QLabel("Idle")
		self.activeFileLabel.setWordWrap(True)
		self.activeFileLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.activeFileLabel.setMinimumWidth(100)
		self.Subtendlayout.addWidget(self.activeFileLabel)
		self.Subtendlayout.addStretch(1)
		self.Boxlayout.addLayout(self.Subtendlayout)
		self.Boxlayout.addStretch(1)
		
		MainWindow.setCentralWidget(self.centralWidget)
		headFont=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)
		
		# #define buttons/widgets
		self.reloadButton = QtGui.QPushButton('Load')
		scalingLabel=QtGui.QLabel("Active axis for scaling")
		scalingLabel.setFont(headFont)
		self.xsButton=QtGui.QRadioButton("x")
		self.ysButton=QtGui.QRadioButton("y")
		self.zsButton=QtGui.QRadioButton("z")
		self.zsButton.setChecked(True)
		self.scalingButtonGroup = QtGui.QButtonGroup()
		self.scalingButtonGroup.addButton(self.xsButton)
		self.scalingButtonGroup.addButton(self.ysButton)
		self.scalingButtonGroup.addButton(self.zsButton)
		self.scalingButtonGroup.setExclusive(True)
		scaleBoxlayout = QtGui.QGridLayout()
		scaleBoxlayout.addWidget(self.xsButton,1,1)
		scaleBoxlayout.addWidget(self.ysButton,1,2)
		scaleBoxlayout.addWidget(self.zsButton,1,3)
		
		horizLine1=QtGui.QFrame()
		horizLine1.setFrameStyle(QtGui.QFrame.HLine)
		mirrorLabel=QtGui.QLabel("Mirroring")
		mirrorLabel.setFont(headFont)
		self.mirrorXbutton = QtGui.QPushButton('ZY')
		self.mirrorYbutton = QtGui.QPushButton('ZX')

		horizLine2=QtGui.QFrame()
		horizLine2.setFrameStyle(QtGui.QFrame.HLine)
		alignLabel=QtGui.QLabel("Alignment")
		alignLabel.setFont(headFont)
		self.transXlabel=QtGui.QLabel("Translate x:")
		self.transX = QtGui.QLineEdit()
		self.transX.setText('0')
		self.transX.setMinimumWidth(50)
		self.transYlabel=QtGui.QLabel("Translate y:")
		self.transY = QtGui.QLineEdit()
		self.transY.setText('0')
		self.transY.setMinimumWidth(50)
		self.transButton=QtGui.QPushButton('Translate floating')
		
		self.alignButtonGroup = QtGui.QButtonGroup()
		self.alignOutlineButton=QtGui.QRadioButton("Outline")
		self.alignPointCloudButton = QtGui.QRadioButton("Point cloud")
		self.alignOutlineButton.setChecked(True)
		self.alignButtonGroup.addButton(self.alignOutlineButton)
		self.alignButtonGroup.addButton(self.alignPointCloudButton)
		self.alignButtonGroup.setExclusive(True)
		self.alignButton = QtGui.QPushButton('Align')
		
		horizLine3=QtGui.QFrame()
		horizLine3.setFrameStyle(QtGui.QFrame.HLine)
		averageLabel=QtGui.QLabel("Averaging")
		averageLabel.setFont(headFont)
		self.averageButton = QtGui.QPushButton('Average')
		
		horizLine4=QtGui.QFrame()
		horizLine4.setFrameStyle(QtGui.QFrame.HLine)
		self.writeButton=QtGui.QPushButton('Write')
		
		horizLine5=QtGui.QFrame()
		horizLine5.setFrameStyle(QtGui.QFrame.HLine)
		self.statusLabel=QtGui.QLabel("Idle")
		self.statusLabel.setWordWrap(True)
		self.statusLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.statusLabel.setMinimumWidth(50)

		# #add to formlayout
		mainUiBox.addRow(self.reloadButton)
		mainUiBox.addRow(scalingLabel)
		mainUiBox.addRow(scaleBoxlayout)
		mainUiBox.addRow(horizLine1)
		mainUiBox.addRow(mirrorLabel)
		mainUiBox.addRow(self.mirrorXbutton,self.mirrorYbutton)
		mainUiBox.addRow(horizLine2)
		mainUiBox.addRow(alignLabel)
		mainUiBox.addRow(self.transXlabel,self.transX)
		mainUiBox.addRow(self.transYlabel,self.transY)
		mainUiBox.addRow(self.transButton)
		mainUiBox.addRow(self.alignOutlineButton,self.alignPointCloudButton)
		mainUiBox.addRow(self.alignButton)
		mainUiBox.addRow(horizLine3)
		mainUiBox.addRow(averageLabel)
		mainUiBox.addRow(self.averageButton)
		mainUiBox.addRow(horizLine4)
		mainUiBox.addRow(self.writeButton)
		mainUiBox.addRow(horizLine5)
		mainUiBox.addRow(self.statusLabel)
		
		self.Boxlayout.addLayout(mainUiBox)
		
	def initialize(self):
		self.vtkWidget.start()

class aa_interactor(QtGui.QMainWindow):
	"""
	Sets up the main VTK window, reads file and sets connections between UI and interactor
	"""
	def __init__(self, parent = None):
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = ali_avg()
		self.ui.setupUi(self)
		self.ren = vtk.vtkRenderer()
		self.ren.SetBackground(0.1, 0.2, 0.4)
		self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
		self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
		style=vtk.vtkInteractorStyleTrackballCamera()
		style.AutoAdjustCameraClippingRangeOn()
		self.iren.SetInteractorStyle(style)
		self.ren.GetActiveCamera().ParallelProjectionOn()
		self.cp=self.ren.GetActiveCamera().GetPosition()
		self.fp=self.ren.GetActiveCamera().GetFocalPoint()
		self.iren.AddObserver("KeyPressEvent", self.keypress)
		
		self.PointSize=2
		self.LineWidth=1
		self.Zaspect=1.0
		self.limits=np.empty(6)
		self.picking=False
		self.Offset=0
		self.mirrored=False
		self.aligned=False
		
		self.ui.reloadButton.clicked.connect(lambda: self.get_input_data(None))
		self.ui.mirrorXbutton.clicked.connect(lambda: self.flipside('x'))
		self.ui.mirrorYbutton.clicked.connect(lambda: self.flipside('y'))
		self.ui.transButton.clicked.connect(lambda: self.shift())
		self.ui.alignButton.clicked.connect(lambda: self.align())
		self.ui.averageButton.clicked.connect(lambda: self.average())
		self.ui.writeButton.clicked.connect(lambda: self.write())
	
	def shift(self):
		if hasattr(self,'fActor'): #then remove this actor and the associated outline actor
			self.ren.RemoveActor(self.fActor)
			self.ren.RemoveActor(self.fOutlineActor)

		#get x and y tranformations
		gx=float(self.ui.transX.text())
		gy=float(self.ui.transY.text())
		transl=np.array([gx, gy, 0]);
		
		#apply operation
		self.flp=self.flp+transl
		self.fO=self.fO+transl
		
		color=(255, 205, 52)
		self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
		self.ren.AddActor(self.fActor)
		self.fOutlineActor, self.fOPC = gen_outline(self.fO,color,self.PointSize)
		self.ren.AddActor(self.fOutlineActor)
		
		#recalculate extents of interactor
		RefMin=np.amin(np.vstack((self.flp,self.rp)),axis=0)
		RefMax=np.amax(np.vstack((self.flp,self.rp)),axis=0)

		extents=RefMax-RefMin #extents
		rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
		self.limits=[RefMin[0]-rl, \
		RefMax[0]+rl, \
		RefMin[1]-rl, \
		RefMax[1]+rl, \
		RefMin[2],RefMax[2]]

		#add axes
		self.add_axis(self.limits,[1,1,1])
		
		s,nl,axs=self.get_scale()

		self.fActor.SetScale(s)
		self.fActor.Modified()

		self.add_axis(nl,axs)
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()	
	
	def write(self):
		mat_contents=sio.loadmat(self.fileo)
		
		new={'transM':self.transM,'aa':self.ap}
		
		mat_contents.update(new) #update the dictionary
			
		sio.savemat(self.fileo,mat_contents)	
		self.ui.statusLabel.setText("Wrote data. Idle.")
	
	def average(self):
		# if not self.aligned:
			# self.ui.statusLabel.setText("Align prior to averaging. Idle.")
			# return
		
		self.ui.statusLabel.setText("Averaging, applying grid . . .")
		QtGui.qApp.processEvents()
		
		#temporarily shift all data such that it appears in the first cartesian quadrant
		tT=np.amin(self.rO,axis=0)
		self.rO, self.fO, self.rp, self.flp=self.rO-tT, self.fO-tT, self.rp-tT, self.flp-tT
		
		#use max to get a 'window' for assessing grid spacing
		RefMax=np.amax(self.rO,axis=0)
		RefMin=np.amin(self.rO,axis=0)
		windowVerts=np.matrix([[0.25*RefMin[0], 0.25*RefMin[1]],
		[0.25*RefMin[0], 0.25*(RefMax[1])],
		[0.25*(RefMax[1]), 0.25*(RefMax[1])],
		[0.25*(RefMax[0]), 0.25*(RefMin[1])]]);
		
		p=path.Path(windowVerts)
		inWindow=p.contains_points(self.rp[:,:2]) #first 2 columns of RefPoints is x and y
		
		windowed=self.rp[inWindow,:2]
		
		gs=squareform(pdist(windowed,'euclidean')) #does the same thing as pdist2
		
		gsize=np.mean(np.sort(gs)[:,1]) #sort the distances, find the mean distance between non self-referencing points
		
		#grid the reference based on gsize, bumping out the grid by 10% in either direction
		grid_x, grid_y = np.meshgrid(
		np.linspace(1.1*RefMin[0],1.1*RefMax[0],int((1.1*RefMax[0]-1.1*RefMin[0])/gsize)),
		np.linspace(1.1*RefMin[1],1.1*RefMax[1],int((1.1*RefMax[1]-1.1*RefMin[1])/gsize)), 
		indexing='xy')
		
		#apply the grid to the reference data
		grid_Ref=griddata(self.rp[:,:2],self.rp[:,-1],(grid_x,grid_y),method='linear')
		
		#apply the grid to the aligned data
		grid_Align=griddata(self.flp[:,:2],self.flp[:,-1],(grid_x,grid_y),method='linear')
		
		self.ui.statusLabel.setText("Averaging using grid . . .")
		QtGui.qApp.processEvents()
		
		#average z values
		grid_Avg=(grid_Ref+grid_Align)/2
		
		#make sure that there isn't anything averaged outside the floating outline
		p=path.Path(self.rO[:,:2])
		inTest=np.hstack((np.ravel(grid_x.T)[np.newaxis].T,np.ravel(grid_y.T)[np.newaxis].T))
		inOutline=p.contains_points(inTest)
		
		#averaged points
		self.ap = np.hstack((inTest[inOutline,:], \
					np.ravel(grid_Avg.T)[np.newaxis].T[inOutline]))
					
		#move everything back to original location
		self.rO, self.fO, self.rp, self.flp, self.ap = \
		self.rO+tT, self.fO+tT, self.rp+tT, self.flp+tT, self.ap+tT
		
		self.ui.statusLabel.setText("Rendering . . .")
		QtGui.qApp.processEvents()
		
		#show it
		color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
		_, self.aActor, _, = gen_point_cloud(self.ap,color,self.PointSize)
		self.ren.AddActor(self.aActor)
		
		s,nl,axs=self.get_scale()

		self.aActor.SetScale(s)
		self.aActor.Modified()
		
		#update
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		self.ui.statusLabel.setText("Averaging complete. Idle.")
		self.averaged=True
	
	def flipside(self,flipDirection):
		self.ui.statusLabel.setText("Starting mirroring . . .")
		self.ui.vtkWidget.update()
		#delete the floating actor
		if hasattr(self,'fActor'): #then remove this actor and the associated outline actor
			self.ren.RemoveActor(self.fActor)
			self.ren.RemoveActor(self.fOutlineActor)
		else:
			print "Need to have data loaded before manipulating . . .\n"

		if flipDirection == "x":
			self.flp[:,0]=-self.flp[:,0]
			self.fO[:,0]=-self.fO[:,0]
		elif flipDirection == "y":
			self.flp[:,1]=-self.flp[:,1]
			self.fO[:,1]=-self.fO[:,1]
		
		color=(255, 205, 52)
		self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
		self.ren.AddActor(self.fActor)
		self.fOutlineActor, self.fOPC = gen_outline(self.fO,color,self.PointSize)
		self.ren.AddActor(self.fOutlineActor)
		
		#recalculate extents of interactor
		RefMin=np.amin(np.vstack((self.flp,self.rp)),axis=0)
		RefMax=np.amax(np.vstack((self.flp,self.rp)),axis=0)

		extents=RefMax-RefMin #extents
		rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
		self.limits=[RefMin[0]-rl, \
		RefMax[0]+rl, \
		RefMin[1]-rl, \
		RefMax[1]+rl, \
		RefMin[2],RefMax[2]]

		#add axes
		self.add_axis(self.limits,[1,1,1])
		
		s,nl,axs=self.get_scale()

		self.fActor.SetScale(s)
		self.fActor.Modified()

		self.add_axis(nl,axs)
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		self.ui.statusLabel.setText("Mirror operation complete. Idle.")
		self.mirrored=True
		
	def align(self):
		self.ui.statusLabel.setText("Starting alignment . . .")
		QtGui.qApp.processEvents()
		
		icp=vtk.vtkIterativeClosestPointTransform()
		if self.ui.alignPointCloudButton.isChecked():
			icp.SetSource(self.fPC)
			icp.SetTarget(self.rPC)
		elif self.ui.alignOutlineButton.isChecked():
			icp.SetSource(self.fOPC)
			icp.SetTarget(self.rOPC)
		
		if hasattr(self,'fActor'): #then remove this actor and the associated outline actor
			self.ren.RemoveActor(self.fActor)
			self.ren.RemoveActor(self.fOutlineActor)
			

		
		icp.SetMaximumNumberOfIterations(200)
		icp.StartByMatchingCentroidsOn()
		icp.Modified()
		icp.Update()
		icp.Inverse()
		
		self.transM=np.zeros(shape=(4,4))
		for i in range(4):
			for j in range(4):
				self.transM[i,j]=icp.GetMatrix().GetElement(i, j)
		self.transM=np.linalg.inv(self.transM)

		#apply operation
		self.flp=np.dot(self.flp,self.transM[0:3,0:3])+self.transM[0:3,-1]
		self.fO=np.dot(self.fO,self.transM[0:3,0:3])+self.transM[0:3,-1]
		
		color=(255, 205, 52)
		self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
		self.ren.AddActor(self.fActor)
		self.fOutlineActor, self.fOPC = gen_outline(self.fO,color,self.PointSize)
		self.ren.AddActor(self.fOutlineActor)
		
		#recalculate extents of interactor
		RefMin=np.amin(np.vstack((self.flp,self.rp)),axis=0)
		RefMax=np.amax(np.vstack((self.flp,self.rp)),axis=0)

		extents=RefMax-RefMin #extents
		rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
		self.limits=[RefMin[0]-rl, \
		RefMax[0]+rl, \
		RefMin[1]-rl, \
		RefMax[1]+rl, \
		RefMin[2],RefMax[2]]

		#add axes
		self.add_axis(self.limits,[1,1,1])
		
		s,nl,axs=self.get_scale()

		self.fActor.SetScale(s)
		self.fActor.Modified()

		self.add_axis(nl,axs)
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		if self.mirrored==False:
			self.ui.statusLabel.setText("WARNING alignment proceeding without a mirror operation. Alignment complete. Idle.")
		else:
			self.ui.statusLabel.setText("Alignment complete. Idle.")
		
		self.aligned = True
		
	def get_grid(RefPoints):
	
		p=path.Path(windowVerts)
		inWindow=p.contains_points(RefPoints[:,:2]) #first 2 columns of RefPoints is x and y

		windowed=RefPoints[inWindow,:2]
		gs=squareform(pdist(windowed,'euclidean')) #does the same thing as pdist2
		gsize=np.mean(np.sort(gs)[:,1]) #sort the distances, find the closest, non self-referencing points

		grid_x, grid_y = np.meshgrid(np.linspace(RefMin[0],RefMax[0],int((RefMax[0]-RefMin[0])/gsize)),
	    np.linspace(RefMin[1],RefMax[1],int((RefMax[1]-RefMin[1])/gsize)),indexing='xy')
		points=RefPoints[:,:2]

		grid_RefVal=griddata(points,RefPoints[:,-1], (grid_x, grid_y), method='linear')
		
	def get_input_data(self,filem):
		"""
		Loads the content of a *.mat file pertaining to this particular step
		"""
		
		if hasattr(self,'rActor'): #then remove everything
			self.ren.RemoveActor(self.rActor)
			self.ren.RemoveActor(self.fActor)
			self.ren.RemoveActor(self.rOutlineActor)
			self.ren.RemoveActor(self.fOutlineActor)
			
		if hasattr(self,'aActor'):
			self.ren.RemoveActor(self.aActor)
		
		if filem == None:
			filem, _, =get_file('*.mat')
		
		if filem: #check variables
			mat_contents = sio.loadmat(filem)
			self.fileo=filem
			try:
				self.rp=mat_contents['ref']['rawPnts'][0][0]
				ind=mat_contents['ref']['mask'][0][0][0]
				self.rO=mat_contents['ref']['x_out'][0][0]
				
				self.rp=self.rp[np.where(ind)]
				
				color=(242, 101, 34)
				self.rPC, self.rActor, _, = gen_point_cloud(self.rp,color,self.PointSize)
				self.ren.AddActor(self.rActor)
				self.rOutlineActor, self.rOPC = gen_outline(self.rO,color,self.PointSize)
				self.ren.AddActor(self.rOutlineActor)
				
				#do other one
				self.flp=mat_contents['float']['rawPnts'][0][0]
				ind=mat_contents['float']['mask'][0][0][0]
				self.fO=mat_contents['float']['x_out'][0][0]
				
				self.flp=self.flp[np.where(ind)]

				color=(255, 205, 52)
				self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
				self.ren.AddActor(self.fActor)
				self.fOutlineActor, self.fOPC = gen_outline(self.fO,color,self.PointSize)
				self.ren.AddActor(self.fOutlineActor)
				
				# update status
				self.ui.activeFileLabel.setText("Current analysis file:%s"%filem)
				
				RefMin=np.amin(np.vstack((self.flp,self.rp)),axis=0)
				RefMax=np.amax(np.vstack((self.flp,self.rp)),axis=0)

				extents=RefMax-RefMin #extents
				rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
				self.limits=[RefMin[0]-rl, \
				RefMax[0]+rl, \
				RefMin[1]-rl, \
				RefMax[1]+rl, \
				RefMin[2],RefMax[2]]

				#add axes
				self.add_axis(self.limits,[1,1,1])

			except:
				print "Couldn't read in both sets of data."
			
		else:
			print 'Invalid *.mat file'
			return
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()

	def keypress(self,obj,event):
		key = obj.GetKeyCode()

		if key =="1":
			xyview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="2":
			yzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="3":
			xzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key=="z":
			self.Zaspect=self.Zaspect*2
			s,nl,axs=self.get_scale()
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(s)
				self.pointActor.Modified()
			if hasattr(self,'rActor'):
				self.rActor.SetScale(s)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				self.fActor.SetScale(s)
				self.fActor.Modified()
			if hasattr(self,'aActor')
				self.aActor.SetScale(s)
				self.aActor.Modified()
			
			self.add_axis(nl,axs)



		elif key=="x":
			self.Zaspect=self.Zaspect*0.5
			s,nl,axs=self.get_scale()
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(s)
			if hasattr(self,'rActor'):
				self.rActor.SetScale(s)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				self.fActor.SetScale(s)
				self.fActor.Modified()
			if hasattr(self,'aActor'):
				self.aActor.SetScale(s)
				self.aActor.Modified()

			self.add_axis(nl,axs)


		elif key=="c":
			self.Zaspect=1.0
			s,_,_,=self.get_scale()
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(s)
			if hasattr(self,'rActor'):
				self.rActor.SetScale(s)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				self.fActor.SetScale(s)
				self.fActor.Modified()
			if hasattr(self,'aActor'):
				# self.fActor.SetScale(1,1,self.Zaspect)
				self.aActor.SetScale(s)
				self.aActor.Modified()
			self.add_axis(self.limits,[1,1,1])
			self.ren.ResetCamera()

		elif key=="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(self.ui.vtkWidget._RenderWindow)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("Avg_aligned.png")
			writer.Write()
			print 'Screen output saved to %s' %os.path.join(os.getcwd(),'Avg_aligned.png')

		elif key=="a":
			if hasattr(self,'ax3D'):
				flip_visible(self.ax3D)
			
		elif key == "o":
			if hasattr(self,'outlineActor'):
				flip_visible(self.outlineActor)
		
		elif key == "f":
			if hasattr(self,'ax3D'):
				flip_colors(self.ren,self.ax3D)
				
		elif key == "r":
			if self.picking == True:
				self.picking =False
				self.show_picking()
			else:
				self.picking =True
				self.show_picking()
		
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()

	def get_scale(self):
		'''
		Returns array for the keypress function based on what radio button is selected.
		'''
		if self.ui.xsButton.isChecked():
			s=np.array([self.Zaspect,1,1])
			nl=np.append([self.limits[0]*self.Zaspect,self.limits[1]*self.Zaspect],self.limits[2:])
			axs=np.array([1/self.Zaspect,1,1])
			
		elif self.ui.ysButton.isChecked():
			s=np.array([1,self.Zaspect,1])
			nl=np.append(self.limits[0:2],([self.limits[2]*self.Zaspect,self.limits[3]*self.Zaspect],self.limits[4:]))
			axs=np.array([1,1/self.Zaspect,1])
		else:
			s=np.array([1,1,self.Zaspect])
			nl=np.append(self.limits[0:4],([self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect]))
			axs=np.array([1,1,1/self.Zaspect])
		return s,nl,axs

	def add_axis(self,limits,scale):
		if hasattr(self,"ax3D"):
			self.ren.RemoveActor(self.ax3D)
		self.ax3D = vtk.vtkCubeAxesActor()
		self.ax3D.ZAxisTickVisibilityOn()
		self.ax3D.SetXTitle('X')
		self.ax3D.SetXUnits('mm')
		self.ax3D.SetYTitle('Y')
		self.ax3D.SetYUnits('mm')
		self.ax3D.SetZTitle('Z')
		self.ax3D.SetZUnits('mm')
		self.ax3D.SetBounds(limits)
		self.ax3D.SetZAxisRange(limits[-2]*scale[2],limits[-1]*scale[2])
		self.ax3D.SetXAxisRange(limits[0]*scale[0],limits[1]*scale[0])
		self.ax3D.SetYAxisRange(limits[2]*scale[1],limits[3]*scale[1])
		self.ax3D.SetCamera(self.ren.GetActiveCamera())
		self.ren.AddActor(self.ax3D)
		if not(self.ren.GetBackground()==(0.1, 0.2, 0.4)):
			flip_colors(self.ren,self.ax3D)
		
if __name__ == '__main__':
	if len(sys.argv)>2:
		ali_avg_interactor(sys.argv[1])
	else:
		ali_avg_interactor()