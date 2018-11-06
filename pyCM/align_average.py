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
1.2 - Refactored to use PyQt interface and eliminated global variables
1.3 - Refactored to use PyQt5, Python 3
'''
__author__ = "M.J. Roy"
__version__ = "1.3"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2018"

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
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from .pyCMcommon import *
from pkg_resources import Requirement, resource_filename



def aa_def(*args,**kwargs):
	"""
	Main function, builds qt interaction
	"""	
	app = QtWidgets.QApplication.instance()
	if app is None:
		app = QtWidgets.QApplication(sys.argv)
	
	spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
	splash_pix = QtGui.QPixmap(spl_fname,'PNG')
	splash = QtWidgets.QSplashScreen(splash_pix)
	splash.setMask(splash_pix.mask())

	splash.show()
	app.processEvents()
	
	window = aa_interactor(None)

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
		MainWindow.setWindowTitle("pyCM - Alignment and averaging tool v%s" %__version__)
		self.centralWidget = QtWidgets.QWidget(MainWindow)
		if hasattr(MainWindow,'setCentralWidget'):
			MainWindow.setCentralWidget(self.centralWidget)
		else:
			self.centralWidget=MainWindow
		self.mainlayout=QtWidgets.QGridLayout(self.centralWidget)
		
		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
		
		mainUiBox = QtWidgets.QGridLayout()
		
		self.vtkWidget.setMinimumSize(QtCore.QSize(1050, 600))
		sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
		sizePolicy.setHorizontalStretch(10)
		sizePolicy.setVerticalStretch(10)
		sizePolicy.setHeightForWidth(self.vtkWidget.sizePolicy().hasHeightForWidth())
		self.vtkWidget.setSizePolicy(sizePolicy)
		
		self.statLabel=QtWidgets.QLabel("Idle")
		self.statLabel.setWordWrap(True)
		self.statLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.statLabel.setMinimumWidth(100)
		sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
		sizePolicy.setHorizontalStretch(0)
		sizePolicy.setVerticalStretch(0)
		sizePolicy.setHeightForWidth(self.statLabel.sizePolicy().hasHeightForWidth())
		self.statLabel.setSizePolicy(sizePolicy)

		headFont=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)

		
		# #define buttons/widgets
		# self.reloadButton = QtWidgets.QPushButton('Load')
		scalingLabel=QtWidgets.QLabel("Active axis for scaling")
		scalingLabel.setFont(headFont)
		self.xsButton=QtWidgets.QRadioButton("x")
		self.ysButton=QtWidgets.QRadioButton("y")
		self.zsButton=QtWidgets.QRadioButton("z")
		self.zsButton.setChecked(True)
		self.scalingButtonGroup = QtWidgets.QButtonGroup()
		self.scalingButtonGroup.addButton(self.xsButton)
		self.scalingButtonGroup.addButton(self.ysButton)
		self.scalingButtonGroup.addButton(self.zsButton)
		self.scalingButtonGroup.setExclusive(True)
		scaleBoxlayout = QtWidgets.QGridLayout()
		scaleBoxlayout.addWidget(self.xsButton,1,1)
		scaleBoxlayout.addWidget(self.ysButton,1,2)
		scaleBoxlayout.addWidget(self.zsButton,1,3)
		
		horizLine1=QtWidgets.QFrame()
		horizLine1.setFrameStyle(QtWidgets.QFrame.HLine)
		mirrorLabel=QtWidgets.QLabel("Mirroring")
		mirrorLabel.setFont(headFont)
		self.mirrorXbutton = QtWidgets.QPushButton('ZY')
		self.mirrorYbutton = QtWidgets.QPushButton('ZX')

		horizLine2=QtWidgets.QFrame()
		horizLine2.setFrameStyle(QtWidgets.QFrame.HLine)
		alignLabel=QtWidgets.QLabel("Alignment")
		alignLabel.setFont(headFont)
		self.transXlabel=QtWidgets.QLabel("Translate x:")
		self.transX = QtWidgets.QLineEdit()
		self.transX.setText('0')
		self.transYlabel=QtWidgets.QLabel("Translate y:")
		self.transY = QtWidgets.QLineEdit()
		self.transY.setText('0')
		self.transButton=QtWidgets.QPushButton('Translate floating')
		
		self.alignButtonGroup = QtWidgets.QButtonGroup()
		self.alignOutlineButton=QtWidgets.QRadioButton("Outline")
		self.alignPointCloudButton = QtWidgets.QRadioButton("Point cloud")
		self.alignOutlineButton.setChecked(True)
		self.alignButtonGroup.addButton(self.alignOutlineButton)
		self.alignButtonGroup.addButton(self.alignPointCloudButton)
		self.alignButtonGroup.setExclusive(True)
		self.alignButton = QtWidgets.QPushButton('Align')
		self.alignButton.setStyleSheet("background-color : None ")
		
		horizLine3=QtWidgets.QFrame()
		horizLine3.setFrameStyle(QtWidgets.QFrame.HLine)
		averageLabel=QtWidgets.QLabel("Averaging")
		averageLabel.setFont(headFont)
		self.averageButton = QtWidgets.QPushButton('Average')
		self.averageButton.setStyleSheet("background-color : None ")
		
		horizLine4=QtWidgets.QFrame()
		horizLine4.setFrameStyle(QtWidgets.QFrame.HLine)
		self.writeButton=QtWidgets.QPushButton('Write')
		
		horizLine5=QtWidgets.QFrame()
		horizLine5.setFrameStyle(QtWidgets.QFrame.HLine)
		# self.statusLabel=QtWidgets.QLabel("Idle")
		# self.statusLabel.setWordWrap(True)
		# self.statusLabel.setFont(QtGui.QFont("Helvetica",italic=True))


		#add widgets to ui
		# mainUiBox.addWidget(self.reloadButton,0,0,1,2)
		mainUiBox.addWidget(scalingLabel,0,0,1,2)
		mainUiBox.addLayout(scaleBoxlayout,1,0,1,2)
		mainUiBox.addWidget(horizLine1,2,0,1,2)
		mainUiBox.addWidget(mirrorLabel,3,0,1,2)
		mainUiBox.addWidget(self.mirrorXbutton,4,0,1,1)
		mainUiBox.addWidget(self.mirrorYbutton,4,1,1,1)
		mainUiBox.addWidget(horizLine2,5,0,1,2)
		mainUiBox.addWidget(alignLabel,6,0,1,2)
		mainUiBox.addWidget(self.transXlabel,7,0,1,1)
		mainUiBox.addWidget(self.transX,7,1,1,1)
		mainUiBox.addWidget(self.transYlabel,8,0,1,1)
		mainUiBox.addWidget(self.transY,8,1,1,1)
		mainUiBox.addWidget(self.transButton,9,0,1,2)
		mainUiBox.addWidget(self.alignOutlineButton,10,0,1,1)
		mainUiBox.addWidget(self.alignPointCloudButton,10,1,1,1)
		mainUiBox.addWidget(self.alignButton,11,0,1,2)
		mainUiBox.addWidget(horizLine3,12,0,1,2)
		mainUiBox.addWidget(averageLabel,13,0,1,2)
		mainUiBox.addWidget(self.averageButton,14,0,1,2)
		mainUiBox.addWidget(horizLine4,15,0,1,2)
		mainUiBox.addWidget(self.writeButton,16,0,1,2)
		mainUiBox.addWidget(horizLine5,17,0,1,2)
		# mainUiBox.addWidget(self.statusLabel,18,0,1,2)

		mainUiBox.setColumnMinimumWidth(0,mainUiBox.columnMinimumWidth(0))
		mainUiBox.setColumnMinimumWidth(1,mainUiBox.columnMinimumWidth(0))
		
		lvLayout=QtWidgets.QVBoxLayout()
		lhLayout=QtWidgets.QHBoxLayout()
		lvLayout.addLayout(mainUiBox)
		lvLayout.addStretch(1)
		lhLayout.addLayout(lvLayout)
		lhLayout.addStretch(2)



		self.mainlayout.addWidget(self.vtkWidget,0,0,1,1)
		self.mainlayout.addLayout(lhLayout,0,1,1,1)
		self.mainlayout.addWidget(self.statLabel,1,0,1,2)
		
	def initialize(self):
		self.vtkWidget.start()

class aa_interactor(QtWidgets.QWidget):
	"""
	Sets up the main VTK window, reads file and sets connections between UI and interactor
	"""
	def __init__(self, parent):
		super(aa_interactor,self).__init__(parent)
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
		
		# self.ui.reloadButton.clicked.connect(lambda: self.get_input_data(None))
		self.ui.mirrorXbutton.clicked.connect(lambda: self.flipside('x'))
		self.ui.mirrorYbutton.clicked.connect(lambda: self.flipside('y'))
		self.ui.transButton.clicked.connect(lambda: self.shift())
		self.ui.alignButton.clicked.connect(lambda: self.align())
		self.ui.averageButton.clicked.connect(lambda: self.average())
		self.ui.writeButton.clicked.connect(lambda: self.write())
	
	def shift(self):
	
		self.unsaved_changes=True
		if hasattr(self,'fActor'): #then remove this actor and the associated outline actor
			self.ren.RemoveActor(self.fActor)
			self.ren.RemoveActor(self.fOutlineActor)

		#get x and y tranformations
		gx=float(self.ui.transX.text())
		gy=float(self.ui.transY.text())
		self.user_transl=np.array([gx, gy, 0]);
		
		#apply operation
		self.flp=self.flp+self.user_transl
		self.fO=self.fO+self.user_transl
		
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
		
		mat_vars=sio.whosmat(self.fileo)
		if not set(['transM', 'aa', 'mirror']).isdisjoint([item for sublist in mat_vars for item in sublist]): #tell the user that they might overwrite their data
			ret=QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
				"There is already data associated with this analysis step saved. Overwrite and invalidate subsequent steps?", \
				QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
			if ret == QtWidgets.QMessageBox.No: #don't overwrite
				return
			else:
				#delete fitting parameters with pyCMcommon helper function, which negates FEA pre-processing as well.
				clear_mat(self.fileo,['x_out','aa_mask','spline_x']) 
					
		mat_contents=sio.loadmat(self.fileo)
		
		new={'transM':self.transM,'aa':self.ap, 'mirror':self.mirror_plane}
		
		mat_contents.update(new) #update the dictionary
			
		sio.savemat(self.fileo,mat_contents)	
		self.ui.statLabel.setText("Wrote data.")
		self.unsaved_changes=False
		
	def average(self):
		
		self.unsaved_changes=True
		
		self.ui.statLabel.setText("Averaging, applying grid . . .")
		QtWidgets.QApplication.processEvents()
		
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
		
		self.ui.statLabel.setText("Averaging using grid . . .")
		QtWidgets.QApplication.processEvents()
		
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
		
		self.ui.statLabel.setText("Rendering . . .")
		QtWidgets.QApplication.processEvents()
		
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
		self.ui.statLabel.setText("Averaging complete.")
		self.averaged=True
		self.ui.averageButton.setStyleSheet("background-color :rgb(77, 209, 97);")
		
	
	def flipside(self,flipDirection):
		self.ui.statLabel.setText("Starting mirroring . . .")
		self.ui.vtkWidget.update()
		#delete the floating actor
		if hasattr(self,'fActor'): #then remove this actor and the associated outline actor
			self.ren.RemoveActor(self.fActor)
			self.ren.RemoveActor(self.fOutlineActor)


		if flipDirection == "x":
			self.flp[:,0]=-self.flp[:,0]
			self.fO[:,0]=-self.fO[:,0]
			self.mirror_plane="x"
			self.ui.mirrorXbutton.setStyleSheet("background-color :rgb(77, 209, 97);")
			self.ui.mirrorYbutton.setStyleSheet("background-color: None")
		elif flipDirection == "y":
			self.flp[:,1]=-self.flp[:,1]
			self.fO[:,1]=-self.fO[:,1]
			self.mirror_plane="y"
			self.ui.mirrorYbutton.setStyleSheet("background-color :rgb(77, 209, 97);")
			self.ui.mirrorXbutton.setStyleSheet("background-color: None")
		
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
		self.ui.statLabel.setText("Mirror operation complete.")
		self.mirrored=True
		
	def align(self):
		self.unchanged_changes=True
		self.ui.statLabel.setText("Starting alignment . . .")
		QtWidgets.QApplication.processEvents()
		
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
		
		# make sure to remove any user-defined translation from transM
		if hasattr(self,'user_transl'):
			self.transM[0:3,-1]=self.transM[0:3,-1]+np.asmatrix(self.user_transl)

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
			self.ui.statLabel.setText("WARNING alignment proceeding without a mirror operation. Alignment complete.")
		else:
			self.ui.statLabel.setText("Alignment complete.")
		
		self.aligned = True
		self.ui.alignButton.setStyleSheet("background-color :rgb(77, 209, 97);")
		
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
			self.user_transl=np.array([0,0,0]) #needs to be reset
			
		if hasattr(self,'aActor'):
			self.ren.RemoveActor(self.aActor)
		
		if filem == None:
			filem, _, =get_file('*.mat')
		
		if filem: #check variables
			mat_contents = sio.loadmat(filem)
			self.fileo=filem
			if 'aa' in mat_contents:

				
				#draw floating and reference datasets
				
				self.rp=mat_contents['ref']['rawPnts'][0][0]
				ind=mat_contents['ref']['mask'][0][0][0]
				self.rO=mat_contents['ref']['x_out'][0][0]
				
				self.rp=self.rp[np.where(ind)]
				
				color=(242, 101, 34)
				self.rPC, self.rActor, _, = gen_point_cloud(self.rp,color,self.PointSize)
				self.ren.AddActor(self.rActor)
				self.rOutlineActor, self.rOPC = gen_outline(self.rO,color,self.PointSize)
				self.ren.AddActor(self.rOutlineActor)
				
				s,nl,axs=self.get_scale()
				
				self.rActor.SetScale(s)
				self.rActor.Modified()
				
				#do other one, but with transformed floating points
				self.flp=mat_contents['float']['rawPnts'][0][0]
				ind=mat_contents['float']['mask'][0][0][0]
				self.fO=mat_contents['float']['x_out'][0][0]
				
				self.flp=self.flp[np.where(ind)]
				self.flipside(mat_contents['mirror'])
				# flipping creates actors, so remove the actors immediately before transforming
				self.ren.RemoveActor(self.fActor)
				self.ren.RemoveActor(self.fOutlineActor)

				self.transM=mat_contents['transM']
				self.flp=np.dot(self.flp,self.transM[0:3,0:3])+self.transM[0:3,-1]
				self.fO=np.dot(self.fO,self.transM[0:3,0:3])+self.transM[0:3,-1]
				
				color=(255, 205, 52)
				self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
				self.ren.AddActor(self.fActor)
				self.fOutlineActor, self.fOPC = gen_outline(self.fO,color,self.PointSize)
				self.ren.AddActor(self.fOutlineActor)
				
				#show aligned and averaged data
				self.ap=mat_contents['aa']
				color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
				_, self.aActor, _, = gen_point_cloud(self.ap,color,self.PointSize)
				self.ren.AddActor(self.aActor)
		
				

				self.aActor.SetScale(s)
				self.aActor.Modified()
				

				
				self.ui.statLabel.setText("This dataset has already been aligned and averaged.")
				self.aligned = True
				self.ui.alignButton.setStyleSheet("background-color :rgb(77, 209, 97);")
				self.mirrored=True
				self.averaged=True
				self.ui.averageButton.setStyleSheet("background-color :rgb(77, 209, 97);")
			else:
				self.ui.statLabel.setText("This dataset has not been previously aligned.")
				# QtWidgets.QApplication.processEvents()

				#'''
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
					print("Couldn't read in both sets of data.")
				
		else:
			print("Invalid *.mat file")
			return
		self.unsaved_changes=False
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
			if hasattr(self,'aActor'):
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
			print("Screen output saved to %s" %os.path.join(os.getcwd(),'Avg_aligned.png'))

		elif key=="a":
			if hasattr(self,'ax3D'):
				flip_visible(self.ax3D)
			
		elif key == "o":
			if hasattr(self,'outlineActor'):
				flip_visible(self.outlineActor)
		
		elif key == "f":
			if hasattr(self,'ax3D'):
				flip_colors(self.ren,self.ax3D)
				
				
		elif key=="l":
			self.get_input_data(None)
		
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
		aa_def(sys.argv[1])
	else:
		aa_def()