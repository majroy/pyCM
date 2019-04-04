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
r - enter/exit picking mode, LMB is used to generate a selection window. Exiting 
	picking mode will highlight selected points.
z - increase aspect ratio
x - decrease aspect ratio
c - return to default aspect
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
a - toggles axes
o - toggles outline (if present)
r - starts picking
-------------------------------------------------------------------------------
1.1 - Fixed array orientation, clipping issue, compass scaling and sped up writing output
      Added ReadMask
1.2 - Fixed window handling, now exits cleanly
1.3 - Modified to run in Python 3.x, uses VTK keyboard interrupts to start picking, Qt button for this function has been commented out.
1.4 - Added the ability to 'level' incoming data based on AFRC input
1.5 - Added SVD analysis/transformations
1.6 - Added ability to read PC-DMIS csv files
1.7 - Added directions to SVD decomposition.
'''
__author__ = "M.J. Roy"
__version__ = "1.7"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2019"

import sys
import os.path
from pkg_resources import Requirement, resource_filename
import numpy as np
import scipy.io as sio
import vtk
import vtk.util.numpy_support as vtk_to_numpy
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from pyCM.pyCMcommon import *

nosio=False

def mask_def(*args,**kwargs):
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
	
	window = pnt_interactor(None)
	if len(args)==2: 
		pnt_interactor.get_input_data(window,args[0],args[1])
	elif len(args)==1: 
		pnt_interactor.get_input_data(window,None,args[0])
	else: 
		pnt_interactor.get_input_data(window,None,None)


	window.show()
	splash.finish(window)
	window.iren.Initialize() # Need this line to actually show the render inside Qt

	ret = app.exec_()
	
	if sys.stdin.isatty() and not hasattr(sys,'ps1'):
		sys.exit(ret)
	else:
		return window

		
class pt_main_window(object):
	"""
	Class to build qt interaction, including VTK widget
	setupUi builds, initialize starts VTK widget
	"""

	def setupUi(self, MainWindow):
		MainWindow.setWindowTitle(("pyCM - Point editor v%s" %__version__))
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
		
		
		#define buttons/widgets
		self.reloadButton = QtWidgets.QPushButton('New profile')
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
		self.levelButton=QtWidgets.QRadioButton("Translate to mean z value")
		svdLabel=QtWidgets.QLabel("Perform SVD reorientation")
		svdLabel.setFont(headFont)
				
		self.rxButton_pos=QtWidgets.QRadioButton("Rx+")
		self.ryButton_pos=QtWidgets.QRadioButton("Ry+")
		self.rxButton_neg=QtWidgets.QRadioButton("Rx-")
		self.ryButton_neg=QtWidgets.QRadioButton("Ry-")
		
		svdButtonGroup = QtWidgets.QButtonGroup()
		svdButtonGroup.addButton(self.rxButton_pos)
		svdButtonGroup.addButton(self.ryButton_pos)
		svdButtonGroup.addButton(self.rxButton_neg)
		svdButtonGroup.addButton(self.ryButton_neg)
		svdButtonGroup.setExclusive(False)
		svdBoxlayout = QtWidgets.QGridLayout()
		svdBoxlayout.addWidget(self.rxButton_pos,1,1)
		svdBoxlayout.addWidget(self.rxButton_neg,1,2)
		svdBoxlayout.addWidget(self.ryButton_pos,1,3)
		svdBoxlayout.addWidget(self.ryButton_neg,1,4)
		
		self.reduce = QtWidgets.QSpinBox()
		self.reduce.setValue(0)
		self.reduceButton = QtWidgets.QPushButton('Reduce')
		self.revertButton = QtWidgets.QPushButton('Undo all/reload')

		reduceBoxlayout= QtWidgets.QGridLayout()
		reduceBoxlayout.addWidget(self.reduce,1,1)
		reduceBoxlayout.addWidget(self.reduceButton,1,2)
		
		horizLine1=QtWidgets.QFrame()
		horizLine1.setFrameStyle(QtWidgets.QFrame.HLine)
		pickLabel=QtWidgets.QLabel("Pick options")
		pickLabel.setFont(headFont)

		self.pickHelpLabel=QtWidgets.QLabel("Press R to activate")
		self.pickActiveLabel=QtWidgets.QLabel("Pick active")
		self.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }")
		self.pickActiveLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.undoLastPickButton=QtWidgets.QPushButton('Undo last pick')
		horizLine2=QtWidgets.QFrame()
		horizLine2.setFrameStyle(QtWidgets.QFrame.HLine)
		outputLabel=QtWidgets.QLabel("Write output")
		outputLabel.setFont(headFont)
		self.refButton=QtWidgets.QRadioButton("Reference")
		
		self.floatButton=QtWidgets.QRadioButton("Floating")
		
		self.refButton.setChecked(True)
		self.writeButtonGroup = QtWidgets.QButtonGroup()
		self.writeButtonGroup.addButton(self.floatButton)
		self.writeButtonGroup.addButton(self.refButton)
		self.writeButtonGroup.setExclusive(True)
		self.writeButton=QtWidgets.QPushButton('Write')
		horizLine3=QtWidgets.QFrame()
		horizLine3.setFrameStyle(QtWidgets.QFrame.HLine)
		showLabel=QtWidgets.QLabel("Load result")
		showLabel.setFont(headFont)
		self.showRefButton=QtWidgets.QRadioButton("Reference")
		self.showRefButton.setChecked(True)
		self.showFloatButton=QtWidgets.QRadioButton("Floating")
		self.showButtonGroup = QtWidgets.QButtonGroup()
		self.showButtonGroup.addButton(self.showFloatButton)
		self.showButtonGroup.addButton(self.showRefButton)
		self.showButtonGroup.setExclusive(True)
		
		self.showButton=QtWidgets.QPushButton("View")
		horizLine4=QtWidgets.QFrame()
		horizLine4.setFrameStyle(QtWidgets.QFrame.HLine)

		#add widgets to ui
		mainUiBox.addWidget(self.reloadButton,0,0,1,2)
		mainUiBox.addWidget(scalingLabel,1,0,1,2)
		mainUiBox.addLayout(scaleBoxlayout,2,0,1,2)
		mainUiBox.addWidget(self.levelButton,3,0,1,2)

		mainUiBox.addWidget(horizLine1,4,0,1,2)
		mainUiBox.addWidget(pickLabel,5,0,1,2)
		mainUiBox.addLayout(reduceBoxlayout,6,0,1,2)
		mainUiBox.addWidget(self.pickHelpLabel,7,0,1,1)
		mainUiBox.addWidget(self.pickActiveLabel,7,1,1,1)
		mainUiBox.addWidget(self.undoLastPickButton,8,0,1,2)
		mainUiBox.addWidget(svdLabel,9,0,1,2)
		mainUiBox.addLayout(svdBoxlayout,10,0,1,2)
		mainUiBox.addWidget(self.revertButton,11,0,1,2)
		mainUiBox.addWidget(horizLine2,12,0,1,2)
		mainUiBox.addWidget(outputLabel,13,0,1,2)
		mainUiBox.addWidget(self.refButton,14,0,1,1)
		mainUiBox.addWidget(self.floatButton,14,1,1,1)
		mainUiBox.addWidget(self.writeButton,15,0,1,2)
		mainUiBox.addWidget(horizLine3,16,0,1,2)
		mainUiBox.addWidget(showLabel,17,0,1,2)
		mainUiBox.addWidget(self.showRefButton,18,0,1,1)
		mainUiBox.addWidget(self.showFloatButton,18,1,1,1)
		mainUiBox.addWidget(self.showButton,19,0,1,2)
		mainUiBox.addWidget(horizLine4,20,0,1,2)

		lvLayout=QtWidgets.QVBoxLayout()
		lvLayout.addLayout(mainUiBox)
		lvLayout.addStretch(1)
		
		self.mainlayout.addWidget(self.vtkWidget,0,0,1,1)
		self.mainlayout.addLayout(lvLayout,0,1,1,1)
		self.mainlayout.addWidget(self.statLabel,1,0,1,2)
		
	def initialize(self):
		self.vtkWidget.start()



class pnt_interactor(QtWidgets.QWidget):

	def __init__(self, parent):
		super(pnt_interactor,self).__init__(parent)
		self.ui = pt_main_window()
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
		self.refWritten = False
		self.floatWritten = False
		
		self.ui.reloadButton.clicked.connect(lambda: self.get_input_data(None,None))
		self.ui.undoLastPickButton.clicked.connect(lambda: self.undo_pick())
		self.ui.writeButton.clicked.connect(lambda: self.write_new())
		self.ui.revertButton.clicked.connect(lambda: self.undo_revert())
		self.ui.reduceButton.clicked.connect(lambda: self.reduce_pnts())
		self.ui.levelButton.clicked.connect(lambda: self.level_pnts())
		self.ui.rxButton_pos.clicked.connect(lambda: self.svd('x',False))
		self.ui.ryButton_pos.clicked.connect(lambda: self.svd('y',False))
		self.ui.rxButton_neg.clicked.connect(lambda: self.svd('x',True))
		self.ui.ryButton_neg.clicked.connect(lambda: self.svd('y',True))
		self.ui.showButton.clicked.connect(lambda: self.load_mat())
	
	def svd(self,dir,reverse):
		'''
		Moves point cloud and outline to the centroid of the point cloud, finds SVD difference between X & Y axes of masked point cloud, and applies transformation, and then moves it back to the starting point.
		'''
		
		color=(70, 171, 176)
		
		self.ren.RemoveActor(self.pointActor)
		self.ren.RemoveActor(self.outlineActor)
		
		#then move all points to have centroid at x,y=0
		#get translation vector
		t=np.mean(self.rawPnts,axis=0)
		RP=self.rawPnts
		RP[:,0]=RP[:,0]-t[0]
		RP[:,1]=RP[:,1]-t[1]
		RP[:,2]=RP[:,2]-t[2]
		
		OP=self.Outline
		OP[:,0]=OP[:,0]-t[0]
		OP[:,1]=OP[:,1]-t[1]
		OP[:,2]=OP[:,2]-t[2]
		
		#debug
		# _,_,vh = np.linalg.svd(RP) #vh is transpose from MATLAB's svd, returns normalised vectors
		# #rows of vh are orthnormal vectors
		# # print('X:',vh[0,:] / np.linalg.norm(vh[0,:]))
		# # print('Y:',vh[1,:] / np.linalg.norm(vh[1,:]))
		# # print('Z:',vh[2,:] / np.linalg.norm(vh[2,:]))
		
		# #handles the case if the dataset is net convex vs. concave
		# if vh[2,-1]<0:
			# c=np.array([0,0,-1])
		# else: 
			# c=np.array([0,0,1])
		
		# vh_y_norm = np.array([vh[2,0],0,vh[2,2]]) / np.linalg.norm(np.array([vh[2,0],0,vh[2,2]])) #xz plane projection
		# vh_x_norm = np.array([0,vh[2,1],vh[2,2]]) / np.linalg.norm(np.array([0,vh[2,1],vh[2,2]])) #yz plane projection
		
		# #solve for angle, update console
		# a_y=np.arccos(np.clip(np.dot(vh_y_norm,c), -1.0, 1.0))
		# a_x=np.arccos(np.clip(np.dot(vh_x_norm,c), -1.0, 1.0))
		# print('SVD difference about X and Y axis in degrees prior to transform:\n'a_x*57.3,a_y*57.3)
		
		# Ry=np.matrix([[np.cos(-a_y),0,np.sin(-a_y)],[0,1,0],[-np.sin(-a_y),0,np.cos(-a_y)]])
		# Rx=np.matrix([[1,0,0],[0,np.cos(-a_x),-np.sin(-a_x)],[0,np.sin(-a_x),np.cos(-a_x)]])
	
		#debug
		# if hasattr(self,'svd_arrow_actor'):
			# self.ren.RemoveActor(self.svd_arrow_actor)
			# self.ren.RemoveActor(self.ref1_arrow_actor)
			# self.ren.RemoveActor(self.ref2_arrow_actor)
		#arrow size is 10% max size of domain
		# asize=np.maximum(self.limits[1]-self.limits[0],self.limits[3]-self.limits[2])*0.10
		# self.svd_arrow_actor=draw_arrow(t,asize,-vh[2,:],self.ren,False,(1,0,0))
		# self.ref1_arrow_actor=draw_arrow(t,asize,-vh[0,:],self.ren,False,(0,1,0)) #xaxis, green
		# self.ref2_arrow_actor=draw_arrow(t,asize,-vh[1,:],self.ren,False,(0,0,3)) #yaxis, blue
		
		#find rotation and pickup which rotation to apply based on masked points
		print('Before SVD:')
		Rx0,Ry0=get_svd_rotation_matrix(RP[self.bool_pnt,:])
		
		if reverse:
			Rx0,Ry0=np.linalg.inv(Rx0),np.linalg.inv(Ry0)
		
		if dir == 'y':
			RP = Ry0*RP.T
			OP = Ry0*OP.T
		else: 
			RP = Rx0*RP.T
			OP = Ry0*OP.T
			
		RP = RP.T		
		OP = OP.T
		
		#check rotation
		print('After SVD:')
		Rx1,Ry1=get_svd_rotation_matrix(RP[self.bool_pnt,:])

		# #add translation back on
		RP[:,0]=RP[:,0]+t[0]
		RP[:,1]=RP[:,1]+t[1]
		RP[:,2]=RP[:,2]+t[2]
		
		OP[:,0]=OP[:,0]+t[0]
		OP[:,1]=OP[:,1]+t[1]
		OP[:,2]=OP[:,2]+t[2]

		#update status UI
		if np.allclose(Rx1,np.eye(3)) and np.allclose(Ry1,np.eye(3)):
			#returned identity matrix and therefore 'aligned'
			self.ui.statLabel.setText("SVD completed. See console for results.")
			
		#update everything
		self.rawPnts = np.asarray(RP)
		self.Outline = np.asarray(OP)
		
		#update both outline and actors
		self.vtkPntsPolyData, \
		self.pointActor, self.colors = \
		gen_point_cloud(self.rawPnts,color,self.PointSize)
		
		#modify point coloration based on mask
		#find points to be painted red
		localind=np.asarray(range(len(self.bool_pnt)))
		localind=localind[np.where(np.logical_not(self.bool_pnt))]
		
		for i in localind:
			#turn them red
			self.colors.SetTuple(i,(255,0,0))

		self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
		self.vtkPntsPolyData.Modified()
		self.ren.AddActor(self.pointActor)
		self.ren.AddActor(self.outlineActor)
		
		
		s,nl,axs=self.get_scale()

		self.pointActor.SetScale(s)
		self.outlineActor.SetScale(s)
		
		self.pointActor.Modified()
		self.outlineActor.Modified()

		self.add_axis(nl,axs)
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
	
		self.unsaved_changes=True
	
	def undo_revert(self):
		'''
		Reloads all data based on filec & filep (if it exists), will re-initialize data read in from results file to be unmasked.
		'''

		try:
			self.get_input_data(self.filep,self.filec)
			self.unsaved_changes=True
		except: #its been loaded from an existing results file
			ret=QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
			"Existing mask of profile will be lost, continue?", \
			QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
			if ret == QtWidgets.QMessageBox.No: #don't overwrite
				return
			else:
				#flip all values in bool_pnt & update color
				localind=np.asarray(range(len(self.bool_pnt)))
				localind=localind[np.where(np.logical_not(self.bool_pnt))]
				
				for i in localind:
					#show them as being unmasked
					self.colors.SetTuple(i,(70, 171, 176))
				self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
				self.vtkPntsPolyData.Modified()
				self.ui.vtkWidget.update()
				self.ui.vtkWidget.setFocus()
				
				#re-initialise the mask
				self.bool_pnt=np.ones(self.bool_pnt.shape,dtype='bool')
				#set flag on ui to show that data has been modified
				self.unsaved_changes=True

	def level_pnts(self):
		'''
		Translates outline and profile by the mean of z so that scaling occurs about 0.
		'''
		color=(70, 171, 176)
		self.ren.RemoveActor(self.pointActor)
		self.ren.RemoveActor(self.outlineActor)
		#adjust to z mean of outline
		self.Outline[:,2]=self.Outline[:,2]-np.mean(self.Outline[:,2])

		#adjust to z mean of point cloud
		self.rawPnts[:,2]=self.rawPnts[:,2]-np.mean(self.rawPnts[:,2])
		self.outlineActor, _ =gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
		
		#get limits
		try:
			RefMin=np.amin(np.hstack((self.Outline,self.rawPnts)),axis=0)
			RefMax=np.amax(np.hstack((self.Outline,self.rawPnts)),axis=0)
		except:
			RefMin=np.amin(self.rawPnts,axis=0)
			RefMax=np.amax(self.rawPnts,axis=0)
		
		extents=RefMax-RefMin #extents
		rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
		self.limits=[RefMin[0]-rl, \
		  RefMax[0]+rl, \
		  RefMin[1]-rl, \
		  RefMax[1]+rl, \
		  RefMin[2],RefMax[2]]

		#add axes
		self.add_axis(self.limits,[1,1,1])
		
		
		self.vtkPntsPolyData, \
		self.pointActor, self.colors = \
		gen_point_cloud(self.rawPnts,color,self.PointSize)
		
		self.ren.AddActor(self.pointActor)
		self.ren.AddActor(self.outlineActor)
		
		self.pointActor.Modified()
		self.outlineActor.Modified()
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		
		
		
				
	def reduce_pnts(self):
		'''
		Reduces the number of points according to the percentage of what's in the spinbox
		0 -> means nothing, 10 means leave 90 percent of the points
		'''
		color=(70, 171, 176)
		red=(100-float(self.ui.reduce.value()))/100
		ind=np.linspace(0, len(self.rawPnts)-1, num=int(red*len(self.rawPnts)))
		self.rawPnts=self.rawPnts[tuple(ind.astype(int)),:]

		self.ren.RemoveActor(self.pointActor)
		self.vtkPntsPolyData, \
		self.pointActor, self.colors = \
		gen_point_cloud(self.rawPnts,color,self.PointSize)
		self.bool_pnt=np.ones(len(self.rawPnts), dtype=bool)
		self.ren.AddActor(self.pointActor)
				
		s,nl,axs=self.get_scale()

		self.pointActor.SetScale(s)
		self.pointActor.Modified()

		self.add_axis(nl,axs)
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
	
	def load_mat(self):
		"""
		Loads the content of a *.mat file pertaining to this particular step
		"""
		
		color=(70, 171, 176)
		
		if self.ui.showRefButton.isChecked():
			str_d='ref'

		if self.ui.showFloatButton.isChecked():
			str_d='float'
			
			
		if hasattr(self,'pointActor'):
			self.ren.RemoveActor(self.pointActor)
		if hasattr(self,'outlineActor'):
			self.ren.RemoveActor(self.outlineActor)
		
		if hasattr(self,'fileo'): #check variables
			if self.fileo == None:
				return
			mat_contents = sio.loadmat(self.fileo)
			#check contents
			if 'ref' in mat_contents: 	
				self.ui.refButton.setStyleSheet("background-color :rgb(77, 209, 97);")
				self.refWritten = True
			if 'float' in mat_contents: 	
				self.ui.floatButton.setStyleSheet("background-color :rgb(77, 209, 97);")
				self.floatWritten = True
			try:
				self.rawPnts=mat_contents[str_d]['rawPnts'][0][0]
				self.bool_pnt=mat_contents[str_d]['mask'][0][0][0]
				self.Outline=mat_contents[str_d]['x_out'][0][0]
				
				self.outlineActor, _ =gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
				self.ren.AddActor(self.outlineActor)
			
				self.vtkPntsPolyData, \
				self.pointActor, self.colors = \
				gen_point_cloud(self.rawPnts,color,self.PointSize)
				
				#find points to be painted red
				localind=np.asarray(range(len(self.bool_pnt)))
				localind=localind[np.where(np.logical_not(self.bool_pnt))]
				
				for i in localind:
					#turn them red
					self.colors.SetTuple(i,(255,0,0))

				self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
				self.vtkPntsPolyData.Modified()
				
				self.ren.AddActor(self.pointActor)
				
			except:
				QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
				"The %s dataset could not be loaded."%(str_d))
			
			#get limits
			try:
				RefMin=np.amin(np.hstack((self.Outline,self.rawPnts)),axis=0)
				RefMax=np.amax(np.hstack((self.Outline,self.rawPnts)),axis=0)
			except:
				RefMin=np.amin(self.rawPnts,axis=0)
				RefMax=np.amax(self.rawPnts,axis=0)
			
			extents=RefMax-RefMin #extents
			rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
			self.limits=[RefMin[0]-rl, \
			  RefMax[0]+rl, \
			  RefMin[1]-rl, \
			  RefMax[1]+rl, \
			  RefMin[2],RefMax[2]]

			#add axes
			self.add_axis(self.limits,[1,1,1])


			

		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		

	def write_new(self):
		
		if self.ui.refButton.isChecked():
			str_d='ref'
			self.refWritten=True
			
		if self.ui.floatButton.isChecked():
			str_d='float'
			self.floatWritten=True

		
		if not hasattr(self,'fileo'):
			self.fileo, _, = get_open_file('*.mat',os.getcwd())
			if self.fileo:
				x_o=self.rawPnts[self.bool_pnt,0]
				y_o=self.rawPnts[self.bool_pnt,1]
				z_o=self.rawPnts[self.bool_pnt,2]
				sio.savemat(self.fileo,{str_d : {'x_out':self.Outline,'rawPnts':self.rawPnts,'mask': self.bool_pnt,'x':x_o,'y':y_o,'z':z_o,'fname':self.filec}})
				if self.ui.refButton.isChecked():
					self.ui.refButton.setStyleSheet("background-color :rgb(77, 209, 97);")
			
				if self.ui.floatButton.isChecked():
					self.ui.floatButton.setStyleSheet("background-color : rgb(77, 209, 97);")
				#reset flag on ui to show that data has been modified
				
		else:
			if not self.fileo:
				self.fileo, _, = get_open_file('*.mat',os.getcwd())

			mat_vars=sio.whosmat(self.fileo)
			if str_d in [item for sublist in mat_vars for item in sublist]: #tell the user that they might overwrite their data
				ret=QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
				"There is already data for this step - doing this will invalidate all further existing analysis steps. Continue?"%(str_d), \
				QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
				if ret == QtWidgets.QMessageBox.No: #don't overwrite
					return
			
			
			mat_contents=sio.loadmat(self.fileo)
			
			x_o=self.rawPnts[self.bool_pnt,0]
			y_o=self.rawPnts[self.bool_pnt,1]
			z_o=self.rawPnts[self.bool_pnt,2]
			
			new={str_d : {'x_out':self.Outline,'rawPnts':self.rawPnts,'mask': self.bool_pnt,'x':x_o,'y':y_o,'z':z_o}}
			
			mat_contents.update(new) #update the dictionary
			if self.ui.refButton.isChecked():
				self.ui.refButton.setStyleSheet("background-color : rgb(77, 209, 97);")
			
			if self.ui.floatButton.isChecked():
				self.ui.floatButton.setStyleSheet("background-color : rgb(77, 209, 97);")
			sio.savemat(self.fileo,mat_contents)	
			#update status
			self.ui.statLabel.setText("Wrote %s data to output file %s."%(str_d,self.fileo))
		
		#check on write
		if self.refWritten==True and self.floatWritten==True:
			self.unsaved_changes=False



			
	def undo_pick(self):
		if hasattr(self,"lastSelectedIds"):
			for i in range(self.lastSelectedIds.GetNumberOfTuples()):
				#turn them from red to starting color
				self.colors.SetTuple(self.lastSelectedIds.GetValue(i),(70, 171, 176))
				self.bool_pnt[self.lastSelectedIds.GetValue(i)]=True
			self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
			self.vtkPntsPolyData.Modified()
			self.ui.vtkWidget.update()
		else:
			self.ui.statLabel.setText("No picked selection to revert.")
			
	def picker_callback(self,obj,event):
		
		extract = vtk.vtkExtractSelectedFrustum()
	
		fPlanes=obj.GetFrustum() #collection of planes based on unscaled display
	
		#scale frustum to account for the zaspect
		scaledPlanes=vtk.vtkPlanes()
		scaledNormals=vtk.vtkDoubleArray()
		scaledNormals.SetNumberOfComponents(3)
		scaledNormals.SetNumberOfTuples(6)
		scaledOrigins=vtk.vtkPoints()
		for j in range(6):
			i=fPlanes.GetPlane(j)
			k=i.GetOrigin()
			q=i.GetNormal()
			scaledOrigins.InsertNextPoint(k[0],k[1],k[2]/float(self.Zaspect))
			scaledNormals.SetTuple(j,(q[0],q[1],q[2]*float(self.Zaspect)))
		scaledPlanes.SetNormals(scaledNormals)
		scaledPlanes.SetPoints(scaledOrigins)
			
		
		extract.SetFrustum(scaledPlanes)
		extract.SetInputData(self.vtkPntsPolyData)
		extract.Update()
		extracted = extract.GetOutput()
		
		ids = vtk.vtkIdTypeArray()
		ids = extracted.GetPointData().GetArray("vtkOriginalPointIds")

		
		if ids:
			#store them in an array for an undo operation
			self.lastSelectedIds=ids
			for i in range(ids.GetNumberOfTuples()):
				#turn them red
				self.colors.SetTuple(ids.GetValue(i),(255,0,0))
				self.bool_pnt[ids.GetValue(i)]=False
		
			self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
			self.vtkPntsPolyData.Modified()
		
		
		self.ui.vtkWidget.update()
		#set flag on ui to show that data has been modified
		self.unsaved_changes=True
			
	def show_picking(self):
		#Updates when the 'r' button is pressed to provide a link between VTK & Qt hooks
		if self.picking == True:
			self.ui.pickActiveLabel.setStyleSheet("QLabel { background-color : red; color : white; }");
		else:
			self.ui.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }");
	
	def start_pick(self):
		#Required to change interactor
		style=vtk.vtkInteractorStyleRubberBandPick()
		self.iren.SetInteractorStyle(style)
		picker = vtk.vtkAreaPicker()
		self.iren.SetPicker(picker)
		picker.AddObserver("EndPickEvent", self.picker_callback)
		

		
	def get_input_data(self,filep,filec):
		color=(70, 171, 176)
		if hasattr(self,'pointActor'):
			self.ren.RemoveActor(self.pointActor)
		if hasattr(self,'outlineActor'):
			self.ren.RemoveActor(self.outlineActor)
		if hasattr(self,'rActor'):
			self.ren.RemoveActor(self.rActor)
		if hasattr(self,'fActor'):
			self.ren.RemoveActor(self.fActor)
		self.ui.levelButton.setChecked(False)
		
		
		if filep is None:
			filep,startdir=get_file('*.txt','Select the *.txt perimeter file:')

		elif not(os.path.isfile(filep)):
			print('Perimeter file invalid.')

		#return focus
		self.ui.vtkWidget.setFocus()
		
		if filec is None:
			filec,startdir=get_file('*.txt') #get filec
		
		#catch if cancel was pressed on file dialog or if a bad path was specified
		if filec == None or (filec != None and not(os.path.isfile(filec))):
			if hasattr(self,'vtkPntsPolyData'):
				print('No file selected, retaining current data.')
			else:
				return


		if filep != None: #because filediag can be cancelled
		
			#get delimiter
			with open(filep) as f:
				first_line = f.readline()
				if ',' in first_line: #NAMRC formatted file
					self.Outline=np.genfromtxt(filep,delimiter=",")
					print('NAMRC outline data type recognised.')
				else:
					self.Outline=np.genfromtxt(filep)
			self.outlineActor, _ =gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
			self.ren.AddActor(self.outlineActor)			
			self.filep=filep	
		_, ext = os.path.splitext(filec)
		
		if ext.lower() == '.txt':		
			self.rawPnts=np.genfromtxt(filec)
		elif ext.lower() == '.mat':
			try:
				self.rawPnts, self.Outline = read_mat(filec)
				self.outlineActor, _ =gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
				self.ren.AddActor(self.outlineActor)
			except:
				print('Could not add actors to visualisation.')
		elif ext.lower() == '.csv':
			self.rawPnts=np.genfromtxt(filec,skip_header=1,delimiter=',',usecols=(0,1,2))
		
		
		self.vtkPntsPolyData, \
		self.pointActor, self.colors = \
		gen_point_cloud(self.rawPnts,color,self.PointSize)
		self.bool_pnt=np.ones(len(self.rawPnts), dtype=bool)
		self.ren.AddActor(self.pointActor)


		self.filec=filec
		
		#get limits
		try:
			RefMin=np.amin(np.hstack((self.Outline,self.rawPnts)),axis=0)
			RefMax=np.amax(np.hstack((self.Outline,self.rawPnts)),axis=0)
		except:
			RefMin=np.amin(self.rawPnts,axis=0)
			RefMax=np.amax(self.rawPnts,axis=0)
		
		extents=RefMax-RefMin #extents
		rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
		self.limits=[RefMin[0]-rl, \
		  RefMax[0]+rl, \
		  RefMin[1]-rl, \
		  RefMax[1]+rl, \
		  RefMin[2],RefMax[2]]

		#add axes
		self.add_axis(self.limits,[1,1,1])
		
		#update status
		self.ui.statLabel.setText("Current perimeter file:%s    Current point cloud file:%s"%(filep,filec))
		
		#update
		self.ren.ResetCamera()
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
				# self.rActor.SetScale(1,1,self.Zaspect)
				self.rActor.SetScale(s)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				# self.fActor.SetScale(1,1,self.Zaspect)
				self.fActor.SetScale(s)
				self.fActor.Modified()
			if hasattr(self,'svd_pointActor'):
				self.svd_pointActor.SetScale(s)
				self.svd_pointActor.Modified()
				
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
			# nl=np.append(self.limits[0:4],[self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect])
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
			self.add_axis(self.limits,[1,1,1])
			self.ren.ResetCamera()

		elif key=="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(self.ui.vtkWidget._RenderWindow)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("point_cloud.png")
			writer.Write()
			print('Screen output saved to %s' %os.path.join(os.getcwd(),'point_cloud.png'))

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
				self.start_pick()
		
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
	

	
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

def get_svd_rotation_matrix(RP):
	'''
	Returns the rotation matrix about the X or Y axis required to take z component of the orthonormal matrix of RP to either 0,0,1 or 0,0,-1 depending on concavity - if reverse is true negative angles are applied. 
	'''
	_,_,vh = np.linalg.svd(RP) #vh is transpose from MATLAB's svd, returns normalised vectors
	#rows of vh are orthnormal vectors
	# print('X:',vh[0,:] / np.linalg.norm(vh[0,:]))
	# print('Y:',vh[1,:] / np.linalg.norm(vh[1,:]))
	# print('Z:',vh[2,:] / np.linalg.norm(vh[2,:]))
	
	#handles the case if the dataset is net convex vs. concave
	if vh[2,-1]<0:
		c=np.array([0,0,-1])
	else: 
		c=np.array([0,0,1])
	
	
	try:
		vh_y_norm = np.array([vh[2,0],0,vh[2,2]]) / np.linalg.norm(np.array([vh[2,0],0,vh[2,2]])) #xz plane projection
		vh_x_norm = np.array([0,vh[2,1],vh[2,2]]) / np.linalg.norm(np.array([0,vh[2,1],vh[2,2]])) #yz plane projection

		#solve for angle, update console
		a_y=np.arccos(np.clip(np.dot(vh_y_norm,c), -1.0, 1.0))
		a_x=np.arccos(np.clip(np.dot(vh_x_norm,c), -1.0, 1.0))

		if np.isnan(a_x) or np.isnan(a_y):
			raise ValueError('Potential division by zero during SVD decomposition.')
		elif a_x==0 and a_y==0:
			raise ValueError('Point cloud normals are aligned near working precision to main axes.')
		else:
			print('Difference about X and Y axis in degrees:\n',a_x*57.3,a_y*57.3)

		Ry=np.matrix([[np.cos(-a_y),0,np.sin(-a_y)],[0,1,0],[-np.sin(-a_y),0,np.cos(-a_y)]])
		Rx=np.matrix([[1,0,0],[0,np.cos(-a_x),-np.sin(-a_x)],[0,np.sin(-a_x),np.cos(-a_x)]])
		return Rx,Ry
	except ValueError as err:
		print(err)
		return np.asmatrix(np.eye(3,3)),np.asmatrix(np.eye(3,3)) #return identity matrices

	


if __name__ == '__main__':
	if len(sys.argv)==3:
		mask_def(sys.argv[1],sys.argv[2])
	elif len(sys.argv)==2:
		mask_def(sys.argv[1])
	else:
		mask_def()
