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
o - togles outline (if present)
-------------------------------------------------------------------------------
ver 1.2 16-11-06
1.1 - Fixed array orientation, clipping issue, compass scaling and sped up writing output
      Added ReadMask
1.2 - Fixed window handling, now exits cleanly
1.3 - Modified to run in Python 3.x, uses VTK keyboard interrupts to start picking, Qt button for this function has been commented out.
'''
__author__ = "M.J. Roy"
__version__ = "1.3"
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
from pyCMcommon import *


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
	
	window = pnt_interactor()
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
		MainWindow.setObjectName("MainWindow")
		MainWindow.setWindowTitle("pyCM - Point editor v%s" %__version__)
		MainWindow.resize(1280, 720)
		self.centralWidget = QtWidgets.QWidget(MainWindow)
		self.Boxlayout = QtWidgets.QHBoxLayout(self.centralWidget)
		self.Subtendlayout=QtWidgets.QVBoxLayout()
		mainUiBox = QtWidgets.QFormLayout()

		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
		self.vtkWidget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
		self.vtkWidget.setMinimumSize(1150, 700); #leave 100 px on the size for i/o

		self.Subtendlayout.addWidget(self.vtkWidget)
		self.statLabel=QtWidgets.QLabel("Idle")
		self.statLabel.setWordWrap(True)
		self.statLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.statLabel.setMinimumWidth(100)
		self.Subtendlayout.addWidget(self.statLabel)
		self.Subtendlayout.addStretch(1)
		self.Boxlayout.addLayout(self.Subtendlayout)
		self.Boxlayout.addStretch(1)
		
		MainWindow.setCentralWidget(self.centralWidget)
		headFont=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)
		
		#define buttons/widgets
		self.reloadButton = QtWidgets.QPushButton('Load')
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
		self.reduce = QtWidgets.QSpinBox()
		self.reduce.setValue(0)
		self.reduceButton = QtWidgets.QPushButton('Reduce')
		self.revertButton = QtWidgets.QPushButton('Reload source')

		reduceBoxlayout= QtWidgets.QGridLayout()
		reduceBoxlayout.addWidget(self.reduce,1,1)
		reduceBoxlayout.addWidget(self.reduceButton,1,2)
		
		
		# self.reloadButton.setMaximumWidth(50)
		horizLine1=QtWidgets.QFrame()
		horizLine1.setFrameStyle(QtWidgets.QFrame.HLine)
		pickLabel=QtWidgets.QLabel("Pick options")
		pickLabel.setFont(headFont)
		# self.pickerButton = QtWidgets.QPushButton('Picker')
		# self.pickerButton.setMaximumWidth(50)
		self.pickHelpLabel=QtWidgets.QLabel("Press R to activate")
		self.pickActiveLabel=QtWidgets.QLabel("Pick active")
		self.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }");
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
		self.releaseButton=QtWidgets.QPushButton('Release')
		self.releaseButton.setEnabled(False)
		horizLine3=QtWidgets.QFrame()
		horizLine3.setFrameStyle(QtWidgets.QFrame.HLine)
		self.loadMatButton=QtWidgets.QPushButton('Load *.mat')

		#add to formlayout
		mainUiBox.addRow(self.reloadButton)
		mainUiBox.addRow(scalingLabel)
		mainUiBox.addRow(scaleBoxlayout)
		mainUiBox.addRow(horizLine1)
		mainUiBox.addRow(pickLabel)
		mainUiBox.addRow(reduceBoxlayout)
		# mainUiBox.addRow(self.pickerButton,self.pickActiveLabel)
		mainUiBox.addRow(self.pickHelpLabel,self.pickActiveLabel)
		mainUiBox.addRow(self.undoLastPickButton)
		
		mainUiBox.addRow(self.revertButton)
		#add formlayout to main ui
		
		mainUiBox.addRow(horizLine2)
		mainUiBox.addRow(outputLabel)
		mainUiBox.addRow(self.refButton,self.floatButton)
		mainUiBox.addRow(self.writeButton,self.releaseButton)
		mainUiBox.addRow(horizLine3)
		mainUiBox.addRow(self.loadMatButton)

		
		self.Boxlayout.addLayout(mainUiBox)
		
	def initialize(self):
		self.vtkWidget.start()

class pnt_interactor(QtWidgets.QMainWindow):
	"""
	Sets up the main VTK window, reads file and sets connections between UI and interactor
	"""
	def __init__(self, parent = None):
		QtWidgets.QMainWindow.__init__(self, parent)
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

		self.ui.reloadButton.clicked.connect(lambda: self.get_input_data(None,None))
		# self.ui.pickerButton.clicked.connect(lambda: self.start_pick())
		self.ui.undoLastPickButton.clicked.connect(lambda: self.undo_pick())
		self.ui.writeButton.clicked.connect(lambda: self.write_new())
		self.ui.releaseButton.clicked.connect(lambda: self.release_mat())
		self.ui.loadMatButton.clicked.connect(lambda: self.load_mat())
		self.ui.revertButton.clicked.connect(lambda: self.undo_revert())
		self.ui.reduceButton.clicked.connect(lambda: self.reduce_pnts())
	
	def undo_revert(self):
		'''
		Reloads all data based on filec & filep (if it exists)
		'''
		if hasattr(self,'filep'):
			self.get_input_data(self.filep,self.filec)
		else:
			self.get_input_data(None,self.filec)
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
		
		if hasattr(self,'pointActor'):
			self.ren.RemoveActor(self.pointActor)
		if hasattr(self,'outlineActor'):
			self.ren.RemoveActor(self.outlineActor)
		
		filem, _, =get_file('*.mat')
		
		if filem: #check variables
			mat_contents = sio.loadmat(filem)
			try:
				rp=mat_contents['ref']['rawPnts'][0][0]
				ind=mat_contents['ref']['mask'][0][0][0]
				
				rp=rp[np.where(ind)]
				
				color=(242, 101, 34)
				_, self.rActor, _, = gen_point_cloud(rp,color,self.PointSize)
				self.ren.AddActor(self.rActor)
				
				#do other one
				fp=mat_contents['float']['rawPnts'][0][0]
				ind=mat_contents['float']['mask'][0][0][0]
				
				fp=fp[np.where(ind)]

				
				color=(255, 205, 52)
				_, self.fActor, _, = gen_point_cloud(fp,color,self.PointSize)
				self.ren.AddActor(self.fActor)
			
				# update status
				self.ui.statLabel.setText("Current file:%s"%filem)
				
				RefMin=np.amin(np.vstack((fp,rp)),axis=0)
				RefMax=np.amax(np.vstack((fp,rp)),axis=0)

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
			print('Invalid *.mat file')
			return
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		

	def write_new(self):
		
		if self.ui.refButton.isChecked():
			str_d='ref'
			
		if self.ui.floatButton.isChecked():
			str_d='float'
		
		if not hasattr(self,'fileo'):
			self.fileo, _, = get_open_file('*.mat',os.getcwd())
			self.ui.releaseButton.setEnabled(True)
			if self.fileo:
				x_o=self.rawPnts[self.bool_pnt,0]
				y_o=self.rawPnts[self.bool_pnt,1]
				z_o=self.rawPnts[self.bool_pnt,2]
				sio.savemat(self.fileo,{str_d : {'x_out':self.Outline,'rawPnts':self.rawPnts,'mask': self.bool_pnt,'x':x_o,'y':y_o,'z':z_o,'fname':self.filec}})
				print('Wrote %s data to %s'%(str_d,self.fileo))
		else:
			mat_vars=sio.whosmat(self.fileo)
			if str_d in mat_vars[0]: #tell the user that they might overwrite their data
				ret=QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
				"The %s dataset has already been written. Overwrite?"%(str_d), \
				QtWidgets.QMessageBox.No,QtWidgets.QMessageBox.Yes,QtWidgets.QMessageBox.NoButton)
				if ret == QtWidgets.QMessageBox.No: #don't overwrite
					return
			
			mat_contents=sio.loadmat(self.fileo)
			
			x_o=self.rawPnts[self.bool_pnt,0]
			y_o=self.rawPnts[self.bool_pnt,1]
			z_o=self.rawPnts[self.bool_pnt,2]
			
			new={str_d : {'x_out':self.Outline,'rawPnts':self.rawPnts,'mask': self.bool_pnt,'x':x_o,'y':y_o,'z':z_o,'fname':self.filec}}
			
			mat_contents.update(new) #update the dictionary
				
			sio.savemat(self.fileo,mat_contents)	
			print('Wrote %s data'%(str_d))
			
	def undo_pick(self):
		if hasattr(self,"lastSelectedIds"):
			for i in range(self.lastSelectedIds.GetNumberOfTuples()):
				#turn them from red to starting color
				self.colors.SetTuple(self.lastSelectedIds.GetValue(i),(70, 171, 176))
				self.bool_pnt[self.lastSelectedIds.GetValue(i)]=True
			self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
			self.vtkPntsPolyData.Modified()
			self.ui.vtkWidget.update()
			
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
		
		
		if filep is None:
			filep,startdir=get_file('*.txt','Select the *.txt perimeter file (optional):')

		elif not(os.path.isfile(filep)):
			print('Perimeter file invalid.')

		
		if filec is None:
			filec,startdir=get_file('*.txt') #get filec
		
		#catch if cancel was pressed on file dialog or if a bad path was specified
		if filec == None or (filec != None and not(os.path.isfile(filec))):
			print(filec)
			if hasattr(self,'vtkPntsPolyData'):
				print('No file selected, retaining current data.')
			else:
				print('No file selected; exiting.')
				exit()
		
		if filep != '': #because filediag can be cancelled
			self.Outline=np.genfromtxt(filep)
			self.outlineActor, _ =gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
			self.ren.AddActor(self.outlineActor)			
			self.filep=filep
			
		_, ext = os.path.splitext(filec)
		

		elif ext == '.txt':
			self.rawPnts=np.genfromtxt(filec)
		elif ext == '.mat':
			try:
				self.rawPnts, self.Outline = read_uom_mat(filec)
				self.outlineActor, _ =gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
				self.ren.AddActor(self.outlineActor)
			except:
				exit()
		
		
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
	
	def release_mat(self):
		if hasattr(self,'fileo'):
			print('Releasing focus on %s . . .'%self.fileo)
			del self.fileo
	
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
		perimFile=sys.argv[1]
		pcloudFile=sys.argv[2]
		mask_def(perimFile,pcloudFile)
	elif len(sys.argv)>1:
		pcloudFile=sys.argv[1]
		mask_def(pcloudFile)
	else:
		mask_def()


