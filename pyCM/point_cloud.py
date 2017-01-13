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
z - increase z-aspect ratio
x - decrease z-aspect ratio
c - return to default z-aspect
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
a - toggles axes
o - togles outline (if present)
-------------------------------------------------------------------------------
ver 1.2 16-11-06
1.1 - Fixed array orientation, clipping issue, compass scaling and sped up writing output
      Added ReadMask
1.2 - Fixed window handling, now exits cleanly
1.3 - Completely refactored to reciprocate packages from later tools
'''
__author__ = "M.J. Roy"
__version__ = "1.2"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2017"

import sys
import os.path
from pkg_resources import Requirement, resource_filename
import numpy as np
import scipy.io as sio
import vtk
import vtk.util.numpy_support as vtk_to_numpy
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt4 import QtCore, QtGui
from pyCMcommon import *


nosio=False

def mask_def(*args,**kwargs):
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
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.Boxlayout = QtGui.QHBoxLayout(self.centralWidget)
		self.Subtendlayout=QtGui.QVBoxLayout()
		mainUiBox = QtGui.QFormLayout()

		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
		self.vtkWidget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		self.vtkWidget.setMinimumSize(1150, 700); #leave 100 px on the size for i/o

		self.Subtendlayout.addWidget(self.vtkWidget)
		self.statLabel=QtGui.QLabel("Idle")
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
		self.reloadButton = QtGui.QPushButton('Load')
		# self.reloadButton.setMaximumWidth(50)
		horizLine1=QtGui.QFrame()
		horizLine1.setFrameStyle(QtGui.QFrame.HLine)
		pickLabel=QtGui.QLabel("Pick options")
		pickLabel.setFont(headFont)
		self.pickerButton = QtGui.QPushButton('Picker')
		self.pickerButton.setMaximumWidth(50)
		self.pickActiveLabel=QtGui.QLabel("Pick active")
		self.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }");
		self.pickActiveLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.undoLastPickButton=QtGui.QPushButton('Undo last pick')
		horizLine2=QtGui.QFrame()
		horizLine2.setFrameStyle(QtGui.QFrame.HLine)
		outputLabel=QtGui.QLabel("Write output")
		outputLabel.setFont(headFont)
		self.refButton=QtGui.QRadioButton("Reference")
		self.floatButton=QtGui.QRadioButton("Floating")
		self.refButton.setChecked(True)
		self.writeButtonGroup = QtGui.QButtonGroup()
		self.writeButtonGroup.addButton(self.floatButton)
		self.writeButtonGroup.addButton(self.refButton)
		self.writeButtonGroup.setExclusive(True)
		self.writeButton=QtGui.QPushButton('Write')
		horizLine3=QtGui.QFrame()
		horizLine3.setFrameStyle(QtGui.QFrame.HLine)
		self.loadMatButton=QtGui.QPushButton('Load *.mat')

		
		#add to formlayout
		mainUiBox.addRow(self.reloadButton)
		mainUiBox.addRow(horizLine1)
		mainUiBox.addRow(pickLabel)
		mainUiBox.addRow(self.pickerButton,self.pickActiveLabel)
		mainUiBox.addRow(self.undoLastPickButton)
		#add formlayout to main ui
		
		mainUiBox.addRow(horizLine2)
		mainUiBox.addRow(outputLabel)
		mainUiBox.addRow(self.refButton,self.floatButton)
		mainUiBox.addRow(self.writeButton)
		mainUiBox.addRow(horizLine3)
		mainUiBox.addRow(self.loadMatButton)

		
		self.Boxlayout.addLayout(mainUiBox)
		
	def initialize(self):
		self.vtkWidget.start()

class pnt_interactor(QtGui.QMainWindow):
	"""
	Sets up the main VTK window, reads file and sets connections between UI and interactor
	"""
	def __init__(self, parent = None):
		QtGui.QMainWindow.__init__(self, parent)
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
		self.ui.pickerButton.clicked.connect(lambda: self.start_pick())
		self.ui.undoLastPickButton.clicked.connect(lambda: self.undo_pick())
		self.ui.writeButton.clicked.connect(lambda: self.write_new())
		self.ui.loadMatButton.clicked.connect(lambda: self.load_mat())
	
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
			self.ui.statLabel.setText("Current point file:%s"%filem)
			
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
			self.AddAxis(self.limits,1)
			
		except:
			print "Couldn't read in both sets of data."
		
		else:
			print 'Invalid *.mat file'
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
			
			if self.fileo:
				x_o=self.rawPnts[self.bool_pnt,0]
				y_o=self.rawPnts[self.bool_pnt,1]
				z_o=self.rawPnts[self.bool_pnt,2]
				sio.savemat(self.fileo,{str_d : {'x_out':self.Outline,'rawPnts':self.rawPnts,'mask': self.bool_pnt,'x':x_o,'y':y_o,'z':z_o,'fname':self.filec}})
				print 'Wrote %s data'%(str_d)
		else:
			mat_vars=sio.whosmat(self.fileo)
			if str_d in mat_vars[0]: #tell the user that they might overwrite their data
				ret=QtGui.QMessageBox.warning(self, "pyCM Warning", \
				"The %s dataset has already been written. Overwrite?"%(str_d), \
				QtGui.QMessageBox.No,QtGui.QMessageBox.Yes,QtGui.QMessageBox.NoButton)
				if ret == QtGui.QMessageBox.No: #don't overwrite
					return
			print self.fileo
			mat_contents=sio.loadmat(self.fileo)
			
			x_o=self.rawPnts[self.bool_pnt,0]
			y_o=self.rawPnts[self.bool_pnt,1]
			z_o=self.rawPnts[self.bool_pnt,2]
			
			new={str_d : {'x_out':self.Outline,'rawPnts':self.rawPnts,'mask': self.bool_pnt,'x':x_o,'y':y_o,'z':z_o,'fname':self.filec}}
			
			mat_contents.update(new) #update the dictionary
				
			sio.savemat(self.fileo,mat_contents)	
			print 'Wrote %s data'%(str_d)
			
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
		for j in xrange(6):
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
		if self.picking == True:
			self.ui.pickActiveLabel.setStyleSheet("QLabel { background-color : red; color : white; }");
		else:
			self.ui.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }");
	
	def start_pick(self):

		style=vtk.vtkInteractorStyleRubberBandPick()		
		self.iren.SetInteractorStyle(style)
		picker = vtk.vtkAreaPicker()
		self.iren.SetPicker(picker)

		picker.AddObserver("EndPickEvent", self.picker_callback)
		self.ui.vtkWidget.setFocus()

		
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
		
		Perim=False #unless otherwise spec'ed
		
		if filep is None:
			filep,startdir=get_file('*.txt','Select the *.txt perimeter file (optional):')
			if filep != None:
				Perim=True
		elif not(os.path.isfile(filep)):
			print 'Perimeter file invalid.'

		
		if filec is None:
			filec,startdir=get_file('*.txt') #get filec
		
		#catch if cancel was pressed on file dialog or if a bad path was specified
		if filec == None or (filec != None and not(os.path.isfile(filec))):
			if hasattr(self,'vtkPntsPolyData'):
				print 'No file selected, retaining current data.'
			else:
				print 'No file selected; exiting.'
				exit()
		
		
		if filep != None: #because filediag can be cancelled
			self.Outline=np.genfromtxt(filep)
			self.outlineActor=gen_outline(self.Outline,tuple(np.array(color)/float(255)),self.PointSize)
			self.ren.AddActor(self.outlineActor)			
		
		_, ext = os.path.splitext(filec)
		
		if ext == '.dat':
			self.rawPnts=np.genfromtxt(filec,skiprows=1)
		elif ext == '.txt':
			self.rawPnts=np.genfromtxt(filec)
		elif ext == '.mat':
			mc= sio.loadmat(filec) #mat contents
			try:
				self.rawPnts=np.hstack((mc['x'],mc['y'],mc['z']))
				self.Outline=mc['x_out']
				Perim=True
			except KeyError:
				print "Couldn't read variables from file."
				return
		
		
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
		self.AddAxis(self.limits,1)
		
		#update status
		self.ui.statLabel.setText("Current point file:%s"%filec)
		
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
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(1,1,self.Zaspect)
				self.pointActor.Modified()
			if hasattr(self,'rActor'):
				self.rActor.SetScale(1,1,self.Zaspect)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				self.fActor.SetScale(1,1,self.Zaspect)
				self.fActor.Modified()
			nl=np.append(self.limits[0:4],[self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect])
			self.AddAxis(nl,1/self.Zaspect)



		elif key=="x":
			self.Zaspect=self.Zaspect*0.5
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'rActor'):
				self.rActor.SetScale(1,1,self.Zaspect)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				self.fActor.SetScale(1,1,self.Zaspect)
				self.fActor.Modified()
			nl=np.append(self.limits[0:4],[self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect])
			self.AddAxis(nl,1/self.Zaspect)


		elif key=="c":
			self.Zaspect=1.0
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'rActor'):
				self.rActor.SetScale(1,1,self.Zaspect)
				self.rActor.Modified()
			if hasattr(self,'fActor'):
				self.fActor.SetScale(1,1,self.Zaspect)
				self.fActor.Modified()
			self.AddAxis(self.limits,1)


		elif key=="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(self.ui.vtkWidget._RenderWindow)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("point_cloud.png")
			writer.Write()
			print 'Screen output saved to %s' %os.path.join(currentdir,'point_cloud.png')

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
	
	def AddAxis(self,limits,scale):
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
		self.ax3D.SetZAxisRange(limits[-2]*scale,limits[-1]*scale)
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


