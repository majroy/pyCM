#!/usr/bin/env python
'''
Uses VTK python to allow for preprocessing FEAs associated with the
 contour method. Full interaction requires a 3-button mouse and keyboard.
-------------------------------------------------------------------------------
Current mapping is as follows:
LMB   - rotate about point cloud centroid.
MMB   - pan
RMB   - zoom
1     - view 1, default, looks down z axis onto xy plane
2     - view 2, looks down x axis onto zy plane
3     - view 3, looks down y axis onto zx plane
z     - increase z-aspect ratio of displacement BC
x     - decrease z-aspect ratio of displacement BC
c     - return to default z-aspect
f     - flip colors from white on dark to dark on white
i     - save output to .png in current working directory
r     - remove/reinstate compass/axes
o     - remove/reinstate outline
LMB+p - The p button with the Left mouse button allow 
        for selecting rigid body boundary conditions. 
        Click first and then press p to select.
e     - allows the user to change their FEA exec location
-------------------------------------------------------------------------------
ver 0.1 17-01-04
1.0 - Refactored for PyQt5 & Python 3.x
'''
__author__ = "M.J. Roy"
__version__ = "1.0"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2018"

import os,sys,time,yaml
import subprocess as sp
from pkg_resources import Requirement, resource_filename
import numpy as np
import scipy.io as sio
from scipy.interpolate import bisplev,interp1d
from scipy.spatial.distance import pdist, squareform
import vtk
from vtk.util.numpy_support import vtk_to_numpy as v2n
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from pyCMcommon import *



def FEAtool(*args, **kwargs):
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
	
	window = msh_interactor()
	if len(args)==2: msh_interactor.getInputData(window,args[0],args[1])
	elif len(args)==1: msh_interactor.getInputData(window,args[0],None)
	else: msh_interactor.getInputData(window,"",None)


	window.show()
	splash.finish(window)
	window.iren.Initialize() # Need this line to actually show the render inside Qt

	ret = app.exec_()
	
	if sys.stdin.isatty() and not hasattr(sys,'ps1'):
		sys.exit(ret)
	else:
		# print window.ren.GetActiveCamera().GetPosition(),window.ren.GetActiveCamera().GetFocalPoint()
		return window
		# http://stackoverflow.com/questions/3394835/args-and-kwargs
class sf_MainWindow(object):
	"""
	Class to build qt interaction, including VTK widget
	setupUi builds, initialize starts VTK widget
	"""
	
	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.setWindowTitle("pyCM - FEA preprocessing v%s" %__version__)
		MainWindow.resize(1280, 720)
		self.centralWidget = QtWidgets.QWidget(MainWindow)
		self.Boxlayout = QtWidgets.QHBoxLayout(self.centralWidget)
		mainUiBox = QtWidgets.QFormLayout()

		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
		self.vtkWidget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
		self.vtkWidget.setMinimumSize(1150, 640); #leave 100 px on the size for i/o

		self.Boxlayout.addWidget(self.vtkWidget)
		self.Boxlayout.addStretch(1)
		MainWindow.setCentralWidget(self.centralWidget)

		horizLine1=QtWidgets.QFrame()
		horizLine1.setFrameStyle(QtWidgets.QFrame.HLine)
		outlineLabel=QtWidgets.QLabel("Outline modification")
		headFont=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)
		outlineLabel.setFont(headFont)
		seedLengthLabel=QtWidgets.QLabel("Spacing")
		self.seedLengthInput = QtWidgets.QLineEdit()
		numSeedLabel=QtWidgets.QLabel("Quantity")
		self.numSeed = QtWidgets.QSpinBox()
		self.numSeed.setMinimum(4)
		self.numSeed.setMaximum(750)
		self.numSeed.setValue(100)

		self.updateOutlineButton = QtWidgets.QPushButton('Update')
		self.spacingButton=QtWidgets.QRadioButton("Spacing")
		self.quantityButton=QtWidgets.QRadioButton("Quantity")
		self.quantityButton.setChecked(True)
		self.outlineButtonGroup = QtWidgets.QButtonGroup()
		self.outlineButtonGroup.addButton(self.spacingButton)
		self.outlineButtonGroup.addButton(self.quantityButton)
		self.outlineButtonGroup.setExclusive(True)
		self.geoButton=QtWidgets.QRadioButton("geo")
		self.dxfButton=QtWidgets.QRadioButton("dxf")
		self.geoButton.setChecked(True)
		self.outlineButtonGroup = QtWidgets.QButtonGroup()
		self.outlineButtonGroup.addButton(self.geoButton)
		self.outlineButtonGroup.addButton(self.dxfButton)
		self.outlineButtonGroup.setExclusive(True)
		self.outlineButton = QtWidgets.QPushButton('Write')
		self.outlineButton.setMinimumWidth(50)
		
		horizLine2=QtWidgets.QFrame()
		horizLine2.setFrameStyle(QtWidgets.QFrame.HLine)
		meshscriptLabel=QtWidgets.QLabel("Generate mesh")
		meshscriptLabel.setFont(headFont)
		lengthLabel=QtWidgets.QLabel("Length")
		self.lengthInput = QtWidgets.QLineEdit()
		numPartLabel=QtWidgets.QLabel("Partitions")
		self.numPart = QtWidgets.QSpinBox()
		self.numPart.setValue(10)
		self.numPart.setMinimum(3)
		self.gmshButton=QtWidgets.QRadioButton("Gmsh")
		self.abaButton=QtWidgets.QRadioButton("Abaqus")
		self.gmshButton.setChecked(True)
		self.codeButtonGroup = QtWidgets.QButtonGroup()
		self.codeButtonGroup.addButton(self.gmshButton)
		self.codeButtonGroup.addButton(self.abaButton)
		self.codeButtonGroup.setExclusive(True)
		
		self.quadButton=QtWidgets.QRadioButton("quads")
		self.tetButton=QtWidgets.QRadioButton("tets")
		self.quadButton.setChecked(True)
		self.meshscriptButton = QtWidgets.QPushButton('Write')
		self.runmeshButton = QtWidgets.QPushButton('Run')
		self.readVTKbutton = QtWidgets.QPushButton('Read VTK mesh')
		self.mtypeButtonGroup = QtWidgets.QButtonGroup()
		self.mtypeButtonGroup.addButton(self.tetButton)
		self.mtypeButtonGroup.addButton(self.quadButton)
		self.mtypeButtonGroup.setExclusive(True)
		
		horizLine3=QtWidgets.QFrame()
		horizLine3.setFrameStyle(QtWidgets.QFrame.HLine)
		bcLabel=QtWidgets.QLabel("Impose BCs & material")
		bcLabel.setFont(headFont)
		# self.rigidBodyButton = QtWidgets.QPushButton('Rigid body')
		self.imposeSpline = QtWidgets.QPushButton('Impose spline fit')
		rBLabel=QtWidgets.QLabel("Rigid body BCs")
		self.rigidBodyButton = QtWidgets.QPushButton('Choose')
		self.rigidBodyUndoButton = QtWidgets.QPushButton('Revert')
		materialLabel=QtWidgets.QLabel("Material properties")
		poissonLabel=QtWidgets.QLabel("Poisson's ratio, v")
		modulusLabel=QtWidgets.QLabel("Modulus, E (MPa)")
		self.poissonInput = QtWidgets.QLineEdit()
		self.poissonInput.setText("%.3f"%(0.3))
		self.modulusInput = QtWidgets.QLineEdit()
		self.modulusInput.setText("%7.0f"%(200000))

		horizLine4=QtWidgets.QFrame()
		horizLine4.setFrameStyle(QtWidgets.QFrame.HLine)
		FEALabel=QtWidgets.QLabel("Compose FEA")
		FEALabel.setFont(headFont)
		self.CalculixButton=QtWidgets.QRadioButton("Calculix")
		self.AbaqusButton=QtWidgets.QRadioButton("Abaqus")
		self.CalculixButton.setChecked(True)
		self.ctypeButtonGroup = QtWidgets.QButtonGroup()
		self.ctypeButtonGroup.addButton(self.AbaqusButton)
		self.ctypeButtonGroup.addButton(self.CalculixButton)
		# self.ctypeButtonGroup.addButton(self.CodeAsterButton)
		self.ctypeButtonGroup.setExclusive(True)
		self.runFEAButton=QtWidgets.QRadioButton("Run")
		self.runFEAButton.setChecked(True)
		self.goButton = QtWidgets.QPushButton('Write')
		
		horizLine5=QtWidgets.QFrame()
		horizLine5.setFrameStyle(QtWidgets.QFrame.HLine)
		self.statLabel=QtWidgets.QLabel("Idle")
		self.statLabel.setWordWrap(True)
		self.statLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.statLabel.setMinimumWidth(50)

		mainUiBox.addRow(horizLine1)
		mainUiBox.addRow(outlineLabel)
		mainUiBox.addRow(self.quantityButton,self.spacingButton)
		mainUiBox.addRow(seedLengthLabel,self.seedLengthInput)
		mainUiBox.addRow(numSeedLabel,self.numSeed)
		mainUiBox.addRow(self.geoButton,self.dxfButton)
		mainUiBox.addRow(self.updateOutlineButton,self.outlineButton)
		mainUiBox.addRow(horizLine2)
		mainUiBox.addRow(meshscriptLabel)
		mainUiBox.addRow(self.gmshButton,self.abaButton)
		mainUiBox.addRow(lengthLabel,self.lengthInput)
		mainUiBox.addRow(numPartLabel,self.numPart)
		mainUiBox.addRow(self.quadButton,self.tetButton)
		mainUiBox.addRow(self.meshscriptButton,self.runmeshButton)
		mainUiBox.addRow(self.readVTKbutton)
		mainUiBox.addRow(horizLine3)
		mainUiBox.addRow(bcLabel)
		mainUiBox.addRow(self.imposeSpline)
		mainUiBox.addRow(rBLabel)
		mainUiBox.addRow(self.rigidBodyButton,self.rigidBodyUndoButton)
		mainUiBox.addRow(materialLabel)
		mainUiBox.addRow(poissonLabel,self.poissonInput)
		mainUiBox.addRow(modulusLabel,self.modulusInput)
		mainUiBox.addRow(horizLine4)
		mainUiBox.addRow(FEALabel)
		mainUiBox.addRow(self.CalculixButton,self.AbaqusButton)
		mainUiBox.addRow(self.goButton,self.runFEAButton)
		mainUiBox.addRow(horizLine5)
		mainUiBox.addRow(self.statLabel)

		self.Boxlayout.addLayout(mainUiBox)

	
	def initialize(self):
		self.vtkWidget.start()

class msh_interactor(QtWidgets.QMainWindow):
	"""
	Sets up the main VTK window, reads file and sets connections between UI and interactor
	"""
	def __init__(self, parent = None):
		QtWidgets.QMainWindow.__init__(self, parent)
		self.ui = sf_MainWindow()
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
		self.iren.AddObserver("KeyPressEvent", self.Keypress)
		
		self.PointSize=2
		self.LineWidth=1
		self.Zaspect=1.0
		self.limits=np.empty(6)

		#check config file
		try:
			self.filec = resource_filename("pyCM","pyCMconfig.yml")#needs to be pyCM
		except: #resource_filename will inform if the directory doesn't exist
			# self.filec=resource_filename("numpy","pyCMconfig.yml")#needs to be pyCM
			print("Did not find config file in the pyCM installation directory.")
		try:
			with open(self.filec,'r') as ymlfile:
				self.cfg = yaml.load(ymlfile)	
		except:
			try:
				print("Trying to get input")
				self.cfg= GetFEAconfig(['','',''],self.filec)
			except:
				sys.exit("Failed to set config file. Quitting.")
				
	
		#Sets all connections between gui and functions

		
		self.ui.updateOutlineButton.clicked.connect(lambda: self.ModOutline())
		self.ui.outlineButton.clicked.connect(lambda: self.WriteOut())
		self.ui.meshscriptButton.clicked.connect(lambda: self.WriteGeo())
		self.ui.runmeshButton.clicked.connect(lambda: self.RunMeshScript())
		self.ui.readVTKbutton.clicked.connect(lambda: self.DisplayMesh())
		self.ui.rigidBodyButton.clicked.connect(lambda: self.ImposeRigidBody())
		self.ui.rigidBodyUndoButton.clicked.connect(lambda: self.UndoRigidBody())
		self.ui.imposeSpline.clicked.connect(lambda: self.ImposeSplineFit())
		self.ui.goButton.clicked.connect(lambda: self.doFEA())

	def getInputData(self,filer,outputd):
		if outputd==None and (not os.path.isfile(filer)):
			self.outputd=None
			self.filer,startdir=get_file("*.mat")
			if self.filer == None:
				#if its called from the commandline
				if sys.stdin.isatty() and not hasattr(sys,'ps1'):
					sys.exit("No file identified")
				else:
					print("No surface file specified, exiting.")
					#return to interactive python
					exit()
		elif not outputd==None and filer==None:
			self.outputd=outputd
			if not os.path.exists(self.outputd): #make the directory if it doesn't exist
				os.makedirs(self.outputd)
		
		elif os.path.isfile(filer):
			self.outputd=os.path.dirname(filer)
			self.filer = filer
		else:
			sys.exit("Arguments not specified correctly. Quitting.")

		#Read in reference data, calculate relevant details
		try:
			mat_contents = sio.loadmat(self.filer)
			try:
				#read in for ImposeSplineFit function
				bsplinerep=mat_contents['spline_x']['tck'][0][0]
				#recast tck as a tuple
				self.tck=tuple()
				for j in range(5):
					self.tck=self.tck+tuple(bsplinerep[0,j])
				#outline for mesh generation
				self.Outline=mat_contents['x_out']
				RefMin=np.amin(self.Outline,axis=0)
				RefMax=np.amax(self.Outline,axis=0)
				self.limits=[RefMin[0],RefMax[0],RefMin[1],RefMax[1],1,1]
				minLength=np.minimum(self.limits[1]-self.limits[0],self.limits[3]-self.limits[2])
				gs=squareform(pdist(self.Outline[:,:2],'euclidean')) #does the same thing as MATLAB's pdist2
				self.Dist=np.mean(np.sort(gs)[:,1]) #get length of first element
				self.ui.lengthInput.setText(str(round(3*minLength)))
				self.ui.seedLengthInput.setText("%4.4f"%self.Dist)
				self.ui.numSeed.setValue(len(self.Outline))

			except KeyError:
				print("Couldn't read variables from file. Quitting.")
				return
			
		except KeyError:
			print("Error reading reference data")
			return
		color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
		self.outlineActor, _, = gen_outline(self.Outline,color,self.PointSize)
		self.ren.AddActor(self.outlineActor)
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		
		# self.outlineActor=self.DisplayOutline(self.Outline,(0.2784,0.6745,0.6941),self.PointSize)

		
	def UndoRigidBody(self):

		if not hasattr(self,"pickedCornerInd"):
			return #ie do nothing
		#remove all actors/data from ImposeRigidBody and revert axis
		for actorInd in self.pickedActorInd:
			self.ren.RemoveActor(self.a[actorInd])
		self.ui.vtkWidget.update()

		del self.a
		del self.pickedCornerInd
		del self.pickedActorInd
		l=self.limits
		self.limits[0:4]=[l[0]+self.asize,l[1]-self.asize,l[2]+self.asize,l[3]-self.asize]
		self.AddAxis(np.append(self.limits[0:4],[self.BClimits[-2]*self.Zaspect,self.BClimits[-1]*self.Zaspect]),1/self.Zaspect)
		
	def ImposeRigidBody(self):
		"""
		Displays and identifies where rigid body BCs can be imposed, calls 
		"""
		
		if hasattr(self,"pickedCornerInd"):
			print("Undo/revert from current BCs first")
			return
		
		if not hasattr(self,"corners"):
			msg=QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Information)
			msg.setText("Impose contour displacement BCs before imposing rigid body BCs")
			msg.setWindowTitle("pyCM Error")
			msg.exec_()
			return
				
		directions=np.array([[0,-1,0],[0,-1,0],[1,0,0],[1,0,0],
                                [0,1,0],[0,1,0],[-1,0,0],[-1,0,0],
					  [0,-1,0],[0,-1,0],[1,0,0],[1,0,0],
                                [0,1,0],[0,1,0],[-1,0,0],[-1,0,0]]) #corresponds to ccw
		if not self.OutlineIsCCW: directions=np.flipud(directions)
		
		#arrow size is 5% max size of domain
		self.asize=np.maximum(self.limits[1]-self.limits[0],self.limits[3]-self.limits[2])*0.05
		
		self.a=[] #arrow actors on 'front' face
		self.aInd=np.empty([8,2]) #index of corners and their 
		for c in range(len(self.corners)):
			if c==0:
				self.a.append(DrawArrow(self.corners[c,:],self.asize,directions[c,:],self.ren))
				self.a.append(DrawArrow(self.corners[c,:],self.asize,directions[-1,:],self.ren))
				self.aInd[c,:]=[c,c+1]
			else:
				self.a.append(DrawArrow(self.corners[c,:],self.asize,directions[c*2-1,:],self.ren))
				self.a.append(DrawArrow(self.corners[c,:],self.asize,directions[c*2,:],self.ren))
				self.aInd[c,:]=[c*2,c*2+1]
				
		
		
		#bump out axis limits
		l=self.limits
		self.limits[0:4]=[l[0]-self.asize,l[1]+self.asize,l[2]-self.asize,l[3]+self.asize]
		#there will be a BCactor, so BClimits will exist
		self.AddAxis(np.append(self.limits[0:4],[self.BClimits[-2]*self.Zaspect,self.BClimits[-1]*self.Zaspect]),1/self.Zaspect)
		self.ui.vtkWidget.update()
		
		self.picks=0
		self.pickedCornerInd=[]
		self.pickedActorInd=[]
		self.iren.AddObserver("EndPickEvent",self.checkPick)
		
	def checkPick(self,object,event):
		"""
		Activates two pick events, one for each discrete rigid body BC. Needs ImposeRigidBody to run.
		"""
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		picker=vtk.vtkPropPicker()
		
		clickPos=self.iren.GetEventPosition()
		
		picker.Pick(clickPos[0],clickPos[1],0,self.ren)
		NewPickedActor = picker.GetActor()
		
		if NewPickedActor:
			
			i=int(NewPickedActor.GetAddressAsString('vtkPolyData')[5:], 16)
			#compare it against the addresses in a
			count=0
			
			for actor in self.a:
				ai=int(actor.GetAddressAsString('vtkPolyData')[5:], 16)
				if ai==i:
					self.pickedCornerInd.append(np.where(self.aInd==count)[0][0])
					if self.picks == 0:
						#highlight both arrows for first pick
						self.a[int(self.aInd[self.pickedCornerInd,0][0])].GetProperty().SetColor(1,0,0)
						self.pickedActorInd.append(int(self.aInd[self.pickedCornerInd,0][0]))
						self.a[int(self.aInd[self.pickedCornerInd,1][0])].GetProperty().SetColor(1,0,0)
						self.pickedActorInd.append(int(self.aInd[self.pickedCornerInd,1][0]))
						self.picks+=1
					else:
						#get the second
						#make sure that it hasn't already been selected
						if count not in self.pickedActorInd and len(self.pickedActorInd)<3:
							self.a[count].GetProperty().SetColor(1,0,0)
							self.pickedActorInd.append(count)
							self.picks+=1
				else:
					count+=1
			

		#check if pick condition has been satisfied
		if self.picks == 2 and len(self.pickedActorInd)==3:
			#delete all other arrow actors and stop the pick observer
			aDel=list(set(range(len(self.a)))-set(self.pickedActorInd))
			for actorInd in aDel:
				self.ren.RemoveActor(self.a[actorInd])
			self.iren.RemoveObservers("EndPickEvent")
			#corners and their coordinates are stored in 
			# print self.corners[self.pickedCornerInd,:]
			#node numbers of the picked corners are stored in
			# print self.cornerInd[self.pickedCornerInd]
			self.ui.statLabel.setText("Rigid body BCs selected . . . Idle")

	def ModOutline(self):
		"""
		Assumes that outline is unordered, finds corners and sorts, returns a new outline with either the node count indicated or an even number of nodes according to the length indicated.
		"""
		
		#if there is already a respaced outline, then remove it from the display
		if hasattr(self,"respacedOutlineActor"):
			self.ren.RemoveActor(self.respacedOutlineActor)
			self.rsOutline=[]
		
		#if it doesn't have corners, then calculate them
		if not hasattr(self,"corners"):
			d=np.array([])
			for j in range(len(self.Outline[:,0])):
				d=np.append(d,
				np.sqrt((self.limits[0]-self.Outline[j,0])**2+(self.limits[2]-self.Outline[j,1])**2)
				)
			ind=np.where(d==np.amin(d))[0][0] #to avoid making ind an array
			

			#reorder the points so that ind is first
			self.Outline=np.vstack((self.Outline[ind::,:],self.Outline[0:ind+1,:]))

			c_target=np.array([
			[self.limits[0],self.limits[3]], #xmin,ymax
			[self.limits[1],self.limits[3]], #xmax,ymax
			[self.limits[1],self.limits[2]] #xmax,ymin
			])
			ind=np.array([])
			for i in c_target:
				d=np.array([])
				for j in range(len(self.Outline[:,0])):
					d=np.append(d,
					np.sqrt((i[0]-self.Outline[j,0])**2+(i[1]-self.Outline[j,1])**2)
						)
				ind=np.append(ind,np.where(d==np.amin(d)))
		
		outlineCornerInd=np.sort(np.append(ind,0)).astype(int)
		
		
		#write routine for either 'spacing' based approach or total number of seeds
		#create empty array to receive respaced points
		
		conv=True
		count=0
		while conv:
			if self.ui.quantityButton.isChecked():
				if 'dist' not in vars():
					#find the perimeter first
					P=respace_equally(self.Outline,1)[1]

					#calculate mean distance between points
					numSeed=float(self.ui.numSeed.value())
					dist=P/float(numSeed) #minus 3 corners
					respacedOutline=np.array([]).reshape(0,2)
				else:
					dist=P/float(numSeed)
					respacedOutline=np.array([]).reshape(0,2)

				#now move through the corners to ID
				for j in range(3):
					nPts=respace_equally(self.Outline[outlineCornerInd[j]:outlineCornerInd[j+1]+1,0:2],dist)[-1]
					X=respace_equally(self.Outline[outlineCornerInd[j]:outlineCornerInd[j+1]+1,0:2],int(nPts+1))[0] #handles the lack of the last point for connectivity
					respacedOutline=np.vstack((respacedOutline,X[0:-1,:]))
				nPts=respace_equally(self.Outline[outlineCornerInd[3]::,0:2],dist)[-1]
				X=respace_equally(self.Outline[outlineCornerInd[3]::,0:2],int(nPts+1))[0]
				respacedOutline=np.vstack((respacedOutline,X[0:-1,:])) #last addition closes profile

			if self.ui.spacingButton.isChecked():
				if 'dist' not in vars():
					dist=float(self.ui.seedLengthInput.text())
				else:
					dist=dist+0.01*dist #increase by 1%

				respacedOutline=np.array([]).reshape(0,2)
				for j in range(3):
					X=respace_equally(self.Outline[outlineCornerInd[j]:outlineCornerInd[j+1]+1,0:2],dist)[0]
					respacedOutline=np.vstack((respacedOutline,X[0:-1,:]))

				X=respace_equally(self.Outline[outlineCornerInd[3]::,0:2],dist)[0]
				respacedOutline=np.vstack((respacedOutline,X[0:-1,:])) #do not close profile

			
			#write warning to the status line if 
			if not np.fmod(len(respacedOutline),2)==0:
				if self.ui.spacingButton.isChecked():
					dist+=0.25
				if self.ui.quantityButton.isChecked():
					self.ui.statLabel.setText("Odd number of outline seeds, retrying . . .")
					self.ui.vtkWidget.update()
					numSeed+=1
			else:
				conv=False
				self.ui.statLabel.setText("Idle") #clears error on recalculation
			count+=1
			
		#Write zeros to the z coordinate of the outline
		respacedOutline=np.hstack((respacedOutline,np.zeros([len(respacedOutline[:,0]),1])))

		#display and get respaced outline actor
		self.respacedOutlineActor, _ =gen_outline(respacedOutline,(1,0,0),self.PointSize+3)
		
		self.ren.AddActor(self.respacedOutlineActor)
		# self.respacedOutlineActor.GetProperty().SetPointSize(12)
		self.ui.vtkWidget.update()
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		self.AddAxis(self.limits,1)
		
		# self.respacedOutlineActor=self.DisplayOutline(respacedOutline,(1,0,0),self.PointSize+3)
		self.rsOutline=respacedOutline
		self.Dist=dist
		#update GUI to report back the number of points
		self.ui.seedLengthInput.setText("%4.4f"%dist)
		self.ui.numSeed.setValue(len(respacedOutline))

	# def DisplayOutline(self,pts,color,size):
		# points=vtk.vtkPoints()
		# for i in pts:
			# points.InsertNextPoint(i)
		# lineseg=vtk.vtkPolygon()
		# lineseg.GetPointIds().SetNumberOfIds(len(pts))
		# for i in range(len(pts)):
			# lineseg.GetPointIds().SetId(i,i)
		# linesegcells=vtk.vtkCellArray()
		# linesegcells.InsertNextCell(lineseg)
		# outline=vtk.vtkPolyData()
		# outline.SetPoints(points)
		# outline.SetVerts(linesegcells)
		# outline.SetLines(linesegcells)
		# self.Omapper=vtk.vtkPolyDataMapper()
		# self.Omapper.SetInputData(outline)
		# outlineActor=vtk.vtkActor()
		# outlineActor.SetMapper(self.Omapper)
		# outlineActor.GetProperty().SetColor(color)
		# outlineActor.GetProperty().SetPointSize(size)
		# self.ren.AddActor(outlineActor)
		# self.AddAxis(self.limits,1)
		# self.ren.ResetCamera()
		# self.ui.vtkWidget.update()
		# self.ren.ResetCamera()
		# return outlineActor

	def WriteOut(self):

		if hasattr(self,"rsOutline"):
			Outline=self.rsOutline
		else:
			Outline=self.Outline
		
		N=len(Outline)
		if self.ui.geoButton.isChecked():
			try:
				filename,self.ofile=get_open_file("*.geo",self.outputd)
				fid=open(filename,'w+')
			except:
				return
			pc=0
			lc=0
			for i in range(N):
				pc+=1
				outputline  = "Point(%i) = {%8.8f,%8.8f,%8.8f};\n" \
								% (pc,self.Outline[i,0],self.Outline[i,1],0)
				fid.write(outputline)
			for i in range(N-1):
				lc+=1
				outputline = "Line(%i) = {%i,%i};\n"%(lc,lc,lc+1)
				fid.write(outputline)
			lc+=1
			fid.write("Line(%i) = {%i,%i};\n"%(lc,lc,lc-N+1)) #last line to enclose

			fid.write("Line Loop(%i) = {%i:%i};\n" %(N+1,1,lc))
			sec=N+1
			#write plane
			fid.write("Plane Surface(%i) = {%i};\n" %(N+2,N+1))


		elif self.ui.dxfButton.isChecked():
			try:
				self.ofile_d,_=get_open_file("*.dxf",self.outputd)
				fid=open(self.ofile_d,'w+')
			except Exception as e:
				print("Didn't write dxf file.")
				return
			
			
			fid.write("0\nSECTION\n2\nENTITIES\n0\n")
			for i in range(N-1):
				outputline = "LINE\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n0\n" \
						%(Outline[i,0],Outline[i,1],0,Outline[i+1,0],Outline[i+1,1],0)
				fid.write(outputline)
			#last line
			fid.write("LINE\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n0\n" \
						%(Outline[-1,0],Outline[-1,1],0,Outline[0,0],Outline[0,1],0))
			fid.write("ENDSEC\n0\nEOF\n")
			fid.close()
	
	def RunMeshScript(self):
	
		if self.ui.gmshButton.isChecked():
			self.ui.statLabel.setText("Running Gmsh script . . .")
			QtWidgets.QApplication.processEvents()
			execStr=(self.cfg['FEA']['gmshExec'])
			self.vtkFile=self.geofile[0:-3]+"vtk"
			try:
				if self.ui.tetButton.isChecked():
					#make sure second order tets are generated
					out=sp.check_output([execStr,"-3","-order","2",self.geofile,"-o",self.vtkFile])
				else:
					out=sp.check_output([execStr,"-3",self.geofile,"-o",self.vtkFile])
				print("Gmsh output log:")
				print("----------------")
				print(out.decode("utf-8"))
				print("----------------")
				self.ui.statLabel.setText("Gmsh VTK file written . . . Idle")
			except sp.CalledProcessError as e:
				print("Gmsh command failed for some reason.")
				print(e.decode("utf-8"))
				self.ui.statLabel.setText("Gmsh call failed . . . Idle")
				
		#Run equivalent Abaqus command chain, with the addition to converting mesh to legacy VTK file
		if self.ui.abaButton.isChecked():
			self.ui.statLabel.setText("Running Abaqus CAE script . . .")
			QtWidgets.QApplication.processEvents()
			execStr=(self.cfg['FEA']['abaqusExec'])#.encode('string-escape')
			currentPath=os.getcwd()
			abaqusCAEfile=os.path.basename(self.abapyfile)
			abaqusrunloc=os.path.dirname(self.abapyfile)
			os.chdir(abaqusrunloc)
			try:
				self.ui.statLabel.setText("Writing Abaqus .inp file . . . ")
				out=sp.check_output([execStr,"cae","noGUI="+abaqusCAEfile],shell=True)
				print("Abaqus CAE log:")
				print("----------------")
				print(out.decode("utf-8"))
				print("----------------")
				self.ui.statLabel.setText("Converting Abaqus .inp file to VTK . . . ")
				#the CAE script will generate an input deck with the same prefix as the dxf file.
				self.vtkFile=self.ofile_d[0:-4]+"_inp.vtk"
				ConvertInptoVTK(self.ofile_d[0:-3]+"inp",self.vtkFile)
				self.ui.statLabel.setText("Abaqus .inp file converted . . . Idle")
				os.chdir(currentPath)
			except sp.CalledProcessError as e:
				print("Abaqus CAE command failed for some reason.")
				print(e.decode("utf-8"))
				self.ui.statLabel.setText("Abaqus CAE call failed . . . Idle")
			
			
	def WriteGeo(self):
		QtWidgets.QApplication.processEvents()
		if self.ui.gmshButton.isChecked():
			try:
				self.geofile,_=get_open_file("*.geo",self.outputd)
				fid=open(self.geofile,'w+')
			except:
				return
		else:
			if not hasattr(self,'ofile_d'):
				msg=QtWidgets.QMessageBox()
				msg.setIcon(QtWidgets.QMessageBox.Information)
				msg.setText("Write dxf outline prior to writing Abaqus Python script.")
				msg.setWindowTitle("pyCM Error")
				msg.exec_()
				return
			try:
				self.abapyfile,_=get_open_file("*.py",self.outputd)
				fid=open(self.abapyfile,'w+')
			except:
				return

		if hasattr(self,"rsOutline"):
			Outline=self.rsOutline[0:-1,0:-1]
		else:
			Outline=self.Outline[0:-1,0:-1]
		N=len(Outline)

		NumNodesDeep=self.ui.numPart.value()
		ExtrudeDepth=float(self.ui.lengthInput.text())
		
		cent=np.mean(Outline,axis=0)
		
		Bias_u=(ExtrudeDepth/float(self.Dist))**(1/float(NumNodesDeep-1)) #upper bound
		
		B_range=np.linspace(Bias_u/2,Bias_u,1000)
		Intersection=self.Dist*(1-np.power(B_range,NumNodesDeep))/(1-B_range)
		b=np.where(Intersection>ExtrudeDepth)
		Bias=B_range[b[0][0]]

		L=np.array([])
		for j in range(1,NumNodesDeep+1):
			L=np.append(L,self.Dist*(1-Bias**j)/(1-float(Bias)))

		L[-1]=ExtrudeDepth
		L=L/float(ExtrudeDepth)

		if hasattr(self,"geofile") and self.ui.gmshButton.isChecked():
			pc=0
			lc=0
			for i in range(N):
				pc+=1
				outputline  = "Point(%i) = {%8.8f,%8.8f,%8.8f};\n" \
								% (pc,Outline[i,0],Outline[i,1],0)
				fid.write(outputline)
			for i in range(N-1):
				lc+=1
				outputline = "Line(%i) = {%i,%i};\n"%(lc,lc,lc+1)
				fid.write(outputline)
			lc+=1
			fid.write("Line(%i) = {%i,%i};\n"%(lc,lc,lc-N+1)) #last line to enclose

			fid.write("Line Loop(%i) = {%i:%i};\n" %(N+1,1,lc))
			sec=N+1
			#write plane
			fid.write("Plane Surface(%i) = {%i};\n" %(N+2,N+1))
			if self.ui.quadButton.isChecked():
				fid.write("Recombine Surface {%i};\n\n" %(N+2)) #for quads, otherwise tets
			sec+=1
			
			fid.write("OutOfPlane[]= Extrude {0, 0, %8.8f} {\n Surface{%i};\n Layers{ {"%(ExtrudeDepth,sec)) 
			for i in range(len(L)-1):
				fid.write("1,")
			fid.write("1}, {")
			for i in range(len(L)-1):
				fid.write("%2.4f,"%L[i])
			if self.ui.quadButton.isChecked(): #quads vs. tets
				fid.write("%2.4f} };\n Recombine;};\n \n//EOF"%L[-1])
			else:
				fid.write("%2.4f} };};\n \n//EOF"%L[-1])
			
			self.ui.statLabel.setText("Gmsh geo file written . . . Idle")
			
		if hasattr(self,"abapyfile") and self.ui.abaButton.isChecked():
			self.ui.statLabel.setText("Writing Abaqus CAE script . . .")

			
			if self.ui.quadButton.isChecked():
				ElemType="C3D8"
			else:
				ElemType="C3D10"
			#move to metadata once dist channels are sorted
			s1="""
# RawInpWriter.py
# Abaqus python script to automatically generate C3D20R element mesh 
# Geometry and mesh is based on *.dxf file
# Produced by pyCM
####################################################################
# 03/12/13 MJR initial script format
import os
from abaqus import *
from abaqusConstants import *
from abaqus import backwardCompatibility
backwardCompatibility.setValues(reportDeprecated=False)

from caeModules import *
from driverUtils import executeOnCaeStartup
from dxf2abq import importdxf


# Parameters employed:\n"""

			s2="""
DXF_file=os.path.normpath(DXF_file)
executeOnCaeStartup()

Mdb()
importdxf(fileName=DXF_file)

s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)

s.retrieveSketch(sketch=mdb.models['Model-1'].sketches[os.path.basename(OutputFname)])

p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=Depth)
s.unsetPrimaryObject()

f = p.faces
faces = f.findAt((CentPoint, ))
p.Set(faces=faces, name='SURFACE')

del mdb.models['Model-1'].sketches['__profile__']

a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)

p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=OFF)


session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
pickedEdges2 = e1.findAt((EdgePoint, ))
a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges2, minSize=MinLength, 
    maxSize=MaxLength, constraint=FINER)

elemType1 = mesh.ElemType(elemCode=ElemType, elemLibrary=STANDARD)
a = mdb.models['Model-1'].rootAssembly
c1 = a.instances['Part-1-1'].cells
cells1 = c1.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells1, )
"""

			s3="""
a.setElementType(regions=pickedRegions, elemTypes=(elemType1,))
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-1-1'], )
a.generateMesh(regions=partInstances)

mdb.models['Model-1'].setValues(noPartsInputFile=ON)

mdb.Job(name=OutputFname, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
mdb.jobs[OutputFname].writeInput(consistencyChecking=OFF)

a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].view.fitView()\n"""

			fid.write("%s"%s1)
			fid.write("DXF_file=r'%s'\n" %self.ofile_d)
			fid.write("Depth=%8.8f\n" %ExtrudeDepth);
			fid.write("NN=%i\n" %NumNodesDeep);
			fid.write("MinLength=%8.8f\n" %(L[0]*ExtrudeDepth));
			if (L[-1]*ExtrudeDepth-L[-2]*ExtrudeDepth)>(L[0]*ExtrudeDepth):
				fid.write("MaxLength=%8.8f\n" %(L[-1]*ExtrudeDepth-L[-2]*ExtrudeDepth));
			else:
				fid.write("MaxLength=%8.8f\n" %(L[0]*ExtrudeDepth));
			fid.write("EdgePoint=(%8.8f,%8.8f,%8.8f)\n" %(Outline[0,0],Outline[0,1],L[0]*ExtrudeDepth));
			fid.write("CentPoint=(%8.8f,%8.8f,%8.8f)\n" %(cent[0],cent[1],0));
			fid.write("OutputFname='%s'\n"
			%(os.path.splitext(os.path.basename(self.ofile_d))[0]))
			fid.write("ElemType=%s\n"%ElemType)
			fid.write("%s"%s2)
			
			if self.ui.quadButton.isChecked():
				fid.write("a.setMeshControls(regions=cells1, algorithm=ADVANCING_FRONT)\n")
			else:
				fid.write("a.setMeshControls(regions=cells1, elemShape=TET, technique=FREE)\n")
			fid.write("%s"%s3)
			self.ui.statLabel.setText("Abaqus CAE script written . . . Idle")
		fid.close()
		
		
	def DisplayMesh(self):
	
		if hasattr(self,"meshActor"):
			self.ren.RemoveActor(self.meshActor)
			del self.meshSource
			del self.mesh
			if hasattr(self,"BCactor"):
				self.ren.RemoveActor(self.BCactor)
				del self.corners
			if hasattr(self,"pickedCornerInd"):
				self.UndoRigidBody()
				
			# if hasattr(self,"labelActor"):
				# self.ren.RemoveActor(self.labelActor) #debug
		if not hasattr(self,"vtkFile"):
			self.vtkFile, startdir = get_file('*.vtk')
		self.ui.statLabel.setText("Reading . . .")
		QtWidgets.QApplication.processEvents()
		self.meshSource=vtk.vtkUnstructuredGridReader()
		self.meshSource.SetFileName(self.vtkFile)
		self.meshSource.Update()
		self.ui.statLabel.setText("Filtering out non volumetric elements . . .")
		QtWidgets.QApplication.processEvents()
		#mesh as-read
		om = self.meshSource.GetOutput()
		#get cell types
		tcs=vtk.vtkCellTypes()
		om.GetCellTypes(tcs)
		#if it's not a 1st order quad or 2nd order tet, gmsh will return non uniform element types.
		if tcs.IsType(12)==1:
			self.mainCellType=12 #1st order quad
		elif tcs.IsType(24)==1:
			self.mainCellType=24 #2nd order tet
		# print "Cells before thresholding:",om.GetNumberOfCells() #debug
		#build int array of types
		cellTypeArray=vtk.vtkIntArray()
		cellTypeArray.SetName("Type")
		cellTypeArray.SetNumberOfComponents(1)
		cellTypeArray.SetNumberOfTuples(om.GetNumberOfCells())
		om.GetCellData().AddArray(cellTypeArray)
		for ind in range(om.GetNumberOfCells()):
			cellTypeArray.SetValue(ind,om.GetCell(ind).GetCellType())
		#generate threshold filter
		t=vtk.vtkThreshold()
		t.SetInputData(om)
		t.ThresholdByUpper(self.mainCellType)
		t.SetInputArrayToProcess(0,0,0,1,"Type")
		t.Update()
		
		self.mesh=t.GetOutput()
		# print "Cells after thresholding:",self.mesh.GetNumberOfCells() #debug

		self.ui.statLabel.setText("Rendering . . .")
		QtWidgets.QApplication.processEvents()
		
		# print "Read VTK mesh file:" #debug
		# print "No. points:",self.mesh.GetNumberOfPoints()
		# print "No. elements:",self.mesh.GetNumberOfCells()
		bounds=self.mesh.GetBounds()
		
		edges=vtk.vtkExtractEdges()
		edges.SetInputConnection(self.meshSource.GetOutputPort())
		edges.Update()

		self.meshMapper=vtk.vtkDataSetMapper()
		self.meshMapper.SetInputData(self.mesh)
		
		self.meshActor = vtk.vtkActor()
		self.meshActor.SetMapper(self.meshMapper)
		
		self.meshActor.GetProperty().SetLineWidth(1)
		self.meshActor.GetProperty().SetColor(0,0.9020,0.9020) #abaqus
		# self.meshActor.GetProperty().SetColor(0,1,0.6039) #gmsh
		self.meshActor.GetProperty().SetEdgeColor([0.8, 0.8, 0.8])
		self.meshActor.GetProperty().EdgeVisibilityOn()
		self.ren.AddActor(self.meshActor)
		#update z extents of interactor limits
		self.limits[4]=bounds[4]
		self.limits[5]=bounds[5]
		self.AddAxis(self.limits,1)
		self.ui.vtkWidget.update()
		self.ui.statLabel.setText("Mesh displayed . . . Idle")
		

	def ImposeSplineFit(self):
		"""
		Draws/identifies BCs from spline object read in from the .mat file. Also identifies candidate corners for rigid body BCs
		"""
		#make sure there's a mesh to work on
		if not hasattr(self,"meshActor"):
			msg=QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Information)
			msg.setText("Need a mesh before imposing BCs")
			msg.setWindowTitle("pyCM Error")
			msg.exec_()
			return
		
		self.ui.statLabel.setText("Locating surface elements . . .")
		QtWidgets.QApplication.processEvents()
		#create a locator from a bounding box for candidate cells.
		#Unless the mesh is highly refined, this locator will a few layers of elements from the z=0 plane
		vil = vtk.vtkIdList()
		locator = vtk.vtkCellTreeLocator()
		locator.SetDataSet(self.mesh)
		locator.BuildLocator()
		locator.FindCellsWithinBounds(self.mesh.GetBounds()[0:4]+(-0.1,self.Dist),vil)
		
		
		#vtk datatypes to hold info from locator filter
		nfaces=vtk.vtkCellArray()
		rptIds=vtk.vtkIdList()
		ptIds=vtk.vtkIdList()
		self.BCelements=np.array([])
		
		#push nodes of cells/elements id'ed into data structures
		count=0
		for i in range(vil.GetNumberOfIds()):
			self.mesh.GetFaceStream(vil.GetId(i),ptIds)
			cellType=self.mesh.GetCellType(vil.GetId(i))
			if cellType == self.mainCellType:
				count+=1
				nfaces.InsertNextCell(ptIds)
				self.BCelements=np.append(self.BCelements,vil.GetId(i))
				
		# for i in xrange(vilrf.GetNumberOfIds()):
			# print self.mesh.GetCellPoints(vil.GetId(i))
		# print(v2n(nfaces.GetData()))
		#convert the vtklist to numpy array, resize accordingly; 1D array consists of number of nodes/element, followed by node/point number
		rawPIds=v2n(nfaces.GetData()) 
		
		#make new matrix to hold node/point number connectivity
		if self.mainCellType == 12: #quads
			SurfPoints=np.resize(rawPIds,(int(len(rawPIds)/float(9)),9))
			BCunit=4
		elif self.mainCellType == 24: #2nd order tets
			SurfPoints=np.resize(rawPIds,(int(len(rawPIds)/float(11)),11))
			BCunit=6
		#remove point count column
		SurfPoints=SurfPoints[:,1::]
		
		#define vtk data structures for BC display
		BCpnts=vtk.vtkPoints()
		self.BCcells=vtk.vtkCellArray()
		
		#create array to append all of the points that are found so their id's can be used later to write specific BC's.
		self.BCindex=np.array([])
		
		#same for node label display
		self.BCnodeLabel=vtk.vtkStringArray()
		self.BCnodeLabel.SetNumberOfComponents(1)
		self.BCnodeLabel.SetName("NodeID")
		
		self.ui.statLabel.setText("Imposing nodal displacements . . .")
		QtWidgets.QApplication.processEvents()
		#build new mesh of shell elements to show the BC surface
		ccount=0
		ppcount=0
		for j in SurfPoints: #rows of surfpoints=elements
			localCell=vtk.vtkPolygon()
			#need to pre-allocate
			localCell.GetPointIds().SetNumberOfIds(BCunit)
			
			pcount=0
			for i in j:
				# count how many points there are with z==0; make sure it's a BCunit's worth
				if self.mesh.GetPoint(i)[2]==0:
					pcount+=1
			if pcount==BCunit: #then this isn't an element with just one node on the surface (tets)
				localcellind=0
				for i in j:
					#restart the process
					p=self.mesh.GetPoint(i)
					if p[2]==0:
						#make sure the pnt hasn't been added to the points
						if i not in self.BCindex:
							#add calculate values for z according to spline
							BCpnts.InsertNextPoint(np.append(p[:2],bisplev(p[0],p[1],self.tck)))
							self.BCindex=np.append(self.BCindex,i)
							localCell.GetPointIds().SetId(localcellind,ppcount)
							ppcount+=1
						else:
							localCell.GetPointIds().SetId(localcellind,
								np.where(self.BCindex==i)[0][0])
						localcellind+=1
						
				self.BCcells.InsertNextCell(localCell)
				ccount+=BCunit

		self.BCnodeLabel.SetNumberOfValues(len(self.BCindex))
		for j in range(len(self.BCindex)):
			self.BCnodeLabel.SetValue(j,str(int(self.BCindex[j])))

		self.ui.statLabel.setText("Rendering . . .")
		QtWidgets.QApplication.processEvents()
		BCPolyData=vtk.vtkPolyData()
		BCPolyData.SetPoints(BCpnts)
		BCPolyData.GetPointData().AddArray(self.BCnodeLabel)
		BCPolyData.SetPolys(self.BCcells)
		
		cBCPolyData=vtk.vtkCleanPolyData() #remove shared edges
		cBCPolyData.SetInputData(BCPolyData)

		BCmapper=vtk.vtkPolyDataMapper()

		BCmapper.SetInputConnection(cBCPolyData.GetOutputPort())
		
		# debug
		# pointLabels=vtk.vtkPointSetToLabelHierarchy()
		# pointLabels.SetInputData(BCPolyData)
		# pointLabels.SetLabelArrayName("NodeID")
		# pointLabels.GetTextProperty().SetColor(1, 0.0, 0.0)
		# pointLabels.GetTextProperty().BoldOn()
		# pointLabels.GetTextProperty().ItalicOn()
		# pointLabels.GetTextProperty().ShadowOn()
		# pointLabels.GetTextProperty().SetFontSize(5)
		# pointLabels.Update()
		# labelMapper = vtk.vtkLabelPlacementMapper()
		# labelMapper.SetInputConnection(pointLabels.GetOutputPort())
		# self.labelActor = vtk.vtkActor2D()
		# self.labelActor.SetMapper(labelMapper)
		# self.ren.AddActor2D(self.labelActor)
		# """
		
		self.BCactor=vtk.vtkActor()
		self.BCactor.SetMapper(BCmapper)

		self.BCactor.GetProperty().SetColor(1,0.804,0.204) #mustard
		self.BCactor.GetProperty().SetInterpolationToFlat()
		# self.BCactor.GetProperty().SetRepresentationToSurface() #??
		self.BCactor.GetProperty().EdgeVisibilityOn()
		self.meshActor.GetProperty().SetOpacity(0.5)
		self.ren.AddActor(self.BCactor)

		#make new ax3D to cover the extents of the BC surface
		self.BClimits=BCpnts.GetBounds()
		self.AddAxis(self.BClimits,1)
		self.ui.vtkWidget.update()

		self.ui.statLabel.setText("Finding corners . . .")
		QtWidgets.QApplication.processEvents()
		#create np matrix to store surface points (and find corners)
		self.BCpnts=v2n(BCpnts.GetData()) #BCsearch will be fast
		allNodes=v2n(self.mesh.GetPoints().GetData())

		c_target=np.array([
		[self.limits[0],self.limits[2]], #xmin,ymin
		[self.limits[0],self.limits[3]], #xmin,ymax
		[self.limits[1],self.limits[3]], #xmax,ymax
		[self.limits[1],self.limits[2]] #xmax,ymin
		])
		
		self.OutlineIsCCW=False #always will be false based on the order of c_target above
		
		ind=np.array([])
		for i in c_target:
			d=np.array([])
			for j in self.BCpnts:
				d=np.append(d,
				np.sqrt((i[0]-j[0])**2+(i[1]-j[1])**2))
			ind=np.append(ind,np.where(d==np.amin(d)))
		
		self.cornerInd=self.BCindex[ind.astype(int)]
		self.corners=self.BCpnts[ind.astype(int),:]
		
		# print(self.corners,self.cornerInd)
		
		#back face
		bfc_target=np.array([
		[self.limits[0],self.limits[2],self.limits[5]], #xmin,ymin,zmax
		[self.limits[0],self.limits[3],self.limits[5]], #xmin,ymax,zmax
		[self.limits[1],self.limits[3],self.limits[5]], #xmax,ymax,zmax
		[self.limits[1],self.limits[2],self.limits[5]] #xmax,ymin,zmax
		])
		
		
		#create point locator for nodes on back face
		backCornerLocator=vtk.vtkPointLocator()
		backCornerLocator.SetDataSet(self.mesh)
		backCornerLocator.AutomaticOn()
		backCornerLocator.BuildLocator()
		
		for i in bfc_target:
			target=backCornerLocator.FindClosestPoint(i)
			self.cornerInd=np.append(self.cornerInd,target)
			self.corners=np.vstack((self.corners,self.mesh.GetPoint(target)))
						
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		
		self.ui.statLabel.setText("Ready for rigid body BCs . . . Idle")
		QtWidgets.QApplication.processEvents()
		
		
	def Keypress(self,obj,event):
		key = obj.GetKeyCode()

		if key =="1":
			xyview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="2":
			yzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="3":
			xzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)

		elif key=="z":
			self.Zaspect=self.Zaspect*2
			# self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'BCactor'):
				self.BCactor.SetScale(1,1,self.Zaspect)
				nl=np.append(self.limits[0:4],[self.BClimits[-2]*self.Zaspect,self.BClimits[-1]*self.Zaspect])
			self.AddAxis(nl,1/self.Zaspect)

		elif key=="x":
			self.Zaspect=self.Zaspect*0.5
			# self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'BCactor'):
				self.BCactor.SetScale(1,1,self.Zaspect)
				nl=np.append(self.limits[0:4],[self.BClimits[-2]*self.Zaspect,self.BClimits[-1]*self.Zaspect])
			self.AddAxis(nl,1/self.Zaspect)

		elif key=="c":
			self.Zaspect=1.0
			# self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'BCactor'):
				self.BCactor.SetScale(1,1,self.Zaspect)
				nl=np.append(self.limits[0:4],self.BClimits[-2::])
				self.AddAxis(nl,1)

		elif key=="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(self.ui.vtkWidget._RenderWindow)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("mesh.png")
			writer.Write()
			print("Screen output saved to %s" %os.path.join(currentdir,'mesh.png'))
		
		elif key=="r":
			flip_visible(self.ax3D)
			
		elif key =="f": #flip color scheme for printing
			flip_colors(self.ren,self.ax3D)
				
		elif key == "o":
			flip_visible(self.outlineActor)

		elif key == "e":
			try:
				with open(self.filec,'r') as ymlfile:
					readcfg = yaml.load(ymlfile)
				l=[readcfg['FEA']['abaqusExec'],readcfg['FEA']['gmshExec'],readcfg['FEA']['ccxExec']]
				self.cfg=GetFEAconfig(l,self.filec)
			except:
				"Couldn't find config file where it normally is." 

		self.ui.vtkWidget.update()

	def AddAxis(self,limits,scale):
		if hasattr(self,"ax3D"):
			self.ren.RemoveActor(self.ax3D)
		self.ax3D = vtk.vtkCubeAxesActor()
		self.ax3D.ZAxisTickVisibilityOn()
		self.ax3D.SetXTitle('X')
		self.ax3D.SetYTitle('Y')
		self.ax3D.SetZTitle('Z')
		self.ax3D.SetBounds(limits)
		self.ax3D.SetZAxisRange(limits[-2]*scale,limits[-1]*scale)
		self.ax3D.SetCamera(self.ren.GetActiveCamera())
		self.ren.AddActor(self.ax3D)
		self.ax3D.SetFlyModeToOuterEdges()

	def doFEA(self):
		"""
		Packages up/writes either Abaqus or Calculix FEA input deck - subtle difference between the two in terms of how integration points are reported, however both share the same extension. Will call RunFEA if it's been specified
		"""
		#so either there isn't a 'pick' attribute or the number of picks is insufficient
		if not hasattr(self,"picks"):
			msg=QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Information)
			msg.setText("Need to impose complete BCs before conducting analysis.")
			msg.setWindowTitle("pyCM Error")
			msg.exec_()
			return
		elif not self.picks==2:
			msg=QtWidgets.QMessageBox()
			msg.setIcon(QtWidgets.QMessageBox.Information)
			msg.setText("Need to specify appropriate rigid body BCs before conducting analysis.")
			msg.setWindowTitle("pyCM Error")
			msg.exec_()
			return

		if self.ui.CalculixButton.isChecked():
			try:
				self.ofile_FEA,_=get_open_file("*.ccx.inp",self.outputd)
				fid=open(self.ofile_FEA,'wb+')
			except:
				return
		elif self.ui.AbaqusButton.isChecked():
			try:
				self.ofile_FEA,_=get_open_file("*.abq.inp",self.outputd)
				fid=open(self.ofile_FEA,'wb+')
			except:
				return
		
		nodes=v2n(self.mesh.GetPoints().GetData())
		nodes=np.column_stack((np.arange(1,len(nodes)+1),nodes+1))
		
		cells=v2n(self.mesh.GetCells().GetData())
		#determine element type based on first entry in cells, 8-C3D8, 10-C3D10
		elType=cells[0]
		cells=np.resize(cells+1,(int(len(cells)/float(elType+1)),elType+1))
		cells=np.column_stack((np.arange(1,len(cells)+1),cells[:,1::]))

		#'top' element set
		nR=len(self.BCelements) % 16 #max no of input entries/line
		
		if not nR==0:
			BCelemsq=np.reshape(self.BCelements[0:-nR],(int(len(self.BCelements[0:-nR])/float(16)),16))
			
		else: #well, the remainder is 0
			BCelemsq=np.reshape(self.BCelements,(int(len(self.BCelements)/float(16)),16))
		BCelemsq=BCelemsq+1 #because elements start numbering at 1
		
		fid.write(str.encode('*HEADING\n'))

		fid.write(str.encode('**pyCM input deck, converted from VTK format\n'))
		fid.write(str.encode('**%s\n'%self.ofile_FEA))

		#dump nodes
		fid.write(str.encode('*NODE\n'))
		np.savetxt(fid,nodes,fmt='%i,%.6f,%.6f,%.6f',delimiter=',')
		#dump 'cells'
		fid.write(str.encode('*ELEMENT, TYPE=C3D%i\n'%elType))
		np.savetxt(fid,cells,fmt='%i',delimiter=',')
		#generate element set to apply material properties
		fid.write(str.encode('*ELSET, ELSET=DOMAIN, GENERATE\n'))
		fid.write(str.encode('%i,%i,%i\n'%(1,len(cells),1)))
		fid.write(str.encode('*ELSET, ELSET=BC\n'))
		np.savetxt(fid,BCelemsq,fmt='%i',delimiter=',')
		if not nR==0:
			for i in self.BCelements[-nR:]:
				if i==self.BCelements[-1]:
					fid.write(str.encode('%i'%(i+1)))
				else:
					fid.write(str.encode('%i,'%(i+1)))
			fid.write(str.encode('\n'))
		#write/apply material properties
		fid.write(str.encode('*SOLID SECTION, ELSET=DOMAIN, MATERIAL=USERSPEC\n'))
		fid.write(str.encode('*MATERIAL, NAME=USERSPEC\n'))
		fid.write(str.encode('*ELASTIC, TYPE=ISO\n'))
		fid.write(str.encode('%7.0f,%.3f\n'%(float(self.ui.modulusInput.text()),float(self.ui.poissonInput.text()))))
		fid.write(str.encode('*STEP, NAME=CONFORM\n'))
		fid.write(str.encode('*STATIC\n'))
		fid.write(str.encode('*BOUNDARY\n'))
		fid.write(str.encode('%i, 1,2, 0\n'%(self.cornerInd[self.pickedCornerInd[0]]+1)))
		fid.write(str.encode('%i, 2, 0\n'%(self.cornerInd[self.pickedCornerInd[1]]+1)))
		fid.write(str.encode('*BOUNDARY\n'))
		for ind in range(len(self.BCindex)):
			fid.write(str.encode('%i, 3,, %6.6f\n'%(self.BCindex[ind]+1,self.BCpnts[ind,2])))
		fid.write(str.encode('*EL FILE\n'))
		fid.write(str.encode('S,E\n'))#get all stresses and strains just to be safe.
		fid.write(str.encode('*EL PRINT, ELSET=BC\n'))
		if self.ui.CalculixButton.isChecked():
			fid.write(str.encode('S\n'))#Coords by default
		elif self.ui.AbaqusButton.isChecked():
			fid.write(str.encode('COORD,S\n'))#have to specify coords
		fid.write(str.encode('*ENDSTEP'))
		
		fid.close()
		
		if self.ui.runFEAButton.isChecked():
			self.RunFEA()
		
	def RunFEA(self):
		'''
		Runs FEA according to specified method & entries in config file
		'''
		if self.ui.CalculixButton.isChecked():
			#get the exe from cfg
			execStr=(self.cfg['FEA']['ccxExec'])
			self.ui.statLabel.setText("Running Calculix . . .")
			QtWidgets.QApplication.processEvents()
			try:
				out=sp.check_output([execStr,"-i",self.ofile_FEA[:-4]])
				
				print("Calculix output log:")
				print("----------------")
				print(out.decode("utf-8"))
				print("----------------")
				self.ui.statLabel.setText("Calculix run completed . . . Idle")
			except sp.CalledProcessError as e:
				print("Calculix command failed for some reason.")
				print(e.decode("utf-8"))
				self.ui.statLabel.setText("Calculix call failed . . . Idle")
				
		if self.ui.AbaqusButton.isChecked():
			execStr=(self.cfg['FEA']['abaqusExec'])
			self.ui.statLabel.setText("Running Abaqus . . .")
			QtWidgets.QApplication.processEvents()
			try:
				currentPath=os.getcwd()
				abaqusrunloc,abaqusinpfile=os.path.split(self.ofile_FEA)
				os.chdir(abaqusrunloc)
				out=sp.check_output(["abaqus","job="+abaqusinpfile[:-4],"int"],shell=True)

				print("Abaqus output log:")
				print("----------------")
				print(out.decode("utf-8"))
				print("----------------")
				self.ui.statLabel.setText("Abaqus run completed . . . Idle")
				os.chdir(currentPath)
			except sp.CalledProcessError as e:
				print("Abaqus command failed for some reason.")
				print(e.decode("utf-8"))
				self.ui.statLabel.setText("Abaqus call failed . . . Idle")


def ConvertInptoVTK(infile,outfile):
	"""
	Converts abaqus inp file into a legacy ASCII vtk file. First order quads (C3D8) and third order tets (C3D10) are supported.
	"""
	fid = open(infile)
	
	#flags for identifying sections of the inp file
	inpKeywords=["*Node", "*Element", "*Nset", "*Elset"]
	
	#map abaqus mesh types to vtk objects
	vtkType={}
	vtkType['C3D8']=12
	vtkType['C3D10']=24

	#create counter for all lines in the inp file, and array to store their location
	i=0
	lineFlag=[];
	#read file and find both where keywords occur as well as the element type used
	while 1:
		lines = fid.readlines(100000)
		if not lines:
			break
		for line in lines:
			i+=1
			for keyword in inpKeywords:
				if line[0:len(keyword)]==keyword:
					lineFlag.append(i)
					if keyword=="*Element":
						line = line.replace("\n", "")
						CellNum=vtkType[line.split("=")[-1]]

	fid.close()
	#use genfromtxt to read between lines id'ed by lineFlag to pull in nodes and elements
	Nodes=np.genfromtxt(infile,skip_header=lineFlag[0],skip_footer=i-lineFlag[1]+1,delimiter=",")
	Elements=np.genfromtxt(infile,skip_header=lineFlag[1],skip_footer=i-lineFlag[2]+1,delimiter=",")
	#Now write it in VTK format to a new file starting with header
	fid=open(outfile,'wb+')
	fid.write(str.encode('# vtk DataFile Version 2.0\n'))
	fid.write(str.encode('%s,created by pyCM\n'%outfile[:-4]))
	fid.write(str.encode('ASCII\n'))
	fid.write(str.encode('DATASET UNSTRUCTURED_GRID\n'))
	fid.write(str.encode('POINTS %i double\n'%len(Nodes)))
	
	#dump nodes
	np.savetxt(fid,Nodes[:,1::],fmt='%.6f')
	fid.write(str.encode('\n'))
	fid.write(str.encode('CELLS %i %i\n'%(len(Elements),len(Elements)*len(Elements[0,:]))))
	#Now elements, stack the number of nodes in the element instead of the element number
	Cells=np.hstack((np.ones([len(Elements[:,0]),1])*len(Elements[0,1::]),Elements[:,1::]-1))
	np.savetxt(fid,Cells,fmt='%i')
	fid.write(str.encode('\n'))

	#Write cell types
	fid.write(str.encode('CELL_TYPES %i\n'%len(Elements)))
	CellType=np.ones([len(Elements[:,0]),1])*CellNum
	np.savetxt(fid,CellType,fmt='%i')

	fid.close()

def respace_equally(X,input):
	distance=np.sqrt(np.sum(np.diff(X,axis=0)**2,axis=1))
	s=np.insert(np.cumsum(distance),0,0)
	Perimeter=np.sum(distance)

	if not isinstance(input,(int)):
		nPts=round(Perimeter/input)
	else:
		nPts=input
	
	sNew=np.linspace(0,s[-1],nPts)
	fx = interp1d(s,X[:,0])
	fy = interp1d(s,X[:,1])
	
	Xnew=fx(sNew)
	Ynew=fy(sNew)
	
	X_new=np.stack((Xnew,Ynew),axis=-1)
	return X_new,Perimeter,nPts

def DrawArrow(startPoint,length,direction,renderer):
	"""
	Draws and scales an arrow with a defined starting point, direction and length, adds to the renderer, returns the actor
	"""
	arrowSource=vtk.vtkArrowSource()
	arrowSource.SetShaftRadius(0.12)
	arrowSource.SetTipRadius(0.35)
	arrowSource.SetTipLength(0.7)
	arrowSource.InvertOn()
	endPoint=startPoint+length*direction
	normalizedX=(endPoint-startPoint)/length

	
	arbitrary=np.array([1,1,1]) #can be replaced with a random vector
	normalizedZ=np.cross(normalizedX,arbitrary/np.linalg.norm(arbitrary))
	normalizedY=np.cross(normalizedZ,normalizedX)
	
	# Create the direction cosine matrix by writing values directly to an identity matrix
	matrix = vtk.vtkMatrix4x4()
	matrix.Identity()
	for i in range(3):
		matrix.SetElement(i, 0, normalizedX[i])
		matrix.SetElement(i, 1, normalizedY[i])
		matrix.SetElement(i, 2, normalizedZ[i])
		
	#Apply transforms
	transform = vtk.vtkTransform()
	transform.Translate(startPoint)
	transform.Concatenate(matrix)
	transform.Scale(length, length, length)
 
	# Transform the polydata
	transformPD = vtk.vtkTransformPolyDataFilter()
	transformPD.SetTransform(transform)
	transformPD.SetInputConnection(arrowSource.GetOutputPort())
	
	#Create mapper and actor
	mapper = vtk.vtkPolyDataMapper()
	mapper.SetInputConnection(transformPD.GetOutputPort())
	actor = vtk.vtkActor()
	actor.SetMapper(mapper)
	renderer.AddActor(actor)
	return actor
	
	
	
def GetFEAconfig(inputlist,filec):
	'''
	Creates a GUI window to let the user specify FEA executable paths and writes them to a config file. Reads configs.
	'''
	getFEAconfigDialog = QtWidgets.QDialog()

	dui = Ui_getFEAconfigDialog()
	dui.setupUi(getFEAconfigDialog)
	dui.abaExec.setText(inputlist[0])
	dui.gmshExec.setText(inputlist[1])
	dui.ccxExec.setText(inputlist[2])
	dui.ConfigFileLoc.setText(filec)

	getFEAconfigDialog.exec_()
	# getFEAconfigDialog.show()

	try:
		with open(filec,'r') as ymlfile:
			return yaml.load(ymlfile)
	except:
		print("Couldn't read config file for some reason.")

	
if __name__ == "__main__":
	# currentdir=os.getcwd()

	if len(sys.argv)>2:
		RefFile=sys.argv[1]
		outDir=sys.argv[2]
		FEAtool(RefFile,outDir)
	elif len(sys.argv)>1:
		RefFile=sys.argv[1]
		FEAtool(RefFile)
	else:
		FEAtool()

