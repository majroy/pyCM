#!/usr/bin/env python
'''
Uses VTK python to allow for fitting an averaged dataset associated with the
 contour method. Full interaction requires a 3-button mouse and keyboard.
-------------------------------------------------------------------------------
Current mapping is as follows:
LMB - rotate about point cloud centroid.
MMB - pan
RMB - zoom
1 - view 1, default, looks down z axis onto xy plane
2 - view 2, looks down x axis onto zy plane
3 - view 3, looks down y axis onto zx plane
z - increase z-aspect ratio
x - decrease z-aspect ratio
c - return to default z-aspect
Shift-z - increase size of points
Shift-x - decrease size of points
Shift-c - return to default point size
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
r - remove/reinstate compass/axes
o - remove/reinstate outline
e - write output
-------------------------------------------------------------------------------
ver 1.1 17-17-03
'''
__author__ = "M.J. Roy"
__version__ = "1.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"

import os,sys,time
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt4 import QtCore, QtGui
from PyQt4.QtGui import QApplication
import numpy as np
import scipy.io as sio
from scipy.interpolate import griddata,bisplrep,bisplev
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import Delaunay
from matplotlib import path
from pkg_resources import Requirement, resource_filename
from pyCMcommon import *


def fit_surface(*args, **kwargs):
	app = QApplication(sys.argv)
	
	spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
	splash_pix = QtGui.QPixmap(spl_fname,'PNG')
	splash = QtGui.QSplashScreen(splash_pix)
	splash.setMask(splash_pix.mask())

	splash.show()
	app.processEvents()
	
	window = sf_interactor()

	if len(args)==1:
		sf_interactor.get_input_data(window,args[0])
	else:
		sf_interactor.get_input_data(window,None)
	
	window.show()
	splash.finish(window)
	window.iren.Initialize() # Need this line to actually show the render inside Qt
	
	ret = app.exec_()
	
	if sys.stdin.isatty() and not hasattr(sys,'ps1'):
		sys.exit(ret)
	else:
		return window

class sf_MainWindow(object):

	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.setWindowTitle("pyCM - surface fitting v%s" %__version__)
		MainWindow.resize(1280, 720)
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.Boxlayout = QtGui.QHBoxLayout(self.centralWidget)
		self.Subtendlayout=QtGui.QVBoxLayout()
		splineBox = QtGui.QFormLayout()
		sectionBox = QtGui.QGridLayout()


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
		
		# self.Boxlayout.addWidget(self.vtkWidget)
		# self.Boxlayout.addStretch()
		# MainWindow.setCentralWidget(self.centralWidget)
		self.reloadButton = QtGui.QPushButton('Load')
		splineLabel=QtGui.QLabel("Bivariate spline fitting")
		splineLabel.setFont(QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold))
		self.numLabel1=QtGui.QLabel("Knot spacing (x)")
		self.numEdit1 = QtGui.QLineEdit()
		self.numEdit1.setText('4')
		self.numEdit1.setMinimumWidth(50)
		
		self.numLabel2=QtGui.QLabel("Knot spacing (y)")
		self.numEdit2 = QtGui.QLineEdit()
		self.numEdit2.setText('4')
		self.numEdit2.setMinimumWidth(50)
		
		self.numLabel3=QtGui.QLabel("Spline order (x)")
		self.numEdit3 = QtGui.QSpinBox()
		self.numEdit3.setValue(3)
		self.numEdit3.setMinimum(1)
		self.numEdit3.setMaximum(5)
		self.numEdit3.setMinimumWidth(50)
		
		self.numLabel4=QtGui.QLabel("Spline order (y)")
		self.numEdit4 = QtGui.QSpinBox()
		self.numEdit4.setValue(3)
		self.numEdit4.setMinimum(1)
		self.numEdit4.setMaximum(5)
		self.numEdit4.setMinimumWidth(50)
		
		self.numLabel5=QtGui.QLabel("Smoothing:")
		self.numEdit5 = QtGui.QLineEdit()
		self.numEdit5.setText('0')
		self.numEdit5.setMinimumWidth(50)
		
		self.numLabel6=QtGui.QLabel("Display resolution:")
		self.drx=QtGui.QLineEdit()
		drx_label=QtGui.QLabel("x")

		self.yMin=QtGui.QLineEdit()
		self.dry=QtGui.QLineEdit()
		dry_label=QtGui.QLabel("y")
		
		
		self.updateButton = QtGui.QPushButton('Update')
		self.updateButton.setMinimumWidth(50)
		
		self.statusLabel=QtGui.QLabel("Idle")
		self.statusLabel.setWordWrap(True)
		self.statusLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.statusLabel.setMinimumWidth(50)
		horizLine1=QtGui.QFrame()
		horizLine1.setFrameStyle(QtGui.QFrame.HLine)
		horizLine2=QtGui.QFrame()
		horizLine2.setFrameStyle(QtGui.QFrame.HLine)
		horizLine3=QtGui.QFrame()
		horizLine3.setFrameStyle(QtGui.QFrame.HLine)
		horizLine4=QtGui.QFrame()
		horizLine4.setFrameStyle(QtGui.QFrame.HLine)		

		sectionLabel=QtGui.QLabel("Data sectioning")
		sectionLabel.setFont(QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold))
		sectionIntList=[]
		self.xMin=QtGui.QLineEdit()
		sectionIntList.append(self.xMin)
		x0_label=QtGui.QLabel("x0")

		self.xMax=QtGui.QLineEdit()
		sectionIntList.append(self.xMax)
		x1_label=QtGui.QLabel("x1")

		self.yMin=QtGui.QLineEdit()
		y0_label=QtGui.QLabel("y0")
		sectionIntList.append(self.yMin)
		
		self.yMax=QtGui.QLineEdit()
		y1_label=QtGui.QLabel("y1")
		sectionIntList.append(self.yMax)
		self.sectionButton = QtGui.QPushButton('Section')
		self.revertButton = QtGui.QPushButton('Revert')
		self.writeButton = QtGui.QPushButton('Write')
		self.writeButton.setMinimumWidth(50)

		#splineBox is the main container, the sectionBox is nested within
		splineBox.addRow(self.reloadButton)
		splineBox.addRow(horizLine1)
		splineBox.addRow(splineLabel)
		splineBox.addRow(self.numLabel1,self.numEdit1)
		splineBox.addRow(self.numLabel2,self.numEdit2)
		splineBox.addRow(self.numLabel3,self.numEdit3)
		splineBox.addRow(self.numLabel4,self.numEdit4)
		splineBox.addRow(self.numLabel5,self.numEdit5)
		splineBox.addRow(self.numLabel6)
		splineBox.addRow(drx_label,self.drx)
		splineBox.addRow(dry_label,self.dry)
		splineBox.addRow(self.updateButton)
		
		splineBox.addRow(horizLine2)
		splineBox.addRow(sectionBox)
		splineBox.addRow(horizLine3)
		splineBox.addRow(self.writeButton)
		splineBox.addRow(horizLine4)
		splineBox.addRow(self.statusLabel)

		
		sectionBox.addWidget(sectionLabel,1,1,1,4)
		sectionBox.addWidget(x0_label,2,1)
		sectionBox.addWidget(sectionIntList[0],2,2)
		sectionBox.addWidget(x1_label,2,3)
		sectionBox.addWidget(sectionIntList[1],2,4)
		sectionBox.addWidget(y0_label,3,1)
		sectionBox.addWidget(sectionIntList[2],3,2)
		sectionBox.addWidget(y1_label,3,3)
		sectionBox.addWidget(sectionIntList[3],3,4)
		sectionBox.addWidget(self.sectionButton,5,1,1,2)
		sectionBox.addWidget(self.revertButton,5,3,1,2)
		self.Boxlayout.addLayout(splineBox)

	
	def initialize(self):
		self.vtkWidget.start()


class sf_interactor(QtGui.QMainWindow):

	def __init__(self, parent = None):
		QtGui.QMainWindow.__init__(self, parent)
		self.ui = sf_MainWindow()
		self.ui.setupUi(self)
		self.ren = vtk.vtkRenderer()
		self.ren.SetBackground(0.1, 0.2, 0.4)
		self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
		self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
		self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
		self.ren.GetActiveCamera().ParallelProjectionOn()
		self.cp=self.ren.GetActiveCamera().GetPosition()
		self.fp=self.ren.GetActiveCamera().GetFocalPoint()
		self.iren.AddObserver("KeyPressEvent", self.keypress)

		self.PointSize=2
		self.LineWidth=1
		self.Zaspect=1.0
		self.limits=np.empty(6)

		self.ui.reloadButton.clicked.connect(lambda: self.get_input_data(None))
		self.ui.updateButton.clicked.connect(lambda: self.onUpdateSpline())
		self.ui.sectionButton.clicked.connect(lambda: self.Cut())
		self.ui.revertButton.clicked.connect(lambda: self.RemoveCut())
		self.ui.writeButton.clicked.connect(lambda: self.WriteOutput())


	def get_input_data(self,filer):
		"""
		Loads the content of a *.mat file pertaining to this particular step
		"""
		
		if hasattr(self,'pointActor'): #then remove everything
			self.ren.RemoveActor(self.pointActor)
			self.ren.RemoveActor(self.outlineActor)
		
		if hasattr(self,'splineActor'):
			self.ren.RemoveActor(self.splineActor)

		
		if filer == None:
			filer, _, =get_file('*.mat')
		
		if filer: #check variables
			mat_contents = sio.loadmat(filer)
			self.fileo=filer
			try:
				pts=mat_contents['aa']

				self.pts=pts[~np.isnan(pts).any(axis=1)] #remove all nans
				self.RefOutline=np.concatenate(mat_contents['ref']['x_out'],axis=0)[0]
				RefMin=np.amin(self.RefOutline,axis=0)
				RefMax=np.amax(self.RefOutline,axis=0)
				self.limits=[RefMin[0],RefMax[0],RefMin[1],RefMax[1],np.amin(self.pts[:,-1],axis=0),np.amax(self.pts[:,-1],axis=0)]
				self.RefMin,self.RefMax=RefMin,RefMax
				self.ui.xMin.setText('%.3f'%self.limits[0])
				self.ui.xMax.setText('%.3f'%self.limits[1])
				self.ui.yMin.setText('%.3f'%self.limits[2])
				self.ui.yMax.setText('%.3f'%self.limits[3])
			
				#Generate actors
				color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
				_, self.pointActor, _, = gen_point_cloud(self.pts,color,self.PointSize)
				self.ren.AddActor(self.pointActor)
				self.outlineActor, _, = gen_outline(self.RefOutline,color,self.PointSize)
				self.ren.AddActor(self.outlineActor)

				#add axes
				self.add_axis(self.limits,[1,1,1])
				
			except: #Exception as e: print str(e)
				
				print "Couldn't read in both sets of data."
			
		else:
			print 'Invalid *.mat file'
			return
		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
		# update status
		self.ui.activeFileLabel.setText("Current analysis file:%s"%filer)
			
	def onUpdateSpline(self):
		p,ro,rmin,rmax=self.pts,self.RefOutline,self.RefMin,self.RefMax
		self.ui.statusLabel.setText("Fitting . . .")
		QtGui.qApp.processEvents()
		#read input parameters
		gx=float(self.ui.numEdit1.text())
		gy=float(self.ui.numEdit2.text())
		kx=self.ui.numEdit3.value()
		ky=self.ui.numEdit4.value()
		s=float(self.ui.numEdit5.text())

		tx=np.linspace(rmin[0],rmax[0],int((rmax[0]-rmin[0])/gx))
		ty=np.linspace(rmin[1],rmax[1],int((rmax[1]-rmin[1])/gy))

		#make sure both x & y have enough values in either direction
		if len(tx)<3 or len(ty)<3:
			self.ui.statusLabel.setText("Grid too large . . .")
			self.ui.updateButton.setEnabled(True)
			return
		tx=np.insert(tx,0,[rmin[0]] * kx) #to make sure knots are repeated at edges
		tx=np.insert(tx,-1,[rmax[0]] * kx)
		ty=np.insert(ty,0,[rmin[1]] * ky) #to make sure knots are repeated at edges
		ty=np.insert(ty,-1,[rmax[1]] * ky)
		
		#try spline fitting
		try:
			self.tck = bisplrep(p[:,0], p[:,1], p[:,2], kx=kx, ky=ky, tx=tx, ty=ty, task=-1) #get spline representation

			#now evaluate for show
			if not hasattr(self,'dryval'): #then it won't have drxval
				rx=np.linspace(rmin[0],rmax[0],int((rmax[0]-rmin[0])/(gx/2)))
				ry=np.linspace(rmin[1],rmax[1],int((rmax[1]-rmin[1])/(gy/2)))
				self.drxval,self.dryval=rx[1]-rx[0],ry[1]-ry[0]
				self.ui.drx.setText('%.3f'%(self.drxval))
				self.ui.dry.setText('%.3f'%(self.dryval))
			else: #read the res directly from UI, apply and update
				self.drxval,self.dryval=float(self.ui.drx.text()),float(self.ui.dry.text())
				rx=np.linspace(rmin[0],rmax[0],int((rmax[0]-rmin[0])/(self.drxval)))
				ry=np.linspace(rmin[1],rmax[1],int((rmax[1]-rmin[1])/(self.dryval)))
				self.drxval,self.dryval=rx[1]-rx[0],ry[1]-ry[0]
				self.ui.drx.setText('%.3f'%(self.drxval))
				self.ui.dry.setText('%.3f'%(self.dryval))
				
			grid_x, grid_y = np.meshgrid(rx,ry,indexing='xy')
			grid_x=np.transpose(grid_x)
			grid_y=np.transpose(grid_y)

			points2D=np.vstack([grid_x.flatten(),grid_y.flatten()]).T
			tri=Delaunay(points2D)
			tri= tri.simplices
			
			znew = bisplev(grid_x[:,0], grid_y[0,:], self.tck)

			points3D=np.vstack([points2D[:,0],points2D[:,1],znew.flatten()]).T
			
			#crop according to what's in the outline
			#temporarily move the outline to the first quadrant
			t_offset=np.array([0.1*rmin[0],0.1*rmin[1]])
			
			for j in range(0,2):
				if t_offset[j]<0:
					points3D[:,j]=points3D[:,j]-t_offset[j]
					ro[:,j]=ro[:,j]-t_offset[j]


			pth=path.Path(ro[:,:2])
			inOutline=pth.contains_points(points3D[:,:2])

			ind=np.array(range(0,len(inOutline)))
			ind_out=ind[~inOutline]


			for i in ind_out:
				tri=tri[np.where(~np.any(tri==i,axis=1))] #pull out triangles which contain points outside
			
			#move things back
			for j in range(0,2):
				if t_offset[j]<0:
					points3D[:,j]=points3D[:,j]+t_offset[j]
					ro[:,j]=ro[:,j]+t_offset[j]

			self.ui.statusLabel.setText("Rendering . . .")
			self.DisplaySplineFit(points3D,tri)

			zeval = np.empty(np.size(p[:,2]))
			for i in range(len(zeval)):
				zeval[i] = bisplev(p[i,0], p[i,1], self.tck)

			a=(zeval-p[:,2])
			RSME=(np.sum(a*a)/len(a))**0.5

			self.ui.statusLabel.setText("RSME: %2.2f micron . . . Idle"%(RSME*1000))

		except ValueError as ve:
			print ve
			splineFail= True
			self.ui.statusLabel.setText("Last fit failed . . . Idle")
		
		self.ui.updateButton.setEnabled(True)

	def DisplaySplineFit(self,p,t):

		if hasattr(self,'splineActor'):
			self.ren.RemoveActor(self.splineActor)
		self.SplinePoints = vtk.vtkPoints()
		triangles = vtk.vtkCellArray()
		
		#load up points
		for i in p:
			self.SplinePoints.InsertNextPoint(i)

		for i in t:
			triangle=vtk.vtkTriangle()
			for j in range(0,3):
				triangle.GetPointIds().SetId(j,i[j])
			triangles.InsertNextCell(triangle)

		trianglePolyData = vtk.vtkPolyData()
		trianglePolyData.SetPoints(self.SplinePoints)
		trianglePolyData.SetPolys(triangles)
		#filter so that edges are shared
		self.cSplinePolyData = vtk.vtkCleanPolyData()
		self.cSplinePolyData.SetInputData(trianglePolyData)

		# Create a mapper and actor for smoothed dataset
		self.Smapper = vtk.vtkPolyDataMapper()
		self.Smapper.SetInputConnection(self.cSplinePolyData.GetOutputPort())

		self.splineActor = vtk.vtkActor()
		self.splineActor.SetMapper(self.Smapper)
		self.splineActor.GetProperty().SetInterpolationToFlat()
		self.splineActor.GetProperty().SetRepresentationToSurface()
		self.splineActor.GetProperty().SetColor(1,0.804,0.204)
		self.splineActor.GetProperty().SetOpacity(0.75)
		self.splineActor.GetProperty().SetLighting(False)
		self.splineActor.SetScale(1,1,self.Zaspect)
		self.ren.AddActor(self.splineActor)
		self.ui.vtkWidget.update()

	def Cut(self):
		pts=np.array([float(self.ui.xMin.text()),float(self.ui.xMax.text()),float(self.ui.yMin.text()),float(self.ui.yMax.text())])

		planex1 = vtk.vtkPlane()
		planex1.SetOrigin(pts[0],0,0)
		planex1.SetNormal(1,0,0)
		
		planex2 = vtk.vtkPlane()
		planex2.SetOrigin(pts[1],0,0)
		planex2.SetNormal(-1,0,0)
		
		planey1 = vtk.vtkPlane()
		planey1.SetOrigin(0,pts[2],0)
		planey1.SetNormal(0,1,0)
		
		planey2 = vtk.vtkPlane()
		planey2.SetOrigin(0,pts[3],0)
		planey2.SetNormal(0,-1,0)
		
		planeCollection = vtk.vtkPlaneCollection()
		planeCollection.AddItem(planex1)
		planeCollection.AddItem(planex2)
		planeCollection.AddItem(planey1)
		planeCollection.AddItem(planey2)
		
		Omapper=self.outlineActor.GetMapper()
		Pmapper=self.pointActor.GetMapper()
		Omapper.SetClippingPlanes(planeCollection)
		Pmapper.SetClippingPlanes(planeCollection)

		if hasattr(self,'Smapper'):
			self.Smapper.SetClippingPlanes(planeCollection)
		nl=np.array(self.ax3D.GetBounds())
		self.ren.RemoveActor(self.ax3D)
		#add axes
		self.add_axis(self.limits,[1,1,1])

		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
	
	def RemoveCut(self):
		self.ui.xMin.setText('%.3f'%self.limits[0])
		self.ui.xMax.setText('%.3f'%self.limits[1])
		self.ui.yMin.setText('%.3f'%self.limits[2])
		self.ui.yMax.setText('%.3f'%self.limits[3])
		Omapper=self.outlineActor.GetMapper()
		Pmapper=self.pointActor.GetMapper()
		Omapper.RemoveAllClippingPlanes()
		Pmapper.RemoveAllClippingPlanes()
		if hasattr(self,'Smapper'):
			self.Smapper.RemoveAllClippingPlanes()
		self.ren.RemoveActor(self.ax3D)
		#add axes
		self.add_axis(self.limits,[1,1,1])

		
		#update
		self.ren.ResetCamera()
		self.ui.vtkWidget.update()
		self.ui.vtkWidget.setFocus()
	
	def DisplayPointCloud(self,pts,min,max):

		self.p = vtk.vtkPoints()
		vertices = vtk.vtkCellArray()


		#load up points
		for i in pts:
			pId= self.p.InsertNextPoint(i)
			vertices.InsertNextCell(1)
			vertices.InsertCellPoint(pId)
			
		pC = vtk.vtkPolyData()
		pC.SetPoints(self.p)
		pC.SetVerts(vertices)
		
		
		self.Pmapper = vtk.vtkDataSetMapper()
		self.Pmapper.SetInputData(pC)

		self.pointCloud = vtk.vtkCleanPolyData()
		self.pointCloud.SetInputData(pC)
		self.Pmapper.SetInputConnection(self.pointCloud.GetOutputPort())
		self.pointActor=vtk.vtkActor()
		self.pointActor.SetMapper(self.Pmapper)

		self.pointActor.GetProperty().SetColor(0.2784,0.6745,0.6941)
		
		self.ren.AddActor(self.pointActor)
		limits=np.array([min[0],max[0],min[1],max[1],min[2],max[2]])
		self.AddAxis(limits,1)
		self.ui.vtkWidget.update()
		self.ren.ResetCamera()

	def DisplayOutline(self,pts):
		points=vtk.vtkPoints()
		for i in pts:
			points.InsertNextPoint(i)
		lineseg=vtk.vtkPolygon()
		lineseg.GetPointIds().SetNumberOfIds(len(pts))
		for i in range(len(pts)):
			lineseg.GetPointIds().SetId(i,i)
		linesegcells=vtk.vtkCellArray()
		linesegcells.InsertNextCell(lineseg)
		outline=vtk.vtkPolyData()
		outline.SetPoints(points)
		outline.SetLines(linesegcells)
		self.Omapper=vtk.vtkPolyDataMapper()
		self.Omapper.SetInputData(outline)
		self.outlineActor=vtk.vtkActor()
		self.outlineActor.SetMapper(self.Omapper)
		self.outlineActor.GetProperty().SetColor(0.2784,0.6745,0.6941)
		self.ren.AddActor(self.outlineActor)
		self.ui.vtkWidget.update()


	def keypress(self,obj,event):
		key = obj.GetKeyCode()

		if key =="1":
			XYView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="2":
			YZView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="3":
			XZView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)

		elif key=="z":
			self.Zaspect=self.Zaspect*2
			s,nl,axs=self.get_scale()
			if hasattr(self,'splineActor'):
				self.splineActor.SetScale(s)
				self.splineActor.Modified()
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(s)
				self.pointActor.Modified()

			self.add_axis(nl,axs)

		elif key=="x":
			self.Zaspect=self.Zaspect*0.5
			s,nl,axs=self.get_scale()
			if hasattr(self,'splineActor'):
				self.splineActor.SetScale(s)
				self.splineActor.Modified()
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(s)
				self.pointActor.Modified()

			self.add_axis(nl,axs)

		elif key=="c":
			self.Zaspect=1.0
			s,_,_,=self.get_scale()
			if hasattr(self,'splineActor'):
				self.splineActor.SetScale(s)
				self.splineActor.Modified()
			if hasattr(self,'pointActor'):
				self.pointActor.SetScale(s)
				self.pointActor.Modified()

			self.add_axis(self.limits,[1,1,1])
			self.ren.ResetCamera()

		elif key=="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(self.ui.vtkWidget._RenderWindow)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("spline_fit.png")
			writer.Write()
			print 'Screen output saved to %s' %os.path.join(currentdir,'spline_fit.png')
		
		elif key=="r":
			FlipVisible(self.ax3D)
			
		elif key =="f": #flip color scheme for printing
			FlipColors(self.ren,self.ax3D,None)
			FlipColors(self.ren,self.pointActor,1)
			if hasattr(self,'splineActor'):
				FlipColors(self.ren,self.splineActor,0)
				
		elif key == "o":
			FlipVisible(self.outlineActor)
			
		elif key == "Z":
			self.PointSize=updatePointSize(self.pointActor,self.PointSize*2)
			if hasattr(self,'splineActor'):
				self.LineWidth=updateLineWidth(self.splineActor,self.LineWidth*2)
			
		elif key == "X":
			self.PointSize=updatePointSize(self.pointActor,self.PointSize*0.5)
			if hasattr(self,'splineActor'):
				self.LineWidth=updateLineWidth(self.splineActor,self.LineWidth*0.5)
			
		elif key == "C":
			self.PointSize=updatePointSize(self.pointActor,1)
			if hasattr(self,'splineActor'):
				self.LineWidth=updateLineWidth(self.splineActor,1)

		elif key == "e":
			self.WriteOutput()
		
		elif key == "q":
			if sys.stdin.isatty():
				sys.exit("Surface fitting complete.")
			else:
				print 'Surface fitting completed.'
				return

		self.ui.vtkWidget.update()

	def WriteOutput(self):
		
		if hasattr(self,'splineActor'): #then spline fitting has been done
			mat_contents=sio.loadmat(self.fileo)
			
			coefs=[np.reshape(self.tck[2],(len(self.tck[0])-self.tck[3]-1,-1))]
			number=np.array([len(self.tck[0]),len(self.tck[1])])
			order=np.array([self.tck[3], self.tck[4]])
			new={'spline_x': {'form': 'B-', 'knots': [self.tck[0], self.tck[1]], 'coefs': coefs, 'number': number, 'order':order, 'dim': 1, 'tck': self.tck},  'x_out':self.RefOutline}
			
			mat_contents.update(new)
			
			sio.savemat(self.fileo,mat_contents)
			
			self.ui.statusLabel.setText("Output written to %s. Idle." %self.fileo)
		else:
			self.ui.statusLabel.setText("Nothing to write. Idle.")

	def get_scale(self):
		'''
		Returns array for the keypress function
		'''
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
		self.ax3D.SetCamera(self.ren.GetActiveCamera())
		self.ren.AddActor(self.ax3D)
		if not(self.ren.GetBackground()==(0.1, 0.2, 0.4)):
			flip_colors(self.ren,self.ax3D)

def GetFile():
	root = tk.Tk()
	root.withdraw()
	filer = askopenfilename(title='Select the AVERAGED data file:',
		initialdir=os.getcwd(),
		filetypes =(("MAT File", "*.mat"),("All Files","*.*")))
	
	if filer == '':
		if sys.stdin.isatty() and not hasattr(sys,'ps1'):
			sys.exit("No file selected; exiting.")
		else:
			print 'No file selected; exiting.'
			filer = None
			startdir = None
		
	else:
		startdir = os.path.dirname(filer)
	return filer, startdir


def XYView(renderer, camera,cp,fp):
	camera.SetPosition(0,0,cp[2]+0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,1,0)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()

def YZView(renderer, camera,cp,fp):
	camera.SetPosition(cp[2]+0,0,0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()


def XZView(renderer,camera,cp,fp):
	vtk.vtkObject.GlobalWarningDisplayOff() #otherwise there's crystal eyes error . . .
	camera.SetPosition(0,cp[2]+0,0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()
	


def FlipVisible(actor):
	if actor.GetVisibility():
		actor.VisibilityOff()
	else:
		actor.VisibilityOn()

def updatePointSize(actor,NewPointSize):
	actor.GetProperty().SetPointSize(NewPointSize)
	actor.Modified()
	return NewPointSize

def updateLineWidth(actor,NewLineWidth):
	actor.GetProperty().SetLineWidth(NewLineWidth)
	actor.Modified()
	return NewLineWidth

def FlipColors(ren,actor,contrast):
	if ren.GetBackground()==(0.1, 0.2, 0.4):
		if hasattr(actor,'GetXAxesLinesProperty'):
			actor.GetTitleTextProperty(0).SetColor(0,0,0)
			actor.GetLabelTextProperty(0).SetColor(0,0,0)
			actor.GetXAxesLinesProperty().SetColor(0,0,0)
			actor.SetXTitle('x') #there's a vtk bug here . . .
			
			actor.GetTitleTextProperty(1).SetColor(0,0,0)
			actor.GetLabelTextProperty(1).SetColor(0,0,0)
			actor.GetYAxesLinesProperty().SetColor(0,0,0)

			actor.GetTitleTextProperty(2).SetColor(0,0,0)
			actor.GetLabelTextProperty(2).SetColor(0,0,0)
			actor.GetZAxesLinesProperty().SetColor(0,0,0)
			ren.SetBackground(1, 1, 1)
		else:
			if contrast == 1:
				actor.GetProperty().SetColor(0.2784,0.6745,0.6941)
			else:
				actor.GetProperty().SetColor(1,0.804,0.204)

	else:
		if hasattr(actor,'GetXAxesLinesProperty'):
			actor.GetTitleTextProperty(0).SetColor(1,1,1)
			actor.GetLabelTextProperty(0).SetColor(1,1,1)
			actor.GetXAxesLinesProperty().SetColor(1,1,1)
			actor.SetXTitle('X')
			
			actor.GetTitleTextProperty(1).SetColor(1,1,1)
			actor.GetLabelTextProperty(1).SetColor(1,1,1)
			actor.GetYAxesLinesProperty().SetColor(1,1,1)
			actor.SetYTitle('Y')
			
			actor.GetTitleTextProperty(2).SetColor(1,1,1)
			actor.GetLabelTextProperty(2).SetColor(1,1,1)
			actor.GetZAxesLinesProperty().SetColor(1,1,1)
			ren.SetBackground(0.1, 0.2, 0.4)
		else:
			if contrast == 1:
				actor.GetProperty().SetColor(0.0353, 0.1922, 0.2706)
			else:
				actor.GetProperty().SetColor(0.8039, 0.3490, 0.2902)

if __name__ == "__main__":
	if len(sys.argv)>2:
		sf_interactor(sys.argv[1])
	else:
		sf_interactor()