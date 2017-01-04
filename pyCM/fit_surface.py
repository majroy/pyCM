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
e - write output and exit
-------------------------------------------------------------------------------
ver 0.1 16-11-26
'''
__author__ = "M.J. Roy"
__version__ = "0.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"


from PyQt4 import QtCore, QtGui
# from PyQt4 import *
from PyQt4.QtGui import QApplication
import vtk
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import os,sys,time
import numpy as np
import scipy.io as sio
from scipy.interpolate import griddata,bisplrep,bisplev
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import Delaunay
from matplotlib import path
import Tkinter as tk
from tkFileDialog import askopenfilename
from tkFileDialog import askdirectory


def fit_surface(*args, **kwargs):
	app = QApplication(sys.argv)
	
	splash_pix = QtGui.QPixmap('meta/pyCM_logo.png')
	splash = QtGui.QSplashScreen(splash_pix,QtCore.Qt.WindowStaysOnTopHint)
	splash.setMask(splash_pix.mask())

	#future placeholder for version control
	# font=QtGui.QFont("Helvetica [Cronyx]",12,weight=QtGui.QFont.Bold)
	# splash.setFont(font)
	# splash.showMessage("Established connections",QtCore.Qt.AlignBottom | QtCore.Qt.AlignVCenter);

	splash.show()
	app.processEvents()
	
	window = sf_Interactor()

	window.show()
	splash.finish(window)
	window.iren.Initialize() # Need this line to actually show the render inside Qt
	if sys.stdin.isatty() and not hasattr(sys,'ps1'):
		sys.exit(app.exec_())
	else:
		return app.exec_()

class sf_MainWindow(object):

	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.setWindowTitle("pyCM - surface fitting v%s" %__version__)
		MainWindow.resize(1280, 720)
		self.centralWidget = QtGui.QWidget(MainWindow)
		self.Boxlayout = QtGui.QHBoxLayout(self.centralWidget)
		splineBox = QtGui.QFormLayout()
		sectionBox = QtGui.QGridLayout()


		self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
		self.vtkWidget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
		self.vtkWidget.setMinimumSize(1100, 640); #leave 100 px on the size for i/o

		self.Boxlayout.addWidget(self.vtkWidget)
		self.Boxlayout.addStretch()
		MainWindow.setCentralWidget(self.centralWidget)

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
		
		self.updateButton = QtGui.QPushButton('Update')
		self.updateButton.setMinimumWidth(50)
		
		self.statLabel=QtGui.QLabel("Idle")
		self.statLabel.setFont(QtGui.QFont("Helvetica",italic=True))
		self.statLabel.setMinimumWidth(50)
		horizLine1=QtGui.QFrame()
		horizLine1.setFrameStyle(QtGui.QFrame.HLine)
		horizLine2=QtGui.QFrame()
		horizLine2.setFrameStyle(QtGui.QFrame.HLine)
		horizLine3=QtGui.QFrame()
		horizLine3.setFrameStyle(QtGui.QFrame.HLine)

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
		self.writeButton = QtGui.QPushButton('Save')
		self.writeButton.setMinimumWidth(50)

		#splineBox is the main container, the sectionBox is nested within
		splineBox.addRow(horizLine1)
		splineBox.addRow(splineLabel)
		splineBox.addRow(self.numLabel1,self.numEdit1)
		splineBox.addRow(self.numLabel2,self.numEdit2)
		splineBox.addRow(self.numLabel3,self.numEdit3)
		splineBox.addRow(self.numLabel4,self.numEdit4)
		splineBox.addRow(self.numLabel5,self.numEdit5)
		splineBox.addRow(self.updateButton)
		splineBox.addRow(self.statLabel)
		splineBox.addRow(horizLine2)
		splineBox.addRow(sectionBox)
		splineBox.addRow(horizLine3)
		splineBox.addRow(self.writeButton)

		
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


class sf_Interactor(QtGui.QMainWindow):

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
		self.iren.AddObserver("KeyPressEvent", self.Keypress)

		self.PointSize=2
		self.LineWidth=1
		self.Zaspect=1.0
		self.limits=np.empty(6)

		if (len(sys.argv)<2) or (not os.path.isfile(sys.argv[1])):
			self.filer,startdir=GetFile()
			if self.filer == None:
				#if its called from the commandline
				if sys.stdin.isatty() and not hasattr(sys,'ps1'):
					sys.exit("No file identified")
				else:
					#hopefully return to interactive python
					return
		elif len(args)==2:
			filer=args[0]
			self.outputd=args[1]
			if not os.path.exists(self.outputd): #make the directory if it doesn't exist
				os.makedirs(self.outputd)
		else:
			sys.exit("Arguments not specified correctly. Quitting.")

		#Read in reference data, calculate relevant details
		try:
			mat_contents = sio.loadmat(self.filer)
			try:
				avg=mat_contents['avg']
				pts=avg['pts'][0]
				pts=np.concatenate(pts,axis=0)
				pts=pts[~np.isnan(pts).any(axis=1)] #remove all nans
				self.RefOutline=np.concatenate(mat_contents['ali']['x_out'],axis=0)[0]
				RefMin=np.amin(self.RefOutline,axis=0)
				RefMax=np.amax(self.RefOutline,axis=0)
				self.limits=[RefMin[0],RefMax[0],RefMin[1],RefMax[1],np.amin(pts[:,-1],axis=0),np.amax(pts[:,-1],axis=0)]
				self.ui.xMin.setText('%.3f'%self.limits[0])
				self.ui.xMax.setText('%.3f'%self.limits[1])
				self.ui.yMin.setText('%.3f'%self.limits[2])
				self.ui.yMax.setText('%.3f'%self.limits[3])
			except KeyError:
				print "Couldn't read variables from file. Quitting."
				return
			
		except KeyError:
			print "Error reading reference data"
			return
		self.DisplayPointCloud(pts,RefMin,RefMax)
		self.DisplayOutline(self.RefOutline)
		self.ui.updateButton.clicked.connect(lambda: self.onUpdateSpline(pts,self.RefOutline,RefMin,RefMax))
		self.ui.sectionButton.clicked.connect(lambda: self.Cut())
		self.ui.revertButton.clicked.connect(lambda: self.RemoveCut())
		self.ui.writeButton.clicked.connect(lambda: self.WriteOutput())
		if self.filer == None:
			return

	def onUpdateSpline(self,p,ro,rmin,rmax):
		
		self.ui.statLabel.setText("Fitting . . .")
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
			self.ui.statLabel.setText("Grid too large . . .")
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
			rx=np.linspace(rmin[0],rmax[0],int((rmax[0]-rmin[0])/(gx/2)))
			ry=np.linspace(rmin[1],rmax[1],int((rmax[1]-rmin[1])/(gy/2)))
			
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

			self.ui.statLabel.setText("Rendering . . .")
			self.DisplaySplineFit(points3D,tri)

			zeval = np.empty(np.size(p[:,2]))
			for i in range(len(zeval)):
				zeval[i] = bisplev(p[i,0], p[i,1], self.tck)

			a=(zeval-p[:,2])
			RSME=(np.sum(a*a)/len(a))**0.5

			self.ui.statLabel.setText("RSME: %2.2f micron . . . Idle"%(RSME*1000))

		except ValueError as ve:
			print ve
			splineFail= True
			self.ui.statLabel.setText("Last fit failed . . . Idle")
		
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
				print j,i[j]
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
		

		self.Omapper.SetClippingPlanes(planeCollection)
		self.Pmapper.SetClippingPlanes(planeCollection)

		if hasattr(self,'Smapper'):
			self.Smapper.SetClippingPlanes(planeCollection)
		nl=np.array(self.ax3D.GetBounds())
		self.ren.RemoveActor(self.ax3D)
		self.AddAxis([pts[0],pts[1],pts[2],pts[3],nl[-2],nl[-1]],1)
		self.ui.vtkWidget.update()
	
	def RemoveCut(self):
		self.ui.xMin.setText('%.3f'%self.limits[0])
		self.ui.xMax.setText('%.3f'%self.limits[1])
		self.ui.yMin.setText('%.3f'%self.limits[2])
		self.ui.yMax.setText('%.3f'%self.limits[3])
		self.Omapper.RemoveAllClippingPlanes()
		self.Pmapper.RemoveAllClippingPlanes()
		if hasattr(self,'Smapper'):
			self.Smapper.RemoveAllClippingPlanes()
		self.ren.RemoveActor(self.ax3D)
		self.AddAxis(self.limits,1)
		self.ui.vtkWidget.update()
	
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


	def Keypress(self,obj,event):
		key = obj.GetKeyCode()

		if key =="1":
			XYView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="2":
			YZView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="3":
			XZView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)

		elif key=="z":
			self.Zaspect=self.Zaspect*2
			self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'splineActor'):
				self.splineActor.SetScale(1,1,self.Zaspect)
			self.ren.RemoveActor(self.ax3D)
			nl=np.append(self.limits[0:4],[self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect])
			self.AddAxis(nl,1/self.Zaspect)

		elif key=="x":
			self.Zaspect=self.Zaspect*0.5
			self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'splineActor'):
				self.splineActor.SetScale(1,1,self.Zaspect)
			self.ren.RemoveActor(self.ax3D)
			nl=np.append(self.limits[0:4],[self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect])
			self.AddAxis(nl,1/self.Zaspect)

		elif key=="c":
			self.Zaspect=1.0
			self.pointActor.SetScale(1,1,self.Zaspect)
			if hasattr(self,'splineActor'):
				self.splineActor.SetScale(1,1,self.Zaspect)
			self.ax3D.SetFlyModeToOuterEdges()
			self.ren.RemoveActor(self.ax3D)
			self.AddAxis(self.limits,1)

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
		QtGui.qApp.processEvents()
		if hasattr(self,'splineActor'): #then spline fitting has been done
			if hasattr(self,'outputd') is False:
				self.outputd = askdirectory(title="Choose output directory.",initialdir=currentdir)
			if self.outputd is '':
				print "No output written."
				pass
			Prefix=os.path.basename(self.filer)
			Prefix=Prefix.split('.')
			fname=os.path.join(self.outputd,Prefix[0][0:-4]+"_sur.mat")
			coefs=[np.reshape(self.tck[2],(len(self.tck[0])-self.tck[3]-1,-1))]
			number=np.array([len(self.tck[0]),len(self.tck[1])])
			order=np.array([self.tck[3], self.tck[4]])
			sio.savemat(fname,{'spline_x': {'form': 'B-', 'knots': [self.tck[0], self.tck[1]], 'coefs': coefs, 'number': number, 'order':order, 'dim': 1, 'tck': self.tck},  'x_out':self.RefOutline})
			print "Output saved to %s" %fname
		else:
			print 'Nothing to save.'

	def AddAxis(self,limits,scale):
		self.ax3D = vtk.vtkCubeAxesActor()
		self.ax3D.ZAxisTickVisibilityOn()
		self.ax3D.SetXTitle('X')
		self.ax3D.SetYTitle('Y')
		self.ax3D.SetZTitle('Z')
		self.ax3D.SetBounds(limits)
		# self.ax3D.SetLabelScaling(True,1,1,2) #ints are powers to raise each by
		self.ax3D.SetZAxisRange(limits[-2]*scale,limits[-1]*scale)
		self.ax3D.SetCamera(self.ren.GetActiveCamera())
		self.ren.AddActor(self.ax3D)
		self.ax3D.SetFlyModeToOuterEdges()

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
	currentdir=os.getcwd()

	if len(sys.argv)>2:
		RefFile=os.path.join(currentdir,sys.argv[1])
		outDir=os.path.join(currentdir,sys.argv[2])
		fit_surface(RefFile,outDir)
	elif len(sys.argv)>1:
		RefFile=os.path.join(currentdir,sys.argv[1])
		fit_surface(RefFile)
	else:
		fit_surface()


