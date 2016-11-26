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
'''
__author__ = "M.J. Roy"
__version__ = "1.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
#############################################
import sys
import os.path

global tk, tkFileDialog, askdirectory, vtk, VN, np, sio, nosio
try:
	import Tkinter as tk
	from tkFileDialog import askopenfilename
	from tkFileDialog import askdirectory
except ImportError:
	print "Error. Tkinter not found; check version of Python."

try:
	import vtk
	import vtk.util.numpy_support as VN
except ImportError:
	print "Error. VTK not found."

try:
	import numpy as np
	import numpy.matlib
except ImportError:
	print "Error. Numpy not found."

nosio=False
try:
	import scipy.io as sio
	from scipy.interpolate import griddata
	from scipy.spatial.distance import pdist, squareform
	from matplotlib import path


except ImportError:
	print "Scipy is not installed or misconfigured. Output will be in text delimited format."
	nosio=True
	

def align_average(*args, **kwargs):
	global Rotating, Panning, Zooming
	global iren, renWin, ren, defaultCameraFocalPoint, defaultCameraPosition, FloatCent
	global Zaspect, pointSize
	global FloatOutline, FloatPoints, OutlineActor2, Outline1, Outline2, qpointActor, transActor
	global outputd, CurrPnts, AlignmentComplete, RefPoints2

	root = tk.Tk()
	root.withdraw()
	
	currentdir=os.getcwd()
	outputd=None
	Perim=True
	
	if len(args)==0:
		filer = askopenfilename(title='Select the REFERENCE data file:',
			initialdir=currentdir,
			filetypes =(("MAT File", "*.mat"),("All Files","*.*")))
		startdir = os.path.dirname(filer)
		if filer == '':
			print 'No file selected; exiting.'
			return
			
		filef = askopenfilename(title='Select the FLOATING data file:',
			initialdir=startdir,
			filetypes =(("MAT File", "*.mat"),("All Files","*.*")))
		startdir = os.path.dirname(filer)
		if filef == '':
			print 'No file selected; exiting.'
			return
	#various argument collections
	elif len(args)==2:
		filer=args[0]
		filef=args[1]
	elif len(args)==3:
		filer=args[0]
		filef=args[1]
		outputd=args[2]
		if not os.path.exists(outputd): #make the directory if it doesn't exist
			os.makedirs(outputd)
	else:
		print 'Arguments not specified correctly. Quitting.'
		return
	

#Read in reference data, calculate relevant details
	try:
		mat_contents = sio.loadmat(filer)
		try:
			RefPoints=np.hstack((mat_contents['x'],mat_contents['y'],mat_contents['z']))
			RefOutline=(mat_contents['x_out'])
		except KeyError:
			print "Couldn't read variables from file. Quitting."
			return
	except KeyError:
		print "Error reading reference data"
		return

	RefCent=np.mean(RefPoints,axis=0)
	#get extents of dataset
	RefMin=np.amin(RefOutline,axis=0)
	RefMax=np.amax(RefOutline,axis=0)
	
	#define the grid spacing for averaging, has to be in 1st quadrant
	windowVerts=np.matrix([[0, 0],
	[0, 0.1*(RefMax[1]-RefMin[1])],
	[0.1*(RefMax[0]-RefMin[0]), 0.1*(RefMax[1]-RefMin[1])],
	[0.1*(RefMax[0]-RefMin[0]), 0]]);


	p=path.Path(windowVerts)
	inWindow=p.contains_points(RefPoints[:,:2]) #first 2 columns of RefPoints is x and y

	windowed=RefPoints[inWindow,:2]
	gs=squareform(pdist(windowed,'euclidean')) #does the same thing as pdist2
	gsize=np.mean(np.sort(gs)[:,1]) #sort the distances, find the closest, non self-referencing points

	grid_x, grid_y = np.meshgrid(np.linspace(RefMin[0],RefMax[0],int((RefMax[0]-RefMin[0])/gsize)),
	    np.linspace(RefMin[1],RefMax[1],int((RefMax[1]-RefMin[1])/gsize)),indexing='xy')
	points=RefPoints[:,:2]

	grid_RefVal=griddata(points,RefPoints[:,-1], (grid_x, grid_y), method='linear')


#Read in reference data, calculate relevant details
	try:
		mat_contents = sio.loadmat(filef)
		try:
			FloatPoints=np.hstack((mat_contents['x'],mat_contents['y'],mat_contents['z']))
			FloatOutline=(mat_contents['x_out'])
		except KeyError:
			print "Couldn't read variables from file. Quitting."
			return
	except KeyError:
		print "Error reading reference data"
		return

	FloatCent=np.mean(FloatPoints,axis=0)
	#get extents of dataset

	pointSize=2

	# Create instances
	pointCloud = VtkPointCloud()
	qpointCloud = VtkPointCloud()

	#calculate the offset between the floating and reference dataset.

	oa=np.array([0,0,40]); #offset of floating data from reference
	
	offset=np.zeros(shape=(2,3))
	
	# offset[0]=np.add(-RefCent);
	offset[0]=-RefCent
	for k in RefPoints:
		pointCloud.addPoint(np.add(k,offset[0]))
	OutlineActor1,Outline1=initializeOutline(np.add(RefOutline,offset[0]),(0.95,0.3961,0.1333))

	# Add points to to the VtkHighlightPointCloud Instance
	offset[1]=np.add(-FloatCent,oa);
	FloatOutline=np.add(FloatOutline,offset[1])
	for k in FloatPoints:
		qpointCloud.addPoint(np.add(k,offset[1]))
	OutlineActor2,Outline2=initializeOutline(FloatOutline,(1,0.804,0.204))

	pointActor = pointCloud.vtkActor
	pointActor.GetProperty().SetColor(0.95,0.3961,0.1333)
	qpointActor = qpointCloud.vtkActor
	qpointActor.GetProperty().SetColor(1,0.804,0.204)

	# Add axes and fix automatic scaling issue
	
	
	axes3D = vtk.vtkCubeAxesActor() #try 2D
	axes3D.ZAxisLabelVisibilityOff()
	axes3D.ZAxisTickVisibilityOff()
	axes3D.SetXTitle('X')
	axes3D.SetYTitle('Y')
	



	# Create the Renderer and assign actors to it. A renderer is like a
	# viewport. It is part or all of a window on the screen and it is
	# responsible for drawing the actors it has.  We also set the
	# background color here.
	ren = vtk.vtkRenderer()
	
	ren.AddActor(pointActor)
	ren.AddActor(qpointActor)

	axes3D.SetBounds(RefMin[0]-RefCent[0],RefMax[0]-RefCent[0],
	RefMin[1]-RefCent[1],RefMax[1]-RefCent[1],
	RefMin[2]-RefCent[2],RefMax[2]-RefCent[2])



	ren.AddActor(OutlineActor1)
	ren.AddActor(OutlineActor2)
	ren.SetBackground(0.1, 0.2, 0.4)

	# Finally we create the render window which will show up on the screen
	# We put our renderer into the render window using AddRenderer. We
	# also set the size to be 720p
	renWin = vtk.vtkRenderWindow()
	renWin.AddRenderer(ren)
	renWin.SetSize(1280, 720)

	# Define custom interaction.
	iren = vtk.vtkRenderWindowInteractor()
	iren.SetInteractorStyle(None)
	iren.SetRenderWindow(renWin)
	AlignmentComplete=False
	AveragingComplete=False


	# Add the observers to watch for particular events. These invoke
	# Python functions.
	Rotating = 0
	Panning = 0
	Zooming = 0
	Zaspect=1

	def Keypress(obj, event):
		global iren, ren, renWin, Zaspect, pointSize, FloatOutline, FloatPoints, FloatCent, OutlineActor2, Outline1, Outline2, transActor, transformPointCloud, AvgActor, AvgPointCloud, AlignmentComplete, AveragingComplete, CurrPnts, inOutline, AvgPoints_grid, AvgPoints
		key = obj.GetKeySym()
		if key == "e":
			if AlignmentComplete and AveragingComplete:
				WriteOutput()
				sys.exit("align_average completed.")
			else:
				sys.exit("Alignment/averaging not performed, not writing any output.")
			return
		elif key =="h":
			#handle outline
			DeleteActor(OutlineActor2)
			FloatOutline=np.add(FloatOutline,-offset[1])
			FloatCent=np.mean(FloatOutline,axis=0)
			FloatOutline=np.add(FloatOutline,-FloatCent)
			FloatOutline[:,0]=-FloatOutline[:,0]
			FloatOutline=np.add(FloatOutline,FloatCent+offset[1])
			OutlineActor2,Outline2=initializeOutline(FloatOutline,(1,0.804,0.204))
			ren.AddActor(OutlineActor2)
			renWin.Render()
			#different function for the point cloud where points are moved, and not recreated
			flipSide('x',qpointCloud,qpointActor,offset[1,-1])

		elif key =="k":
			#handle outline
			DeleteActor(OutlineActor2)
			FloatOutline=np.add(FloatOutline,-offset[1])
			FloatCent=np.mean(FloatOutline,axis=0)
			FloatOutline=np.add(FloatOutline,-FloatCent)
			FloatOutline[:,1]=-FloatOutline[:,1]
			FloatOutline=np.add(FloatOutline,FloatCent+offset[1])
			OutlineActor2,Outline2=initializeOutline(FloatOutline,(1,0.804,0.204))
			ren.AddActor(OutlineActor2)
			renWin.Render()
			#different function for the point cloud where points are moved, and not recreated
			flipSide('y',qpointCloud,qpointActor,offset[1,-1])
		elif key=="a":
			if not AlignmentComplete:
				print "Alignment started . . ."
				print "Removing scaling . . ."
				oldZ=Zaspect;
				#clear any scaling
				updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),1,offset[1][-1])
				Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),1,offset[0][-1]) #because each call will update Zaspect

				#Align the floating points with the reference points using vtk's native icp filter
				icp=vtk.vtkIterativeClosestPointTransform()

				icp.SetSource(qpointCloud.geometry)
				icp.SetTarget(pointCloud.geometry)
				# icp.GetLandmarkTransform().SetModeToRigidBody()
				icp.SetMaximumNumberOfIterations(200)
				icp.StartByMatchingCentroidsOn()
				icp.Modified()
				icp.Update()
				icp.Inverse()

				transM=np.zeros(shape=(4,4))
				for i in range(4):
					for j in range(4):
						transM[i,j]=icp.GetMatrix().GetElement(i, j)
				transM=np.linalg.inv(transM)
				print "Transformation matrix for the floating data onto the reference:"
				print(transM)
				
				icpTransformFilter = vtk.vtkTransformPolyDataFilter()
				if vtk.VTK_MAJOR_VERSION <= 5:
					icpTransformFilter.SetInput(qpointCloud.geometry)
				else:
					icpTransformFilter.SetInputData(qpointCloud.geometry)

				icpTransformFilter.SetTransform(icp)
				icpTransformFilter.Update()
				transformedSource = icpTransformFilter.GetOutput()
				transformPointCloud=VtkPointCloud()

				#Apply transformation matrix to the floating point cloud
				CurrPnts=np.empty(FloatPoints.shape)
				for k in range(int(qpointCloud.points.GetNumberOfPoints())):
					dummy=np.dot(np.asarray(qpointCloud.points.GetPoint(k)),transM[0:3,0:3])+transM[0:3,-1] 
					transformPointCloud.addPoint(dummy)
					CurrPnts[k,:]=dummy

				# Or apply points directly from the filter
				# CurrPnts=np.empty(FloatPoints.shape)
				# for k in range(int(transformedSource.GetNumberOfPoints())):
					# dummy=np.asarray(transformedSource.GetPoint(k))
					# transformPointCloud.addPoint(dummy)
					# CurrPnts[k,:]=dummy

				transActor = transformPointCloud.vtkActor
				transActor.GetProperty().SetColor(1,0.9098,0.6863)


				ren.AddActor(transActor)
				renWin.Render()

				print "Re-applying scaling . . ."
				updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),oldZ,offset[0][-1])

				updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),oldZ,offset[1][-1])
			
				Zaspect=updateZaspect(transformPointCloud,transActor,transformPointCloud.points.GetNumberOfPoints(),oldZ,offset[0][-1])
				
				
				print "Alignment complete."
				AlignmentComplete=True
			else:
				print "Restart routine to re-try alignment."

		elif key=="A":
			if not AlignmentComplete:
				print "Alignment started using perimeter . . ."

				#Align the floating points with the reference points using vtk's native icp filter
				icp=vtk.vtkIterativeClosestPointTransform()

				icp.SetSource(Outline2)
				icp.SetTarget(Outline1)
				icp.GetLandmarkTransform().SetModeToRigidBody()
				icp.SetMaximumNumberOfIterations(20)
				icp.StartByMatchingCentroidsOn()
				icp.Modified()
				icp.Update()
				icp.Inverse()
				
				transM=np.zeros(shape=(4,4))
				for i in range(4):
					for j in range(4):
						transM[i,j]=icp.GetMatrix().GetElement(i, j)
				transM=np.linalg.inv(transM)
				print "Transformation matrix for the floating data onto the reference:"
				print(transM)
				
				icpTransformFilter = vtk.vtkTransformPolyDataFilter()
				if vtk.VTK_MAJOR_VERSION <= 5:
					icpTransformFilter.SetInput(Outline2)
				else:
					icpTransformFilter.SetInputData(Outline2)

				icpTransformFilter.SetTransform(icp)
				icpTransformFilter.Update()

				transformedSource = icpTransformFilter.GetOutput()
				transformPointCloud=VtkPointCloud()

				#Apply transformation matrix to the floating point cloud
				CurrPnts=np.empty(FloatPoints.shape)
				for k in range(int(qpointCloud.points.GetNumberOfPoints())):
					dummy=np.dot(np.asarray(qpointCloud.points.GetPoint(k)),transM[0:3,0:3])+transM[0:3,-1] 
					transformPointCloud.addPoint(dummy)
					CurrPnts[k,:]=dummy

				transActor = transformPointCloud.vtkActor
				transActor.GetProperty().SetColor(1,0.9098,0.6863)


				ren.AddActor(transActor)
				renWin.Render()
				
				print "Alignment using perimeter complete."
				AlignmentComplete=True
			else:
				print "Restart routine to re-try alignment."
		elif key =="v":
			if AlignmentComplete:
				print "Averaging started . . ."
				#clear any scaling
				print "Removing scaling . . ."
				oldZ=Zaspect;
				updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),1,offset[1][-1])
				updateZaspect(transformPointCloud,transActor,transformPointCloud.points.GetNumberOfPoints(),1,offset[0][-1])
				Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),1,offset[0][-1]) #because each call will update Zaspect
				#grid CurrPoints using grid_x and grid_y, clear scaling as CurrPnts is coming from a filter
				points=np.add(CurrPnts,-offset[0])

				grid_AlignVal=griddata(points[:,:2],points[:,-1], (grid_x, grid_y), method='linear')
				#average them
				AvgPoints_grid=(grid_AlignVal+grid_RefVal)/2
				#find points in the outline
				test_outline=RefOutline[:,:2]
				#do column-wise raveling, replicate Matlab handling
				p=path.Path(test_outline)
				inTest=np.hstack((np.ravel(grid_x.T)[np.newaxis].T,np.ravel(grid_y.T)[np.newaxis].T))

				inOutline=p.contains_points(inTest) 
				#build a 3D list of points containing those points 'inside' the outline
				AvgPoints=np.hstack((inTest[inOutline,:],
					np.ravel(AvgPoints_grid.T)[np.newaxis].T[inOutline]))

				#create new actor, etc
				AvgPointCloud=VtkPointCloud()
				for k in AvgPoints:
					AvgPointCloud.addPoint(np.add(k,offset[0]))
				AvgActor = AvgPointCloud.vtkActor

				AvgActor.GetProperty().SetColor(0.2784,0.6745,0.6941)

				ren.AddActor(AvgActor)
				renWin.Render()
				print "Re-applying scaling . . ."
				updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),oldZ,offset[0][-1])

				updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),oldZ,offset[1][-1])
			
				updateZaspect(transformPointCloud,transActor,transformPointCloud.points.GetNumberOfPoints(),oldZ,offset[0][-1])
				
				Zaspect=updateZaspect(AvgPointCloud,AvgActor,AvgPointCloud.points.GetNumberOfPoints(),oldZ,offset[0][-1])
				
				print "Averaging complete."
				AveragingComplete=True
			else:
				print "Align before averaging . . ."
		elif key =="d":
			FlipVisible(qpointActor)
			FlipVisible(OutlineActor2)
			FlipVisible(pointActor)
			FlipVisible(OutlineActor1)
		elif key =="1":
			XYView(ren, ren.GetActiveCamera())
		elif key =="2":
			YZView(ren, ren.GetActiveCamera())
		elif key =="3":
			XZView(ren, ren.GetActiveCamera())
		#change z aspect ratio
		elif key=="z":
			Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),Zaspect*2,offset[0][-1])

			updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),Zaspect*2,offset[1][-1])
			try:
				updateZaspect(transformPointCloud,transActor,transformPointCloud.points.GetNumberOfPoints(),Zaspect*2,offset[0][-1])
			except: #the alignment might not have taken place.
				pass
			try:
				updateZaspect(AvgPointCloud,AvgActor,AvgPointCloud.points.GetNumberOfPoints(),Zaspect*2,offset[0][-1])
			except: #averaging might not have taken place.
				pass


		elif key=="x":
			Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),Zaspect*0.5,offset[0][-1])
			updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),Zaspect*0.5,offset[1][-1])
			try:
				updateZaspect(transformPointCloud,transActor,transformPointCloud.points.GetNumberOfPoints(),Zaspect*0.5,offset[0][-1])
			except: #the alignment might not have taken place.
				pass
			try:
				updateZaspect(AvgPointCloud,AvgActor,AvgPointCloud.points.GetNumberOfPoints(),Zaspect*0.5,offset[0][-1])
			except: #averaging might not have taken place.
				pass

		elif key=="c":
			#make sure the old scaling factor is used, not the updated one.
			updateZaspect(qpointCloud,qpointActor,qpointCloud.points.GetNumberOfPoints(),1,offset[1][-1])
			try:
				updateZaspect(transformPointCloud,transActor,transformPointCloud.points.GetNumberOfPoints(),1,offset[0][-1])
			except: #the alignment might not have taken place.
				pass
			try:
				updateZaspect(AvgPointCloud,AvgActor,AvgPointCloud.points.GetNumberOfPoints(),1,offset[0][-1])
			except: #averaging might not have taken place.
				pass
			Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),1,offset[0][-1]) #because each call will update Zaspect, so do this one last
			renWin.Render()

		elif key =="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(renWin)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("Avg.png")
			writer.Write()
			print 'Screen output saved to %s' %os.path.join(currentdir,'Avg.png')
		elif key =="r":
			FlipVisible(axes3D)



		elif key =="f": #flip color scheme for printing
			if ren.GetBackground()==(0.1, 0.2, 0.4):
				axes3D.GetTitleTextProperty(0).SetColor(0,0,0)
				axes3D.GetLabelTextProperty(0).SetColor(0,0,0)
				axes3D.GetXAxesLinesProperty().SetColor(0,0,0)
				axes3D.SetXTitle('x') #there's a vtk bug here . . .
				
				axes3D.GetTitleTextProperty(1).SetColor(0,0,0)
				axes3D.GetLabelTextProperty(1).SetColor(0,0,0)
				axes3D.GetYAxesLinesProperty().SetColor(0,0,0)

				axes3D.GetTitleTextProperty(2).SetColor(0,0,0)
				axes3D.GetLabelTextProperty(2).SetColor(0,0,0)

				ren.SetBackground(1, 1, 1)
			else:
				axes3D.GetTitleTextProperty(0).SetColor(1,1,1)
				axes3D.GetLabelTextProperty(0).SetColor(1,1,1)
				axes3D.GetXAxesLinesProperty().SetColor(1,1,1)
				axes3D.SetXTitle('X')
				
				axes3D.GetTitleTextProperty(1).SetColor(1,1,1)
				axes3D.GetLabelTextProperty(1).SetColor(1,1,1)
				axes3D.GetYAxesLinesProperty().SetColor(1,1,1)
				axes3D.SetYTitle('Y')
				
				axes3D.GetTitleTextProperty(2).SetColor(1,1,1)
				axes3D.GetLabelTextProperty(2).SetColor(1,1,1)

				ren.SetBackground(0.1, 0.2, 0.4)

			renWin.Render()



	def WriteOutput():
	
		global outputd, iren, renWin
		
		# dtype1={'names':['f{}'.format(i) for i in range(ReadPoints.shape[1])],
			# 'formats':ReadPoints.shape[1] * [ReadPoints.dtype]} #make each row an entry in a set
		# dtype2={'names':['f{}'.format(i) for i in range(CurrPnts.shape[1])],
			# 'formats':CurrPnts.shape[1] * [CurrPoints.dtype]} #make each row an entry in a set

		##################################
		currentdir=os.getcwd()
		if outputd is None:
			outputd = askdirectory(title="Choose output directory.",
							initialdir=currentdir)
		if outputd is '':
			print "No output written."
			close_window(iren)
			del renWin, iren
			return
		else:
			Prefix=os.path.basename(filer)
			Prefix=Prefix.split('.')
			print"Now writing output . . ."

			sio.savemat(os.path.join(outputd,Prefix[0]+"_avg.mat"),
							{'ali': {'xf':CurrPnts[:,0]-offset[0][0],
							'yf':CurrPnts[:,1]-offset[0][1],
							'zf':CurrPnts[:,2]-offset[0][2],
							'xr':RefPoints[:,0],'yr':RefPoints[:,1],'zr':RefPoints[:,2],
							'x_out':RefOutline},
							'avg' : { 'in' : inOutline, 'xi':grid_x , 'yi':grid_y  , 'zi': AvgPoints_grid, 'pts': AvgPoints}})
			# print np.transpose(np.asmatrix(RefPoints2[:,0])), np.transpose(np.asmatrix(RefPoints2[:,1])), np.transpose(np.asmatrix(RefPoints2[:,2]))
			print "Output saved to %s" %outputd
			close_window(iren)
			del renWin, iren

	iren.AddObserver("LeftButtonPressEvent", ButtonEvent)
	iren.AddObserver("LeftButtonReleaseEvent", ButtonEvent)
	iren.AddObserver("MiddleButtonPressEvent", ButtonEvent)
	iren.AddObserver("MiddleButtonReleaseEvent", ButtonEvent)
	iren.AddObserver("RightButtonPressEvent", ButtonEvent)
	iren.AddObserver("RightButtonReleaseEvent", ButtonEvent)
	iren.AddObserver("MouseMoveEvent", MouseMove)
	iren.AddObserver("KeyPressEvent", Keypress)

	iren.Initialize()
	renWin.Render()
	#get default camera position as inputs for 'views'
	defaultCameraPosition=ren.GetActiveCamera().GetPosition()
	defaultCameraFocalPoint=ren.GetActiveCamera().GetFocalPoint()
	
	axes3D.SetCamera(ren.GetActiveCamera())
	ren.AddActor(axes3D)

	renWin.SetWindowName("UoM Contour Method - Alignment and Averaging v%s" %__version__)
	iren.Start()

	#if the user otherwise closes the interactor window
	if 'iren' in globals():
		print "Interactor closed, no output written."
		close_window(iren)
		del renWin, iren

def flipSide(flipDirection,instance,actor,Offset):
	NumPoints=instance.points.GetNumberOfPoints()
	for i in range(NumPoints):
		xyz=instance.points.GetPoint(i)
		lp=(xyz[0], xyz[1], xyz[2]-Offset)
		if flipDirection == "x": #then mirror in XZ plane
			mlp=((-lp[0], lp[1], lp[2]+Offset))
		else: #mirror in ZY
			mlp=(lp[0], -lp[1], lp[2]+Offset)
		instance.points.SetPoint(i,mlp[0],mlp[1],mlp[2])
	instance.geometry.Modified()
	actor.Modified()
	renWin.Render()

def initializeOutline(perimeter,color):
	points = vtk.vtkPoints()
	for j in perimeter:
		points.InsertNextPoint(j)

	polygon=vtk.vtkPolygon()
	polygon.GetPointIds().SetNumberOfIds(len(perimeter))
	for i in range(len(perimeter)):
		polygon.GetPointIds().SetId(i,i)
	
	polygons=vtk.vtkCellArray()
	polygons.InsertNextCell(polygon)
	polyData=vtk.vtkPolyData()
	polyData.SetPoints(points)
	polyData.SetLines(polygons) #change to SetPoly for filled polygon
	polygonMapper=vtk.vtkPolyDataMapper()
	if vtk.VTK_MAJOR_VERSION<= 5:
		polygonMapper.SetInput(polyData)
	else:
		polygonMapper.SetInputData(polyData)
	perimActor=vtk.vtkActor()
	perimActor.SetMapper(polygonMapper)
	perimActor.GetProperty().SetColor(color)
	return perimActor,polyData

##Common functions
# Routines that translate the events into camera motions.
# Handle the mouse button events.
def ButtonEvent(obj, event):
	global Rotating, Panning, Zooming
	global iren, renWin, ren
	if event == "LeftButtonPressEvent":
		Rotating = 1
	elif event == "LeftButtonReleaseEvent":
		Rotating = 0
	elif event == "MiddleButtonPressEvent":
		Panning = 1
	elif event == "MiddleButtonReleaseEvent":
		Panning = 0
	elif event == "RightButtonPressEvent":
		Zooming = 1
	elif event == "RightButtonReleaseEvent":
		Zooming = 0
# This one is associated with the left mouse button. It translates x
# and y relative motions into camera azimuth and elevation commands.
def Rotate(renderer, camera, x, y, lastX, lastY, centerX, centerY):
	camera.Azimuth(lastX-x)
	camera.Elevation(lastY-y)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCameraClippingRange() #needed to remove clipping artefact
	renWin.Render()

def XYView(renderer, camera):
	camera.SetPosition(0,0,defaultCameraPosition[2]+0)
	camera.SetFocalPoint(defaultCameraFocalPoint)
	camera.SetViewUp(0,1,0)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renWin.Render()

def YZView(renderer, camera):
	camera.SetPosition(defaultCameraPosition[2]+0,0,0)
	camera.SetFocalPoint(defaultCameraFocalPoint)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renWin.Render()

def XZView(renderer, camera):
	camera.SetPosition(0,defaultCameraPosition[2]+0,0)
	camera.SetFocalPoint(defaultCameraFocalPoint)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renWin.Render()

# Pan translates x-y motion into translation of the focal point and
# position.
def Pan(renderer, camera, x, y, lastX, lastY, centerX, centerY):
	FPoint = camera.GetFocalPoint()
	FPoint0 = FPoint[0]
	FPoint1 = FPoint[1]
	FPoint2 = FPoint[2]

	PPoint = camera.GetPosition()
	PPoint0 = PPoint[0]
	PPoint1 = PPoint[1]
	PPoint2 = PPoint[2]

	renderer.SetWorldPoint(FPoint0, FPoint1, FPoint2, 1.0)
	renderer.WorldToDisplay()
	DPoint = renderer.GetDisplayPoint()
	focalDepth = DPoint[2]

	APoint0 = centerX+(x-lastX)
	APoint1 = centerY+(y-lastY)

	renderer.SetDisplayPoint(APoint0, APoint1, focalDepth)
	renderer.DisplayToWorld()
	RPoint = renderer.GetWorldPoint()
	RPoint0 = RPoint[0]
	RPoint1 = RPoint[1]
	RPoint2 = RPoint[2]
	RPoint3 = RPoint[3]

	if RPoint3 != 0.0:
		RPoint0 = RPoint0/RPoint3
		RPoint1 = RPoint1/RPoint3
		RPoint2 = RPoint2/RPoint3

	camera.SetFocalPoint( (FPoint0-RPoint0)/2.0 + FPoint0,
						  (FPoint1-RPoint1)/2.0 + FPoint1,
						  (FPoint2-RPoint2)/2.0 + FPoint2)
	camera.SetPosition( (FPoint0-RPoint0)/2.0 + PPoint0,
						(FPoint1-RPoint1)/2.0 + PPoint1,
						(FPoint2-RPoint2)/2.0 + PPoint2)
	camera.ParallelProjectionOn()
	renWin.Render()

# Dolly converts y-motion into a camera dolly commands.
def Dolly(renderer, camera, x, y, lastX, lastY, centerX, centerY):
	dollyFactor = pow(1.02,(0.5*(y-lastY)))
	if camera.GetParallelProjection():
		parallelScale = camera.GetParallelScale()*dollyFactor
		camera.SetParallelScale(parallelScale)
	else:
		camera.Dolly(dollyFactor)
		renderer.ResetCameraClippingRange()
	camera.ParallelProjectionOn()
	renWin.Render()

def DeleteActor(actor):
	global ren, renWin, iren
	ren.RemoveActor(actor)
	renWin.Render()
	

def FlipVisible(actor):
	global ren, renWin, iren
	if actor.GetVisibility():
		actor.VisibilityOff()
	else:
		actor.VisibilityOn()
	renWin.Render()
	
# General high-level logic
def MouseMove(obj, event):
	global Rotating, Panning, Zooming
	global iren, renWin, ren
	lastXYpos = iren.GetLastEventPosition()
	lastX = lastXYpos[0]
	lastY = lastXYpos[1]
	
	xypos = iren.GetEventPosition()
	x = xypos[0]
	y = xypos[1]
	
	center = renWin.GetSize()
	centerX = center[0]/2.0
	centerY = center[1]/2.0

	if Rotating:
		Rotate(ren, ren.GetActiveCamera(), x, y, lastX, lastY,
			   centerX, centerY)
	elif Panning:
		Pan(ren, ren.GetActiveCamera(), x, y, lastX, lastY, centerX,
			centerY)
	elif Zooming:
		Dolly(ren, ren.GetActiveCamera(), x, y, lastX, lastY,
			  centerX, centerY)

def updateZaspect(instance,actor,NumPoints,NewZaspect,Offset):
	if NewZaspect==1:
		scale=1.0/Zaspect
	else:
		scale=NewZaspect/Zaspect
	for i in range(NumPoints):
		xyz=instance.points.GetPoint(i)

		#because tuples are immutable
		lp=(xyz[0], xyz[1], xyz[2]-Offset)
		
		mlp=(lp[0], lp[1], lp[2]*scale)
		instance.points.SetPoint(i,mlp[0],mlp[1],mlp[2]+Offset)

	instance.geometry.Modified()
	actor.Modified()
	renWin.Render()
	
	return NewZaspect


#Create a class to generate a point cloud
class VtkPointCloud:

	def __init__(self, zMin=0, zMax=1, maxNumPoints=1e3):
		self.maxNumPoints = maxNumPoints
		self.points = vtk.vtkPoints()
		self.vertices = vtk.vtkCellArray()
		self.geometry=vtk.vtkPolyData()
		self.geometry.SetPoints(self.points)
		self.geometry.SetVerts(self.vertices)
		mapper = vtk.vtkPolyDataMapper()
		if vtk.VTK_MAJOR_VERSION<= 5:
			mapper.SetInput(self.geometry)
		else:
			mapper.SetInputData(self.geometry)
		mapper.SetColorModeToDefault()
		mapper.SetScalarRange(zMin, zMax)
		mapper.SetScalarVisibility(1)
		self.vtkActor = vtk.vtkActor()
		self.vtkActor.SetMapper(mapper)
		self.vtkActor.GetProperty().SetPointSize(pointSize)
		

	def addPoint(self, point):
		pointId = self.points.InsertNextPoint(point[:])
		self.vertices.InsertNextCell(1)
		self.vertices.InsertCellPoint(pointId)
		self.vertices.Modified()
		self.points.Modified()
		
	def clearPoints(self):
		self.points = vtk.vtkPoints()
		self.vertices = vtk.vtkCellArray()
		self.geometry.SetPoints(self.points)
		self.geometry.SetVerts(self.vertices)

def close_window(iren):
	render_window=iren.GetRenderWindow()
	render_window.Finalize()
	iren.TerminateApp()

if __name__ == '__main__':
	currentdir=os.getcwd()
	# checkImport()
	if len(sys.argv)>3:
		RefFile=os.path.join(currentdir,sys.argv[1])
		FloatFile=os.path.join(currentdir,sys.argv[2])
		outDir=os.path.join(currentdir,sys.argv[3])
		align_average(RefFile,FloatFile,outDir)
	elif len(sys.argv)>2:
		RefFile=os.path.join(currentdir,sys.argv[1])
		FloatFile=os.path.join(currentdir,sys.argv[2])
		align_average(RefFile,FloatFile)
	else:
		align_average()


