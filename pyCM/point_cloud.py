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
p - enter/exit picking mode, LMB is used to generate a selection window. Exiting 
	picking mode will highlight selected points.
u - update. Will remove highlighted points from the point cloud.
Shift-u - soft undo, will deselect any selected points
z - increase z-aspect ratio
x - decrease z-aspect ratio
c - return to default z-aspect
Shift-z - increase size of points
Shift-x - decrease size of points
Shift-c - return to default point size
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
r - remove compass
e - Update output file, write output and exit
-------------------------------------------------------------------------------
ver 1.2 16-11-06
1.1 - Fixed array orientation, clipping issue, compass scaling and sped up writing output
      Added ReadMask
1.2 - Fixed window handling, now exits cleanly
'''
__author__ = "M.J. Roy"
__version__ = "1.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
#############################################
import sys
import os.path

# def checkImport():
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
except ImportError:
	print "Error. Numpy not found."

nosio=False
try:
	import scipy.io as sio
except ImportError:
	print "Scipy is not installed or misconfigured. Output will be in text delimited format."
	nosio=True
	

def MaskDef(*args, **kwargs):
	global Rotating, Panning, Zooming, Selection, Selecting, rectActor
	global iren, renWin, ren, defaultCameraFocalPoint, defaultCameraPosition, centroid
	global xyOffClickPosition, xyOnClickPosition
	global Zaspect, pointSize, pointsToRemove, CurrPnts
	global outputd

	root = tk.Tk()
	root.withdraw()
	
	currentdir=os.getcwd()
	outputd=None
	Perim=True
	
	if len(args)==0:
		filep = askopenfilename(title='Select the perimeter file. (Optional)',
								initialdir=currentdir)
		startdir = os.path.dirname(filep)
		if filep == '':
			print 'No perimeter selected.'
			Perim=False
			startdir=os.getcwd()
			pass
			
		filec = askopenfilename(title='Select the point cloud file. (Mandatory)',
								initialdir=startdir)
		startdir = os.path.dirname(filec)
		if filec == '':
			print 'No file selected; exiting.'
			return
	#various argument collections
	elif len(args)==1:
		Perim=False
		filec=args[0]
	elif len(args)==2:
		filep=args[0]
		filec=args[1]
	elif len(args)==3:
		filep=args[0]
		filec=args[1]
		outputd=args[2]
		if not os.path.exists(outputd): #make the directory if it doesn't exist
			os.makedirs(outputd)
	else:
		print 'Arguments not specified correctly. Quitting.'
		return
	
	if Perim:
		ReadPerimeter=np.genfromtxt(filep)
		ReadPerimeter = np.append(ReadPerimeter,[ReadPerimeter[0,:]],axis=0) #for closed contour
	#Exception for *.dat files
	Extension=os.path.splitext(filec)[1]
	if Extension == '.dat':
		ReadPoints=np.genfromtxt(filec,skiprows=1)
	elif Extension == '.mat':
		if nosio:
			print 'Error. Requires Scipy to read file.'
			return
		else:
			mat_contents = sio.loadmat(filec)
			try:
				ReadPoints=np.hstack((mat_contents['x'],mat_contents['y'],mat_contents['z']))
				#ReadPoints=np.transpose(ReadPoints)
				ReadPerimeter=(mat_contents['x_out'])
				Perim=True
			except KeyError:
				print "Couldn't read variables from file. Quitting."
				return
	else:
		ReadPoints=np.genfromtxt(filec)
	CurrPnts=ReadPoints #Periodically updated
	OutputMask=[] #Updated at completion
	pointsToRemove=[] #Periodically updated
	centroid=np.mean(ReadPoints,axis=0)
	#get extents of dataset
	min=np.amin(ReadPoints,axis=0)
	max=np.amax(ReadPoints,axis=0)
	pointSize=2


			

	# Create instances
	pointCloud = VtkPointCloud()
	qPointCloud = VtkHighlightPointCloud()


	# Add points to to the VtkPointCloud Instance
	for k in ReadPoints:
		pointCloud.addPoint(np.add(k,-centroid))
		
	NumPoints=len(ReadPoints)

	#test/create perimeter
	if Perim:
		perimActor=initializePerimeter(ReadPerimeter,centroid)
	else: perimActor=None

	pointActor = pointCloud.vtkActor
	qPointActor = qPointCloud.vtkActor

	# Add axes and fix automatic scaling issue
	axes = vtk.vtkAxesActor()

	axes.GetYAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
	axes.GetXAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
	axes.GetZAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
	axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(7*pointSize)
	axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(7*pointSize)
	axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(7*pointSize)
	axScale=np.amax(np.fabs(ReadPoints))*0.1 #scale to dataset
	axes.SetTotalLength(axScale,axScale,axScale)



	# Create the Renderer and assign actors to it. A renderer is like a
	# viewport. It is part or all of a window on the screen and it is
	# responsible for drawing the actors it has.  We also set the
	# background color here.
	ren = vtk.vtkRenderer()
	
	ren.AddActor(pointActor)
	ren.AddActor(qPointActor)
	ren.AddActor(axes)

	if Perim:
		ren.AddActor(perimActor)
	ren.SetBackground(0.1, 0.2, 0.4)

	# Finally we create the render window which will show up on the screen
	# We put our renderer into the render window using AddRenderer. We
	# also set the size to be 800 pixels by 640. Because it's 1995.
	renWin = vtk.vtkRenderWindow()
	renWin.AddRenderer(ren)
	renWin.SetSize(800, 640)
	# renWin.SetWindowName("UoM Contour Method - Point Editor")

	# Define custom interaction.
	iren = vtk.vtkRenderWindowInteractor()
	iren.SetInteractorStyle(None)
	iren.SetRenderWindow(renWin)



	# Add the observers to watch for particular events. These invoke
	# Python functions.
	Rotating = 0
	Panning = 0
	Zooming = 0
	Selecting = 0
	Selection = False
	rectActor=None
	Zaspect=1

	def MySelection(xmin,ymin,xmax,ymax):
		global pointsToRemove
		hsel = vtk.vtkHardwareSelector()
		hsel.SetRenderer(ren)
		hsel.SetArea(xmin,ymin,xmax,ymax)
		result=hsel.Select()
		numNodes = result.GetNumberOfNodes()
		if (numNodes < 1):
			# print("No visible cells") #debug
			pass
		else:
			sel_node = result.GetNode(0)
			q=sel_node.GetSelectionList()
			NowPicked=VN.vtk_to_numpy(sel_node.GetSelectionList()).tolist()
			pointsToRemove=np.append(pointsToRemove,NowPicked)
			pointsToRemove=np.unique(pointsToRemove)
			
			#for view only
			for i in NowPicked:
				q = pointCloud.points.GetPoint(i)
				qPointCloud.addPoint(q)
		renWin.Render()

	def Keypress(obj, event):
		global Selection, iren, ren, renWin, rectActor, Zaspect, pointSize, pointsToRemove, CurrPnts
		key = obj.GetKeySym()
		if key == "e":
			WriteOutput()
		elif key =="1":
			XYView(ren, ren.GetActiveCamera())
		elif key =="2":
			YZView(ren, ren.GetActiveCamera())
		elif key =="3":
			XZView(ren, ren.GetActiveCamera())
		#change z aspect ratio
		elif key=="z":
			Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),Zaspect*2,0)
			updateZaspect(qPointCloud,qPointActor,qPointCloud.points.GetNumberOfPoints(),Zaspect*2,0)
		elif key=="Z":
			pointSize=updatePointSize(pointActor,pointSize+1)
		elif key=="x":
			Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),Zaspect*0.5,0)
			updateZaspect(qPointCloud,qPointActor,qPointCloud.points.GetNumberOfPoints(),Zaspect*0.5,0)
		elif key=="X":
			pointSize=updatePointSize(pointActor,pointSize-1)
		elif key=="c":
			#clear out all scaling.
			oldZaspect=Zaspect
			Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),1,0) #because each call will update Zaspect
			#make sure the old scaling factor is used, not the updated one.
			updateZaspect(qPointCloud,qPointActor,qPointCloud.points.GetNumberOfPoints(),Zaspect*(1/oldZaspect),0)
			Zaspect=1
		elif key=="C":
			pointSize=updatePointSize(pointActor,2)
		elif key =="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(renWin)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("PointCloud.png")
			writer.Write()
			print 'Screen output saved to %s' %os.path.join(currentdir,'PointCloud.png')
		elif key =="r":
			if axes.GetVisibility():
				axes.VisibilityOff()
			else:
				axes.VisibilityOn()
			renWin.Render()
		elif key =="p":
			if Selection==True:
				DeleteActor(rectActor)
				rectActor=None
				Selection=False
				#sort the on/off selection points
				if xyOnClickPosition[0]>=xyOffClickPosition[0]:
					xmax=xyOnClickPosition[0]
					xmin=xyOffClickPosition[0]
				else:
					xmax=xyOffClickPosition[0]
					xmin=xyOnClickPosition[0]
				if xyOnClickPosition[1]>=xyOffClickPosition[1]:
					ymax=xyOnClickPosition[1]
					ymin=xyOffClickPosition[1]
				else:
					ymax=xyOffClickPosition[1]
					ymin=xyOnClickPosition[1]
				MySelection(xmin,ymin,xmax,ymax)
			else:
				Selection=True

		#Update by removing all points from CurrPnts contained in pointsToRemove
		elif key =="u":
			localPoints=pointCloud.points.GetNumberOfPoints()
			pointCloud.clearPoints()
			qPointCloud.clearPoints()
			mask=np.in1d(range(localPoints),pointsToRemove)
			mask=np.invert(mask)
			maskedPoints=CurrPnts[mask]		
			for k in maskedPoints:
				pointCloud.addPoint(np.add(k,-centroid))
			renWin.Render()
			CurrPnts=maskedPoints
			pointsToRemove=[]

		#update by removing the current selection
		elif key =="U": #simple undo
			pointsToRemove=[]
			qPointCloud.clearPoints()
			renWin.Render()
			
		elif key =="f": #flip color scheme for printing
			if ren.GetBackground()==(0.1, 0.2, 0.4):
				axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,0)
				axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,0)
				axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,0)
				ren.SetBackground(1, 1, 1)
				pointActor.GetProperty().SetColor(0,0,0)
				if perimActor is not None:
					perimActor.GetProperty().SetColor(0,0,0)
			else:
				axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,1,1)
				axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,1,1)
				axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,1,1)
				pointActor.GetProperty().SetColor(1,1,1)
				ren.SetBackground(0.1, 0.2, 0.4)
				if perimActor is not None:
					perimActor.GetProperty().SetColor(1,1,1)
			renWin.Render()

	#Find the difference between CurrPnts and ReadPoints, store as mask. Write output depending on what's available.
	#Follow the same syntax as used previously: separate variables for x, y and z, outline as x_out
	def WriteOutput():
		global outputd, OutputMask, iren, renWin
		ncols=ReadPoints.shape[1] #number of columns
		dtype={'names':['f{}'.format(i) for i in range(ncols)],
			'formats':ncols * [ReadPoints.dtype]} #make each row an entry in a set
		OutputMask=np.in1d(ReadPoints.view(dtype), CurrPnts.view(dtype))#boolean
		#output mask is passed as integers b/c logical conversion to Matlab are buggy
		OutputMask=np.int8(OutputMask)
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
			Prefix=os.path.basename(filec)
			Prefix=Prefix.split('.')
			print"Now writing output . . ."
			if nosio:
				np.savetxt(os.path.join(outputd,Prefix[0]+"_mod.csv"),
								np.column_stack((ReadPoints,OutputMask)),
								delimiter=',',fmt='%f %f %f %i')
			else:
				if Perim:
					sio.savemat(os.path.join(outputd,Prefix[0]+"_mod.mat"), 
								dict(
								x=np.transpose(np.asmatrix(CurrPnts[:,0])),
								y=np.transpose(np.asmatrix(CurrPnts[:,1])),
								z=np.transpose(np.asmatrix(CurrPnts[:,2])),
								rawPnts=ReadPoints,mask=OutputMask,x_out=ReadPerimeter))
				else:
					sio.savemat(os.path.join(outputd,Prefix[0]+"_mod.mat"), 
								dict(x=CurrPnts[:,0],y=CurrPnts[:,1],z=CurrPnts[:,2],
								rawPnts=ReadPoints,mask=OutputMask))
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
	# print defaultCameraPosition, defaultCameraFocalPoint #debug

	renWin.SetWindowName("UoM Contour Method - Point Cloud Mask Definition v%s" %__version__)
	iren.Start()

	#if the user otherwise closes the interactor window
	if 'iren' in globals():
		print "Interactor closed, no output written."
		close_window(iren)
		del renWin, iren

def ReadMask(*args, **kwargs):
	global Rotating, Panning, Zooming, Selecting, Selection, defaultCameraPosition, defaultCameraFocalPoint
	global iren, renWin, ren
	global Zaspect, pointSize, centroid

	root = tk.Tk()
	root.withdraw()
	
	currentdir=os.getcwd()
	outputd=None
	Perim=True
	Show=True
	
	if len(args)==0:
		filec = askopenfilename(title='Select the data file (.mat)',
								initialdir=currentdir)
		startdir = os.path.dirname(filec)
		if filec == '':
			print 'No file selected; exiting.'
			return
	#various argument collections
	elif len(args)==1:
		filec=args[0]
	elif len(args)==2:
		filec=args[0]
		Show=False
	else:
		print 'Arguments not specified correctly. Quitting.'
		return

	if nosio:
		print 'Error. Requires Scipy to read file.'
		return
	else:
		mat_contents = sio.loadmat(filec)
		try:
			ReadPoints=np.hstack((mat_contents['x'],mat_contents['y'],mat_contents['z']))
			try:
				ReadMask=(mat_contents['mask'])
				RawPnts=(mat_contents['rawPnts'])
				ReadMask=np.bool_(ReadMask[0])
				maskedPnts=RawPnts[~ReadMask,:] #display points that are 'off'
			except KeyError:
				print "No mask found. Run MaskDef first."
				return
			try:
				ReadPerimeter=(mat_contents['x_out'])
			except KeyError:
				Perim=false
		except KeyError:
			print "Couldn't read variables from file. Quitting."
			return

	if Show:
		centroid=np.mean(ReadPoints,axis=0)
		#get extents of dataset
		min=np.amin(ReadPoints,axis=0)
		max=np.amax(ReadPoints,axis=0)
		pointSize=2



		# Create instances
		pointCloud = VtkPointCloud()
		existingMaskedCloud=VtkPointCloud()
		# Add points to to the VtkPointCloud Instance
		for k in ReadPoints:
			pointCloud.addPoint(np.add(k,-centroid))
		for k in maskedPnts:
			existingMaskedCloud.addPoint(np.add(k,-centroid))
		existingMaskedCloud.vtkActor.GetProperty().SetColor(1,0.549,0) #make the existing masked cloud dark orange

		#test/create perimeter
		if Perim:
			perimActor=initializePerimeter(ReadPerimeter,centroid)
		else: perimActor=None

		pointActor = pointCloud.vtkActor
		existingActor=existingMaskedCloud.vtkActor

		# Add axes and fix automatic scaling issue
		axes = vtk.vtkAxesActor()

		axes.GetYAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
		axes.GetXAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
		axes.GetZAxisCaptionActor2D().GetTextActor().SetTextScaleModeToNone()
		axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(7*pointSize)
		axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(7*pointSize)
		axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetFontSize(7*pointSize)
		axScale=np.amax(np.fabs(ReadPoints))*0.1 #scale to dataset
		axes.SetTotalLength(axScale,axScale,axScale)


		# Create the Renderer and assign actors to it. A renderer is like a
		# viewport. It is part or all of a window on the screen and it is
		# responsible for drawing the actors it has.  We also set the
		# background color here.
		ren = vtk.vtkRenderer()
		ren.AddActor(pointActor)
		ren.AddActor(existingActor)
		ren.AddActor(axes)
		if Perim:
			ren.AddActor(perimActor)
		ren.SetBackground(0.1, 0.2, 0.4)

		# Finally we create the render window which will show up on the screen
		# We put our renderer into the render window using AddRenderer. We
		# also set the size to be 1280 pixels by 720 (720p).
		renWin = vtk.vtkRenderWindow()
		renWin.AddRenderer(ren)
		renWin.SetSize(1280, 720)
		# renWin.SetWindowName("UoM Contour Method - Point Editor")

		# Define custom interaction.
		iren = vtk.vtkRenderWindowInteractor()
		iren.SetInteractorStyle(None)
		iren.SetRenderWindow(renWin)

		# Add the observers to watch for particular events. These invoke
		# Python functions.
		Rotating = 0
		Panning = 0
		Zooming = 0
		Selecting = 0
		Selection = False
		rectActor=None
		Zaspect=1



		iren.AddObserver("LeftButtonPressEvent", ButtonEvent)
		iren.AddObserver("LeftButtonReleaseEvent", ButtonEvent)
		iren.AddObserver("MiddleButtonPressEvent", ButtonEvent)
		iren.AddObserver("MiddleButtonReleaseEvent", ButtonEvent)
		iren.AddObserver("RightButtonPressEvent", ButtonEvent)
		iren.AddObserver("RightButtonReleaseEvent", ButtonEvent)
		iren.AddObserver("MouseMoveEvent", MouseMove)
		def Keypress(obj, event):
			global Selection, ren, renWin, Zaspect, pointSize
			key = obj.GetKeySym()
			if key =="1":
				XYView(ren, ren.GetActiveCamera())
			elif key =="2":
				YZView(ren, ren.GetActiveCamera())
			elif key =="3":
				XZView(ren, ren.GetActiveCamera())
			#change z aspect ratio
			elif key=="z":
				Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),Zaspect*2,0)
				updateZaspect(existingMaskedCloud,existingActor,existingMaskedCloud.points.GetNumberOfPoints(),Zaspect*2,0)
			elif key=="Z":
				pointSize=updatePointSize(pointActor,pointSize+1)
			elif key=="x":
				Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),Zaspect*0.5,0)
				updateZaspect(existingMaskedCloud,existingActor,existingMaskedCloud.points.GetNumberOfPoints(),Zaspect*0.5,0)
			elif key=="X":
				pointSize=updatePointSize(pointActor,pointSize-1)
			elif key=="c":
				#clear out all scaling.
				oldZaspect=Zaspect
				Zaspect=updateZaspect(pointCloud,pointActor,pointCloud.points.GetNumberOfPoints(),1,0) #because each call will update Zaspect
				#make sure the old scaling factor is used, not the updated one.
				updateZaspect(existingMaskedCloud,existingActor,existingMaskedCloud.points.GetNumberOfPoints(),Zaspect*(1/oldZaspect),0)
				Zaspect=1
			elif key=="C":
				pointSize=updatePointSize(pointActor,2)
			elif key =="i":
				im = vtk.vtkWindowToImageFilter()
				writer = vtk.vtkPNGWriter()
				im.SetInput(renWin)
				im.Update()
				writer.SetInputConnection(im.GetOutputPort())
				writer.SetFileName("MaskedPointCloud.png")
				writer.Write()
				print 'Screen output saved to %s' %os.path.join(currentdir,'MaskedPointCloud.png')
			elif key =="r":
				if axes.GetVisibility():
					axes.VisibilityOff()
				else:
					axes.VisibilityOn()
				renWin.Render()

			elif key =="f": #flip color scheme for printing
				if ren.GetBackground()==(0.1, 0.2, 0.4):
					axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,0)
					axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,0)
					axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(0,0,0)
					ren.SetBackground(1, 1, 1)
					pointActor.GetProperty().SetColor(0,0,0)
					if perimActor is not None:
						perimActor.GetProperty().SetColor(0,0,0)
				else:
					axes.GetXAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,1,1)
					axes.GetYAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,1,1)
					axes.GetZAxisCaptionActor2D().GetCaptionTextProperty().SetColor(1,1,1)
					pointActor.GetProperty().SetColor(1,1,1)
					ren.SetBackground(0.1, 0.2, 0.4)
					if perimActor is not None:
						perimActor.GetProperty().SetColor(1,1,1)
				renWin.Render()
		iren.AddObserver("KeyPressEvent", Keypress)

		iren.Initialize()
		renWin.Render()
		#get default camera position as inputs for 'views'
		defaultCameraPosition=ren.GetActiveCamera().GetPosition()
		defaultCameraFocalPoint=ren.GetActiveCamera().GetFocalPoint()
		# print defaultCameraPosition, defaultCameraFocalPoint #debug

		renWin.SetWindowName("UoM Contour Method - Masked Data Viewer v%s" %__version__)
		iren.Start()

		del renWin, iren
	return ReadPoints,ReadMask,RawPnts



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

# A duplicate of the VtkPointCloud class, but makes smaller red points.
# Ideally, this should really be a filter, but it's easier to duplicate a class.
# Minimal overhead since HighlightPoints<<Points
class VtkHighlightPointCloud:

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
		self.vtkActor.GetProperty().SetPointSize(pointSize-2)
		self.vtkActor.GetProperty().SetColor(1,0,0)

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

def initializePerimeter(perimeter,centroid):
	points = vtk.vtkPoints()
	for j in perimeter:
		points.InsertNextPoint(np.add(j,-centroid))

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
	return perimActor

##Common functions
# Routines that translate the events into camera motions.
# Handle the mouse button events.
def ButtonEvent(obj, event):
	global Rotating, Panning, Zooming, Selection, Selecting, rectActor
	global iren, renWin, ren
	global xyOffClickPosition, xyOnClickPosition
	if event == "LeftButtonPressEvent" and not Selection:
		Rotating = 1
	if event == "LeftButtonPressEvent" and Selection:
		Selecting = 1
		xyOnClickPosition=iren.GetEventPosition()
	elif event == "LeftButtonReleaseEvent" and not Selection:
		Rotating = 0
	if event == "LeftButtonReleaseEvent" and Selection:
		xyOffClickPosition=iren.GetEventPosition()
		Selecting = 0
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
	camera.SetPosition(defaultCameraPosition)
	camera.SetFocalPoint(defaultCameraFocalPoint)
	camera.SetViewUp(0,1,0)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renWin.Render()

def YZView(renderer, camera):
	camera.SetPosition(defaultCameraPosition[2]+centroid[0],centroid[1],centroid[2])
	camera.SetFocalPoint(defaultCameraFocalPoint)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renWin.Render()

def XZView(renderer, camera):
	camera.SetPosition(centroid[0],defaultCameraPosition[2]+centroid[1],centroid[2])
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

# General high-level logic
def MouseMove(obj, event):
	global Rotating, Panning, Zooming, Selecting
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
	elif Selecting:
		DrawSelectionBox(ren,xyOnClickPosition[0], y, x, xyOnClickPosition[1])

def DrawSelectionBox(renderer,xmin,ymin,xmax,ymax):
	global rectActor, sel, selectRect
	if rectActor is None:
		sel = vtk.vtkPoints()
		sel.InsertPoint(0, xmin, ymin, 0)
		sel.InsertPoint(1, xmax, ymin, 0)
		sel.InsertPoint(2, xmax, ymax, 0)
		sel.InsertPoint(3, xmin, ymax, 0)
		rect = vtk.vtkCellArray()
		rect.InsertNextCell(5)
		rect.InsertCellPoint(0)
		rect.InsertCellPoint(1)
		rect.InsertCellPoint(2)
		rect.InsertCellPoint(3)
		rect.InsertCellPoint(0)
		selectRect = vtk.vtkPolyData()
		selectRect.SetPoints(sel)
		selectRect.SetLines(rect)
		rectMapper = vtk.vtkPolyDataMapper2D()
		if vtk.VTK_MAJOR_VERSION<= 5:
			rectMapper.SetInput(selectRect)
		else:
			rectMapper.SetInputData(selectRect)
		rectActor = vtk.vtkActor2D()
		rectActor.SetMapper(rectMapper)
		ren.AddActor2D(rectActor)
		rectActor.GetProperty().SetColor(0.5,0.5,0.5)
	else:
		#update everything
		sel.SetPoint(0,xmin,ymin,0)
		sel.SetPoint(1,xmax,ymin,0)
		sel.SetPoint(2,xmax,ymax,0)
		sel.SetPoint(3,xmin,ymax,0)
		selectRect.Modified()
		rectActor.Modified()
	renWin.Render()

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

# Updates the point size for a given point actor
def updatePointSize(actor,NewPointSize):
	actor.GetProperty().SetPointSize(NewPointSize)
	actor.Modified()
	renWin.Render()
	return NewPointSize

def close_window(iren):
	render_window=iren.GetRenderWindow()
	render_window.Finalize()
	iren.TerminateApp()

if __name__ == '__main__':
	currentdir=os.getcwd()
	# checkImport()
	if len(sys.argv)>3:
		perimFile=os.path.join(currentdir,sys.argv[1])
		pcloudFile=os.path.join(currentdir,sys.argv[2])
		outDir=os.path.join(currentdir,sys.argv[3])
		MaskDef(perimFile,pcloudFile,outDir)
	elif len(sys.argv)>2:
		perimFile=os.path.join(currentdir,sys.argv[1])
		pcloudFile=os.path.join(currentdir,sys.argv[2])
		MaskDef(perimFile,pcloudFile)
	elif len(sys.argv)>1:
		pcloudFile=os.path.join(currentdir,sys.argv[1])
		MaskDef(pcloudFile)
	else:
		MaskDef()


