#!/usr/bin/env python
'''
Common objects and functions for the pyCM package
'''
__author__ = "M.J. Roy"
__version__ = "0.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2017"

import os,sys,yaml,math
import vtk
import vtk.util.numpy_support as vtk_to_numpy
import numpy as np
import scipy.io as sio
import h5py
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *


def get_file(*args):
	'''
	Returns absolute path to filename and the directory it is located in from a PyQt5 filedialog. First value is file extension, second is a string which overwrites the window message.
	'''
	ext=args[0]
	ftypeName={}
	ftypeName['*.txt']=["Select the point cloud data file:", "*.txt", "TXT File"]
	ftypeName['*.mat']=["Select MAT data file:", "*.mat", "MAT File"]
	ftypeName['*.vtk']=["Select the legacy VTK file:", "*.vtk", "VTK File"]
	ftypeName['*.dat']=["Select FEA results file:", "*.dat", "DAT File"]
	ftypeName['*.inp']=["Select FEA input file:", "*.inp", "INP File"]
	ftypeName['*.vtu']=["Select postprocessing file to view:", "*.vtu", "VTK XML unstructured grid file"]

	if len(args)==2:
		ftypeName[ext][0] = args[1]
		
	lapp = QApplication.instance()
	if lapp is None:
		lapp = QtWidgetsQApplication([])
	if ext=='*.txt':
		filer = QFileDialog.getOpenFileName(None, ftypeName[ext][0], 
         os.getcwd(),(ftypeName[ext][2]+' ('+ftypeName[ext][1]+');;'+ftypeName['*.mat'][2]+' ('+ftypeName['*.mat'][1]+');;All Files (*.*)'))
	else:
		filer = QFileDialog.getOpenFileName(None, ftypeName[ext][0], 
         os.getcwd(),(ftypeName[ext][2]+' ('+ftypeName[ext][1]+');;All Files (*.*)'))

	if filer[0] == '':
		filer = None
		startdir = None
		return filer, startdir
		
	else:
		return filer[0], os.path.dirname(filer[0])

				
	
def get_open_file(ext,outputd):
	'''
	Returns a the complete path to the file name with ext, starting in outputd. Checks extensions and if an extension is not imposed, it will write the appropriate extension based on ext.
	'''
	ftypeName={}
	ftypeName['*.csv']='Comma delimited pyCM output file'
	ftypeName['*.mat']='pyCM output file'
	ftypeName['*.geo']='Gmsh geometry file'
	ftypeName['*.dxf']='Drawing eXchange Format'
	ftypeName['*.py']='Abaqus Python script'
	ftypeName['*.ccx.inp']='Calculix input file'
	ftypeName['*.abq.inp']='Abaqus input file'
	ftypeName['*out.vtk'] = 'VTK legacy file'
	ftypeName['*.vtu'] = 'VTK XML-based post processing file'
	
	
	
	if outputd==None: id=str(os.getcwd())
	else: id=outputd
	lapp = QApplication.instance()
	if lapp is None:
		lapp = QApplication([])

	filer = QFileDialog.getSaveFileName(None, "Save as:", "",str(ftypeName[ext]+' ('+ext+')'))
	
	if filer == '':
		filer = None
		startdir = None
		return filer, startdir

	if filer:
		filer=str(filer[0])
		if filer.endswith(ext[-4:]) and not filer.endswith(ext[-7:]):
			filer=filer[:-4]+ext.split('*')[-1]

		startdir=os.path.split(filer)
	
		return filer, os.path.split(filer)

def extract_from_mat(targetfile,matfile,field):
	'''
	Function to pull ascii data from .mat file and write to a given location which FEA packages need - called in the instance if the user deletes various auxiliary files, this function will recover them from a mat file.
	'''
	mat_contents = sio.loadmat(matfile)
	fid=open(targetfile,'w+')
	fid.write('%s'%mat_contents[field][0])
	fid.close()

def gen_point_cloud(pts,color,size):
	'''
	Returns vtk objects and actor for a point cloud having size points and color mapped to z limits
	'''

	vtkPnts = vtk.vtkPoints()
	vtkVerts = vtk.vtkCellArray()
	
	if color[0]<=1:
		color=(int(color[0]*255),int(color[1]*255),int(color[2]*255))


	colors=vtk.vtkUnsignedCharArray()
	colors.SetNumberOfComponents(3)
	colors.SetName("color")
	
	
	#load up points
	for i in pts:
		pId= vtkPnts.InsertNextPoint(i)
		vtkVerts.InsertNextCell(1)
		vtkVerts.InsertCellPoint(pId)
		colors.InsertNextTuple(color)

		
	
	pC = vtk.vtkPolyData()
	pC.SetPoints(vtkPnts)
	pC.SetVerts(vtkVerts)
	pC.GetPointData().SetScalars(colors)
	
	vtkPntMapper = vtk.vtkDataSetMapper()
	vtkPntMapper.SetInputData(pC)

	actor=vtk.vtkActor()
	actor.SetMapper(vtkPntMapper)

	actor.GetProperty().SetPointSize(size)
	return pC, actor, colors
	
def gen_outline(pts,color,size):
	'''
	Returns an outline actor with specified pts, color and size. Incoming pnts should be ordered.
	'''
	if color[0]<=1 and color != None:
		color=(int(color[0]*255),int(color[1]*255),int(color[2]*255))
	if color[0]>=1 and color != None:
		color=(color[0]/float(255),color[1]/float(255),color[2]/float(255))
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
	outline.SetVerts(linesegcells)
	outline.SetLines(linesegcells)
	Omapper=vtk.vtkPolyDataMapper()
	Omapper.SetInputData(outline)
	outlineActor=vtk.vtkActor()
	outlineActor.SetMapper(Omapper)
	outlineActor.GetProperty().SetColor(color)
	outlineActor.GetProperty().SetPointSize(size)
	return outlineActor,outline

def xyview(renderer, camera,cp,fp):
	camera.SetPosition(0,0,cp[2]+0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,1,0)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()
	
def xyview_post(renderer, camera,cp,fp):
	# print('cp',cp,'fp',fp)
	camera.SetPosition(0,0,-cp[2]+0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,-1,0)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()

def yzview(renderer, camera,cp,fp):
	camera.SetPosition(cp[2]+0,0,0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()


def xzview(renderer,camera,cp,fp):
	vtk.vtkObject.GlobalWarningDisplayOff() #otherwise there's crystal eyes error . . .
	camera.SetPosition(0,cp[2]+0,0)
	camera.SetFocalPoint(fp)
	camera.SetViewUp(0,0,1)
	camera.OrthogonalizeViewUp()
	camera.ParallelProjectionOn()
	renderer.ResetCamera()
	
def flip_visible(actor):
	'''
	Convenience function for changing the visibility of actors
	'''
	if actor.GetVisibility():
		actor.VisibilityOff()
	else:
		actor.VisibilityOn()	

		
def flip_colors(ren,actor):
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
			
			
def update_point_size(actor,NewPointSize):
	'''
	Updates the incoming actor's pointsize and returns the new one.
	'''
	actor.GetProperty().SetPointSize(NewPointSize)
	actor.Modified()
	renWin.Render()
	return NewPointSize
	
def read_mat(file):
	'''
	Reads a .mat file from UoM's MATLAB routine outliner.m which preprocesses NanoFocus .dat file (Origin format). Returns points and an outline in numpy arrays.
	'''
	try:
		mat_contents = sio.loadmat(file)
		rawPnts=np.hstack((mat_contents['x'],mat_contents['y'],mat_contents['z']))
		outline=(mat_contents['x_out'])
		return rawPnts, outline
	except:
		print("Couldn't data from %s."%file)
		return

def read_abq_dat(file_name):

	i=0
	lineFlag=[]
	keyStrings=["    ELEMENT  PT"," MAXIMUM"]

	fid = open(file_name)
	while 1:
		lines = fid.readlines(100000)
		if not lines:
			break
		for line in lines:
			i+=1
			for keyString in keyStrings:
				if line[0:len(keyString)]==keyString:
					lineFlag.append(i)
				
	fid.close()
	# extract quadrature data for
	# node id, quadrature id, x coord, y coord, z coord, S33
	# No idea why there's a full element's worth of an offset
	return np.genfromtxt(file_name, skip_header=lineFlag[0]+2, skip_footer=i-lineFlag[1]-8, \
							usecols=(0, 2, 3, 4, 5, 6, 7), autostrip=True)
	

def read_ccx_dat(file_name):
	read_res = np.genfromtxt(file_name, skip_header=3, skip_footer=0, \
							usecols=(0, 2, 3, 4), autostrip=True)
	#pad with zeros to match what's being read by read_abq_dat
	padding = np.zeros([len(read_res),3])
	return np.concatenate((np.reshape(read_res[:,0],(len(read_res[:,0]),1)),padding, \
		read_res[:,1::]), axis=1)
	
	
def set_size_policy(target_widget):
	sizePolicy=QSizePolicy(QSizePolicy.MinimumExpanding, QSizePolicy.MinimumExpanding)
	sizePolicy.setHorizontalStretch(0)
	sizePolicy.setVerticalStretch(0)
	sizePolicy.setHeightForWidth(target_widget.sizePolicy().hasHeightForWidth())
	target_widget.setSizePolicy(sizePolicy)
	target_widget.setMinimumSize(1080, 600); #leave 100 px on the size for i/o

def clear_mat(matfile,fields):
	'''
	Remove specified fields from a .mat file, called in conjunction to a 'get_input_data' function for each respective step, keeps mat file current
	'''

	mat_contents=sio.loadmat(matfile) #dictionary of variable names

	for field in fields:
		if field in mat_contents:
			del mat_contents[field]
		
	sio.savemat(matfile,mat_contents)



class Ui_getFEAconfigDialog(object):
	def setupUi(self, getFEAconfigDialog):
		getFEAconfigDialog.setObjectName("getFEAconfigDialog")
		getFEAconfigDialog.resize(630, 201)
		self.verticalLayout = QVBoxLayout(getFEAconfigDialog)
		self.verticalLayout.setObjectName("verticalLayout")
		
		self.horizontalLayout_3 = QHBoxLayout()
		self.horizontalLayout_3.setObjectName("horizontalLayout_3")
		self.dialoglabel = QLabel(getFEAconfigDialog)
		self.horizontalLayout_3.addWidget(self.dialoglabel)
		self.verticalLayout.addLayout(self.horizontalLayout_3)
		
		self.horizontalLayout_2 = QHBoxLayout()
		self.abaExec = QLineEdit(getFEAconfigDialog)
		self.abaLabel = QLabel(getFEAconfigDialog)
		self.horizontalLayout_2.addWidget(self.abaLabel)
		self.horizontalLayout_2.addWidget(self.abaExec)
		self.verticalLayout.addLayout(self.horizontalLayout_2)

		self.horizontalLayout_1 = QHBoxLayout()
		self.gmshExec = QLineEdit(getFEAconfigDialog)
		self.gmshLabel = QLabel(getFEAconfigDialog)
		self.horizontalLayout_1.addWidget(self.gmshLabel)
		self.horizontalLayout_1.addWidget(self.gmshExec)
		self.verticalLayout.addLayout(self.horizontalLayout_1)
		
		self.horizontalLayout_0 = QHBoxLayout()
		self.ccxExec = QLineEdit(getFEAconfigDialog)
		self.ccxLabel = QLabel(getFEAconfigDialog)
		self.horizontalLayout_0.addWidget(self.ccxLabel)
		self.horizontalLayout_0.addWidget(self.ccxExec)
		self.verticalLayout.addLayout(self.horizontalLayout_0)
		
		
		self.horizontalLayout = QHBoxLayout()
		self.ConfigFileLoc = QLabel(getFEAconfigDialog)
		self.ConfigFileLoc.setFont(QFont("Helvetica",italic=True))
		self.ConfigFileLabel = QLabel(getFEAconfigDialog)
		self.pushButton = QPushButton(getFEAconfigDialog)
		self.horizontalLayout.addWidget(self.ConfigFileLabel)
		self.horizontalLayout.addWidget(self.ConfigFileLoc)
		self.horizontalLayout.addWidget(self.pushButton)
		self.verticalLayout.addLayout(self.horizontalLayout)

		self.retranslateUi(getFEAconfigDialog)
		QMetaObject.connectSlotsByName(getFEAconfigDialog)
		
		# connect the two functions
		self.pushButton.clicked.connect(lambda: self.makeConfigChange(getFEAconfigDialog))
		# self.pushButton_2.clicked.connect(self.return_cancel)

	def retranslateUi(self, getFEAconfigDialog):
		_translate = QCoreApplication.translate
		getFEAconfigDialog.setWindowTitle(_translate("getFEAconfigDialog", "pyCM FEA configuration settings"))
		self.dialoglabel.setText(_translate("getFEAconfigDialog", "Please specify the full path to relevant executable, or leave blank if not relevant."))
		self.abaLabel.setText(_translate("getFEAconfigDialog", "Abaqus executable"))
		self.gmshLabel.setText(_translate("getFEAconfigDialog", "Gmsh executable"))
		self.ccxLabel.setText(_translate("getFEAconfigDialog", "Calculix executable"))
		# self.pushButton_2.setText(_translate("getFEAconfigDialog", "Acept"))
		self.ConfigFileLabel.setText(_translate("getFEAconfigDialog", "Current config file:",None))
		self.ConfigFileLoc.setText(_translate("getFEAconfigDialog", "Undefined", None))
		self.pushButton.setText(_translate("getFEAconfigDialog", "Set"))


	def makeConfigChange(self, getFEAconfigDialog):
		try:
			data= dict(FEA = 
			dict(
			abaqusExec = str(self.abaExec.text()),
			gmshExec = str(self.gmshExec.text()),
			ccxExec = str(self.ccxExec.text()),
			)
			)
			with open(str(self.ConfigFileLoc.text()), 'w') as outfile:
				yaml.dump(data, outfile, default_flow_style=False)
		except:
			print("Configuration change failed.")
			

		getFEAconfigDialog.close()
