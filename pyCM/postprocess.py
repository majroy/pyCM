"""
Uses VTK python to allow for postprocessing FEAs associated with the
contour method. Full interaction requires a 3-button mouse and keyboard.
-------------------------------------------------------------------------------
Current mapping is as follows:
LMB   - rotate about point cloud centroid.
MMB   - pan
RMB   - zoom
-------------------------------------------------------------------------------
ver 0.1 21 April 2018
"""
__author__ = "N. Stoyanov"
__version__ = "0.1"
__email__ = "nikola.stoyanov@postgrad.manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, N. Stoyanov 2014-2018"

import sys, time
import vtk
import pandas
import numpy as np
import scipy.io as sio
from shutil import copyfile
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtGui, QtWidgets, QtCore
from vtk.util.numpy_support import vtk_to_numpy as v2n
from pkg_resources import Requirement, resource_filename
from pyCMcommon import *


def post_process_tool(*args, **kwargs):
	"""
	Build QT interaction
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

	window = pp_interactor(None)
	window.show()
	splash.finish(window)
	window.iren.Initialize()

	ret = app.exec_()

	if sys.stdin.isatty() and not hasattr(sys, 'ps1'):
		sys.exit(ret)
	else:
		return window


class post_main_window(object):
	"""
	Setup and handle VTK and QT interaction
	"""

	def setupUi(self, MainWindow):
		"""
		Performs the Qt setup of GUI
		Creates a windows and adds button interaction
		"""

		MainWindow.setWindowTitle("pyCM - FEA postprocessing v%s" %__version__)
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

		self.horiz_line1 = QtWidgets.QFrame()
		self.horiz_line1.setFrameStyle(QtWidgets.QFrame.HLine)
		self.head_font=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)


		self.quadButton = QtWidgets.QRadioButton("quads")
		self.tetButton = QtWidgets.QRadioButton("tets")
		self.quadButton.setChecked(True)
		self.mtypeButtonGroup = QtWidgets.QButtonGroup()
		self.mtypeButtonGroup.addButton(self.tetButton)
		self.mtypeButtonGroup.addButton(self.quadButton)
		self.mtypeButtonGroup.setExclusive(True)

		self.read_sim_data_button = QtWidgets.QPushButton('Extract S33')
		self.display_vtk_button = QtWidgets.QPushButton('Display S33')

		self.horiz_line3 = QtWidgets.QFrame()
		self.horiz_line3.setFrameStyle(QtWidgets.QFrame.HLine)
		self.stress_label = QtWidgets.QLabel("Contour levels")
		self.stress_label.setWordWrap(True)
		self.stress_label.setFont(QtGui.QFont("Helvetica"))
		self.stress_label.setMinimumWidth(50)

		# enter minimum stress to update the display
		inp_min_stress_label = QtWidgets.QLabel("Min S33:")
		self.inp_min_stress = QtWidgets.QLineEdit()
		self.inp_min_stress.setMinimumWidth(50)
		
		num_contours_label = QtWidgets.QLabel("Contours:")
		self.numContour = QtWidgets.QSpinBox()
		self.numContour.setMinimum(3)
		self.numContour.setMaximum(20)
		self.numContour.setValue(5)

		# enter maximum stress to update the display
		inp_max_stress_label = QtWidgets.QLabel("Max S33:")
		self.inp_max_stress = QtWidgets.QLineEdit()
		self.inp_max_stress.setMinimumWidth(50)
		

		# stress update button
		self.updateButton = QtWidgets.QPushButton('Update')
		self.updateButton.setMinimumWidth(50)

		self.horiz_line2 = QtWidgets.QFrame()
		self.horiz_line2.setFrameStyle(QtWidgets.QFrame.HLine)
		self.statLabel = QtWidgets.QLabel("Idle")
		self.statLabel.setWordWrap(True)
		self.statLabel.setFont(QtGui.QFont("Helvetica", italic=True))
		self.statLabel.setMinimumWidth(50)

		mainUiBox.addWidget(self.quadButton,0,0,1,1)
		mainUiBox.addWidget(self.tetButton,0,1,1,1)
		mainUiBox.addWidget(self.read_sim_data_button,1,0,1,1)
		mainUiBox.addWidget(self.display_vtk_button,1,1,1,1)
		mainUiBox.addWidget(self.horiz_line3,2,0,1,2)
		mainUiBox.addWidget(self.stress_label,3,0,1,2)
		mainUiBox.addWidget(inp_min_stress_label,4,0,1,1)
		mainUiBox.addWidget(self.inp_min_stress,4,1,1,1)
		mainUiBox.addWidget(inp_max_stress_label,5,0,1,1)
		mainUiBox.addWidget(self.inp_max_stress,5,1,1,1)
		mainUiBox.addWidget(num_contours_label,6,0,1,1)
		mainUiBox.addWidget(self.numContour,6,1,1,1)
		mainUiBox.addWidget(self.inp_max_stress,7,1,1,1)
		mainUiBox.addWidget(self.updateButton,7,0,1,2)
		mainUiBox.addWidget(self.horiz_line2,8,0,1,2)

		
		lvLayout=QtWidgets.QVBoxLayout()
		lvLayout.addLayout(mainUiBox)
		lvLayout.addStretch(1)
		
		self.mainlayout.addWidget(self.vtkWidget,0,0,1,1)
		self.mainlayout.addLayout(lvLayout,0,1,1,1)
		self.mainlayout.addWidget(self.statLabel,1,0,1,2)
		

		
		def initialize(self):
			self.vtkWidget.start()

class pp_interactor(QtWidgets.QWidget):
	"""
	Sets up the main VTK window
	reads .mat file and sets connection between UI and interactor
	"""

	def __init__(self,parent):
		super(pp_interactor,self).__init__(parent)
		self.ui = post_main_window()
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

		self.vtk_file = None

		self.ui.read_sim_data_button.clicked.connect(lambda: self.load_integration_points())
		self.ui.display_vtk_button.clicked.connect(lambda: self.load_vtk_legacy_file())
		self.ui.updateButton.clicked.connect(lambda: self.update_stress_display())

	def get_input_data(self,filem):
		if filem == None:
			filem,_,=get_file('*.mat')
		
		if filem:
			mat_contents = sio.loadmat(filem)
			self.fileo=filem

		
		try:
			#get dat file
			FEAbasename=mat_contents['FEA_filename'][0]
			filename, _ = os.path.splitext(FEAbasename)
			self.dat_file=filename+'.dat'
			#check if the dat file is valid
			if not os.path.exists(self.dat_file):
				ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
				"Could not find a *.dat file associated with this analysis. Would you like to search for it yourself?", \
				QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
				if ret == QtWidgets.QMessageBox.Yes:
					self.dat_file,_ = get_file("*.dat")
				else:
					return
		except Exception as e:
			print('Could not lock on dat file associated with analysis')
			print(e)
			
		try:
			#get vtk output file - the one written on execution of the FEA
			vtkbasename=mat_contents['vtk_filename'][0]
			filename, _ = os.path.splitext(vtkbasename)
			self.vtk_file=filename[0:-3]+'out.vtk'
			if not os.path.exists(self.vtk_file):
				#attempt to write it from the mat file
				try:
					extract_from_mat(self.vtk_file,self.fileo,'vtk_out')
					print('Re-wrote VTK file.')
				except:
					ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
				"The host directory specified in the originating analysis cannot be written to. Would you like to specify a new location?", \
				QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
					if ret == QtWidgets.QMessageBox.Yes:
						self.vtk_file,_ = get_open_file("*_out.vtk")
						extract_from_mat(self.vtk_file,self.fileo,'vtk_out')
						print('Wrote *_out.vtk file in a different location than specified in results file.')
					else:
						return
		except:
			print('Could not lock on dat file associated with analysis')
			
	
	def update_stress_display(self):
		"""
		Updates the display to reflect the input of min/max stress
		"""
		QtWidgets.QApplication.processEvents()

		min_stress = float(self.ui.inp_min_stress.text())
		max_stress = float(self.ui.inp_max_stress.text())

		self.mesh_mapper.SetScalarRange(min_stress, max_stress)
		self.sbActor.SetNumberOfLabels(self.ui.numContour.value())
		self.ui.vtkWidget.update()
		self.ui.statLabel.setText("Updated stress intervals on contour. Idle.")

		
		
	def load_vtk_legacy_file(self,file):
		"""
		Loads the vtk mesh and displays the scalar data in a color map.
		Allows further postprocessing to be done, such as grayscale and contour plots.
		"""

		if hasattr(self, "mesh_actor"):
			self.ren.RemoveActor(self.mesh_actor)
			self.ren.RemoveActor(self.sbActor)
		
		if file is None:
			file,_ = get_file("*.vtk")
		
		
		self.ui.statLabel.setText("Reading %s for mesh . . ."%file)
		mesh_source = vtk.vtkUnstructuredGridReader()
		mesh_source.SetFileName(file)

		# read scalar to vtk
		mesh_source.SetScalarsName("S33")
		mesh_source.Update()
		mesh_reader_output = mesh_source.GetOutput()

		
		# bounds for axis
		bounds = mesh_reader_output.GetBounds()

		# show element edges
		edges = vtk.vtkExtractEdges()
		edges.SetInputConnection(mesh_source.GetOutputPort())
		edges.Update()

		# lookup table and scalar range for a vtk file
		mesh_lookup_table = vtk.vtkLookupTable()

		# make scalar red = max; blue = min
		self.ui.statLabel.setText("Building lookup table . . .")
		self.draw_color_range(mesh_lookup_table)
		mesh_lookup_table.Build()
		scalar_range = mesh_reader_output.GetScalarRange()

		# mesh data set
		self.mesh_mapper = vtk.vtkDataSetMapper()
		self.mesh_mapper.SetInputData(mesh_reader_output)
		self.mesh_mapper.SetScalarRange(scalar_range)
		self.mesh_mapper.SetLookupTable(mesh_lookup_table)

		#define actors
		self.mesh_actor = vtk.vtkActor()
		# self.sbActor = vtk.vtkScalarBarActor()
		


		# #the scalar bar widget is associated with the qt interactor box
		scalar_bar_widget = vtk.vtkScalarBarWidget()
		scalar_bar_widget.SetInteractor(self.iren)
		# scalar_bar_widget.SetCurrentRenderer(self.ren.GetRenderer())
		# scalar_bar_widget.SetDefaultRenderer(self.ren.GetRenderer())
		scalar_bar_widget.SetEnabled(True)
		scalar_bar_widget.RepositionableOn()
		scalar_bar_widget.On()
		
		# define scalar bar actor
		self.sbActor=scalar_bar_widget.GetScalarBarActor()
		# self.sbActor.SetOrientationToVertical()
		self.sbActor.SetLookupTable(mesh_lookup_table)
		self.sbActor.SetTitle("S33")


		scalarBarRep = scalar_bar_widget.GetRepresentation()
		scalarBarRep.GetPositionCoordinate().SetValue(0.01,0.01)
		scalarBarRep.GetPosition2Coordinate().SetValue(0.09,0.9)
		
		#attempt to change scalebar properties [ineffective]
		propT = vtk.vtkTextProperty()
		propL = vtk.vtkTextProperty()
		propT.SetFontFamilyToArial()
		# propT.ItalicOff()
		propT.BoldOn()
		propL.BoldOff()
		propL.SetFontSize(1)
		propT.SetFontSize(2)
		self.sbActor.SetTitleTextProperty(propT);
		self.sbActor.SetLabelTextProperty(propL);
		self.sbActor.GetLabelTextProperty().SetFontSize(7)
		self.sbActor.GetTitleTextProperty().SetFontSize(7)
		self.sbActor.SetLabelFormat("%.1f")




		#define the mesh actor properties
		self.mesh_actor.SetMapper(self.mesh_mapper)
		self.mesh_actor.GetProperty().SetLineWidth(1)
		self.mesh_actor.GetProperty().EdgeVisibilityOn()

		#display the actors
		self.ren.AddActor(self.mesh_actor)
		self.ren.AddActor(self.sbActor)
		
		#get boundary of mesh
		self.limits = mesh_reader_output.GetBounds()
		
		
		self.ui.vtkWidget.setFocus()
		self.AddAxis(self.limits,1)
		xyview_post(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp) #sorts out the camera issue
		self.ui.vtkWidget.update()
		self.ui.statLabel.setText("Loaded results. Idle.")
		QtWidgets.QApplication.processEvents()
		
		
	def load_integration_points(self):
		"""
		Writes the FE output data - gaussian integration points
		to the legacy vtk file. The FE data is stored in the .dat file for ABAQUS.
		"""
		QtWidgets.QApplication.processEvents()

		# load input files
		# use the original *.inp file, not *.abq.inp
		# dat_file,_ = get_file("*.dat")
		# vtk_file,_ = get_file("*.vtk")

		self.ui.statLabel.setText("Reading results from %s . . ."%self.dat_file)
		try: quadrature_data = read_abq_dat(self.dat_file)
		except Exception as e: 
			print('Could not read specified .dat file.')
			print(e)
			return

		self.ui.statLabel.setText("Retrieving mesh information from %s . . ."%self.vtk_file)
		# # node_data, element_data = self.get_node_data_abq(inp_file)
		node_data, element_data = self.get_node_data_n(self.vtk_file)
		
		#debug
		# print('vtk node',node_data)
		# print('vtk elem',element_data)
		# node_data, element_data = self.get_node_data_abq(inp_file)
		# print('abq node',node_data)
		# print('abq elem',element_data)
		

		# obtain a matrix of node number, x, y, z, stress
		# check if the discretisation was done with brick (C3D8) or tetrahedra (C3D10) elements
		if self.ui.quadButton.isChecked():
			stress_array = self.calculate_quadrature_stress_C3D8(quadrature_data, element_data, node_data)
		else:
			stress_array = self.calculate_quadrature_stress_C3D10(quadrature_data, element_data, node_data)

		# nodes will duplicate in elements
		# sum contributions from nodes
		stress_data_frame = pandas.DataFrame(stress_array)
		stress_data_frame[5] = stress_data_frame.groupby([0])[4].transform(np.mean)

		# remove the individual stress at nodes
		del stress_data_frame[4]

		# drop the duplicates
		stress_data_frame = pandas.DataFrame(stress_data_frame).drop_duplicates(subset=[0])
		stress_data_frame = pandas.DataFrame(stress_data_frame).values
		stress_data_frame = pandas.DataFrame(stress_data_frame)

		# get *.vtk file
		# vtk_file, _ = get_file("*.vtk")
		# copy vtk file
		filename, _ = os.path.splitext(self.vtk_file)
		new_vtk_file = filename+'_post.vtk'
		copyfile(self.vtk_file,new_vtk_file)
		

		# append at the end of the vtk file
		self.append_postprocess_data(new_vtk_file, stress_data_frame)
		
		#show the result
		self.load_vtk_legacy_file(new_vtk_file)

	def append_postprocess_data(self, vtk_file, stress_data_frame):
		"""
		Append the postprocessed data to the end of the *.vtk file
		"""

		# sort the dataframe since the nodes are by connectivity from elements
		# vtk requires them by numbering from 0
		stress_data_frame = stress_data_frame.sort_values(by=[0])

		# open the *.vtk file in append mode
		file_handle = open(vtk_file, 'ab')

		# add vtk headers
		#vtk_header_text = "SCALARS S33 \nLOOKUP_ TABLE default\n"
		file_handle.write(str.encode('POINT_DATA %i\nSCALARS S33 DOUBLE\nLOOKUP_TABLE default\n' % \
									stress_data_frame[4].count()))

		# add the dataframe
		np.savetxt(file_handle, stress_data_frame[4].values, fmt='%.6f')
		file_handle.close()

	def calculate_quadrature_stress_C3D8(self, quadrature_data, element_data, node_data):
		"""
		Calculate the stress values from quadrature and element data
		for element C3D8.
		"""

		# default step for the C3D8 element -> number of nodes
		element_step = 8

		# define the element counter
		element_index = 0

		# define the counter for the shape matrix
		shape_matrix_index = 0

		# define nodal coordinates and stress storage
		stress_array = np.zeros(shape=(len(quadrature_data), 5))

		for row_index in range(0, len(quadrature_data)-1, element_step):
			# extract the stresses at the quadrature points
			quadrature_point_1 = quadrature_data[row_index, 4]
			quadrature_point_2 = quadrature_data[row_index + 1, 4]
			quadrature_point_3 = quadrature_data[row_index + 2, 4]
			quadrature_point_4 = quadrature_data[row_index + 3, 4]
			quadrature_point_5 = quadrature_data[row_index + 4, 4]
			quadrature_point_6 = quadrature_data[row_index + 5, 4]
			quadrature_point_7 = quadrature_data[row_index + 6, 4]
			quadrature_point_8 = quadrature_data[row_index + 7, 4]

			# construct the quadrature stress matrix
			quadrature_stress = np.array([quadrature_point_1,
										quadrature_point_2,
										quadrature_point_3,
										quadrature_point_4,
										quadrature_point_5,
										quadrature_point_6,
										quadrature_point_7,
										quadrature_point_8])

			# extract element row
			element_row = element_data[element_index, :]
			element_index = element_index + 1

			# find the nodal points in the element
			# the int conversion could be done so much better
			# in the future i have to extract this as a structured array
			# and set as int
			nodal_point_1 = node_data[int(element_row[1]) - 1, :]
			nodal_point_2 = node_data[int(element_row[2]) - 1, :]
			nodal_point_3 = node_data[int(element_row[3]) - 1, :]
			nodal_point_4 = node_data[int(element_row[4]) - 1, :]
			nodal_point_5 = node_data[int(element_row[5]) - 1, :]
			nodal_point_6 = node_data[int(element_row[6]) - 1, :]
			nodal_point_7 = node_data[int(element_row[7]) - 1, :]
			nodal_point_8 = node_data[int(element_row[8]) - 1, :]

			# create the square shape function matrix
			#shape_function_matrix = np.zeros(shape=(8,8))

			# obtain the natural coordinates of the gauss points
			# and apply the shape functions
			shape_function_matrix = self.C3D8_quadrature_points()

			# extrapolate from quadrature points to nodal points
			nodal_stress = shape_function_matrix.dot(quadrature_stress)

			# create an array with nodal coordinates and stress
			nodal_data1 = np.array([[nodal_point_1[0], nodal_point_1[1], nodal_point_1[2], nodal_point_1[3], nodal_stress[0]]])
			nodal_data2 = np.array([[nodal_point_2[0], nodal_point_2[1], nodal_point_2[2], nodal_point_2[3], nodal_stress[1]]])
			nodal_data3 = np.array([[nodal_point_3[0], nodal_point_3[1], nodal_point_3[2], nodal_point_3[3], nodal_stress[2]]])
			nodal_data4 = np.array([[nodal_point_4[0], nodal_point_4[1], nodal_point_4[2], nodal_point_4[3], nodal_stress[3]]])
			nodal_data5 = np.array([[nodal_point_5[0], nodal_point_5[1], nodal_point_5[2], nodal_point_5[3], nodal_stress[4]]])
			nodal_data6 = np.array([[nodal_point_6[0], nodal_point_6[1], nodal_point_6[2], nodal_point_6[3], nodal_stress[5]]])
			nodal_data7 = np.array([[nodal_point_7[0], nodal_point_7[1], nodal_point_7[2], nodal_point_7[3], nodal_stress[6]]])
			nodal_data8 = np.array([[nodal_point_8[0], nodal_point_8[1], nodal_point_8[2], nodal_point_8[3], nodal_stress[7]]])

			# collate the data from all nodes
			stress_array[row_index, :] = nodal_data1
			stress_array[row_index + 1, :] = nodal_data2
			stress_array[row_index + 2, :] = nodal_data3
			stress_array[row_index + 3, :] = nodal_data4
			stress_array[row_index + 4, :] = nodal_data5
			stress_array[row_index + 5, :] = nodal_data6
			stress_array[row_index + 6, :] = nodal_data7
			stress_array[row_index + 7, :] = nodal_data8

		return stress_array

	def calculate_quadrature_stress_C3D10(self, quadrature_data, element_data, node_data):
		"""
		Calculate the stress values from quadrature and element data
		for element C3D10
		"""

		# default step for the C3D10 element -> number of nodes
		element_step = 4

		# define the element counter
		element_index = 0

		# define the counter for the shape matrix
		shape_matrix_index = 0

		# define nodal coordinates and stress storage
		stress_array = np.zeros(shape=(int(len(quadrature_data) * 2.5), 5))

		stress_array_row = 0

		for row_index in range(0, len(quadrature_data)-1, element_step):
			# extract the stresses at the quadrature points
			quadrature_point_1 = quadrature_data[row_index, 4]
			quadrature_point_2 = quadrature_data[row_index + 1, 4]
			quadrature_point_3 = quadrature_data[row_index + 2, 4]
			quadrature_point_4 = quadrature_data[row_index + 3, 4]

			# construct the quadrature stress matrix
			quadrature_stress = np.array([quadrature_point_1,
										quadrature_point_2,
										quadrature_point_3,
										quadrature_point_4])

			# extract element row
			element_row = element_data[element_index, :]
			element_index = element_index + 1

			# find the nodal points in the element
			# the int conversion could be done so much better
			# in the future i have to extract this as a structured array
			# and set as int
			nodal_point_1 = node_data[int(element_row[1]) - 1, :]
			nodal_point_2 = node_data[int(element_row[2]) - 1, :]
			nodal_point_3 = node_data[int(element_row[3]) - 1, :]
			nodal_point_4 = node_data[int(element_row[4]) - 1, :]
			nodal_point_5 = node_data[int(element_row[5]) - 1, :]
			nodal_point_6 = node_data[int(element_row[6]) - 1, :]
			nodal_point_7 = node_data[int(element_row[7]) - 1, :]
			nodal_point_8 = node_data[int(element_row[8]) - 1, :]
			nodal_point_9 = node_data[int(element_row[9]) - 1, :]
			nodal_point_10 = node_data[int(element_row[10]) - 1, :]

			# obtain the natural coordinates of the gauss points
			# and apply the shape functions
			shape_function_matrix = self.C3D10_quadrature_points()

			# extrapolate from quadrature points to nodal points
			shape_function_matrix = shape_function_matrix.T
			nodal_stress = shape_function_matrix.dot(quadrature_stress)

			# create an array with nodal coordinates and stress
			nodal_data1 = np.array([[nodal_point_1[0], nodal_point_1[1], nodal_point_1[2], nodal_point_1[3], nodal_stress[0]]])
			nodal_data2 = np.array([[nodal_point_2[0], nodal_point_2[1], nodal_point_2[2], nodal_point_2[3], nodal_stress[1]]])
			nodal_data3 = np.array([[nodal_point_3[0], nodal_point_3[1], nodal_point_3[2], nodal_point_3[3], nodal_stress[2]]])
			nodal_data4 = np.array([[nodal_point_4[0], nodal_point_4[1], nodal_point_4[2], nodal_point_4[3], nodal_stress[3]]])
			nodal_data5 = np.array([[nodal_point_5[0], nodal_point_5[1], nodal_point_5[2], nodal_point_5[3], nodal_stress[4]]])
			nodal_data6 = np.array([[nodal_point_6[0], nodal_point_6[1], nodal_point_6[2], nodal_point_6[3], nodal_stress[5]]])
			nodal_data7 = np.array([[nodal_point_7[0], nodal_point_7[1], nodal_point_7[2], nodal_point_7[3], nodal_stress[6]]])
			nodal_data8 = np.array([[nodal_point_8[0], nodal_point_8[1], nodal_point_8[2], nodal_point_8[3], nodal_stress[7]]])
			nodal_data9 = np.array([[nodal_point_9[0], nodal_point_9[1], nodal_point_9[2], nodal_point_9[3], nodal_stress[8]]])
			nodal_data10 = np.array([[nodal_point_10[0], nodal_point_10[1], nodal_point_10[2], nodal_point_10[3], nodal_stress[9]]])

			# collate the data from all nodes
			stress_array[stress_array_row, :] = nodal_data1
			stress_array[stress_array_row + 1, :] = nodal_data2
			stress_array[stress_array_row + 2, :] = nodal_data3
			stress_array[stress_array_row + 3, :] = nodal_data4
			stress_array[stress_array_row + 4, :] = nodal_data5
			stress_array[stress_array_row + 5, :] = nodal_data6
			stress_array[stress_array_row + 6, :] = nodal_data7
			stress_array[stress_array_row + 7, :] = nodal_data8
			stress_array[stress_array_row + 8, :] = nodal_data9
			stress_array[stress_array_row + 9, :] = nodal_data10
			stress_array_row = stress_array_row + 10

		return stress_array
		
	def get_node_data_n(self, file_name):
		"""
		Reads mesh information. Returns a numpy arrays of nodes & elements in 'abaqus' format.
		"""
		meshSource=vtk.vtkUnstructuredGridReader()
		meshSource.SetFileName(file_name)
		meshSource.Update()

		#get nodes & elements returned to numpy arrays
		nread=v2n(meshSource.GetOutput().GetPoints().GetData())
		#allocate for extra node numbers to be input
		node_data=np.zeros((np.shape(nread)[0],np.shape(nread)[1]+1))
		#reshape according to 'abaqus' standards (node num, coord1, 2 ...)
		node_data[:,0]=np.arange(np.shape(nread)[0])+1
		node_data[:,1:]=nread
		# print([np.shape(np.arange(np.shape(n)[0])].transpose()))
		e=v2n(meshSource.GetOutput().GetCells().GetData())

		#get cell types & reshape element array as needed.
		tcs=vtk.vtkCellTypes()
		meshSource.GetOutput().GetCellTypes(tcs)

		#reshape according to 'abaqus' standards (elem number, connectivity1 2 ...)
		if tcs.IsType(12)==1:
			mainCellType=12 #1st order quad
			element_data=np.resize(e,(int(len(e)/float(9)),9))
			element_data[:,0]=np.arange(np.shape(element_data)[0])+1
			print('Quaaad') #debug
			self.ui.quadButton.setChecked(True)
		elif tcs.IsType(24)==1:
			mainCellType=24 #2nd order tet
			element_data=np.resize(e,(int(len(e)/float(11)),11))
			element_data[:,0]=np.arange(np.shape(element_data)[0])+1
			print('Tet tet tet') #debug
			self.ui.tetButton.setChecked(True)
		# print(np.shape(node_data),np.shape(element_data)) #debug
		return node_data, element_data+1 #add one to the element number to match abaqus format

	def get_node_data(self, file_name):
		"""
		Reads the nodal point coordinates. Returns a numpy array.
		"""

		# set parser keywords
		INP_FILE_NODE_LOOKUP_STR = "*Node"
		INP_FILE_ELEM_LOOKUP_STR = "*Element"
		INP_FILE_ELEM_END_LOOKUP_STR = "*Nset, nset=Part-1-1_SURFACE"

		# initialize
		curr_line = 0
		node_start = 0
		node_end = 0
		elem_start = 0
		elem_end = 0

		with open(file_name) as inp_file:
			p_lines = inp_file.readlines()
			for line in p_lines:
				curr_line = curr_line + 1
				if line.find(INP_FILE_NODE_LOOKUP_STR) >= 0:
					node_start = curr_line
				if line.find(INP_FILE_ELEM_LOOKUP_STR) >= 0:
					node_end = curr_line - 1
					elem_start = curr_line
				if line.find(INP_FILE_ELEM_END_LOOKUP_STR) >= 0:
					elem_end = curr_line - 1

		node_end = curr_line - node_end
		elem_end = curr_line - elem_end

		inp_file.close()

		# extract nodal point data for
		# node id, x coord, y coord, z coord, S33
		node_data = np.genfromtxt(file_name, skip_header=node_start, skip_footer=node_end, \
									delimiter=',')

		element_data = np.genfromtxt(file_name, skip_header=elem_start, skip_footer=elem_end, \
									delimiter=',')

		node_data = node_data.view().reshape(len(node_data), -1)
		element_data = element_data.view().reshape(len(element_data), -1)
		prin( node_data,element_data)
		return node_data, element_data
		
		
	def get_node_data_abq(self, file_name):
		"""
		Reads the nodal point coordinates. Returns a numpy array.
		"""
		fid = open(file_name)

		#flags for identifying sections of the inp file
		inpKeywords=["*Node", "*Element", "*Nset", "*Elset", "*NODE", "*ELEMENT", "*NSET", "*ELSET"]


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


		fid.close()
		#use genfromtxt to read between lines id'ed by lineFlag to pull in nodes and elements
		Nodes=np.genfromtxt(file_name,skip_header=lineFlag[0],skip_footer=i-lineFlag[1]+1,delimiter=",")
		Elements=np.genfromtxt(file_name,skip_header=lineFlag[1],skip_footer=i-lineFlag[2]+1,delimiter=",")
		node_data = Nodes.view().reshape(len(Nodes), -1)
		element_data = Elements.view().reshape(len(Elements), -1)
		return node_data, element_data

	def C3D8_quadrature_points(self):
		"""
		Define the natural coordinates of the quadrature points for C3D8.
		The element is a full integration brick with 8 nodes and 8 quadrature points.
		"""

		# create the square shape function matrix
		shape_function_matrix = np.zeros(shape=(8,8))

		# natural coordinates of the quadrature points
		nat_coord_quadrature_points = np.array([[-1/3**(0.5), -1/3**(0.5), -1/3**(0.5)], \
												[-1/3**(0.5), -1/3**(0.5), 1/3**(0.5)], \
												[-1/3**(0.5), 1/3**(0.5), -1/3**(0.5)], \
												[-1/3**(0.5), 1/3**(0.5), 1/3**(0.5)], \
												[1/3**(0.5), -1/3**(0.5), -1/3**(0.5)], \
												[1/3**(0.5), -1/3**(0.5), 1/3**(0.5)], \
												[1/3**(0.5), 1/3**(0.5), -1/3**(0.5)], \
												[1/3**(0.5), 1/3**(0.5), 1/3**(0.5)]])

		# apply the shape functions to the natural coordinates
		# and built the matrix
		for shape_matrix_index in range(0, 8):
			shape_function_matrix[shape_matrix_index, 0] = 0.125 * (1 - nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 1] = 0.125 * (1 + nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 2] = 0.125 * (1 + nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 3] = 0.125 * (1 - nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 4] = 0.125 * (1 - nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 5] = 0.125 * (1 + nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 - nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 6] = 0.125 * (1 + nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 2])
			shape_function_matrix[shape_matrix_index, 7] = 0.125 * (1 - nat_coord_quadrature_points[shape_matrix_index, 0]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 1]) \
																* (1 + nat_coord_quadrature_points[shape_matrix_index, 2])

		return shape_function_matrix

	def C3D10_quadrature_points(self):
		"""
		Define the natural coordinates of the quadrature points for C3D10.
		The element has 10 nodes and 4 quadrature points.
		"""

		# create the square shape function matrix
		shape_function_matrix = np.zeros(shape=(4, 10))

		# natural coordinates of the quadrature points
		nat_coord_quadrature_points = np.array([[(5+3*(5)**(0.5)), (5-(5)**(0.5)), (5-(5)**(0.5)), (5-(5)**(0.5))], \
												[(5-(5)**(0.5)), (5+3*(5)**(0.5)), (5-(5)**(0.5)), (5-(5)**(0.5))], \
												[(5-(5)**(0.5)), (5-(5)**(0.5)), (5+3*(5)**(0.5)), (5-(5)**(0.5))], \
												[(5-(5)**(0.5)), (5-(5)**(0.5)), (5-(5)**(0.5)), (5+3*(5)**(0.5))]])

		nat_coord_quadrature_points = nat_coord_quadrature_points / 20

		# apply the shape functions to the natural coordinates
		# and built the matrix
		for shape_matrix_index in range(0, 4):
			shape_function_matrix[shape_matrix_index, 0] = nat_coord_quadrature_points[shape_matrix_index, 0] \
														* (2 * nat_coord_quadrature_points[shape_matrix_index, 0] - 1)
			shape_function_matrix[shape_matrix_index, 1] = nat_coord_quadrature_points[shape_matrix_index, 1] \
														* (2 * nat_coord_quadrature_points[shape_matrix_index, 1] - 1)
			shape_function_matrix[shape_matrix_index, 2] = nat_coord_quadrature_points[shape_matrix_index, 2] \
														* (2 * nat_coord_quadrature_points[shape_matrix_index, 2] - 1)
			shape_function_matrix[shape_matrix_index, 3] = nat_coord_quadrature_points[shape_matrix_index, 3] \
														* (2 * nat_coord_quadrature_points[shape_matrix_index, 3] - 1)
			shape_function_matrix[shape_matrix_index, 4] = 4 * nat_coord_quadrature_points[shape_matrix_index, 0] \
														* nat_coord_quadrature_points[shape_matrix_index, 1]
			shape_function_matrix[shape_matrix_index, 5] = 4 * nat_coord_quadrature_points[shape_matrix_index, 1] \
														* nat_coord_quadrature_points[shape_matrix_index, 2]
			shape_function_matrix[shape_matrix_index, 6] = 4 * nat_coord_quadrature_points[shape_matrix_index, 2] \
														* nat_coord_quadrature_points[shape_matrix_index, 0]
			shape_function_matrix[shape_matrix_index, 7] = 4 * nat_coord_quadrature_points[shape_matrix_index, 0] \
														* nat_coord_quadrature_points[shape_matrix_index, 3]
			shape_function_matrix[shape_matrix_index, 8] = 4 * nat_coord_quadrature_points[shape_matrix_index, 1] \
														* nat_coord_quadrature_points[shape_matrix_index, 3]
			shape_function_matrix[shape_matrix_index, 9] = 4 * nat_coord_quadrature_points[shape_matrix_index, 2] \
														* nat_coord_quadrature_points[shape_matrix_index, 3]

		return shape_function_matrix

	def Keypress(self,obj, event):
		key = obj.GetKeySym()
		
		if key =="l": #load with keyboard shortcut
			self.get_input_data(None)
		if key =="1":
			xyview_post(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="2":
			yzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key =="3":
			xzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
		elif key == "Up":
			self.ren.GetActiveCamera().Roll(90)
			self.ren.ResetCamera()

		elif key=="i":
			im = vtk.vtkWindowToImageFilter()
			writer = vtk.vtkPNGWriter()
			im.SetInput(self.ui.vtkWidget._RenderWindow)
			im.Update()
			writer.SetInputConnection(im.GetOutputPort())
			writer.SetFileName("mesh.png")
			writer.Write()
			self.ui.statLabel.setText("Screen output saved to %s" %os.path.join(currentdir,'post.png'))
		
		elif key=="r":
			flip_visible(self.ax3D)
			
		elif key =="f": #flip color scheme for printing
			flip_colors(self.ren,self.ax3D)
				

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
	
	def draw_color_range(self, mesh_lookup_table):
		"""
		Draw the scalar range so that red is max, blue is min
		"""

		mesh_lookup_table.SetHueRange(0.667, 0)

	def draw_iso_surface(self):
		"""
		Draws contours on the object by the use of filters
		http://web.cs.wpi.edu/~matt/courses/cs563/talks/vtk/visualization.html
		"""

if __name__ == "__main__":
	post_process_tool()