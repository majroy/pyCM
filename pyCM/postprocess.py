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

import sys
import vtk
import pandas
import numpy as np
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtGui, QtWidgets
from pkg_resources import Requirement, resource_filename
from .pyCMcommon import *

__author__ = "N. Stoyanov"
__version__ = "0.1"
__email__ = "nikola.stoyanov@postgrad.manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2018"

DAT_FILE_LOOKUP_STR = "E L E M E N T   O U T P U T"
INP_FILE_NODE_LOOKUP_STR = "*Node"
INP_FILE_ELEM_LOOKUP_STR = "*Element, type=C3D8"
INP_FILE_ELEM_END_LOOKUP_STR = "*Nset, nset=Part-1-1_SURFACE, generate"

def post_process_tool():
    """
    Build QT interaction
    """
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    #spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
    #splash_pix = QtGui.QPixmap(spl_fname,'PNG')
    #splash = QtWidgets.QSplashScreen(splash_pix)
    #splash.setMask(splash_pix.mask())

    #splash.show()
    app.processEvents()

    window = MeshInteractor()
    window.show()
    #splash.finish(window)
    window.interactor_ui_box.Initialize()

    ret = app.exec_()

    if sys.stdin.isatty() and not hasattr(sys, 'ps1'):
        sys.exit(ret)
    else:
        return window


class QtMainWindow(object):
    """
    Setup and handle VTK and QT interaction
    """

    def __init__(self, main_window):
        """
        Performs the Qt setup of GUI

        Creates a windows and adds button interaction
        """

        main_window.setObjectName("Main Window")
        main_window.setWindowTitle("pyCM - FEA postprocessing v%s" %__version__)
        main_window.resize(1280, 720)

        self.central_widget = QtWidgets.QWidget(main_window)
        self.box_layout = QtWidgets.QHBoxLayout(self.central_widget)
        self.main_ui_box = QtWidgets.QFormLayout()

        self.vtk_widget = QVTKRenderWindowInteractor(self.central_widget)
        self.vtk_widget.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.vtk_widget.setMinimumSize(1150, 640)

        self.box_layout.addWidget(self.vtk_widget)
        self.box_layout.addStretch(1)
        main_window.setCentralWidget(self.central_widget)

        self.horiz_line1 = QtWidgets.QFrame()
        self.horiz_line1.setFrameStyle(QtWidgets.QFrame.HLine)
        self.outline_label = QtWidgets.QLabel("Outline modification")
        self.head_font=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)
        self.outline_label.setFont(self.head_font)

        self.read_sim_data_button = QtWidgets.QPushButton('Extract S33')
        self.display_vtk_button = QtWidgets.QPushButton('Display S33')

        self.horiz_line3 = QtWidgets.QFrame()
        self.horiz_line3.setFrameStyle(QtWidgets.QFrame.HLine)
        self.stress_label = QtWidgets.QLabel("Input min/max S33")
        self.stress_label.setWordWrap(True)
        self.stress_label.setFont(QtGui.QFont("Helvetica"))
        self.stress_label.setMinimumWidth(50)

        # enter minimum stress to update the display
        self.inp_min_stress = QtWidgets.QLabel("Min S33:")
        self.inp_min_stress = QtWidgets.QLineEdit()
        self.inp_min_stress.setMinimumWidth(50)

        # enter maximum stress to update the display
        self.inp_max_stress = QtWidgets.QLabel("Max S33:")
        self.inp_max_stress = QtWidgets.QLineEdit()
        self.inp_max_stress.setMinimumWidth(50)

        # stress update button
        self.updateButton = QtWidgets.QPushButton('Update')
        self.updateButton.setMinimumWidth(50)

        self.horiz_line2 = QtWidgets.QFrame()
        self.horiz_line2.setFrameStyle(QtWidgets.QFrame.HLine)
        self.stat_label = QtWidgets.QLabel("Idle")
        self.stat_label.setWordWrap(True)
        self.stat_label.setFont(QtGui.QFont("Helvetica", italic=True))
        self.stat_label.setMinimumWidth(50)

        self.main_ui_box.addRow(self.read_sim_data_button, self.display_vtk_button)
        self.main_ui_box.addRow(self.horiz_line3)
        self.main_ui_box.addRow(self.stress_label)
        self.main_ui_box.addRow(self.inp_min_stress)
        self.main_ui_box.addRow(self.inp_max_stress)
        self.main_ui_box.addRow(self.updateButton)
        self.main_ui_box.addRow(self.horiz_line2)
        self.main_ui_box.addRow(self.stat_label)
        self.box_layout.addLayout(self.main_ui_box)

class MeshInteractor(QtWidgets.QMainWindow):
    """
    Sets up the main VTK window
    reads .mat file and sets connection between UI and interactor
    """

    def __init__(self, parent = None):
        #read config file and setup code to be inserted here
        QtWidgets.QMainWindow.__init__(self, parent)
        self.main_ui_box = QtMainWindow(self)

        self.main_renderer = vtk.vtkRenderer()
        self.main_renderer.SetBackground(0.1, 0.2, 0.4)

        self.main_ui_box.vtk_widget.GetRenderWindow().AddRenderer(self.main_renderer)
        self.interactor_ui_box = self.main_ui_box.vtk_widget.GetRenderWindow().GetInteractor()

        self.style = vtk.vtkInteractorStyleTrackballCamera()
        self.style.AutoAdjustCameraClippingRangeOn()

        self.interactor_ui_box.SetInteractorStyle(self.style)
        self.main_renderer.GetActiveCamera().ParallelProjectionOn()
        self.renderer_position = self.main_renderer.GetActiveCamera().GetPosition()
        self.renderer_focal_point = self.main_renderer.GetActiveCamera().GetFocalPoint()

        self.mesh_actor = vtk.vtkActor()
        self.scalar_bar_actor = vtk.vtkScalarBarActor()

        self.vtk_file = None
        self.mesh_mapper = vtk.vtkDataSetMapper()

        self.main_ui_box.read_sim_data_button.clicked.connect(lambda: self.load_integration_points())
        self.main_ui_box.display_vtk_button.clicked.connect(lambda: self.load_vtk_legacy_file())
        self.main_ui_box.updateButton.clicked.connect(lambda: self.update_stress_display())

    def update_stress_display(self):
        """
        Updates the display to reflect the input of min/max stress
        """
        QtWidgets.QApplication.processEvents()

        min_stress = float(self.main_ui_box.inp_min_stress.text())
        max_stress = float(self.main_ui_box.inp_max_stress.text())

        self.mesh_mapper.SetScalarRange(min_stress, max_stress)
        self.main_ui_box.vtk_widget.update()
        self.stat_label = QtWidgets.QLabel("Updated!")

    def load_vtk_legacy_file(self):
        """
        Loads the vtk mesh and displays the scalar data in a color map.
        Allows further postprocessing to be done, such as grayscale and contour plots.
        """

        QtWidgets.QApplication.processEvents()

        if hasattr(self, "mesh_actor"):
            self.main_renderer.RemoveActor(self.mesh_actor)
            self.main_renderer.RemoveActor(self.scalar_bar_actor)

        # load vtk file
        self.vtk_file,_ = get_file("*.vtk")
        mesh_source = vtk.vtkUnstructuredGridReader()
        mesh_source.SetFileName(self.vtk_file)

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
        self.draw_color_range(mesh_lookup_table)

        # additional post processing work to be added
        # grayscale to amplify structural detail
        #self.draw_grayscale(mesh_lookup_table)
        # draw contours
        #self.draw_iso_surface()

        mesh_lookup_table.Build()

        scalar_range = mesh_reader_output.GetScalarRange()

        # mesh data set
        self.mesh_mapper.SetInputData(mesh_reader_output)
        self.mesh_mapper.SetScalarRange(scalar_range)
        self.mesh_mapper.SetLookupTable(mesh_lookup_table)

        # define scalar bar actor
        self.scalar_bar_actor.SetOrientationToVertical()
        self.scalar_bar_actor.SetLookupTable(mesh_lookup_table)

        #the scalar bar widget is associated with the qt interactor box
        scalar_bar_widget = vtk.vtkScalarBarWidget()
        scalar_bar_widget.SetInteractor(self.interactor_ui_box)
        scalar_bar_widget.SetScalarBarActor(self.scalar_bar_actor)
        scalar_bar_widget.On()

        #define the mesh actor properties
        self.mesh_actor.SetMapper(self.mesh_mapper)
        self.mesh_actor.GetProperty().SetLineWidth(1)
        self.mesh_actor.GetProperty().EdgeVisibilityOn()

        #display the actors
        self.main_renderer.AddActor(self.mesh_actor)
        self.main_renderer.AddActor(self.scalar_bar_actor)

        self.main_ui_box.vtk_widget.update()

    def load_integration_points(self):
        """
        Writes the FE output data - gaussian integration points
        to the legacy vtk file. The FE data is stored in the .dat file for ABAQUS.
        """
        QtWidgets.QApplication.processEvents()

        # load input files
        # use the original *.inp file, not *.abq.inp
        dat_file,_ = get_file("*.dat")
        inp_file,_ = get_file("*.inp")

        # default ABAQUS C3D8 elements only for now
        #quadrature_data = self.get_quadrature_data(dat_file)
        quadrature_data = read_abq_dat(dat_file)
        node_data, element_data = self.get_node_data(inp_file)

        # obtain a matrix of node number, x, y, z, stress
        stress_array = self.calculate_quadrature_stress_C3D8(quadrature_data, element_data, node_data)

        # nodes will duplicate in elements
        # sum contributions from nodes
        stress_data_frame = pandas.DataFrame(stress_array)
        stress_data_frame[5] = stress_data_frame.groupby([0])[4].transform('sum')

        # remove the individual stress at nodes
        del stress_data_frame[4]

        # drop the duplicates
        stress_data_frame = pandas.DataFrame(stress_data_frame).drop_duplicates(subset=[0])
        stress_data_frame = pandas.DataFrame(stress_data_frame).values
        stress_data_frame = pandas.DataFrame(stress_data_frame)

        # get *.vtk file
        vtk_file, _ = get_file("*.vtk")

        # append at the end of the vtk file
        self.append_postprocess_data(vtk_file, stress_data_frame)

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

            # create the square shape function matrix for C3D8
            shape_function_matrix = np.zeros(shape=(8,8))

            # obtain the natural coordinates of the gauss points
            C3D8_qp_natural_coord = self.C3D8_quadrature_points()

            for shape_matrix_index in range(0, 8):
                shape_function_matrix[shape_matrix_index, 0] = self.C3D8_shape_function1( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 1] = self.C3D8_shape_function2( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 2] = self.C3D8_shape_function3( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 3] = self.C3D8_shape_function4( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 4] = self.C3D8_shape_function5( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 5] = self.C3D8_shape_function6( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 6] = self.C3D8_shape_function7( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])
                shape_function_matrix[shape_matrix_index, 7] = self.C3D8_shape_function8( \
                                                                C3D8_qp_natural_coord[shape_matrix_index, :])

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

    def get_node_data(self, file_name):
        """
        Reads the nodal point coordinates. Returns a numpy array.
        """

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
        return node_data, element_data

    def C3D8_quadrature_points(self):
        """
        Define the natural coordinates of the quadrature points for C3D8.
        The element is a full integration brick with 8 nodes and 8 quadrature points.
        """

        # natural coordinates of the quadrature points
        nat_coord_quadrature_points = np.array([[-1/3**(0.5), -1/3**(0.5), -1/3**(0.5)], \
                                                [-1/3**(0.5), -1/3**(0.5), 1/3**(0.5)], \
                                                [-1/3**(0.5), 1/3**(0.5), -1/3**(0.5)], \
                                                [-1/3**(0.5), 1/3**(0.5), 1/3**(0.5)], \
                                                [1/3**(0.5), -1/3**(0.5), -1/3**(0.5)], \
                                                [1/3**(0.5), -1/3**(0.5), 1/3**(0.5)], \
                                                [1/3**(0.5), 1/3**(0.5), -1/3**(0.5)], \
                                                [1/3**(0.5), 1/3**(0.5), 1/3**(0.5)]])
        return nat_coord_quadrature_points

    def C3D8_shape_function1(self, coords):
        """
        Calculate the shape function for the first point
        """
        return 0.125 * (1 - coords[0]) * (1 - coords[1]) * (1 - coords[2])

    def C3D8_shape_function2(self, coords):
        """
        Calculate the shape function for the second point
        """
        return 0.125 * (1 + coords[0]) * (1 - coords[1]) * (1 - coords[2])

    def C3D8_shape_function3(self, coords):
        """
        Calculate the shape function for the third point
        """
        return 0.125 * (1 + coords[0]) * (1 + coords[1]) * (1 - coords[2])

    def C3D8_shape_function4(self, coords):
        """
        Calculate the shape function for the fourth point
        """
        return 0.125 * (1 - coords[0]) * (1 + coords[1]) * (1 - coords[2])
    def C3D8_shape_function5(self, coords):
        """
        Calculate the shape function for the fifth point
        """
        return 0.125 * (1 - coords[0]) * (1 - coords[1]) * (1 + coords[2])

    def C3D8_shape_function6(self, coords):
        """
        Calculate the shape function for the sixth point
        """
        return 0.125 * (1 + coords[0]) * (1 - coords[1]) * (1 + coords[2])

    def C3D8_shape_function7(self, coords):
        """
        Calculate the shape function for the seventh point
        """
        return 0.125 * (1 + coords[0]) * (1 + coords[1]) * (1 + coords[2])

    def C3D8_shape_function8(self, coords):
        """
        Calculate the shape function for the eight point
        """
        return 0.125 * (1 - coords[0]) * (1 + coords[1]) * (1 + coords[2])

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
