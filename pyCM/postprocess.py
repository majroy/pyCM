"""
Uses VTK python to allow for postprocessing FEAs associated with the
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
ver 0.1 17-10-20
"""

import sys
import numpy as np
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtGui, QtWidgets
from pkg_resources import Requirement, resource_filename
from .pyCMcommon import *

__author__ = "M.J. Roy"
__version__ = "0.1"
__email__ = "matthew.roy@manchester.ac.uk"
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


class QtMainWindow:
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

        self.horiz_line2 = QtWidgets.QFrame()
        self.horiz_line2.setFrameStyle(QtWidgets.QFrame.HLine)
        self.stat_label = QtWidgets.QLabel("Idle")
        self.stat_label.setWordWrap(True)
        self.stat_label.setFont(QtGui.QFont("Helvetica", italic=True))
        self.stat_label.setMinimumWidth(50)

        self.main_ui_box.addRow(self.read_sim_data_button, self.display_vtk_button)
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

        #self.interactor_ui_box.AddObserver("KeyPressEvent", self.Keypress)

        # load from config file
        self.point_size = 2
        self.line_width = 1
        self.z_aspect = 1.0
        self.limits = np.empty(6)

        self.mesh_actor = vtk.vtkActor()
        self.scalar_bar_actor = vtk.vtkScalarBarActor()

        self.vtk_file = None

        #self.interactor_ui_box.AddObserver("KeyPressEvent", self.Keypress)

        self.main_ui_box.read_sim_data_button.clicked.connect(lambda: self.load_integration_points())
        self.main_ui_box.display_vtk_button.clicked.connect(lambda: self.load_vtk_legacy_file())

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
        self.vtk_file,_=get_file("*.vtk")
        mesh_source = vtk.vtkUnstructuredGridReader()
        mesh_source.SetFileName(self.vtk_file)

        # vtk will only read the first scalar
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

        # grayscale to amplify structural detail
        #self.draw_grayscale(mesh_lookup_table)

        # draw contours
        #self.draw_iso_surface()

        mesh_lookup_table.Build()

        scalar_range = mesh_reader_output.GetScalarRange()

        #mesh data set
        mesh_mapper = vtk.vtkDataSetMapper()
        mesh_mapper.SetInputData(mesh_reader_output)
        mesh_mapper.SetScalarRange(scalar_range)
        mesh_mapper.SetLookupTable(mesh_lookup_table)

        #define scalar bar actor
        self.scalar_bar_actor.SetOrientationToVertical()
        self.scalar_bar_actor.SetLookupTable(mesh_lookup_table)

        #the scalar bar widget is associated with the qt interactor box
        scalar_bar_widget = vtk.vtkScalarBarWidget()
        scalar_bar_widget.SetInteractor(self.interactor_ui_box)
        scalar_bar_widget.SetScalarBarActor(self.scalar_bar_actor)
        scalar_bar_widget.On()

        #define the mesh actor properties
        self.mesh_actor.SetMapper(mesh_mapper)
        self.mesh_actor.GetProperty().SetLineWidth(1)
        self.mesh_actor.GetProperty().SetColor(0, 0.9020, 0.9020)
        self.mesh_actor.GetProperty().SetEdgeColor([0.8, 0.8, 0.8])
        self.mesh_actor.GetProperty().EdgeVisibilityOn()

        #display the actors
        self.main_renderer.AddActor(self.mesh_actor)
        self.main_renderer.AddActor(self.scalar_bar_actor)

        #self.limits[4] = bounds[4]
        #self.limits[5] = bounds[5]

        #add axis code - causes a crash???!
        #self.AddAxis(self.limits, 1)

        self.main_ui_box.vtk_widget.update()

    def load_integration_points(self):
        """
        Writes the the FE output data - gaussian integration points
        to the legacy vtk file. The FE data is stored in the .dat file for ABAQUS.
        """
        QtWidgets.QApplication.processEvents()

        dat_file,_=get_file("*.dat")
        inp_file,_=get_file("*.inp")

        # we default to vtk element 12 only for now
        quadrature_data = self.get_quadrature_data(dat_file)
        node_data, element_data = self.get_node_data(inp_file)

        stress_array = self.calculate_quadrature_stress_C3D8(quadrature_data, element_data, node_data)

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

        #for row_index in range(0, len(quadrature_data), element_step):
        for row_index in range(0, 1, element_step):
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

            # calculate the nodal stresses
            nodal_stress = shape_function_matrix.dot(quadrature_stress)
        return 1
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

    def get_quadrature_data(self, file_name):
        """
        Reads the quadrature point coordinates and stress values. Returns a numpy array.
        """

        # initialize
        curr_line = 0
        read_flag = False
        feed_flag = False
        num_steps = 0
        row_start = 0
        row_end = 1

        # locate the start and end line of the quadrature block
        with open(file_name) as dat_file:
            p_lines = dat_file.readlines()
            for line in p_lines:
                num_steps = num_steps + 1
                if line.rstrip() and not feed_flag:
                    line = line.strip()
                    # find the element output section
                    if line.find(DAT_FILE_LOOKUP_STR) >= 0:
                        read_flag = True
                    # we can read the quadrature point data
                    # but we need to skip 3 lines
                    if read_flag:
                        curr_line = curr_line + 1
                        # coordinate and stress data reached
                        if curr_line == 5:
                            feed_flag = True
                            row_start = num_steps - 1
                            curr_line = 0
                else:
                    if not line.rstrip():
                        if read_flag is True and curr_line is 0:
                            # count the number of lines at the end of the .dat file
                            row_end = row_end + 1
        dat_file.close()

        # extract quadrature data for
        # node id, quadrature id, x coord, y coord, z coord, S33
        quadrature_data = np.genfromtxt(file_name, skip_header=row_start, skip_footer=11, \
                                    usecols=(0, 2, 3, 4, 7), autostrip=True)
        return quadrature_data

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

    def C3D8_nodal_points(self):
        """
        Define the natural coordinates for the nodal points for C3D8.
        """

        # natural coordinates of the nodal points
        nat_coord_nodal_points = np.array[[-1, -1, -1], \
                                        [1, -1, -1], \
                                        [1, 1, -1], \
                                        [-1, 1, -1], \
                                        [-1, -1, 1], \
                                        [1, -1, 1], \
                                        [1, 1, 1], \
                                        [-1, 1, 1]]
        return nat_coord_nodal_points
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

    def load_scalar_bar(self, vtk_mesh):
        """
        Load the field data in the renderer

        code must be moved to this method
        """

    def load_visualization(self):
        """
        Load the vtk file with the scalar data
        """

    def draw_color_range(self, mesh_lookup_table):
        """
        Draw the scalar range so that red is max, blue is min
        """

        mesh_lookup_table.SetHueRange(0.667, 0)

    def draw_grayscale(self, mesh_lookup_table):
        """
        Draw the scalar range grayscale - to amplify structural detail
        """

        mesh_lookup_table.SetHueRange(0, 0)
        mesh_lookup_table.SetSaturationRange(0, 0)
        mesh_lookup_table.SetValueRange(0.2, 1.0)

    def draw_iso_surface(self):
        """
        Draws contours on the object by the use of filters
        http://web.cs.wpi.edu/~matt/courses/cs563/talks/vtk/visualization.html
        """

    def AddAxis(self, limits, scale):
        """
        Add axis
        """

        if hasattr(self, "ax3D"):
            self.render.RemoveActor(self.ax3D)

        self.ax3D = vtk.vtkCubeAxesActor()
        self.ax3D.ZAxisTickVisibilityOn()

        self.ax3D.SetXTitle('X')
        self.ax3D.SetYTitle('Y')
        self.ax3D.SetZTitle('Z')

        #self.ax3D.SetXUnits('mm')
        #self.ax3D.SetYUnits('mm')
        #self.ax3D.SetZUnits('mm')

        self.ax3D.SetBounds(limits)
        self.ax3D.SetZAxisRange(limits[-2] * scale, limits[-1] * scale)
        self.ax3D.SetCamera(self.main_renderer.GetActiveCamera())
        self.main_renderer.AddActor(self.ax3D)

        self.ax3D.SetFlyModeToOuterEdges()

if __name__ == "__main__":
    post_process_tool()
