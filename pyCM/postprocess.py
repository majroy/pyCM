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
from vtk.Qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtGui
from PyQt5.QtGui import QApplication
import pyCMcommon

__author__ = "M.J. Roy"
__version__ = "0.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2017"


def post_process_tool():
    """
    Build QT interaction
    """
    app = QtGui.QApplication.instance()

    if app is None:
        app = QApplication(sys.argv)

    app.processEvents()
    window = MeshInteractor()
    window.show()
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

        self.central_widget = QtGui.QWidget(main_window)
        self.box_layout = QtGui.QHBoxLayout(self.central_widget)
        self.main_ui_box = QtGui.QFormLayout()

        self.vtk_widget = QVTKRenderWindowInteractor(self.central_widget)
        self.vtk_widget.setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.vtk_widget.setMinimumSize(1150, 640)

        self.box_layout.addWidget(self.vtk_widget)
        self.box_layout.addStretch(1)
        main_window.setCentralWidget(self.central_widget)

        self.horiz_line1 = QtGui.QFrame()
        self.horiz_line1.setFrameStyle(QtGui.QFrame.HLine)
        self.outline_label = QtGui.QLabel("Outline modification")
        self.head_font = QtGui.QFont("Helvetica [Cronyx]", weight=QtGui.QFont.Bold)
        self.outline_label.setFont(self.head_font)

        self.read_mesh_button = QtGui.QPushButton('Read Mesh')
        self.read_data_button = QtGui.QPushButton('Read Dat')

        self.main_ui_box.addRow(self.read_mesh_button, self.read_data_button)
        self.box_layout.addLayout(self.main_ui_box)

class MeshInteractor(QtGui.QMainWindow):
    """
    Sets up the main VTK window
    reads .mat file and sets connection between UI and interactor
    """

    def __init__(self, parent = None):
        #read config file and setup code to be inserted here
        QtGui.QMainWindow.__init__(self, parent)
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

        self.main_ui_box.read_mesh_button.clicked.connect(lambda: self.load_vtk_legacy_file())
        self.main_ui_box.read_data_button.clicked.connect(lambda: self.load_integration_points())
        self.main_ui_box.read_data_button.clicked.connect(lambda: self.load_visualization())

    def load_vtk_xml_file(self, vtkFileAcquired = False):
        """
        Loads the vtk mesh and displays the scalar data in a color map.
        Allows further postprocessing to be done, such as grayscale and contour plots.
        """

        # to remove existing actors
        if hasattr(self, "mesh_actor"):
            self.main_renderer.RemoveActor(self.mesh_actor)
            self.main_renderer.RemoveActor(self.scalar_bar_actor)

        if vtkFileAcquired == False:
            self.vtk_file, startdir = pyCMcommon.get_file('*.vtk')

        # only xml files in vtk support integration points
        # no more legacy
        mesh_source = vtk.vtkXMLUnstructuredGridReader()

        mesh_source.SetFileName(self.vtk_file)

        # vtk will only read the first scalar. If there is more than one we need to specify it
        # in case the field is empty vtk will not display a scalar field
        # this is convevient since we can call the same function after adding the scalar data
        #mesh_source.SetScalarsName("T")
        #mesh_source.Update()
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

    def load_vtk_legacy_file(self, vtkFileAcquired = False):
        """
        Loads the vtk mesh and displays the scalar data in a color map.
        Allows further postprocessing to be done, such as grayscale and contour plots.
        """
        #ConvertInptoVTK("")
        # to remove existing actors
        if hasattr(self, "mesh_actor"):
            self.main_renderer.RemoveActor(self.mesh_actor)
            self.main_renderer.RemoveActor(self.scalar_bar_actor)

        if vtkFileAcquired == False:
            self.vtk_file, startdir = pyCMcommon.get_file('*.vtk')

        mesh_source = vtk.vtkUnstructuredGridReader()

        mesh_source.SetFileName(self.vtk_file)

        # vtk will only read the first scalar. If there is more than one we need to specify it
        # in case the field is empty vtk will not display a scalar field
        # this is convevient since we can call the same function after adding the scalar data
        mesh_source.SetScalarsName("Stress2")
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

        The following steps are performed
         1. For 'n' 0 to MAX Elements;
         2. Get x, y, z of each point on element 'n' in .vtk;
         3. Locate all integration point coordinates that are in the element boundary in .dat;
         4. Average S33 to obtain an elemental stress; - THIS IS WRONG (currently for testing)
         5. Print averaged S33 to CELL_DATA table in .vtk file; - THIS IS WRONG (currently for testing)
         6. 'n'++ and go to step 1.;

         Note: steps 4 and 5 are just exploratory. They need to be rewritten to:
         4. Use element shape functions to obtain S33 at nodes from integration points
            according to the ABAQUS numbering system;
         5. Average all node S33 from neighbouring cells and print to Point data in .vtk file;

         - Pandas is used for the data processing. This will allow statistical manipulations to be made
           with ease if so desired in the future.
         - MMap is used to inspect the files and locate the start/end of element scope (.vtk) and integration points (.dat).
           IF rewriting to Python 3 MMap works in a different way!!! Care must be taken!
        """

        # define start/end positions for points and cells for clarity in the code
        vtk_legacy_file_point_start = 6
        vtk_legacy_file_point_end = None
        vtk_legacy_file_cell_start = None
        vtk_legacy_file_cell_end = None

        # define start/end positions for FE (only ABAQUS for now) integration points for clarity
        fe_file_integration_point_start = None
        fe_file_integration_point_end = None

        # read vtk legacy file
        vtk_file_data = open(self.vtk_file)
        mmap_vtk_file_data = mmap.mmap(vtk_file_data.fileno(), 0, access = mmap.ACCESS_READ)

        # determine the start and end positions of the points and cells
        print (vtk_file_data)
        print (mmap_vtk_file_data.find("CELLS"))


        #vtk_element

        # add the scalar data and call the load_mesh
        self.load_vtk_legacy_file(vtkFileAcquired = True)

    def ConvertInptoVTK(infile, outfile):
        """
        Converts abaqus inp file into a legacy ASCII vtk file. First order quads (C3D8) and third order tets (C3D10) are supported.
        """
        fid = open(infile)

        #flags for identifying sections of the inp file
        inpKeywords = ["*Node", "*Element", "*Nset", "*Elset"]

        #map abaqus mesh types to vtk objects
        vtkType = {}
        vtkType['C3D8'] = 12
        vtkType['C3D10'] = 24

        #create counter for all lines in the inp file, and array to store their location
        i = 0
        lineFlag = [];

        #read file and find both where keywords occur as well as the element type used
        while 1:
            lines = fid.readlines(100000)
            if not lines:
                break
            for line in lines:
                i += 1
                for keyword in inpKeywords:
                    if line[0:len(keyword)] == keyword:
                        lineFlag.append(i)
                        if keyword == "*Element":
                            line = line.replace("\n", "")
                            CellNum = vtkType[line.split("=")[-1]]
        fid.close()

        #use genfromtxt to read between lines id'ed by lineFlag to pull in nodes and elements
        Nodes = np.genfromtxt(infile, skip_header = lineFlag[0], skip_footer = i - lineFlag[1] + 1, delimiter = ",")
        Elements = np.genfromtxt(infile, skip_header = lineFlag[1], skip_footer = i - lineFlag[2] + 1, delimiter = ",")

        #Now write it in VTK format to a new file starting with header
        fid = open(outfile,'w+')
        fid.write('# vtk DataFile Version 2.0\n')
        fid.write('%s,created by pyCM\n' %outfile[:-4])
        fid.write('ASCII\n')
        fid.write('DATASET UNSTRUCTURED_GRID\n')
        fid.write('POINTS %i double\n' %len(Nodes))

        #dump nodes
        np.savetxt(fid, Nodes[:,1::], fmt='%.6f')
        fid.write('\n')
        fid.write('CELLS %i %i\n'%(len(Elements), len(Elements) * len(Elements[0,:])))
        #Now elements, stack the number of nodes in the element instead of the element number
        Cells = np.hstack((np.ones([len(Elements[:, 0]), 1]) * len(Elements[0, 1 :: ]), Elements[:, 1 :: ] - 1))
        np.savetxt(fid,Cells, fmt = '%i')
        fid.write('\n')

        #Write cell types
        fid.write('CELL_TYPES %i\n' % len(Elements))
        CellType = np.ones([len(Elements[ : , 0]), 1]) * CellNum
        np.savetxt(fid, CellType, fmt = '%i')

        fid.close()

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
