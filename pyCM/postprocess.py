#!/usr/bin/env python
'''
Uses VTK python to allow for postprocessing FEAs associated with the
contour method. Full interaction requires a 3-button mouse and keyboard.
1.0 - Heavily refactored for overall version 2.
'''

__author__ = "N. Stoyanov, M.J. Roy"
__version__ = "1.0"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, N. Stoyanov 2014--"

import numpy as np
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib import rc
from registration import read_file as reg_read_file
from pyCM.pyCMcommon import *

def launch(*args, **kwargs):
    '''
    Start Qt/VTK interaction if started independently
    '''
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)

    app.processEvents()
    
    window = interactor(None) #otherwise specify parent widget
    window.show()
    
    # if a data file is specified at launch
    if len(args) == 1:
        window.file = args[0]
        interactor.get_data(window)
    
    ret = app.exec_()
    
    if sys.stdin.isatty() and not hasattr(sys, 'ps1'):
        sys.exit(ret)
    else:
        return window

class post_main_window(object):
    '''
    Class that builds, populates and initializes aspects of Qt interaction
    '''

    def setup(self, MainWindow):
        '''
        Sets up Qt interactor
        '''
        #if called as a script, treat as a Qt mainwindow, otherwise a generic widget for embedding elsewhere
        if hasattr(MainWindow,'setCentralWidget'):
            MainWindow.setCentralWidget(self.centralWidget)
        else:
            self.centralWidget=MainWindow
            MainWindow.setWindowTitle("pyCM - Postprocessing tool v%s" %__version__)
        
        #create new layout to hold both VTK and Qt interactors
        mainlayout=QtWidgets.QHBoxLayout(self.centralWidget)

        #create VTK widget
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)

        #create VTK widget
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(100)
        sizePolicy.setVerticalStretch(100)
        self.vtkWidget.setSizePolicy(sizePolicy)
        
        self.vtkWidget.setMinimumSize(QtCore.QSize(800, 600))

        #make display layout
        display_box = QtWidgets.QGroupBox('Display')
        display_box_layout = QtWidgets.QGridLayout()
        display_box.setLayout(display_box_layout)

        #make combo box for components
        self.component_cb = QtWidgets.QComboBox()
        self.component_cb.setToolTip('Change stress component displayed')
        # self.component_cb.addItems(['\u03C311', '\u03C322', '\u03C333'])
        self.component_cb.setEnabled(False)
        self.mesh_display=QtWidgets.QPushButton("Edges off")
        self.mesh_display.setToolTip('Turn mesh/edges on and off')
        self.mesh_display.setCheckable(True)
        self.mesh_display.setChecked(False)
        self.mesh_display.setEnabled(False)
        self.draw_directions_button = QtWidgets.QPushButton('Show cut')
        self.draw_directions_button.setToolTip('Show cutting direction(s)')
        self.draw_directions_button.setEnabled(False)
        self.draw_directions_button.setCheckable(True)
        display_box_layout.addWidget(self.component_cb,0,0,1,1)
        display_box_layout.addWidget(self.mesh_display,0,1,1,1)
        display_box_layout.addWidget(self.draw_directions_button,0,2,1,1)

        #make contour layout
        contour_layout = QtWidgets.QGridLayout()
        contour_box = QtWidgets.QGroupBox('Contours')
        contour_box.setLayout(contour_layout)
        min_contour_label = QtWidgets.QLabel("Min:")
        self.min_contour = QtWidgets.QDoubleSpinBox()
        self.min_contour.setMinimum(-100000)
        self.min_contour.setMaximum(100000)
        max_contour_label = QtWidgets.QLabel("Max:")
        self.max_contour = QtWidgets.QDoubleSpinBox()
        self.max_contour.setMinimum(-100000)
        self.max_contour.setMaximum(100000)
        num_contour_label = QtWidgets.QLabel("Interval:")
        self.num_contour = QtWidgets.QSpinBox()
        self.num_contour.setToolTip('Number of entries shown on colorbar')
        self.num_contour.setMinimum(3)
        self.num_contour.setMaximum(20)
        self.num_contour.setValue(11)
        self.update_contours_button = QtWidgets.QPushButton('Update')
        self.update_contours_button.setToolTip('Update the contour limits and interval')
        self.update_contours_button.setEnabled(False)
        contour_layout.addWidget(min_contour_label,1,0,1,1)
        contour_layout.addWidget(self.min_contour,1,1,1,1)
        contour_layout.addWidget(max_contour_label,1,2,1,1)
        contour_layout.addWidget(self.max_contour,1,3,1,1)
        contour_layout.addWidget(num_contour_label,1,4,1,1)
        contour_layout.addWidget(self.num_contour,1,5,1,1)
        contour_layout.addWidget(self.update_contours_button,1,6,1,1)

        # line extraction from surface
        extract_layout = QtWidgets.QGridLayout()
        self.extract_box = QtWidgets.QGroupBox('Extract')
        # x, y, z of first point
        start_label = QtWidgets.QLabel("Start")
        start_label.setToolTip('Start point of line trace')
        self.point1_x_coord = QtWidgets.QDoubleSpinBox()
        self.point1_x_coord.setPrefix('x = ')
        self.point1_x_coord.setMinimum(-100000)
        self.point1_x_coord.setMaximum(100000)
        self.point1_y_coord = QtWidgets.QDoubleSpinBox()
        self.point1_y_coord.setPrefix('y = ')
        self.point1_y_coord.setMinimum(-100000)
        self.point1_y_coord.setMaximum(100000)

        # x, y, z of second point
        end_label = QtWidgets.QLabel("End")
        end_label.setToolTip('End point of line trace')
        self.point2_x_coord = QtWidgets.QDoubleSpinBox()
        self.point2_x_coord.setPrefix('x = ')
        self.point2_x_coord.setMinimum(-100000)
        self.point2_x_coord.setMaximum(100000)
        self.point2_y_coord = QtWidgets.QDoubleSpinBox()
        self.point2_y_coord.setPrefix('y = ')
        self.point2_y_coord.setMinimum(-100000)
        self.point2_y_coord.setMaximum(100000)

        #clip settings
        self.clip_active_button=QtWidgets.QPushButton("Clip")
        self.clip_active_button.setToolTip('Show/update clipped model')
        self.clip_active_button.setEnabled(False)

        interval_label=QtWidgets.QLabel("Line interval:")
        self.extract_interval=QtWidgets.QSpinBox()
        self.extract_interval.setToolTip('Number of points to extract along line trace')
        self.extract_interval.setValue(50)
        self.extract_interval.setMinimum(3)
        self.extract_interval.setMaximum(1000)

        self.extract_button = QtWidgets.QPushButton('Update line')
        self.extract_button.setToolTip('Show/update line trace')
        self.extract_button.setEnabled(False)
        self.export_line_status_label = QtWidgets.QLabel("Ready")
        self.export_line_status_label.setWordWrap(True)
        self.export_line_button = QtWidgets.QPushButton('Export line')
        self.export_line_button.setEnabled(False)
        self.export_line_button.setToolTip('Export line trace to file')



        #create figure canvas etc
        #initialize plot
        self.figure = plt.figure(figsize=(4,4))
        plt.rc('font', size = 9)
        plt.text(0.5, 0.5, "'Update' for plot", ha='center', style='italic', fontweight = 'bold', color='lightgray', size= 12)
        plt.axis('off')
        #changes the background of the plot, otherwise white
        # self.figure.patch.set_facecolor((242/255,242/255,242/255))
        self.canvas = FigureCanvas(self.figure)
        # self.canvas.setMinimumSize(QtCore.QSize(100, 200))
        self.canvas.setMinimumWidth(150)

        #populate extact_box
        extract_layout.addWidget(start_label,0,0,1,1)
        extract_layout.addWidget(self.point1_x_coord,0,1,1,1)
        extract_layout.addWidget(self.point1_y_coord,0,2,1,1)
        extract_layout.addWidget(end_label,1,0,1,1)
        extract_layout.addWidget(self.point2_x_coord,1,1,1,1)
        extract_layout.addWidget(self.point2_y_coord,1,2,1,1)
        extract_layout.addWidget(self.extract_interval,2,0,1,1)
        extract_layout.addWidget(self.clip_active_button,2,1,1,1)
        extract_layout.addWidget(self.extract_button,2,2,1,1)
        extract_layout.addWidget(self.canvas,3,0,1,3)
        
        evlayout=QtWidgets.QVBoxLayout()
        evbutton_layout = QtWidgets.QHBoxLayout()
        verticalSpacer = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        evbutton_layout.addWidget(self.export_line_status_label)
        evbutton_layout.addItem(verticalSpacer)
        evbutton_layout.addWidget(self.export_line_button)
        
        evlayout.addLayout(extract_layout)
        evlayout.addLayout(evbutton_layout)
        self.extract_box.setLayout(evlayout)

        lvlayout=QtWidgets.QVBoxLayout()
        lvlayout.addWidget(display_box)
        lvlayout.addWidget(contour_box)
        lvlayout.addWidget(self.extract_box)

        lvlayout.addStretch(1)
        
        mainlayout.addWidget(self.vtkWidget)
        mainlayout.addStretch(1)
        mainlayout.addLayout(lvlayout)

class interactor(QtWidgets.QWidget):
    '''
    Inherits most properties from a generic QWidget - see interactor docstring elsewhere in this package.
    '''
    
    def __init__(self,parent):
        super(interactor, self).__init__(parent)
        self.ui = post_main_window()
        self.ui.setup(self)
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(vtk.vtkNamedColors().GetColor3d("slategray"))
        self.ren.GradientBackgroundOn()

        self.file = None
        self.active_dir = os.getcwd()

        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style=vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        # self.iren.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ui.vtkWidget.Initialize()

        self.ui.mesh_display.clicked.connect(self.toggle_edges)
        self.ui.draw_directions_button.clicked.connect(self.draw_cut_directions)
        self.ui.clip_active_button.clicked.connect(self.clip)
        self.ui.update_contours_button.clicked.connect(self.update_scale_bar)

        self.ui.extract_button.clicked.connect(self.extract)
        self.ui.export_line_button.clicked.connect(self.export_line)

    def get_data(self):
        if self.file is None:
            self.file, self.active_dir = get_file("*.pyCM",self.active_dir)
        
        #make sure valid file was selected
        if self.file is None or not(os.path.isfile(self.file)):
            return
        
        #clear renderer
        self.ren.RemoveAllViewProps()
        
        self.cut_attr = reg_read_file(self.file)[-1]
        
        self.mesh = read_file_for_fea(self.file)[-1]
        if self.mesh is None:
            return
        
        #get component names
        components = []
        for index in range(self.mesh.GetPointData().GetNumberOfArrays()):
            components.append(self.mesh.GetPointData().GetArrayName(index))
        if not components: #then its a pre-processed mesh that hasn't been run with results
            return
        self.ui.component_cb.addItems(components)
        #find the index of S33
        self.ui.component_cb.setCurrentIndex(components.index("S33"))
        
        #call draw mesh
        self.draw_mesh()
        
        self.ui.component_cb.currentIndexChanged.connect(self.draw_mesh)
        
        self.ui.component_cb.setEnabled(True)
        self.ui.mesh_display.setEnabled(True)
        self.ui.draw_directions_button.setEnabled(True)
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
    
    def draw_mesh(self):
        
        # #Logic which bypasses the initial call from the combobox on loading
        # if self.ui.component_cb.currentText() != '':
            # self.component = self.ui.component_cb.currentText()
        self.component = self.ui.component_cb.currentText()
        
        #clear the renderer completely
        self.ren.RemoveAllViewProps()
        
        edges = True #default
        if self.ui.mesh_display.isChecked():
            edges = False
        
        self.mesh_actor, self.mesh_mapper, self.mesh_lut, mesh_range = generate_ug_actor(self.mesh, self.component, edges)
        self.ui.update_contours_button.setEnabled(True)
        self.ui.extract_button.setEnabled(True)
        self.ui.clip_active_button.setEnabled(True)
        self.ui.export_line_button.setEnabled(True)
    
        #update contour limits
        self.ui.min_contour.setValue(mesh_range[0])
        self.ui.max_contour.setValue(mesh_range[1])
    
        sb_widget = gen_scalar_bar(title = "MPa", num_contours = self.ui.num_contour.value(), side = 'left')
        sb_widget.SetInteractor(self.iren)
        sb_widget.On()
        self.sb_actor = sb_widget.GetScalarBarActor()
        self.sb_actor.SetLabelFormat("%.1f")
        self.sb_actor.SetLookupTable(self.mesh_lut)
        
        if hasattr(self,'axis_actor'):
            self.ren.RemoveActor(self.axis_actor)
        nl = self.mesh_actor.GetBounds()
        self.axis_actor = get_axis(self.ren, nl, 1, z = True)
        self.ren.AddActor(self.axis_actor)
        
        self.ren.AddActor(self.mesh_actor)
        self.ren.AddActor(self.sb_actor)
        
        if self.ui.draw_directions_button.isChecked():
            self.draw_cut_directions()
        
        self.ui.vtkWidget.update()
    
    def update_scale_bar(self):
        '''
        updates the active scale bar with limits and number of intervals from ui
        '''
        r = (self.ui.min_contour.value(),self.ui.max_contour.value())
        self.mesh_mapper.SetScalarRange(r[0], r[1])
        if hasattr(self,'clip_mapper'):
            self.clip_mapper.SetScalarRange(r[0], r[1])
        self.sb_actor.SetNumberOfLabels(self.ui.num_contour.value())
        self.ui.vtkWidget.update()

    def draw_cut_directions(self):
        '''
        Draws direction actors based on the contents of cut_attr. Called on load/display reset
        '''
        
        if not self.cut_attr['ref']:
            return
        
        cut_attr = self.cut_attr.copy()
        
        if hasattr(self,'ref_cut_orient_actor'):
            self.ren.RemoveActor(self.ref_cut_orient_actor)
        
        if not self.ui.draw_directions_button.isChecked():
            self.ui.vtkWidget.update()
            return
        
        a_limits = self.mesh.GetBounds()
        a_bump = np.maximum(a_limits[1]-a_limits[0],a_limits[3]-a_limits[2])*0.025
        limits = np.array([
        a_limits[0]-a_bump, a_limits[1]+a_bump,
        a_limits[2]-a_bump, a_limits[3]+a_bump])
        
        if cut_attr['ref']:
            ld = cut_attr['ref']
            if ld['cut_path'].any():
                self.ref_cut_orient_actor = gen_cutting_orientation_actor(\
                limits,\
                ld['cut_dir'],\
                ld['cut_path'])
            else:
                self.ref_cut_orient_actor = gen_cutting_orientation_actor(\
                limits,\
                ld['cut_dir'])

            self.ren.AddActor(self.ref_cut_orient_actor)

        self.ui.vtkWidget.update()

    def clip(self):
        '''
        Activates clipping by hiding mesh_actor and replacing it with a clipped actor based on the points set in the text box. Clipping plane is specified by the plane defined by 'Start','End' and 'Clip'.
        Clipping is removed by specifying zeros for the third point, and by virtue avoids a divide by zero error when calculating the clipping plane normal.
        '''
        
        if hasattr(self,'clipped_actor'):
            self.ren.RemoveActor(self.clipped_actor)
            
        #read points for plane
        p1 = np.array([self.ui.point1_x_coord.value(), self.ui.point1_y_coord.value(), 0])
        p2 = np.array([self.ui.point2_x_coord.value(), self.ui.point2_y_coord.value(), 0])
        p3 = np.array([(p1[0]+p2[0])/2,(p1[1]+p2[1])/2,np.sum(self.mesh_actor.GetBounds()[-2:])/2])
        
        c = p3 == np.zeros(3)
        if c.all() and self.mesh_actor.GetVisibility() == 0: #no clipping plane (p3 = 0,0,0) is specified & mesh is hidden
            flip_visible(self.mesh_actor)

        elif not c.all():
            clipPlane = vtk.vtkPlane()
            clipPlane.SetOrigin(((p1+p2)/2).tolist())
            #solve cross product between p1,p2 and p2,p3
            xnorm = np.cross((p2-p1),(p3-p2))
            xnorm = xnorm / np.sqrt(np.sum(xnorm**2))
            clipPlane.SetNormal(xnorm.tolist())
            
            clipper = vtk.vtkTableBasedClipDataSet() #needs to be table based, otherwise the grid is interpolated
            clipper.SetClipFunction(clipPlane)
            clipper.SetInputData(self.mesh) #needs to remain vtk object
            clipper.GenerateClippedOutputOn()
            clipper.Update()

            self.clip_mapper = vtk.vtkDataSetMapper()
            self.clip_mapper.SetInputData(clipper.GetClippedOutput())
            self.clip_mapper.SetLookupTable(self.mesh_lut)
            self.clip_mapper.SetScalarRange(self.mesh.GetScalarRange())

            self.clipped_actor = vtk.vtkActor()
            self.clipped_actor.SetMapper(self.clip_mapper)
            if self.ui.mesh_display.isChecked():
                self.clipped_actor.GetProperty().EdgeVisibilityOff()
            else:
                self.clipped_actor.GetProperty().EdgeVisibilityOn()
            if self.mesh_actor.GetVisibility() == 1:
                flip_visible(self.mesh_actor)
            self.ren.AddActor(self.clipped_actor)
            if hasattr(self,'axis_actor'):
                self.ren.RemoveActor(self.axis_actor)
            nl = self.clipped_actor.GetBounds()
            self.axis_actor = get_axis(self.ren, nl, 1, z = True)
            self.ren.AddActor(self.axis_actor)
        
        self.ui.vtkWidget.update()
    
    def extract(self):
        p1 = [self.ui.point1_x_coord.value(), self.ui.point1_y_coord.value(), 0]
        p2 = [self.ui.point2_x_coord.value(), self.ui.point2_y_coord.value(), 0]
        self.q = line_query(self.mesh,p1,p2,self.ui.extract_interval.value(),self.component)
        self.x = range(self.q.shape[0])
        self.ui.figure.clear()
        
        # print(self.x,self.q)
        ax = self.ui.figure.add_subplot(111)
        ax.scatter(self.x,self.q[:,-1])
        ax.set_ylabel("%s (MPa)"%self.component)
        ax.set_xlabel("Point number")
        ax.grid(visible=True, which='major', color='#666666', linestyle='-')
        ax.minorticks_on()
        ax.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        self.ui.figure.tight_layout()
        self.ui.canvas.draw()
        
        a_limits = self.mesh.GetBounds()
        rad = np.maximum(a_limits[1]-a_limits[0],a_limits[3]-a_limits[2])*0.005
        
        #remove any line actor currently present
        if hasattr(self,'line_actor'):
            self.ren.RemoveActor(self.line_actor)
            
        self.line_actor = gen_line_actor(p1,p2,50, None, rad)
        self.ren.AddActor(self.line_actor)
        self.ui.export_line_button.setEnabled(True)
        self.ui.vtkWidget.update()
        
    def export_line(self):
        """
        Collects data from ui, writes to a valid file
        """
        
        fileo, _ = get_save_file('*.csv')
        if fileo is None:
            return

        np.savetxt(fileo,
        np.column_stack((self.x,self.q)), 
        delimiter=',',
        header = "%s\nPoint number, x, y, z, %s (MPa)"%(self.file,self.component),
        fmt='%i, %.3f, %.3f, %.3f, %.3f')
        
        self.ui.export_line_status_label.setText('Exported to %s'%fileo)
    
    def toggle_edges(self):
        '''
        changes the visibility of edges on the active mesh_actor
        '''
        if hasattr(self,'mesh_actor'):
            if self.ui.mesh_display.isChecked():
                self.mesh_actor.GetProperty().EdgeVisibilityOff()
            else:
                self.mesh_actor.GetProperty().EdgeVisibilityOn()
        
        if hasattr(self,'clipped_actor'):
            if self.ui.mesh_display.isChecked():
                self.clipped_actor.GetProperty().EdgeVisibilityOff()
            else:
                self.clipped_actor.GetProperty().EdgeVisibilityOn()
        
        self.ui.vtkWidget.update()

    def gen_report(self):
        pass

    def keypress(self,obj,event):
        '''
        VTK Interactor specific keypress binding
        '''
        key = obj.GetKeyCode()
        if key == "1":
            xyview(self.ren)
        elif key == "2":
            yzview(self.ren)
        elif key == "3":
            xzview(self.ren)
        elif key == "Up":
            self.ren.GetActiveCamera().Roll(30)
        elif key == "Down":
            self.ren.GetActiveCamera().Roll(-30)
        elif key == "a":
            if hasattr(self,'axis_actor'):
                flip_visible(self.axis_actor)
        elif key == "l":
            self.file = None
            self.get_data()
            
        elif key=="i":
            if hasattr(self,'mesh'):
                im = vtk.vtkWindowToImageFilter()
                writer = vtk.vtkPNGWriter()
                colors = vtk.vtkNamedColors()
                self.ren.SetBackground(colors.GetColor3d("white"))
                self.ren.GradientBackgroundOff()
                
                self.sb_actor.GetTitleTextProperty().SetColor(0,0,0)
                self.sb_actor.GetLabelTextProperty().SetColor(0,0,0)
                self.axis_actor.GetTitleTextProperty(0).SetColor(0,0,0)
                self.axis_actor.GetLabelTextProperty(0).SetColor(0,0,0)
                self.axis_actor.GetXAxesLinesProperty().SetColor(0,0,0)
                self.axis_actor.GetTitleTextProperty(1).SetColor(0,0,0)
                self.axis_actor.GetLabelTextProperty(1).SetColor(0,0,0)
                self.axis_actor.GetYAxesLinesProperty().SetColor(0,0,0)
                self.axis_actor.GetTitleTextProperty(2).SetColor(0,0,0)
                self.axis_actor.GetLabelTextProperty(2).SetColor(0,0,0)
                self.axis_actor.GetZAxesLinesProperty().SetColor(0,0,0)
                
                im.SetInput(self.ui.vtkWidget._RenderWindow)
                im.Update()
                writer.SetInputConnection(im.GetOutputPort())
                writer.SetFileName("pyCM_capture.png")
                writer.Write()
                
                self.ren.SetBackground(colors.GetColor3d("slategray"))
                self.ren.GradientBackgroundOn()
                self.sb_actor.GetTitleTextProperty().SetColor(1,1,1)
                self.sb_actor.GetLabelTextProperty().SetColor(1,1,1)
                self.axis_actor.GetTitleTextProperty(0).SetColor(1,1,1)
                self.axis_actor.GetLabelTextProperty(0).SetColor(1,1,1)
                self.axis_actor.GetXAxesLinesProperty().SetColor(1,1,1)
                self.axis_actor.GetTitleTextProperty(1).SetColor(1,1,1)
                self.axis_actor.GetLabelTextProperty(1).SetColor(1,1,1)
                self.axis_actor.GetYAxesLinesProperty().SetColor(1,1,1)
                self.axis_actor.GetTitleTextProperty(2).SetColor(1,1,1)
                self.axis_actor.GetLabelTextProperty(2).SetColor(1,1,1)
                self.axis_actor.GetZAxesLinesProperty().SetColor(1,1,1)

        self.ui.vtkWidget.update()

def generate_ug_actor(ug, component, edges):
    '''
    Return an actor with a look up table and range of component selected. Edges are on if true.
    '''
    ug.GetPointData().SetActiveScalars(component)

    #build lookup table according to field
    lut = vtk.vtkLookupTable()
    lut.SetHueRange(0.667, 0)
    lut.Build()
    ug_range = ug.GetScalarRange()

    # map data set
    mesh_mapper = vtk.vtkDataSetMapper()
    mesh_mapper.SetInputData(ug)
    mesh_mapper.SetScalarRange(ug_range)
    mesh_mapper.SetLookupTable(lut)

    actor = vtk.vtkActor()
    actor.SetMapper(mesh_mapper)
    if edges:
        actor.GetProperty().EdgeVisibilityOn()
    else:
        actor.GetProperty().EdgeVisibilityOff()
    actor.GetProperty().SetLineWidth(0)

    return actor, mesh_mapper, lut, ug_range

def line_query(polydata, q1, q2, numPoints,component = None):
    """
    Interpolate the data in polydata over q1 to q2 (list of x,y,z)
    """
    query_point = [q1,q2]
    line = vtk.vtkLineSource()
    line.SetResolution(numPoints)
    line.SetPoint1(q1)
    line.SetPoint2(q2)
    line.Update()
    
    probe = vtk.vtkProbeFilter()
    probe.SetInputConnection(line.GetOutputPort())
    probe.SetSourceData(polydata)
    
    probe.Update() 
    
    #initialize numpy array - number of points in probe potentially != numPoints
    line_pts = np.empty((probe.GetOutput().GetNumberOfPoints(),3)) #x,y,z
    
    #get all points: could also iterate over probe.GetOutput().GetNumberOfPoints()
    for i in range(numPoints):
        line_pts[i,:] = probe.GetOutput().GetPoint(i)
    if component is not None:
        #stack probe value on end column of line_pts
        line_pts = np.hstack((line_pts, \
        np.array([v2n.vtk_to_numpy(probe.GetOutput().GetPointData().GetArray(component))]).T))
    
    return line_pts

if __name__ == "__main__":
    if len(sys.argv)>1:
        launch(sys.argv[1])
    else:
        launch()