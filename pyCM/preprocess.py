#!/usr/bin/env python
'''
Uses VTK Python to allow for fitting an averaged dataset associated with the
 contour method. Full interaction requires a 3-button mouse and keyboard.
1.7 - Updated for overall version 2.
'''

__author__ = "M.J. Roy"
__version__ = "1.7"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014--"

import os,io,sys,yaml
import subprocess as sp
from pkg_resources import Requirement, resource_filename
import numpy as np
from scipy.interpolate import bisplev
import vtk
import vtk.util.numpy_support as v2n
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from pyCM.pyCMcommon import *
from pyCM.icp import *
from pyCM.registration import read_file as reg_read_file
from pyCM.fit_surface import read_file as fs_read_file
from pyCM.extrude_widget import extrude_widget
from pyCM.fea_widget import fea_widget

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

class pre_main_window(object):
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
            MainWindow.setWindowTitle("pyCM - FEA preprocessing v%s" %__version__)
        
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

        #make outline box
        self.outline_box = QtWidgets.QGroupBox('Outline modification')
        outline_layout = QtWidgets.QGridLayout()
        self.outline_box.setLayout(outline_layout)
        
        #make active profile display
        entry_spec_layout = QtWidgets.QHBoxLayout()
        self.entry_spec = QtWidgets.QComboBox()
        self.entry_spec.setToolTip('Set focus on this entry')
        self.caption_rb = QtWidgets.QCheckBox('Caption entries')
        self.caption_rb.setToolTip('Label entries in viewport')
        self.caption_rb.setChecked(False)
        self.show_corners_rb = QtWidgets.QCheckBox('Caption corners')
        self.show_corners_rb.setToolTip("Numbers corners detected in outlines")
        self.show_corners_rb.setChecked(True)
        entry_spec_layout.addWidget(self.entry_spec)
        entry_spec_layout.addWidget(self.caption_rb)
        entry_spec_layout.addWidget(self.show_corners_rb)
        
        self.spacing_rb=QtWidgets.QRadioButton("Spacing")
        self.quantity_rb=QtWidgets.QRadioButton("Quantity")
        self.quantity_rb.setChecked(True)
        self.outline_rb_group = QtWidgets.QButtonGroup()
        self.outline_rb_group.addButton(self.spacing_rb)
        self.outline_rb_group.addButton(self.quantity_rb)
        self.outline_rb_group.setExclusive(True)
        self.seed_by_length_sb = QtWidgets.QDoubleSpinBox()
        self.seed_by_length_sb.setSuffix(' mm')
        self.seed_by_length_sb.setToolTip('Average line segment length describing active outline')
        self.seed_by_length_sb.setMinimum(0.001)
        self.seed_by_length_sb.setMaximum(1000)
        self.seed_by_length_sb.setDecimals(3)
        self.seed_by_num_sb = QtWidgets.QSpinBox()
        self.seed_by_num_sb.setPrefix('N = ')
        self.seed_by_num_sb.setMinimum(4)
        self.seed_by_num_sb.setMaximum(10000)
        self.seed_by_num_sb.setValue(100)
        self.seed_by_num_sb.setToolTip('Number of points on active outline')
        self.update_outline_button = QtWidgets.QPushButton('Update')
        self.update_outline_button.setToolTip('Update spacing on active outline entry')
        self.reset_outline_button = QtWidgets.QPushButton('Reset')
        self.reset_outline_button.setToolTip('Reset outline spacing')
        self.export_outline_button = QtWidgets.QPushButton('Export')
        self.export_outline_button.setToolTip('Export active outline to .dxf file')
        
        #populate outline box
        outline_layout.addLayout(entry_spec_layout,0,0,1,2)
        outline_layout.addWidget(self.spacing_rb, 1, 0, 1, 1)
        outline_layout.addWidget(self.quantity_rb, 1, 1, 1, 1)
        outline_layout.addWidget(self.seed_by_length_sb,2,0,1,1)
        outline_layout.addWidget(self.seed_by_num_sb,2,1,1,1)
        outline_button_layout = QtWidgets.QHBoxLayout()
        outline_button_layout.addWidget(self.update_outline_button)
        outline_button_layout.addWidget(self.reset_outline_button)
        outline_button_layout.addWidget(self.export_outline_button)
        outline_layout.addLayout(outline_button_layout,3,0,1,2)
        
        self.outline_box.setEnabled(False)
        
        #mesh interaction tools
        self.mesh_box = QtWidgets.QGroupBox('Mesh')
        mesh_layout = QtWidgets.QGridLayout()
        self.mesh_box.setLayout(mesh_layout)
        
        self.extrude_mesh_button = QtWidgets.QPushButton('Extrude')
        self.extrude_mesh_button.setToolTip('Generate extruded mesh')
        self.import_mesh_button = QtWidgets.QPushButton('Import')
        self.import_mesh_button.setToolTip('Import *.vtk mesh file')

        
        trans_x_label=QtWidgets.QLabel("Translate x:")
        self.trans_x = QtWidgets.QDoubleSpinBox()
        self.trans_x.setValue(0)
        self.trans_x.setMaximum(300)
        self.trans_x.setMinimum(-300)
        trans_y_label=QtWidgets.QLabel("Translate y:")
        self.trans_y = QtWidgets.QDoubleSpinBox()
        self.trans_y.setValue(0)
        self.trans_y.setMaximum(300)
        self.trans_y.setMinimum(-300)
        rotate_label = QtWidgets.QLabel("Rotation about z:")
        self.rotate_z= QtWidgets.QDoubleSpinBox()
        self.rotate_z.setToolTip('Positive is clockwise')
        self.rotate_z.setValue(0)
        self.rotate_z.setSuffix("\u00b0")
        self.rotate_z.setMaximum(180)
        self.rotate_z.setMinimum(-180)
        self.move_button = QtWidgets.QPushButton("Apply")
        self.move_button.setEnabled(False)
        self.move_button.setSizePolicy(sizePolicy)
        
        trans_algo_layout = QtWidgets.QHBoxLayout()
        self.mirror_button = QtWidgets.QPushButton("Mirror")
        self.mirror_button.setToolTip('Mirror mesh on z = 0 plane')
        self.mirror_button.setEnabled(False)
        trans_algo_layout.addWidget(self.mirror_button)

        node_selection_label = QtWidgets.QLabel('Node selection:')
        self.active_picking_indicator = QtWidgets.QLabel('Active')
        self.active_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
        self.active_picking_indicator.setAlignment(QtCore.Qt.AlignCenter)
        # self.active_picking_indicator.setFixedSize(75, 20)
        self.active_picking_indicator.setToolTip('Press N with interactor in focus to activate/deactivate node selection. When enabled, press A to accept selections and D to deselect.')
        self.reset_node_selection_button = QtWidgets.QPushButton("Reset")
        self.reset_node_selection_button.setToolTip("Reset all selected nodes")
        self.reset_node_selection_button.setEnabled(False)
        
        align_algo_layout = QtWidgets.QHBoxLayout()
        self.corner_best_fit_button = QtWidgets.QPushButton("Corner SVD")
        self.corner_best_fit_button.setToolTip('Solve and apply the best-fit 3D transform to outline corners with single value decomposition')
        self.corner_best_fit_button.setEnabled(False)
        self.vtk_icp_button = QtWidgets.QPushButton("VTK ICP")
        self.vtk_icp_button.setToolTip('Generate and apply a VTK 3D transformation based on corners')
        self.vtk_icp_button.setEnabled(False)
        align_algo_layout.addWidget(self.corner_best_fit_button)
        align_algo_layout.addWidget(self.vtk_icp_button)
        
        mesh_layout.addWidget(self.extrude_mesh_button,0,0,1,1)
        mesh_layout.addWidget(self.import_mesh_button,0,1,1,1)
        mesh_layout.addWidget(self.mirror_button,0,2,1,1)
        mesh_layout.addWidget(trans_x_label,1,0,1,1)
        mesh_layout.addWidget(self.trans_x,1,1,1,1)
        mesh_layout.addWidget(trans_y_label,2,0,1,1)
        mesh_layout.addWidget(self.trans_y,2,1,1,1)
        mesh_layout.addWidget(rotate_label,3,0,1,1)
        mesh_layout.addWidget(self.rotate_z,3,1,1,1)
        mesh_layout.addWidget(self.move_button,1,2,3,1)
        mesh_layout.addLayout(trans_algo_layout,4,0,1,3)
        mesh_layout.addWidget(node_selection_label,4,0,1,1)
        mesh_layout.addWidget(self.active_picking_indicator,4,1,1,1)
        mesh_layout.addWidget(self.reset_node_selection_button,4,2,1,1)
        mesh_layout.addLayout(align_algo_layout,5,0,1,3)
        mesh_layout.setRowStretch(mesh_layout.rowCount(), 1)
        
        self.mesh_box.setEnabled(False)
        
        #BC interaction tools
        self.bc_box = QtWidgets.QGroupBox('Impose boundary conditions and material')
        bc_layout = QtWidgets.QGridLayout()
        self.bc_box.setLayout(bc_layout)
        bc_selection_label = QtWidgets.QLabel('Rigid body BCs:')
        self.bc_picking_indicator = QtWidgets.QLabel('Active')
        self.bc_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
        self.bc_picking_indicator.setAlignment(QtCore.Qt.AlignCenter)
        # self.bc_picking_indicator.setFixedSize(75, 20)
        self.bc_picking_indicator.setToolTip('Press B with interactor in focus to activate/deactivate node selection. When enabled, press A to accept selections and D to deselect.')
        self.reset_bc_selection_button = QtWidgets.QPushButton("Reset")
        self.reset_bc_selection_button.setToolTip("Reset all rigid body boundary conditions")

        impose_fit_layout = QtWidgets.QHBoxLayout()
        self.impose_fit_button = QtWidgets.QPushButton("Impose surface BCs")
        self.impose_fit_button.setToolTip('Impose boundary conditions from fitted surface(s).')
        self.z_aspect_sb = QtWidgets.QSpinBox()
        self.z_aspect_sb.setPrefix("Scaling: ")
        self.z_aspect_sb.setToolTip("Scaling factor applied to Z axis")
        self.z_aspect_sb.setValue(1)
        self.z_aspect_sb.setMinimum(1)
        self.z_aspect_sb.setMaximum(1000)
        impose_fit_layout.addWidget(self.impose_fit_button)
        impose_fit_layout.addWidget(self.z_aspect_sb)
        
        material_props_label = QtWidgets.QLabel('Material properties:')
        self.modulus_sb = QtWidgets.QDoubleSpinBox()
        self.modulus_sb.setPrefix("E = ")
        self.modulus_sb.setSuffix(" MPa")
        self.modulus_sb.setMinimum(1)
        self.modulus_sb.setDecimals(0)
        self.modulus_sb.setMaximum(1000000)
        self.modulus_sb.setValue(200000)
        
        self.poisson_sb = QtWidgets.QDoubleSpinBox()
        self.poisson_sb.setPrefix("\u03bd = ")
        self.poisson_sb.setMinimum(0.1)
        self.poisson_sb.setMaximum(0.5)
        self.poisson_sb.setValue(0.30)
        self.poisson_sb.setDecimals(3)
        
        bc_layout.addWidget(bc_selection_label,0,0,1,1)
        bc_layout.addWidget(self.bc_picking_indicator,0,1,1,1)
        bc_layout.addWidget(self.reset_bc_selection_button,0,2,1,1)
        bc_layout.addLayout(impose_fit_layout,1,0,1,3)
        bc_layout.addWidget(material_props_label,2,0,1,1)
        bc_layout.addWidget(self.modulus_sb,2,1,1,1)
        bc_layout.addWidget(self.poisson_sb,2,2,1,1)
        self.bc_box.setEnabled(False)
        
        #make save box
        self.save_box = QtWidgets.QGroupBox('Write current')
        save_layout = QtWidgets.QHBoxLayout()
        self.save_box.setLayout(save_layout)
        #make save buttons
        self.save_button = QtWidgets.QPushButton("Save")
        self.save_button.setToolTip('Save data to hdf5-formatted results file')
        self.save_label = QtWidgets.QLabel('Ready')
        self.save_label.setWordWrap(True)
        self.run_FEA_rb = QtWidgets.QCheckBox('Run FEA on save')
        self.run_FEA_rb.setToolTip("Launch FEA widget after saving to the data file.")
        self.run_FEA_rb.setChecked(True)
        save_layout.addWidget(self.save_label)
        verticalSpacer = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        save_layout.addItem(verticalSpacer)
        save_layout.addWidget(self.run_FEA_rb)
        save_layout.addWidget(self.save_button)
        self.save_box.setEnabled(False)
        
        lvlayout=QtWidgets.QVBoxLayout()
        lvlayout.addWidget(self.outline_box)
        lvlayout.addWidget(self.mesh_box)
        lvlayout.addWidget(self.bc_box)
        lvlayout.addWidget(self.save_box)
        lvlayout.minimumSize()
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
        self.ui = pre_main_window()
        self.ui.setup(self)
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(vtk.vtkNamedColors().GetColor3d("slategray"))
        self.ren.GradientBackgroundOn()

        self.file = None
        self.point_size = 2
        self.active_dir = os.getcwd()
        self.picking = False
        self.bc_picking = False
        self.trans = np.eye(4) #default
        self.c_trans = np.eye(4) #initialise cumulative transformation
        
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style=vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.iren.AddObserver("KeyPressEvent", self.keypress)

        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ui.vtkWidget.Initialize()

        self.ui.entry_spec.activated.connect(self.focus_outline_entry)
        self.ui.caption_rb.toggled.connect(self.hide_caption)
        self.ui.show_corners_rb.toggled.connect(self.show_corners)
        self.ui.update_outline_button.clicked.connect(self.update_outlines)
        self.ui.reset_outline_button.clicked.connect(self.reset_outlines)
        self.ui.export_outline_button.clicked.connect(self.export_outline)
        self.ui.extrude_mesh_button.clicked.connect(self.generate_extrusion)
        self.ui.import_mesh_button.clicked.connect(self.import_mesh)
        self.ui.mirror_button.clicked.connect(self.mirror_mesh)
        self.ui.move_button.clicked.connect(self.trans_from_ui)
        self.ui.reset_node_selection_button.clicked.connect(self.reset_node_pick)
        self.ui.corner_best_fit_button.clicked.connect(self.trans_from_corners)
        self.ui.vtk_icp_button.clicked.connect(self.trans_from_vtk_icp)
        self.ui.reset_bc_selection_button.clicked.connect(self.reset_bc_pick)
        self.ui.impose_fit_button.clicked.connect(self.impose_surface_bc)
        self.ui.z_aspect_sb.valueChanged.connect(self.update_z_aspect)
        self.ui.save_button.clicked.connect(self.write_output)

    def make_orientation_widget(self):

        axes = vtk.vtkAxesActor()
        self.orientation_widget = vtk.vtkOrientationMarkerWidget()
        self.orientation_widget.SetOutlineColor(1,1,1)
        self.orientation_widget.SetOrientationMarker(axes)
        self.orientation_widget.SetInteractor(self.iren)
        self.orientation_widget.SetViewport(0.0, 0.0, 0.2, 0.2)

        self.orientation_widget.EnabledOn()
        self.orientation_widget.InteractiveOn()

    def get_data(self):
        if self.file is None:
            self.file, self.active_dir = get_file("*.pyCM",self.active_dir)
        
        #make sure valid file was selected
        if self.file is None or not(os.path.isfile(self.file)):
            return
        
        #clear the renderer completely
        self.ren.RemoveAllViewProps()

        #try reading data fit surface
        _, _, _, _, self.bv_tck_list = fs_read_file(self.file)
        _, self.outlines, _, _, _, = reg_read_file(self.file)
        #populate outline box
        self.ui.entry_spec.clear()
        for i in range(len(self.outlines)):
            self.ui.entry_spec.insertItem(i,'Entry %d'%i)
        
        if not self.bv_tck_list or None in self.bv_tck_list:
            #if there aren't fitted surfaces for all outlines
            return
        
        #try reading data from this step
        these_outlines,\
        _,\
        _,\
        these_bc_nodes,\
        trans,\
        mod, pr,\
        self.mesh = read_file_for_fea(self.file)

        #initialize rs_outlines
        self.rs_outline_corners = [None]*len(self.outlines)
        if not these_outlines:
            self.rs_outlines = self.outlines.copy()
        else:
            self.rs_outlines = these_outlines
        
        #test other incoming data
        if trans is not None:
            self.c_trans = trans
        
        if mod is not None:
            self.ui.modulus_sb.setValue(mod)
            self.ui.poisson_sb.setValue(pr)
        
        if self.mesh is not None:
            
            self.ui.mirror_button.setEnabled(True)
            self.ui.move_button.setEnabled(True)
            self.ui.reset_node_selection_button.setEnabled(True)
            self.ui.corner_best_fit_button.setEnabled(True)
            self.ui.vtk_icp_button.setEnabled(True)
            self.ui.bc_box.setEnabled(True)
            self.ui.save_box.setEnabled(True)
            
            self.update_mesh_display()
            self.impose_surface_bc()
            
        if not len(these_bc_nodes) == 0:
            self.bc_nodes = these_bc_nodes
            self.redraw_bc_pick()
        
        #populate ui with 0th entry of outline
        _, perimeter, _, = respace_equally(self.rs_outlines[0],1)
        self.ui.seed_by_length_sb.setValue(perimeter/len(self.rs_outlines[0]))
        self.ui.seed_by_num_sb.setValue(len(self.rs_outlines[0]))
        
        self.update_outlines()
        self.make_orientation_widget()

        self.ui.mesh_box.setEnabled(True)
        self.ui.outline_box.setEnabled(True)
    
    def reset_outlines(self):
    
        target = self.ui.entry_spec.currentIndex()

        self.rs_outlines[target] = self.outlines[target].copy()

        _, perimeter, _, = respace_equally(self.rs_outlines[target],1)
        self.ui.seed_by_length_sb.setValue(perimeter/len(self.rs_outlines[target]))
        self.ui.seed_by_num_sb.setValue(len(self.rs_outlines[target]))

        self.update_outlines()
        
    def update_outlines(self):
        '''
        Calls mod_outline with input from ui
        '''
        
        target = self.ui.entry_spec.currentIndex()

        if self.ui.quantity_rb.isChecked():
            self.rs_outlines[target], self.rs_outline_corners[target] = mod_outline(self.outlines[target], self.ui.seed_by_num_sb.value(), None)
        else:
            self.rs_outlines[target], self.rs_outline_corners[target] = mod_outline(self.outlines[target], None, self.ui.seed_by_length_sb.value())
        
        self.focus_outline_entry()
        self.draw_outlines()
        self.draw_corners()


    def hide_caption(self):
        for actor in self.caption_actor_list:
            flip_visible(actor)
        self.ui.vtkWidget.update()


    def export_outline(self):
        '''
        Gets a file and calls write dxf for active outline
        '''
        
        fileo, _ = get_save_file("*.dxf", self.active_dir)
        
        #make sure valid file was selected
        if fileo is None:
            return
        
        target = self.ui.entry_spec.currentIndex()
        write_dxf(fileo, self.rs_outlines[target])

    def draw_outlines(self):
    
        #remove all avg actors
        if hasattr(self,'rs_outline_actor_list'):
            for actor in self.rs_outline_actor_list:
                self.ren.RemoveActor(actor)
        
        self.rs_outline_actor_list = []
        self.caption_actor_list = []
        
        for i in range(len(self.outlines)):
            o_actor, _ = gen_outline(self.rs_outlines[i],\
            (255, 0, 0),\
            self.point_size)
            o_actor.SetPickable(0)
            self.ren.AddActor(o_actor)
            self.rs_outline_actor_list.append(o_actor)
            
            # debug
            # for idx in range(len(self.rs_outline_corners[i])):
                # print(self.rs_outline_corners[i][idx])
                
            cap_actor = gen_caption_actor('%s'%i, o_actor)
            self.ren.AddActor(cap_actor)
            self.caption_actor_list.append(cap_actor)

        if not self.ui.caption_rb.isChecked():
            self.hide_caption()

        self.ren.ResetCamera()
        self.ui.vtkWidget.update()

    def draw_corners(self):
        '''
        Draws a highlighted 'corner' on each outline and provides a number caption of that
        '''
        
        if hasattr(self,'corner_actor_list'):
            for actor in self.corner_actor_list:
                self.ren.RemoveActor(actor)
            for actor in self.corner_caption_actor_list:
                self.ren.RemoveActor(actor)
        
        self.ball_radius = self.ui.seed_by_length_sb.value()*0.8
        
        self.corner_actor_list = []
        self.corner_caption_actor_list = []
        self.corner_pts = []
        c = 0
        for i in range(len(self.rs_outlines)):
            corner_ind, self.rs_outlines[i] = get_corner_ind(self.rs_outlines[i])
            for ind in corner_ind:
                sphere = vtk.vtkSphereSource()
                sphere.SetPhiResolution(24)
                sphere.SetThetaResolution(24)
                sphere.SetRadius(self.ball_radius)
                sphere.SetCenter(self.rs_outlines[i][ind,:])
                # debug
                # if ind == 0:
                    # print('Corner 0:', outline[ind,:])
                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputConnection(sphere.GetOutputPort())
                local_actor = vtk.vtkActor()
                local_actor.SetMapper(mapper)
                self.ren.AddActor(local_actor)
                self.corner_actor_list.append(local_actor)
                
                cap_actor = gen_caption_actor('%s'%c, local_actor, (1,1,1))
                self.ren.AddActor(cap_actor)
                self.corner_caption_actor_list.append(cap_actor)
                self.corner_pts.append(self.rs_outlines[i][ind,:])
                c+=1

        if not self.ui.show_corners_rb.isChecked():
            #change visibility of actors
            self.show_corners()

    def show_corners(self):
        for actor in self.corner_actor_list:
            flip_visible(actor)
        for actor in self.corner_caption_actor_list:
            flip_visible(actor)
        self.ui.vtkWidget.update()

    def focus_outline_entry(self):
        
        target = self.ui.entry_spec.currentIndex()
        _, perimeter, _, = respace_equally(self.rs_outlines[target],1)
        self.ui.seed_by_length_sb.setValue(perimeter/len(self.rs_outlines[target]))
        self.ui.seed_by_num_sb.setValue(len(self.rs_outlines[target]))
        

    def generate_extrusion(self):
        '''
        Generates & populates the extrusion widget dialog box. Does nothing if the extrusion widget does not return a file on closing, otherwise, calls add_vtk on what is returned.
        '''
        
        target = self.ui.entry_spec.currentIndex()
        outline = self.rs_outlines[target]
        self.focus_outline_entry()
        dist = self.ui.seed_by_length_sb.value()
        limits = get_limits(outline,0)
        min_length=np.minimum(limits[1]-limits[0],limits[3]-limits[2])
        length = round(3*min_length)
        ew = extrude_widget(self, outline, dist, 11, length)
        ew.exec_()
        if ew.vtk_file is not None:
            self.add_vtk(file = ew.vtk_file)
        else:
            return
        
    def import_mesh(self):
        '''
        Imports a mesh either in vtk format and read directly, or inp format, which is converted to an unstructured grid with a call to convert_inp_to_vtk
        '''
        f,_ = get_file("mesh")
        if f is None or not(os.path.isfile(f)):
            return
        if f.endswith('vtk'):
            self.add_vtk(file = f)
        elif f.endswith('inp'):
            convert_inp_to_vtk(f,os.path.splitext(f)[0]+'.vtk')
            self.add_vtk(os.path.splitext(f)[0]+'.vtk')

    def add_vtk(self, file = None):
        if file is None:
            return

        mesh_source = vtk.vtkUnstructuredGridReader()
        mesh_source.SetFileName(file)
        mesh_source.Update()
        
        self.mesh = mesh_source.GetOutput()

        self.ui.mirror_button.setEnabled(True)
        self.ui.move_button.setEnabled(True)
        self.ui.reset_node_selection_button.setEnabled(True)
        self.ui.corner_best_fit_button.setEnabled(True)
        self.ui.vtk_icp_button.setEnabled(True)
        self.ui.bc_box.setEnabled(True)
        self.ui.save_box.setEnabled(True)
        
        self.reset_bc_pick()
        self.reset_node_pick()
        self.reset_fit_display()
        
        self.update_mesh_display()
        self.ren.ResetCamera()


    def update_mesh_display(self):
    
        if hasattr(self,'mesh_actor'):
            self.ren.RemoveActor(self.mesh_actor)
    
        self.mesh_mapper=vtk.vtkDataSetMapper()
        self.mesh_mapper.SetInputData(self.mesh)

        self.mesh_actor = vtk.vtkActor()
        self.mesh_actor.SetMapper(self.mesh_mapper)

        self.mesh_actor.GetProperty().SetLineWidth(1)
        self.mesh_actor.GetProperty().SetColor(0,0.9020,0.9020) #abaqus style
        # self.mesh_actor.GetProperty().SetColor(0,1,0.6039) #gmsh
        self.mesh_actor.GetProperty().SetEdgeColor([0.8, 0.8, 0.8])
        self.mesh_actor.GetProperty().EdgeVisibilityOn()
        self.ren.AddActor(self.mesh_actor)
        self.redraw_node_pick()
        self.redraw_bc_pick()
        self.ui.vtkWidget.update()

    def redraw_node_pick(self):
        '''
        Redraws picked nodes after a transformation or on load directly from selected_nodes
        '''
        if not hasattr(self,'selected_nodes'):
            return
        selected_nodes = self.selected_nodes
        
        num_redraw = sum(x is not None for x in selected_nodes)
        if num_redraw == 0:
            return

        #clear all selected nodes
        self.reset_node_pick()
        
        #redraw and update lists
        color = tuple(vtk.vtkNamedColors().GetColor3d("wheat"))
        
        for ind in range(num_redraw):
            sphere = vtk.vtkSphereSource()
            sphere.SetPhiResolution(24)
            sphere.SetThetaResolution(24)
            sphere.SetRadius(self.ball_radius)
            sphere.SetCenter(self.mesh.GetPoint(selected_nodes[ind]))
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            this_actor = vtk.vtkActor()
            this_actor.SetMapper(mapper)
            this_actor.GetProperty().SetColor(color)
            self.ren.AddActor(this_actor)
            
            #add to lists and generate caption
            self.picked_node_actor_list[ind] = this_actor
            cap_actor = gen_caption_actor('%s'%ind, this_actor, color)
            self.ren.AddActor(cap_actor)
            self.picked_node_caption_actor_list[ind] = cap_actor
            self.selected_nodes[ind] = selected_nodes[ind]

    def actuate_node_pick(self):
        '''
        Starts picking and handles ui button display
        '''
        
        if self.bc_picking:
            self.actuate_bc_pick()
        
        if hasattr(self,'selected_actor'):
            self.ren.RemoveActor(self.selected_actor)
        
        if self.picking:
            #Remove picking observer and re-initialise
            self.iren.RemoveObservers('LeftButtonPressEvent')
            self.iren.AddObserver('LeftButtonPressEvent',self.default_left_button)
            QtWidgets.QApplication.processEvents()
            self.picking = False
            self.ui.active_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
            self.ui.bc_box.setEnabled(True)
        else:
            if not hasattr(self,'selected_nodes'):
                self.reset_node_pick()
            self.iren.AddObserver('LeftButtonPressEvent', self.picker_callback)
            self.picking = True
            self.ui.active_picking_indicator.setStyleSheet("background-color :rgb(77, 209, 97);")
            self.ui.bc_box.setEnabled(False)
    
    def default_left_button(self, obj, event):
        #forward standard events according to the default style`
        self.iren.GetInteractorStyle().OnLeftButtonDown()

    def reset_node_pick(self):
        
        if hasattr(self,'picked_node_actor_list'):
            for actor in self.picked_node_actor_list:
                self.ren.RemoveActor(actor)
            for actor in self.picked_node_caption_actor_list:
                self.ren.RemoveActor(actor)
        
        
        self.picked_node_actor_list = [None]*len(self.corner_actor_list)
        self.picked_node_caption_actor_list = [None]*len(self.corner_actor_list)
        self.selected_nodes = [None]*len(self.corner_actor_list)
        
        self.ui.vtkWidget.update()
        
    def accept_node_pick(self):
    
    
        colors = vtk.vtkNamedColors()
        this_color = tuple(colors.GetColor3d("wheat"))
        this_polydata = vtk.vtkPolyData()
        this_polydata.DeepCopy(self.selected_actor.GetMapper().GetInput())
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(this_polydata)
        this_actor = vtk.vtkActor()
        this_actor.SetMapper(mapper)
        this_actor.GetProperty().SetColor(this_color)
        self.ren.AddActor(this_actor)
        
        #condition if the last actor available is selected
        if self.selected_nodes[-1] is not None:
            self.ren.RemoveActor(self.picked_node_actor_list[-1])
            self.ren.RemoveActor(self.picked_node_caption_actor_list[-1])
            self.picked_node_actor_list[-1] = None
            self.picked_node_caption_actor_list[-1] = None
            self.selected_nodes[-1] = None
        
        #get index of first 'None' entry in self.selected_nodes
        for ind in range(len(self.selected_nodes)):
            if self.selected_nodes[ind] is None:
                break

        self.picked_node_actor_list[ind] = this_actor
        cap_actor = gen_caption_actor('%s'%ind, this_actor, this_color)
        self.ren.AddActor(cap_actor)
        self.picked_node_caption_actor_list[ind] = cap_actor
        self.selected_nodes[ind] = self.selected_node
        
    def delete_node_pick(self):
        
        #get index of first not 'None' entry in self.selected_nodes
        for i in range(len(self.selected_nodes)):
            if self.selected_nodes[i] == None:
                if i > 0: ind = i-1
                else: i = 0
                break
            else:
                ind = i

        try:
            self.ren.RemoveActor(self.picked_node_actor_list[ind])
            self.ren.RemoveActor(self.picked_node_caption_actor_list[ind])
            self.picked_node_actor_list[ind] = None
            self.picked_node_caption_actor_list[ind] = None
            self.selected_nodes[ind] = None
        except:
            return


    def actuate_bc_pick(self):
        '''
        Starts picking and handles ui button display
        '''

        if self.picking:
            self.actuate_node_pick()

        if hasattr(self,'selected_actor'):
            self.ren.RemoveActor(self.selected_actor)
        
        if self.bc_picking:
            self.iren.RemoveObservers('LeftButtonPressEvent')
            self.iren.AddObserver('LeftButtonPressEvent',self.default_left_button)
            QtWidgets.QApplication.processEvents()
            
            QtWidgets.QApplication.processEvents()
            self.bc_picking = False
            self.ui.bc_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
            self.ui.mesh_box.setEnabled(True)
        else:
            if not hasattr(self,'bc_nodes'):
                self.reset_bc_pick()
            self.iren.AddObserver('LeftButtonPressEvent', self.picker_callback)
            self.bc_picking = True
            self.ui.bc_picking_indicator.setStyleSheet("background-color :rgb(77, 209, 97);")
            self.ui.mesh_box.setEnabled(False)

    def reset_bc_pick(self):
        
        if hasattr(self,'bc_actor_list'):
            for actor in self.bc_actor_list:
                self.ren.RemoveActor(actor)
        
        self.bc_actor_list = [None]*2
        self.bc_nodes = [None]*2
        
        self.ui.vtkWidget.update()

    def accept_bc_pick(self):
        
        #condition if the last actor available is selected
        if self.bc_nodes[-1] is not None:
            self.ren.RemoveActor(self.bc_actor_list[-1])
            self.bc_actor_list[-1] = None
            self.bc_nodes[-1] = None
        
        #get index of first 'None' entry in self.selected_nodes
        for ind in range(len(self.bc_nodes)):
            if self.bc_nodes[ind] is None:
                break
        
        self.bc_nodes[ind] = self.selected_node
        self.draw_bc(ind)
    
    def redraw_bc_pick(self):
        if not hasattr(self,'bc_nodes'):
            return
        
        bc_nodes = self.bc_nodes #can't do list comprehensions over objects
        num_redraw = sum(x is not None for x in bc_nodes)
        if num_redraw == 0:
            return
        
        self.reset_bc_pick()
        
        for ind in range(num_redraw):
            self.bc_nodes[ind] = bc_nodes[ind]
            self.draw_bc(ind)
    
    def draw_bc(self,ind):
        '''
        Draws the boundary condition and updates lists according to ind
        '''
    
        pnt = self.mesh.GetPoint(self.bc_nodes[ind])
        
        limits = self.mesh_actor.GetBounds()
        a_len = np.maximum(limits[1]-limits[0],limits[3]-limits[2])*0.025

        c_target=np.array([
        [limits[0]-1,limits[2]-1], #xmin,ymin
        [limits[0]-1,limits[3]+1], #xmin,ymax
        [limits[1]+1,limits[3]+1], #xmax,ymax
        [limits[1]+1,limits[2]-1] #xmax,ymin
        ])

        #identify closest corner based on minimum euclidean distance
        d = np.sqrt((c_target[:,0] - pnt[0])**2 + (c_target[:,1] - pnt[1])**2)
        closest_corner = c_target[np.argmin(d),:]
        #get normalized x & y components of vector from pnt to closest_corner
        v = closest_corner - pnt[0:2]
        
        mapper = vtk.vtkDataSetMapper()
        
        #generate actor based on ind value
        if ind == 0: #x&y orientation
            arrow_directions = \
            np.array([[v[0]/np.abs(v[0]),0,0],[0,v[1]/np.abs(v[1]),0]])
        
            _, a1_pd = draw_arrow(pnt, a_len, arrow_directions[0], None)
            _, a2_pd = draw_arrow(pnt, a_len, arrow_directions[1], None)
            
            append_filter = vtk.vtkAppendFilter()
            append_filter.AddInputData(a1_pd)
            append_filter.AddInputData(a2_pd)
            append_filter.Update()
            
            mapper.SetInputData(append_filter.GetOutput())
        
        else: #just y orientation
            arrow_directions = np.array([[0,v[1]/np.abs(v[1]),0]])
        
            _, a1_pd = draw_arrow(pnt, a_len, arrow_directions[0], None)
            
            mapper.SetInputData(a1_pd)

        arrow_color = tuple(vtk.vtkNamedColors().GetColor3d("red"))
        arrow_actor = vtk.vtkActor()
        arrow_actor.GetProperty().SetColor(arrow_color)
        arrow_actor.SetMapper(mapper)
        self.ren.AddActor(arrow_actor)
        self.bc_actor_list[ind] = arrow_actor
        

        
    def delete_bc_pick(self):
        
        #get index of first not 'None' entry in self.selected_nodes
        for i in range(len(self.bc_nodes)):
            if self.bc_nodes[i] == None:
                if i > 0: ind = i-1
                else: i = 0
                break
            else:
                ind = i

        try:
            self.ren.RemoveActor(self.bc_actor_list[ind])
            self.bc_actor_list[ind] = None
            self.bc_nodes[ind] = None
        except:
            return

    def picker_callback(self, obj, event):
    
        colors = vtk.vtkNamedColors()
        
        picker = vtk.vtkPointPicker()
        picker.SetTolerance(0.005)
        
        pos = self.iren.GetEventPosition()
        
        picker.Pick(pos[0], pos[1], 0, self.ren)
        
        if picker.GetPointId() != -1:
        
            ids = vtk.vtkIdTypeArray()
            ids.SetNumberOfComponents(1)
            ids.InsertNextValue(picker.GetPointId())
            #debug
            # print(picker.GetPointId())
            # print(self.mesh.GetPoint(picker.GetPointId()))
            # val = self.mesh.GetPoint(picker.GetPointId())

            if hasattr(self,'selected_actor'):
                self.ren.RemoveActor(self.selected_actor)
            self.selected_node = picker.GetPointId()
            
            sphere = vtk.vtkSphereSource()
            sphere.SetPhiResolution(24)
            sphere.SetThetaResolution(24)
            sphere.SetRadius(self.ball_radius)
            sphere.SetCenter(self.mesh.GetPoint(picker.GetPointId()))
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            self.selected_actor = vtk.vtkActor()
            self.selected_actor.SetMapper(mapper)
            self.selected_actor.GetProperty().SetColor(colors.GetColor3d("tomato"))
            self.ren.AddActor(self.selected_actor)
        
    def mirror_mesh(self):
        
        self.trans[2,2] = -1
        self.apply_transformation()
        #nodes will have moved, but the connectivity has to be updated otherwise there will be inside-out elements
        
        cell_type = vtk.vtkCellTypes()
        self.mesh.GetCellTypes(cell_type)
    
        if cell_type.IsType(24) == 1: #if it contains 2nd order tets
            vtk_type = 24
            n_nodes_per_element = 10
            local_index = np.array([2,1,0,3,5,4,6,9,8,7])
        elif cell_type.IsType(12) == 1: #if it contains linear quads
            vtk_type = 12
            n_nodes_per_element = 8
            local_index = np.array([4,5,6,7,0,1,2,3])

        mirror_ugrid = vtk.vtkUnstructuredGrid()
        mirror_ugrid.SetPoints(self.mesh.GetPoints())

        cells = v2n.vtk_to_numpy(self.mesh.GetCells().GetData())
        cells.shape = (-1,cells[0]+1)
        
        for row in range(len(cells)):
            this_cell = cells[row,1::]
            mirror_ugrid.InsertNextCell(\
            vtk_type, n_nodes_per_element, this_cell[local_index])
        
        self.mesh = mirror_ugrid
        self.mesh.Modified()

    def trans_from_ui(self):
    
        self.trans[0,-1] = self.ui.trans_x.value()
        self.trans[1,-1] = self.ui.trans_y.value()
        a = np.deg2rad(self.ui.rotate_z.value()) #negative for counterclockwise
        self.trans[0:2,0:2]=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        self.apply_transformation()

    def trans_from_corners(self):
        '''
        Uses corners from outline (B) versus user selected mesh corners (A) and moves the mesh to suit.
        '''
        
        if self.picking:
            self.actuate_node_pick()
        
        #number of not-none entries in selected_nodes
        selected_nodes = self.selected_nodes
        num_compare = sum(x is not None for x in selected_nodes)
        
        if num_compare < 3:
            info_msg('At least three points on the outline and mesh need to be selected. Currently there are %d.'%num_compare)
            return
        
        B = self.corner_pts[:num_compare]
        
        A = []
        for i in range(num_compare):
            pnt = self.selected_nodes[i]
            A.append(np.asarray(self.mesh.GetPoint(pnt)))
        
        _, R, t = best_fit_transform(np.vstack(A), np.vstack(B))
        #pack relevant values into self.trans
        self.trans[:3,-1] = t
        self.trans[0:3,0:3] = R
        self.apply_transformation()

    def trans_from_vtk_icp(self):
        '''
        Get icp transformation matrix from outline polydata
        '''
        
        if self.picking:
            self.actuate_node_pick()
        
        # Get the number of nodes to compare
        selected_nodes = self.selected_nodes
        num_compare = sum(x is not None for x in selected_nodes)
        
        if num_compare < 3:
            info_msg('At least three points on the outline and mesh need to be selected. Currently there are %d.'%num_compare)
            return

        A = []
        for i in range(num_compare):
            pnt = self.selected_nodes[i]
            A.append(np.asarray(self.mesh.GetPoint(pnt)))
        A = np.asarray(A)
        B = np.asarray(self.corner_pts[:num_compare])

        #generate polydata objects of A and B points for vtk filter
        A_pd = gen_point_cloud(A)[1]
        B_pd = gen_point_cloud(B)[1]
        
        self.trans = vtk_icp(A_pd,B_pd)
        self.apply_transformation()

    def apply_transformation(self):
    
        T = self.trans.copy()

        np_pts = v2n.vtk_to_numpy(self.mesh.GetPoints().GetData())
        np_pts = do_transform(np_pts,T)
        self.c_trans = T @ self.c_trans
        self.trans = np.eye(4)
        
        self.mesh.GetPoints().SetData(v2n.numpy_to_vtk(np_pts))
        self.mesh.Modified()
        self.update_mesh_display()
        #if there's already a displacement actor present
        if hasattr(self,'bc_disp_actor_list'):
            self.impose_surface_bc()
    
    def reset_fit_display(self):
        if hasattr(self,'bc_disp_actor_list'):
            for actor in self.bc_disp_actor_list:
                self.ren.RemoveActor(actor)
            del self.bc_disp_nodes, self.bc_disp_val, self.bc_disp_actor_list

    def impose_surface_bc(self):
        
        self.reset_fit_display()
        dist = self.ui.seed_by_length_sb.value()
        
        #find out what side of the xy plane the mesh is
        pos_position = False
        if np.sum(self.mesh_actor.GetBounds()[-2:])/2 > 0:
            pos_position = True
        
        #initialize lists for bc's
        self.bc_disp_nodes = [None] * len(self.outlines)
        self.bc_disp_val = [None] * len(self.outlines)
        self.bc_disp_actor_list = []
        
        for index in range(len(self.outlines)):
            these_nodes, these_disp_val, this_actor = get_surf_bc(\
            self.mesh, \
            self.outlines[index], \
            self.bv_tck_list[index][:5], \
            dist, \
            pos_position)
            #if bc_nodes == None, then nothing was returned from the underlying locator
            if these_nodes is None:
                pass
            else:
                self.bc_disp_nodes[index] = these_nodes
                self.bc_disp_val[index] = these_disp_val
                self.bc_disp_actor_list.append(this_actor)
                self.ren.AddActor(this_actor)

        self.mesh_actor.GetProperty().SetOpacity(0.5)
        self.update_z_aspect()
        self.ui.vtkWidget.update()

    def check_save_state(self,id):
        
        #check if there's something in mesh/point_data
        with h5py.File(self.file, 'r') as f:
            if 'mesh/point_data' in f:
                if len(f['mesh/point_data'].keys()) != 0:
                    overwrite = warning_msg(self, \
                        'Saving will overwrite the existing FEA results requiring running the calculation again. Continue?'
                        )
                    if not overwrite: #fail check
                        return False

        #check itinerary
        if not hasattr(self,'bc_disp_nodes'):
            info_msg('Surface displacement boundary conditions need to be imposed.')
            return False
        if not hasattr(self,'bc_nodes'):
            info_msg('Rigid body boundary conditions need to be imposed.')
            return False
        num_rigid_body_nodes =  sum(x is not None for x in list(self.bc_nodes))
        if num_rigid_body_nodes < 2:
            info_msg('Two rigid body boundary conditions need to be imposed. There are currently %d.'%num_rigid_body_nodes)
            return False
        
        #flag if len(bc_disp_nodes) having valid entries is less than the number of outlines
        test_list = self.bc_disp_nodes
        if any(entry is None for entry in test_list):
            response = warning_msg(self, \
                'There are %d entries corresponding to outlines, and only %d set(s) of surface boundary conditions applied to the mesh. Continue?'%(len(self.rs_outlines), sum(x is not None for x in self.bc_disp_nodes))
                )
            if not response:
                return False
        
        with h5py.File(self.file, 'r') as f:
            #check the id's entry keys
            if id in f.keys():
                g = f['%s'%id]
            else:
                return True
            if 'rigid_body_nodes' in list(g.keys()):
                overwrite = warning_msg(self, \
                'There is existing data in the specified file for the target field, overwrite?'
                )
                if not overwrite: #fail check
                    return False
                else:
                    return True
            else:
                #there isn't an entry in the id
                return True
        

    def write_output(self):
        '''
        Checks whether sufficient input has been received to run the FEA, and loads h5 file up with relevant details.
        '''
        passed = self.check_save_state('bc_prop')
        
        if not passed:
            self.ui.save_label.setText('Ready')
            return
        
        #write data - if it doesn't exist, create entries, otherwise assign values to existing
        self.ui.save_label.setText('Saving . . .')
        QtWidgets.QApplication.processEvents()
        
        #per entry details
        entry_dict = {'outlines': self.rs_outlines, 'surface_nodes': self.bc_disp_nodes, 'bc_disp_val': self.bc_disp_val}
        
        with h5py.File(self.file, 'r+') as f:
            if 'bc_prop' in list(f.keys()):
                del f['bc_prop']
            g = f.create_group('bc_prop')
            for i in range(len(self.rs_outlines)):
                gg = g.create_group(str(i))
                gg.create_dataset('outline', data = self.rs_outlines[i])
                if self.bc_disp_nodes[i] is not None:
                    gg.create_dataset('surface_nodes', data = self.bc_disp_nodes[i])
                    gg.create_dataset('bc_disp_val', data = self.bc_disp_val[i])
            #per mesh entries
            g.create_dataset('rigid_body_nodes', data = np.array(self.bc_nodes))
            g.create_dataset('transform', data = self.c_trans)
            g.attrs['modulus'] = self.ui.modulus_sb.value()
            g.attrs['poisson_ratio'] = self.ui.poisson_sb.value()

        #now write mesh; date_modfied is applied by vtkug_writer
        pt = vtk.vtkPassThrough() #need a filter object to set up pipeline
        pt.SetInputData(self.mesh)
        w = vtkug_writer()
        w.SetInputConnection(pt.GetOutputPort())
        w.SetFileName(self.file)
        w.Update()
        

        self.ui.save_label.setText('Saved.')
        
        if self.ui.run_FEA_rb.isChecked():
            fw = fea_widget(self,self.file)
            fw.exec_()
        

    def update_z_aspect(self):
        '''
        Updates z_aspect and redraws points displayed based on new z_aspect
        '''
        
        if not hasattr(self,'bc_disp_actor_list'):
            return
        
        z_aspect = self.ui.z_aspect_sb.value()
        #scale points, not outlines, so skipping every other entry in the actor list
        for actor in (self.bc_disp_actor_list):
            actor.SetScale((1,1,z_aspect))
            actor.Modified()
        
        self.ui.vtkWidget.update()
        
    def keypress(self,obj,event):
        '''
        VTK Interactor specific keypress binding
        '''
        key = obj.GetKeyCode()
        if key == "z":
            self.ui.z_aspect_sb.setValue(self.ui.z_aspect_sb.value()*2)
        elif key == "x":
            self.ui.z_aspect_sb.setValue(int(self.ui.z_aspect_sb.value()*0.5))
        elif key == "c":
            self.ui.z_aspect_sb.setValue(1)
        elif key == "1":
            xyview(self.ren)
        elif key == "2":
            yzview(self.ren)
        elif key == "3":
            xzview(self.ren)
        elif key == "Up":
            self.ren.GetActiveCamera().Roll(30)
            self.ren.ResetCamera()
        elif key == "Down":
            self.ren.GetActiveCamera().Roll(-30)
            self.ren.ResetCamera()
        elif key == "n":
            self.actuate_node_pick()
        elif key == "b":
            self.actuate_bc_pick()
        elif key == "a":
            if self.picking:
                self.accept_node_pick()
            if self.bc_picking:
                self.accept_bc_pick()
        elif key == "d":
            if self.picking:
                self.delete_node_pick()
            if self.bc_picking:
                self.delete_bc_pick()
        if key == "l":
            self.file = None
            self.get_data()

        self.ui.vtkWidget.update()


def get_surf_bc(mesh, outline, fit, dist, pos):
    '''
    Main method that takes an unstructured grid, applies a locator of outline x dist in volume, oriented either in the positive z side of the xy plane (pos = true) or negative (pos = false). Finds cells & their nodes within this locator and applies fit to the z coordinate. Returns numpy arrays of bc_nodes (nodes of the incoming mesh where bcs are applied), and the displacement bc, along with an actor representing the boundary condition.
    '''
    
    #get cell type and set parameters of the incoming mesh
    cell_type = vtk.vtkCellTypes()
    mesh.GetCellTypes(cell_type)

    if cell_type.IsType(24) == 1:
        n_nodes_per_element = 10
        bc_unit = 6
        f_lookup = np.array([\
        [0,4,1,5,2,6],
        [1,8,3,9,2,5],
        [3,7,0,6,2,9],
        [0,4,1,8,3,7]])
    elif cell_type.IsType(12) == 1:
        n_nodes_per_element = 8
        bc_unit = 4
        f_lookup = np.array([\
        [0,1,2,3],
        [4,5,6,7],
        [0,1,5,4],
        [1,2,6,5],
        [2,3,7,6],
        [3,0,4,7]])
    
    #get bounding box for cell locator
    outline_limits = get_limits(outline)
    if pos: #mesh nodes z coord all > 0
        bounds = tuple(outline_limits[0:4])+(-0.1,dist)
    else:
        bounds = tuple(outline_limits[0:4])+(-dist, 0.1)
    
    #create cell locator that looks on either side of the xy plane for cells.
    surf_candidate_id_list = vtk.vtkIdList()
    locator = vtk.vtkCellTreeLocator()
    locator.SetDataSet(mesh)
    locator.BuildLocator()
    surf_candidate_id_list = vtk.vtkIdList()
    locator.FindCellsWithinBounds(bounds,surf_candidate_id_list)
    locator.Update()
    
    num_cells = surf_candidate_id_list.GetNumberOfIds()
    if num_cells == 0:
        return None, None, None #return nones because no cells were found
    
    cells  = vtk.vtkCellArray()
    pid_list = vtk.vtkIdList()

    for i in range(num_cells):
        mesh.GetFaceStream(surf_candidate_id_list.GetId(i), pid_list)
        cells.InsertNextCell(pid_list)
        
    raw_point_ids = v2n.vtk_to_numpy(cells.GetData()) #1d list of cell_type & connectivity
    
    #filter out any non volumetric elements  using the same basis as gen_filtered_ugrid in pyCMcommon
    cell_offsets = v2n.vtk_to_numpy(cells.GetOffsetsArray())
    cell_points = []#np.array([])
    
    ind = 0 #index of raw_points
    for i in range(1,len(cell_offsets)):
        local = raw_point_ids[ind:cell_offsets[i]+i]
        if len(local) == n_nodes_per_element+1: #if there are low order elements
            cell_points.append(local)
        ind = cell_offsets[i]+i

    #cell_points now contains node number, connectivity1 . . . connectivity_n
    cell_points = np.array(cell_points)
    #remove node number from cell_points
    cell_points = cell_points[:,1::]

    bc_pnts = vtk.vtkPoints()
    bc_cells = vtk.vtkCellArray()
    bc_nodes = []
    bc_val = []
    
    #main loop that finds faces of elements and populates a new polydata object
    p_count = 0
    for element in cell_points:
        local_cell = vtk.vtkPolygon()
        local_cell.GetPointIds().SetNumberOfIds(bc_unit)
       
        found_face = False
        for face in range(len(f_lookup)):
            z_value = np.array([mesh.GetPoint(node)[-1] for node in element[f_lookup[face]]])
            if not np.any(z_value):
                face_ind = face
                found_face = True
                break
        if found_face:
            #build new cells and add points to new vtkPolygon object
            local_cell_ind = 0
            for node_index in element[f_lookup[face_ind]]:
                if node_index not in bc_nodes:
                    p = mesh.GetPoint(node_index)
                    bc_nodes.append(node_index)
                    #assign values of z based on what side of the z axis the mesh is on
                    if pos:
                        bc_z_val = bisplev(p[0],p[1],fit)
                    else:
                        bc_z_val = -bisplev(p[0],p[1],fit)
                    bc_val.append(bc_z_val)
                    bc_pnts.InsertNextPoint(np.append(p[:2],bc_z_val))
                    local_cell.GetPointIds().SetId(local_cell_ind,p_count)
                    p_count += 1
                else:
                    local_cell.GetPointIds().SetId(local_cell_ind,bc_nodes.index(node_index))
                local_cell_ind += 1
            bc_cells.InsertNextCell(local_cell)

    bc_pd = vtk.vtkPolyData()
    bc_pd.SetPoints(bc_pnts)
    bc_pd.SetPolys(bc_cells)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(bc_pd)
    bc_actor = vtk.vtkActor()
    bc_actor.SetMapper(mapper)
    bc_actor.GetProperty().SetColor(1,0.804,0.204) #mustard
    if bc_unit == 4:
        bc_actor.GetProperty().EdgeVisibilityOn()
    elif bc_unit == 6:
        bc_actor.GetProperty().EdgeVisibilityOff()
    bc_actor.GetProperty().SetRepresentationToSurface()
    bc_actor.SetPickable(0)

    return np.array(bc_nodes, dtype = int), np.array(bc_val, dtype = float), bc_actor

def mod_outline(outline,num,length):
    '''
    Based on entries in ui, converge on a respaced outline while preserving corners
    '''
    
    #get corners and re-order outline to match
    corner_ind, indexed_outline = get_corner_ind(outline)
    outline = indexed_outline[:,:2]
    
    #get perimeter
    _, perimeter, _, = respace_equally(outline,1)
    
    
    while 1:
        respaced_outline = np.array([]).reshape(0,2)
        if num != None and length == None:
            if 'dist' not in vars():
                #then quantity, need to get dist
                dist = perimeter / float(num)
                
            for j in range(3):
                #solve for first three segments
                segment = outline[corner_ind[j]:corner_ind[j+1] + 1]
                _, _, npts = respace_equally(segment,dist)
                new_segment, _, _, = respace_equally(segment,int(npts+1))
                respaced_outline = np.vstack((respaced_outline,new_segment[0:-1,:])) #dropping last point
            #do final segment
            segment = outline[corner_ind[3]::]
            _, _, npts = respace_equally(segment,dist)
            new_segment, _, _, = respace_equally(segment,int(npts+1))
            respaced_outline = np.vstack((respaced_outline,new_segment[0:-1,:]))
            
        if num == None and length != None:
            if 'dist' not in vars():
                dist = length
            for j in range(3):
                segment = outline[corner_ind[j]:corner_ind[j+1]+1]
                new_segment, _, _, = respace_equally(segment,dist)
                respaced_outline = np.vstack((respaced_outline,new_segment[0:-1,:])) #dropping last point
            #do final segment
            segment = outline[corner_ind[3]::]
            new_segment, _, _, = respace_equally(segment,dist)
            respaced_outline = np.vstack((respaced_outline,new_segment[0:-1,:]))
        
        if not np.fmod(len(respaced_outline),2) == 0:
            dist += 0.005
        else:
            break
        
    #pad out with 0s
    final_respaced_outline = np.hstack((respaced_outline,np.zeros([len(respaced_outline[:,0]),1])))
        
    return final_respaced_outline, corner_ind

def gen_arrow_polydata(start,length,direction,invert = True):
    '''
    Returns the polydata of a an arrow, suitable for appending together for a single actor.
    '''
    
    source = vtk.vtkArrowSource()
    source.SetShaftRadius(0.024)
    source.SetTipRadius(0.07)
    source.SetTipLength(0.14)
    source.SetTipResolution(25)
    source.SetShaftResolution(25)
    if invert:
        source.InvertOn()
    else: source.InvertOff()

    end = start + (length*direction)
    norm_x =(end - start)/length

    arbitrary=np.array([1,1,1]) #can be replaced with a random vector
    norm_z=np.cross(norm_x,arbitrary/np.linalg.norm(arbitrary))
    norm_y=np.cross(norm_z,norm_x)
    
    # Create the direction cosine matrix by writing values directly to an identity matrix
    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    for i in range(3):
        matrix.SetElement(i, 0, norm_x[i])
        matrix.SetElement(i, 1, norm_y[i])
        matrix.SetElement(i, 2, norm_z[i])
        
    #Apply transforms
    transform = vtk.vtkTransform()
    transform.Translate(start)
    transform.Concatenate(matrix)
    transform.Scale(length, length, length)
 
    # Transform the polydata
    transform_pd = vtk.vtkTransformPolyDataFilter()
    transform_pd.SetTransform(transform)
    transform_pd.SetInputConnection(source.GetOutputPort())
    transform_pd.Update()
    return transform_pd.GetOutput()

if __name__ == "__main__":
    if len(sys.argv)>1:
        launch(sys.argv[1])
    else:
        launch()