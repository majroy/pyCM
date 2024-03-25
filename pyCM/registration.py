#!/usr/bin/env python
'''
Uses VTK python to allow for editing point clouds associated with the contour method. Full interaction requires a 3-button mouse and keyboard. See documentation for keyboard/mouse bindings.
1.8 - Updated for overall version 2.
'''

__author__ = "M.J. Roy"
__version__ = "1.8"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014--"

import sys, os, datetime
from pkg_resources import Requirement, resource_filename
import numpy as np
from scipy.spatial import Delaunay
import vtk
import vtk.util.numpy_support as v2n
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from pyCM.reg_preview import registration_viewer
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
    
    #if a data file is specified at launch
    if len(args) == 1:
        window.output_filename = args[0]
        window.ui.tab_widget.setCurrentIndex(1)

    ret = app.exec_()
    
    if sys.stdin.isatty() and not hasattr(sys, 'ps1'):
        sys.exit(ret)
    else:
        return window

class reg_main_window(object):
    '''
    Class to setup Qt interactor
    '''
    def setup(self, main_window):
        '''
        Sets up Qt interactor
        '''
        
        #if called as a script, treat as a Qt mainwindow, otherwise a generic widget for embedding elsewhere
        if hasattr(main_window,'setCentralWidget'):
            main_window.setCentralWidget(self.centralWidget)
        else:
            self.centralWidget = main_window
            main_window.setWindowTitle("pyCM - registration v%s" %__version__)
        
        self.tab_widget = QtWidgets.QTabWidget(self.centralWidget)
        self.editor_tab = QtWidgets.QWidget(self.tab_widget)
        self.tab_widget.addTab(self.editor_tab,"Editor")
        self.preview_tab = QtWidgets.QWidget(self.tab_widget)
        self.tab_widget.addTab(self.preview_tab, "Preview")
        
        tab_layout = QtWidgets.QHBoxLayout(self.centralWidget)
        tab_layout.addWidget(self.tab_widget)
        #create new layout to hold both VTK and Qt interactors
        mainlayout = QtWidgets.QHBoxLayout(self.editor_tab)
        previewlayout = QtWidgets.QHBoxLayout(self.preview_tab)
        
        #create VTK widgets
        self.vtkWidget = QVTKRenderWindowInteractor(self.tab_widget)
        self.preview_widget = registration_viewer(self.tab_widget)
        
        #create VTK widget
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(100)
        sizePolicy.setVerticalStretch(100)
        self.vtkWidget.setSizePolicy(sizePolicy)
        self.preview_widget.setSizePolicy(sizePolicy)
        
        self.vtkWidget.setMinimumSize(QtCore.QSize(800, 600))
        
        #Controls for editor tab
        #make import box
        import_box = QtWidgets.QGroupBox('Import')
        import_layout = QtWidgets.QGridLayout()
        import_box.setLayout(import_layout)
        #make combo box to specify what is being loaded
        self.load_spec = QtWidgets.QComboBox()
        self.load_spec.setToolTip('Specify what is to be loaded')
        self.load_spec.addItems(['Outline/perimeter','Outlined surface','Unregistered surface'])
        self.entry_spec = QtWidgets.QComboBox()
        self.entry_spec.addItems(['Entry 0','New entry', 'Clear entries'])
        self.load_label = QtWidgets.QLabel("Ready")
        self.load_label.setWordWrap(True)
        self.load_label.setSizePolicy(sizePolicy)
        import_layout.addWidget(self.load_spec,0,1,1,1)
        import_layout.addWidget(self.entry_spec,0,0,1,1)
        import_layout.addWidget(self.load_label,1,0,1,2)
        
        #make display box
        self.display_box = QtWidgets.QGroupBox('Point cloud display')
        display_layout = QtWidgets.QGridLayout()
        self.display_box.setLayout(display_layout)
        self.z_aspect_sb = QtWidgets.QSpinBox()
        self.z_aspect_sb.setPrefix("Scaling: ")
        self.z_aspect_sb.setToolTip("Scaling factor applied to Z axis")
        self.z_aspect_sb.setValue(1)
        self.z_aspect_sb.setMinimum(1)
        self.z_aspect_sb.setMaximum(1000)
        self.color_on_height=QtWidgets.QPushButton("Colour")
        self.color_on_height.setToolTip('Active points coloured according to height (selected), single colour (deselected)')
        self.color_on_height.setCheckable(True)
        self.color_on_height.setChecked(True)
        
        self.reset_point_display = QtWidgets.QPushButton("Reset")
        self.reset_point_display.setToolTip('Reset contour and redraw active points')
        
        self.hide_display = QtWidgets.QPushButton("Hide")
        self.hide_display.setToolTip('Hide/show active points')
        self.hide_display.setCheckable(True)
        
        self.contour_tools = QtWidgets.QGroupBox('Change contour limits')
        contour_layout = QtWidgets.QGridLayout()
        self.contour_tools.setLayout(contour_layout)
        self.min_contour = QtWidgets.QDoubleSpinBox()
        self.min_contour.setPrefix("Min: ")
        self.min_contour.setMinimum(-100000)
        self.min_contour.setMaximum(100000)
        self.min_contour.setDecimals(4)
        self.max_contour = QtWidgets.QDoubleSpinBox()
        self.max_contour.setPrefix("Max: ")
        self.max_contour.setMinimum(-100000)
        self.max_contour.setMaximum(100000)
        self.max_contour.setDecimals(4)
        self.num_contour = QtWidgets.QSpinBox()
        self.num_contour.setPrefix("Interval: ")
        self.num_contour.setToolTip('Number of entries shown on colorbar')
        self.num_contour.setMinimum(3)
        self.num_contour.setMaximum(20)
        self.num_contour.setValue(13)
        self.update_contours_button = QtWidgets.QPushButton('Update')
        self.update_contours_button.setToolTip('Update the contour limits and interval')

        contour_layout.addWidget(self.min_contour,0,0,1,1)
        contour_layout.addWidget(self.max_contour,0,1,1,1)
        contour_layout.addWidget(self.num_contour,1,0,1,1)
        contour_layout.addWidget(self.update_contours_button,1,1,1,1)

        display_layout.addWidget(self.z_aspect_sb,0,0,1,1)
        display_layout.addWidget(self.color_on_height,0,2,1,1)
        display_layout.addWidget(self.hide_display,1,2,1,1)
        display_layout.addWidget(self.reset_point_display,2,2,1,1)
        display_layout.addWidget(self.contour_tools,1,0,2,2)

        self.display_box.setEnabled(False)

        #make point cloud tools
        self.point_processing_box = QtWidgets.QGroupBox('Decimation')
        point_process_layout = QtWidgets.QGridLayout()
        self.point_processing_box.setLayout(point_process_layout)
        
        #make decimate box - point cloud tools
        self.decimate_box = QtWidgets.QGroupBox('Filter')
        decimate_layout = QtWidgets.QGridLayout()
        self.decimate_box.setLayout(decimate_layout)
        self.num_active_points = QtWidgets.QSpinBox()
        self.num_active_points.setToolTip('Number of active points')
        self.num_active_points.setPrefix("N = ")
        self.num_active_points.setMinimum(0)
        self.num_active_points.setMaximum(10000000)
        self.num_active_points.setValue(0)
        self.percent_active_points = QtWidgets.QDoubleSpinBox()
        self.percent_active_points.setSuffix("%")
        self.percent_active_points.setMinimum(0)
        self.percent_active_points.setMaximum(100)
        self.percent_active_points.setValue(0)
        self.percent_active_points.setToolTip('Percentage of total points that are active')
        z_active_points_label = QtWidgets.QLabel('Z threshold:')
        self.z_active_points = QtWidgets.QDoubleSpinBox()
        self.z_active_points.setMinimum(-100000)
        self.z_active_points.setMaximum(100000)
        self.z_active_points.setDecimals(4)
        self.z_active_points.setToolTip('Points with Z value greater than this are active')
        outline_offset_label = QtWidgets.QLabel('Outline offset:')
        self.outline_offset = QtWidgets.QDoubleSpinBox()
        self.outline_offset.setMinimum(-100000)
        self.outline_offset.setMaximum(100000)
        self.outline_offset.setDecimals(4)
        self.outline_offset.setToolTip('Points within an outline offset by this amount are active')
        
        z_norm_points_label = QtWidgets.QLabel('Z norm cutoff:')
        self.z_norm_active_points = QtWidgets.QDoubleSpinBox()
        self.z_norm_active_points.setDecimals(4)
        self.z_norm_active_points.setMinimum(0.9)
        self.z_norm_active_points.setMaximum(0.9999)
        self.z_norm_active_points.setToolTip('Points corresponding to a triangulation with a normal greater than this are active')
        self.z_norm_active_points.setEnabled(False)
        self.triangulate_button = QtWidgets.QPushButton('Triangulate')
        self.triangulate_button.setToolTip('Perform Delaunay triangulation - required for Z norm cutoff filter and outline processing')
        self.triangulated_indicator = QtWidgets.QLabel('Triangulated')
        self.triangulated_indicator.setStyleSheet("background-color : gray; color : darkGray;")
        self.triangulated_indicator.setAlignment(QtCore.Qt.AlignCenter)
        self.triangulate_button.setSizePolicy(sizePolicy)
        self.triangulated_indicator.setEnabled(False)
        
        self.filter_status_label = QtWidgets.QLabel('Ready')
        self.local_decimate=QtWidgets.QRadioButton("Preserve")
        self.local_decimate.setToolTip('On - deactivate points for local processing, Off - permanently discard points')
        self.local_decimate.setChecked(True)
        
        decimate_z_layout = QtWidgets.QGridLayout()
        decimate_z_layout.addWidget(z_active_points_label,0,0,1,1)
        decimate_z_layout.addWidget(self.z_active_points,0,1,1,1)
        decimate_z_layout.addWidget(outline_offset_label,0,2,1,1)
        decimate_z_layout.addWidget(self.outline_offset,0,3,1,1)
        decimate_z_layout.addWidget(z_norm_points_label,1,0,1,1)
        decimate_z_layout.addWidget(self.z_norm_active_points,1,1,1,1)
        decimate_z_layout.addWidget(self.triangulated_indicator,1,2,1,1)
        decimate_z_layout.addWidget(self.triangulate_button,1,3,1,1)
        
        self.decimate_by_number=QtWidgets.QRadioButton("By number")
        self.decimate_by_percent=QtWidgets.QRadioButton("By percentage")
        self.decimate_by_percent.setChecked(True)
        self.decimate_button_group = QtWidgets.QButtonGroup()
        self.decimate_button_group.addButton(self.decimate_by_number)
        self.decimate_button_group.addButton(self.decimate_by_percent)
        self.decimate_button_group.setExclusive(True)
        decimate_local_button_layout = QtWidgets.QHBoxLayout()
        decimate_local_button_layout.addWidget(self.decimate_by_number)
        decimate_local_button_layout.addWidget(self.num_active_points)
        decimate_local_button_layout.addWidget(self.decimate_by_percent)
        decimate_local_button_layout.addWidget(self.percent_active_points)

        self.update_reduce_points_button = QtWidgets.QPushButton('Preview')
        self.update_reduce_points_button.setToolTip('Preview points that will be discarded.')
        self.undo_reduce_points_button = QtWidgets.QPushButton('Reset')
        self.undo_reduce_points_button.setToolTip('Restore all points, including manually set for decimation.')
        self.apply_reduce_points_button = QtWidgets.QPushButton('Apply')
        self.apply_reduce_points_button.setToolTip('Apply previewed filter.')
        self.apply_reduce_points_button.setEnabled(False)

        decimate_layout.addLayout(decimate_local_button_layout,0,0,1,3)
        decimate_layout.addLayout(decimate_z_layout,1,0,1,3)
        decimate_layout.addWidget(self.update_reduce_points_button,2,0,1,1)
        decimate_layout.addWidget(self.apply_reduce_points_button,2,1,1,1)
        decimate_layout.addWidget(self.undo_reduce_points_button,2,2,1,1)
        decimate_layout.addWidget(self.filter_status_label,3,0,2,1)
        decimate_layout.addWidget(self.local_decimate,3,2,1,1)
        
        #manual picking
        self.decimate_manual_box=QtWidgets.QGroupBox('Selection by area')
        decimate_manual_layout = QtWidgets.QGridLayout()
        self.decimate_manual_box.setLayout(decimate_manual_layout)
        self.active_picking_indicator = QtWidgets.QLabel('Active')
        self.active_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
        self.active_picking_indicator.setAlignment(QtCore.Qt.AlignCenter)
        # self.active_picking_indicator.setFixedSize(50, 20)
        self.active_picking_indicator.setToolTip('Press R with interactor in focus to activate/deactivate manual point selection')
        self.active_picking_indicator.setEnabled(False)
        self.undo_last_pick_button=QtWidgets.QPushButton('Undo last')
        self.undo_last_pick_button.setToolTip('Undo last manual selection')
        self.apply_pick_button=QtWidgets.QPushButton('Apply')
        self.apply_pick_button.setToolTip('Remove selected points and update')
        
        decimate_manual_layout.addWidget(self.active_picking_indicator,0,0,1,1)
        decimate_manual_layout.addWidget(self.undo_last_pick_button,0,1,1,1)
        decimate_manual_layout.addWidget(self.apply_pick_button,0,2,1,1)
        
        #populate point_process_layout
        point_process_layout.addWidget(self.decimate_box)
        point_process_layout.addWidget(self.decimate_manual_box)
        
        self.point_processing_box.setEnabled(False)#default
        
        #make reorientation tools box - point cloud tools
        self.reorient_box = QtWidgets.QGroupBox('Reorientation')
        reorient_layout = QtWidgets.QVBoxLayout()
        self.reorient_box.setLayout(reorient_layout)
        
        self.level_button = QtWidgets.QPushButton("Move to mean Z")
        self.level_button.setToolTip('Translate to mean Z value of active points')
        self.centroid_move_button = QtWidgets.QPushButton("Move to centroid")
        self.centroid_move_button.setToolTip('Translate to centroid of active points')
        self.activate_svd_button=QtWidgets.QPushButton('Apply SVD')
        self.activate_svd_button.setToolTip('Align to Z axis by Single Value Decomposition of up to 10,000 active points')
        
        self.rotate_z_box = QtWidgets.QGroupBox('Rotate about Z')
        z_rotate_layout = QtWidgets.QGridLayout()
        self.rotate_z_box.setLayout(z_rotate_layout)
        self.rotate_z= QtWidgets.QDoubleSpinBox()
        self.rotate_z.setToolTip('Positive is clockwise')
        self.rotate_z.setValue(0)
        self.rotate_z.setMaximum(180)
        self.rotate_z.setMinimum(-180)
        self.rotate_z.setSuffix("\u00b0")
        self.impose_rotation_button = QtWidgets.QPushButton('Apply')
        self.impose_rotation_button.setToolTip('Manually impose rotation about z axis')
        self.auto_rotate_button = QtWidgets.QPushButton('Auto')
        self.auto_rotate_button.setToolTip('Attempt to rotate to align longest side to x axis by rotating about z.')
        z_rotate_layout.addWidget(self.rotate_z,1,1)
        z_rotate_layout.addWidget(self.impose_rotation_button,1,2)
        z_rotate_layout.addWidget(self.auto_rotate_button,1,3)
        
        self.reorient_box.setEnabled(False)
        
        #populate reorientation box
        upper_reorient_layout = QtWidgets.QHBoxLayout()
        upper_reorient_layout.addWidget(self.level_button)
        upper_reorient_layout.addWidget(self.centroid_move_button)
        upper_reorient_layout.addWidget(self.activate_svd_button)
        reorient_layout.addLayout(upper_reorient_layout)
        reorient_layout.addWidget(self.rotate_z_box)
        
        self.point_processing_box.setEnabled(False)#default
        
        #make outline processing box
        self.outliner_box = QtWidgets.QGroupBox('Outline processing')
        outliner_layout = QtWidgets.QGridLayout()
        self.outliner_box.setLayout(outliner_layout)
        self.alpha_cutoff = QtWidgets.QDoubleSpinBox()
        self.alpha_cutoff.setPrefix("Alpha: ")
        self.alpha_cutoff.setMinimum(0.001)
        self.alpha_cutoff.setMaximum(10000)
        self.alpha_cutoff.setDecimals(3)
        self.alpha_cutoff.setValue(1)
        self.alpha_cutoff.setToolTip('Alpha shape semiperimeter threshold for outline')
        self.find_alpha_outline_button = QtWidgets.QPushButton('Preview')
        self.find_alpha_outline_button.setToolTip('Generate outline using semiperimeter threshold')
        self.apply_outline_button = QtWidgets.QPushButton('Apply')
        self.apply_outline_button.setToolTip('Apply processed outline to active dataset and register')
        self.reset_outline_button = QtWidgets.QPushButton('Reset')
        self.reset_outline_button.setToolTip('Reset and deregister outline from current dataset')
        self.outline_status_label = QtWidgets.QLabel('Ready')
        
        #populate outline processing box
        outliner_layout.addWidget(self.alpha_cutoff,0,0,1,1)
        outliner_layout.addWidget(self.find_alpha_outline_button,0,1,1,1)
        outliner_layout.addWidget(self.apply_outline_button,0,2,1,1)
        outliner_layout.addWidget(self.reset_outline_button,0,3,1,1)
        outliner_layout.addWidget(self.outline_status_label,1,0,1,4)
        
        self.outliner_box.setEnabled(False) #default

        #make perimeter/outline box
        self.cut_orientation_box = QtWidgets.QGroupBox('Define cutting orientations')
        orientation_layout = QtWidgets.QGridLayout()
        self.cut_orientation_box.setLayout(orientation_layout)
        self.choose_cut_path_button = QtWidgets.QPushButton('Cut path')
        self.choose_cut_path_button.setToolTip('Specify orientation of the path the cut took through the outline')
        self.choose_cut_dir_button = QtWidgets.QPushButton('Cut direction')
        self.choose_cut_dir_button.setToolTip('Specify cutting direction')
        self.apply_cut_directions_button = QtWidgets.QPushButton('Apply')
        self.apply_cut_directions_button.setToolTip('Update viewer with cutting available cutting orientations')
        self.undo_cut_orient_button = QtWidgets.QPushButton('Reset')
        self.undo_cut_orient_button.setToolTip('Reset all cutting specifications')
        orientation_layout.addWidget(self.choose_cut_dir_button,0,0,1,1)
        orientation_layout.addWidget(self.choose_cut_path_button,0,1,1,1)
        orientation_layout.addWidget(self.apply_cut_directions_button,0,2,1,1)
        orientation_layout.addWidget(self.undo_cut_orient_button,0,3,1,1)
        
        self.cut_orientation_box.setEnabled(False) #default
        
        #make save box
        self.save_box = QtWidgets.QGroupBox('Write current')
        save_layout = QtWidgets.QGridLayout()
        self.save_box.setLayout(save_layout)
        #make save buttons
        self.save_button = QtWidgets.QPushButton("Save")
        self.save_button.setToolTip('Save data to hdf5-formatted results file')
        self.save_reference = QtWidgets.QRadioButton('Reference')
        self.save_reference.setChecked(True)
        self.save_reference.setToolTip('Specify active data to be written as the reference dataset')
        self.save_floating = QtWidgets.QRadioButton('Floating')
        self.save_floating.setToolTip('Specify active data to be written as the floating dataset')
        save_button_group = QtWidgets.QButtonGroup()
        save_button_group.addButton(self.save_reference)
        save_button_group.addButton(self.save_floating)
        save_button_group.setExclusive(True)
        self.save_label = QtWidgets.QLabel('Ready')
        self.save_label.setWordWrap(True)
        
        #populate save box
        save_layout.addWidget(self.save_reference,0,0,1,1)
        save_layout.addWidget(self.save_floating,0,1,1,1)
        save_layout.addWidget(self.save_button,0,2,1,1)
        save_layout.addWidget(self.save_label,1,0,1,3)
        
        self.save_box.setEnabled(False)#default
        
        #make vertical layout for interaction groupbox collection
        lvlayout=QtWidgets.QVBoxLayout()
        lvlayout.minimumSize()
        lvlayout.addWidget(import_box)
        lvlayout.addWidget(self.display_box)
        lvlayout.addWidget(self.point_processing_box)
        lvlayout.addWidget(self.reorient_box)
        lvlayout.addWidget(self.outliner_box)
        lvlayout.addWidget(self.cut_orientation_box)
        lvlayout.addWidget(self.save_box)
        lvlayout.addStretch(1)
        mainlayout.addWidget(self.vtkWidget)
        mainlayout.addStretch(1)
        mainlayout.addLayout(lvlayout)
        previewlayout.addWidget(self.preview_widget)
    

class interactor(QtWidgets.QWidget):
    '''
    Inherits most properties from Qwidget, but primes the VTK window, and ties functions and methods to interactors defined in main_window
    '''
    def __init__(self,parent):
        super(interactor, self).__init__(parent)
        self.ui = reg_main_window()
        self.ui.setup(self)
        
        #main (editor) vtk interactor settings
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(vtk.vtkNamedColors().GetColor3d("slategray"))
        self.ren.GradientBackgroundOn()
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ui.vtkWidget.Initialize()

        #default/initial values that are needed to initialize functions at first launch
        self.point_size=2
        self.limits=np.empty(6)
        self.zrange = None
        self.picking = False
        self.registered = False
        self.active_dir = os.getcwd()
        self.output_filename = None
        self.trans = np.eye(4) #default
        self.c_trans = np.eye(4) #initialise cumulative transformation
        
        
        #connect logic from Qt entities
        self.ui.load_spec.activated.connect(self.get_data)
        self.ui.entry_spec.activated.connect(self.change_entries)
        self.ui.color_on_height.clicked.connect(self.redraw)
        self.ui.z_aspect_sb.valueChanged.connect(self.update_z_aspect)
        self.ui.hide_display.clicked.connect(self.hide_point_cloud)
        self.ui.update_contours_button.clicked.connect(self.update_scalar_bar)
        self.ui.reset_point_display.clicked.connect(self.reset_display)
        self.ui.update_reduce_points_button.clicked.connect(self.decimate)
        self.ui.triangulate_button.clicked.connect(self.actuate_triangulate)
        self.ui.apply_reduce_points_button.clicked.connect(self.apply_decimate)
        self.ui.undo_reduce_points_button.clicked.connect(self.reset_decimate)
        self.ui.undo_last_pick_button.clicked.connect(self.undo_pick)
        self.ui.apply_pick_button.clicked.connect(self.apply_decimate)
        self.ui.level_button.clicked.connect(self.trans_mean_z)
        self.ui.centroid_move_button.clicked.connect(self.trans_centroid)
        self.ui.activate_svd_button.clicked.connect(self.svd)
        self.ui.impose_rotation_button.clicked.connect(self.rotate)
        self.ui.auto_rotate_button.clicked.connect(self.auto_rotate)
        self.ui.find_alpha_outline_button.clicked.connect(self.preview_outline)
        self.ui.apply_outline_button.clicked.connect(self.apply_outline)
        self.ui.reset_outline_button.clicked.connect(self.reset_outline)
        self.ui.choose_cut_path_button.clicked.connect(self.choose_path)
        self.ui.choose_cut_dir_button.clicked.connect(self.chose_dir)
        self.ui.apply_cut_directions_button.clicked.connect(self.draw_active_directions)
        self.ui.undo_cut_orient_button.clicked.connect(self.reset_dir)
        self.ui.save_button.clicked.connect(self.write_output)
        self.ui.tab_widget.currentChanged[int].connect(self.update_preview)

    def closeEvent(self,event):
        '''
        Finalize all VTK widgets to negate OpenGL messages/errors
        '''
        self.ui.vtkWidget.close()
        self.ui.preview_widget.vtkWidget.close()

    def change_entries(self):
        '''
        Manages entries for each side
        '''
        entry_count = self.ui.entry_spec.count()
        if self.ui.entry_spec.currentText() == 'New entry':
            self.insert_entries(entry_count + 1)
            self.ui.entry_spec.setCurrentIndex(entry_count-2)
        if self.ui.entry_spec.currentText() == 'Clear entries':
            self.ui.entry_spec.clear()
            self.ui.entry_spec.addItems(['Entry 0','New entry', 'Clear entries'])
            self.ui.entry_spec.setCurrentIndex(0)
        if self.ui.entry_spec.currentIndex() == 0:
            self.ui.reorient_box.setEnabled(True)
            self.ui.cut_orientation_box.setEnabled(True)
        else:
            self.ui.reorient_box.setEnabled(False)
            self.ui.cut_orientation_box.setEnabled(False)

    def insert_entries(self,entry):
        #inserts additional entries, up to and including the 'entry' value if the entry value exceeds the number of values in entry_spec
        
        entry_count = self.ui.entry_spec.count()
        
        req_entries = entry + 3 #0th, new entry, clear entries
        if (req_entries - entry_count) > 0:
            for new_entry in range(req_entries - entry_count):
                self.ui.entry_spec.insertItem(entry_count+new_entry-2,"Entry %d"%(entry_count+new_entry-2))
        

    def get_data(self):
        '''
        Calls get_file, loads that data and then calls draw on the basis of what was read, and activates parts of the ui that are available
        '''
        #call get_file from current working directory. If user selects different directory, change the current working directory to be there.
        file, self.active_dir = get_file('*.txt',self.active_dir)
        
        #check to make sure that a valid file was selected:
        if file is None or not(os.path.isfile(file)):
            self.ui.load_label.setText('Nothing loaded.')
            return
            
        self.ui.load_label.setText('Loading . . .')
        QtWidgets.QApplication.processEvents()
        
        data = read_data(file)
        
        if self.ui.entry_spec.currentIndex() == 0:
            #initialize transformation for this dataset
            if not np.array_equal(self.c_trans, np.eye(4)):
                reset_transform = warning_msg(self, \
                'Re-initialize transformations for this entry?')
                if not reset_transform:
                    self.ui.load_label.setText('Nothing loaded.')
                    return
            self.trans = np.eye(4) #default
            self.c_trans = np.eye(4) #initialise
        
        self.ui.outliner_box.setEnabled(False)
        
        #check what the user is trying to load and call draw_pc on what was read with the appropriate call to draw_pc.
        ls = self.ui.load_spec.currentIndex()
        if ls == 0:
            self.outline = data
            self.outline = do_transform(self.outline,self.c_trans)
            self.draw(self.outline,ls)
            self.registered = True
            self.ui.outline_offset.setEnabled(True)
            self.ui.reorient_box.setEnabled(False)
        if ls == 2:
            self.reset_outline()
            self.ui.outliner_box.setEnabled(True)
        if ls > 0:
            if self.ui.entry_spec.currentIndex() == 0:
                self.ui.reorient_box.setEnabled(True)
                self.ui.cut_orientation_box.setEnabled(True)
            self.raw_pts = data
            self.raw_pts = do_transform(self.raw_pts,self.c_trans)
            self.bool_pnt = np.ones(len(data), dtype=bool)
            self.active_pnt = np.arange(0, len(self.bool_pnt), 1, dtype=int)
            #turn on ui groups
            self.ui.display_box.setEnabled(True)
            self.ui.point_processing_box.setEnabled(True)
            self.ui.save_box.setEnabled(True)
            
            #reset/initialise ui
            self.reset_display() #calls draw with defaults on
            self.reset_decimate()
            self.reset_dir()
        
        self.ui.load_label.setText('Loaded %s'%os.path.basename(file))
        
    def reset_display(self):
        '''
        Resets all display associated ui according to defaults
        '''
        self.limits = get_limits(self.raw_pts[self.active_pnt])
        r = (self.limits[-2],self.limits[-1]) #min/max
        self.ui.min_contour.setValue(r[0])
        self.ui.max_contour.setValue(r[1])
        self.zrange = r
        self.redraw()
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
        
    def redraw(self):
        '''
        Checks the status of the triggering button, and calls draw with the appropriate option/mode. Updates bool_pnt with status of active_pnt.
        '''
        if self.ui.color_on_height.isChecked():
            draw_mode = 0
            self.ui.contour_tools.setEnabled(True)
        else:
            draw_mode = 1
            self.ui.contour_tools.setEnabled(False)
        
        if self.registered:
            draw_option = 1
        else:
            draw_option = 2
        
        self.draw(self.raw_pts[self.active_pnt], draw_option, draw_mode)
        #hide point cloud if button is selected
        if self.ui.hide_display.isChecked():
            self.hide_point_cloud()
        if self.registered:
            self.draw(self.outline, 0)
    
    def decimate(self):
        '''
        Decimates active points according what radio buttons have been selected
        '''
        self.ui.filter_status_label.setText('Applying filters . . .')
        QtWidgets.QApplication.processEvents()
        #Call reduce_pnts with conditions of ui & active points
        active_pnts = self.raw_pts[self.active_pnt]
        
        # ind will return points that remain active
        if self.ui.decimate_by_number.isChecked():
            try:
                ind = reduce_pnts(active_pnts,\
                self.ui.num_active_points.value(),\
                self.ui.z_active_points.value()\
                ,0)
            except: pass
        if self.ui.decimate_by_percent.isChecked():
            try:
                ind = reduce_pnts(active_pnts,\
                self.ui.percent_active_points.value(),\
                self.ui.z_active_points.value(),\
                1)
            except: pass
        
        #if the Z norm cutoff is active
        if self.ui.z_norm_active_points.isEnabled():
            z_norm_ind = self.decimate_by_tri_normal()
            #get intersection of in and the z_norm_ind
            ind = np.intersect1d(ind,z_norm_ind)
        
        if hasattr(self,'outline'):
            if self.outline is not None:
                offset_outline = offset_poly(self.outline,self.ui.outline_offset.value())
                in_poly_bool = in_poly(offset_outline,active_pnts)
                in_poly_ind = np.arange(0, len(active_pnts), 1, dtype=int)
                ind = np.intersect1d(ind,in_poly_ind[in_poly_bool])
        
        #update ui based on what is returned by reduce_pnts
        self.ui.num_active_points.setValue(len(ind))
        self.ui.percent_active_points.setValue(len(ind)/len(self.raw_pts)*100)
        
        #ind will be an index of points to keep, so update self.bool_pnt
        self.bool_pnt[ind] = False #turn those returned false
        self.bool_pnt = np.logical_not(self.bool_pnt) #negate, make true
        
        #update display color of points
        self.update_mask()
        self.ui.apply_reduce_points_button.setEnabled(True)
        self.ui.display_box.setEnabled(False)
        
        self.ui.filter_status_label.setText('Generated preview.')
        
    def decimate_by_tri_normal(self):
        '''
        Removes points whose triangulated normal falls outside of a threshold
        '''
        def tri_normal_z(points,tri):
            '''
            Returns an array of normalized z components of each triangle in a 2D triangulation, tri - an index of points
            '''
            triangles = points[tri.vertices]
            V = triangles[:,2,:] - triangles[:,0,:]
            U = triangles[:,1,:] - triangles[:,0,:]
            n = np.cross(U,V)
            
            mag = np.sqrt((n ** 2).sum(-1)) #inner1d can make faster
            z_norm = np.abs(n[:,-1])/mag
            
            return z_norm
        
        #an array of the z component of each normal for each triangle in self.tri:
        normals = tri_normal_z(self.raw_pts[self.active_pnt,:],self.tri)
        #find triangles that are above the cutoff:
        cutoff = self.ui.z_norm_active_points.value()
        norm_cutoff_ind = normals > cutoff
        #find the unique indices that are in the triangle index to keep:
        ind = np.unique(self.tri.vertices[norm_cutoff_ind,:].copy().flatten())
        return ind
        
    def actuate_triangulate(self):
        '''
        Starts triangulation and handles ui display, or deactivates if it's currently set
        '''
        if self.ui.z_norm_active_points.isEnabled():
            self.ui.triangulated_indicator.setStyleSheet("background-color : gray; color : darkGray;")
            self.ui.z_norm_active_points.setEnabled(False)
            del self.tri
        else:
            self.triangulate()
            self.ui.triangulated_indicator.setStyleSheet("background-color :rgb(77, 209, 97);")
            self.ui.z_norm_active_points.setEnabled(True)
    
    def triangulate(self):
        '''
        Performs a 2D Delaunay triangulation on active points and returns it to self.tri.
        '''
        self.ui.filter_status_label.setText('Triangulating')
        QtWidgets.QApplication.processEvents()
        self.tri = Delaunay(self.raw_pts[self.active_pnt,0:2])
        self.ui.filter_status_label.setText('Triangulation complete.')
        
    def reset_decimate(self):
        '''
        Reset all filters such that all raw points are active
        '''
        self.ui.filter_status_label.setText('Resetting all filters . . .')
        QtWidgets.QApplication.processEvents()
        self.ui.num_active_points.setValue(len(self.raw_pts))
        self.ui.percent_active_points.setValue(100)
        #re-initialize self.bool_pnt and active_pnt
        self.bool_pnt=np.ones(len(self.raw_pts), dtype=bool)
        self.active_pnt = np.arange(0, len(self.bool_pnt), 1, dtype=int)
        
        self.redraw()
        self.reset_display()
        self.update_decimate()

        self.ui.apply_reduce_points_button.setEnabled(False)
        self.ui.update_reduce_points_button.setEnabled(True)
        self.ui.display_box.setEnabled(True)

        self.ui.filter_status_label.setText('Ready')
        
    def update_decimate(self):
        # update the current colour and decimation filter
        self.ui.num_active_points.setValue(len(self.raw_pts[self.active_pnt]))
        self.ui.percent_active_points.setValue(\
        len(self.raw_pts[self.active_pnt]) / \
        len(self.raw_pts)*100)
        if self.ui.z_norm_active_points.isEnabled():
            #handle triangulation
            self.actuate_triangulate()
        self.ui.z_active_points.setValue(self.zrange[0]) #zrange updated by reset_display
        
    def apply_decimate(self):
        '''
        calls redraw, permanently change both raw_pts and bool_pnt if specified
        '''
        self.ui.filter_status_label.setText('Applying filters . . .')
        QtWidgets.QApplication.processEvents()
        #update active_pnt for redraw
        self.active_pnt = self.active_pnt[self.bool_pnt]

        if not self.ui.local_decimate.isChecked():
            #permanently remove points
            ind = self.bool_pnt
            self.raw_pts = self.raw_pts[ind,:]
            self.bool_pnt = self.bool_pnt[ind]
        self.ui.display_box.setEnabled(True)
        # reset bool_pnt
        self.bool_pnt = np.ones(len(self.active_pnt), dtype=bool) #needs to stay the same length as self.color
        self.reset_display()
        self.update_decimate()
        self.ui.apply_reduce_points_button.setEnabled(False) #until preview activates it again
        self.ui.filter_status_label.setText('Ready')
    
    
    def actuate_pick(self):
        '''
        Starts picking and handles ui button display
        '''
        if self.picking:
            self.picking = False
            self.ui.active_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
        else:
            self.picking = True
            self.ui.active_picking_indicator.setStyleSheet("background-color :rgb(77, 209, 97);")
            self.start_pick()

    def start_pick(self):
        '''
        Activates picking and adds the picker observer to the interactor
        '''
        style=vtk.vtkInteractorStyleRubberBandPick()
        self.iren.SetInteractorStyle(style)
        picker = vtk.vtkAreaPicker()
        self.iren.SetPicker(picker)
        picker.AddObserver("EndPickEvent", self.picker_callback)

    def undo_pick(self):
        '''
        Undoes the last pick registered
        '''
        if hasattr(self,"last_selected_ids"):
            #update colours
            for i in range(self.last_selected_ids.GetNumberOfTuples()):
                self.colors.SetTuple(self.last_selected_ids.GetValue(i),self.last_selected_colors[i])
                self.bool_pnt[self.last_selected_ids.GetValue(i)] = True
            
            self.pnt_polydata.GetPointData().SetScalars(self.colors)
            self.pnt_polydata.Modified()
            self.point_actor.Modified()
            self.ui.vtkWidget.update()
            del self.last_selected_ids, self.last_selected_colors
        else:
            return

    def picker_callback(self,obj,event):
        '''
        Manual picking callback function
        '''
        extract = vtk.vtkExtractSelectedFrustum()
        f_planes=obj.GetFrustum() #collection of planes based on unscaled display
        z_aspect = self.ui.z_aspect_sb.value()
        #scale frustum to account for the z_aspect
        scaled_planes=vtk.vtkPlanes()
        scaled_normals=vtk.vtkDoubleArray()
        scaled_normals.SetNumberOfComponents(3)
        scaled_normals.SetNumberOfTuples(6)
        scaled_origins=vtk.vtkPoints()
        for j in range(6):
            i=f_planes.GetPlane(j)
            k=i.GetOrigin()
            q=i.GetNormal()
            scaled_origins.InsertNextPoint(k[0],k[1],k[2]/float(z_aspect))
            scaled_normals.SetTuple(j,(q[0],q[1],q[2]*float(z_aspect)))
        scaled_planes.SetNormals(scaled_normals)
        scaled_planes.SetPoints(scaled_origins)

        extract.SetFrustum(scaled_planes)
        extract.SetInputData(self.pnt_polydata)
        extract.Update()
        extracted = extract.GetOutput()

        ids = vtk.vtkIdTypeArray()
        ids = extracted.GetPointData().GetArray("vtkOriginalPointIds")

        if ids:
            self.last_selected_colors = []
            for i in range(ids.GetNumberOfTuples()):
                self.last_selected_colors.append(self.colors.GetTuple(i))
                self.bool_pnt[ids.GetValue(i)] = False
            #store them in an array for an undo operation, and update self.bool_pnt
            self.last_selected_ids=ids
            self.update_mask()

    def hide_point_cloud(self):
        '''
        Changes the visibility of anything related to point clouds on the display
        '''
        flip_visible(self.point_actor)
        if hasattr(self,'sb_actor'):
            flip_visible(self.sb_actor)
        self.ui.vtkWidget.update()
        if self.ui.contour_tools.isEnabled() and not self.ui.color_on_height.isChecked():
            self.ui.contour_tools.setEnabled(False)

    def draw(self, data, option, mode=0):
        '''
        draws a point cloud-type of interactor using supplied data with one of the following options:
        Option 0 - draws registered perimeter
        Option 1 - draws registered point cloud
        Option 2 - draws unregistered point cloud, clears the renderer first
        mode 0 - colour points based on height
        mode 1 - solid colour
        '''
        if mode == 0:
            # color = (70, 171, 176) #will be overwritten
            color = (255, 255, 255)
        else:
            color = (153,204,255)#default colour
        if option == 0:
            if hasattr(self,'outline_actor'):
                self.ren.RemoveActor(self.outline_actor)
            self.outline_actor, _ = gen_outline(data,color,self.point_size)
            self.ren.AddActor(self.outline_actor)
        if option == 2:
            self.ren.RemoveAllViewProps()
        if option > 0: #not an outline; draw a point cloud
            #if there's already a point cloud actor and scalebar
            if hasattr(self,'point_actor'):
                self.ren.RemoveActor(self.point_actor)
                self.ren.RemoveActor(self.axis_actor)
                if hasattr(self,'sb_actor'):
                    self.ren.RemoveActor(self.sb_actor)
            if mode == 0:
                if self.zrange is None:
                    self.point_actor, \
                    self.pnt_polydata, \
                    colors, lut = \
                    gen_point_cloud(data)
                    self.ui.min_contour.setValue(self.limits[-2])
                    self.ui.max_contour.setValue(self.limits[-1])
                    self.zrange = (self.limits[-2],self.limits[-1])
                else:
                    self.point_actor, \
                    self.pnt_polydata, \
                    self.colors, lut = \
                    gen_point_cloud(data,None,self.zrange)
                #handle scalebar
                sb_widget = gen_scalar_bar()
                sb_widget.SetInteractor(self.iren)
                sb_widget.On()
                self.sb_actor = sb_widget.GetScalarBarActor()
                self.sb_actor.SetLookupTable(lut)
                self.ren.AddActor(self.sb_actor)
            elif mode == 1:
                self.point_actor, \
                self.pnt_polydata, \
                self.colors, lut = \
                gen_point_cloud(data,color)

            self.update_z_aspect()
            self.ren.AddActor(self.point_actor)

        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
    
    def update_mask(self):
        '''
        Updates the coloration of poly data on the basis of what's been masked.
        '''
        localind = np.arange(0, len(self.bool_pnt), 1, dtype=int)
        localind = localind[np.where(np.logical_not(self.bool_pnt))]
        for i in localind:
            self.colors.SetTuple(i,(255,255,255))

        self.pnt_polydata.GetPointData().SetScalars(self.colors)
        self.pnt_polydata.Modified()
        self.point_actor.Modified()
        self.ui.vtkWidget.update()

    def update_scalar_bar(self):
        '''
        calls redraw with a scalarbar range
        '''
        r = (self.ui.min_contour.value(),self.ui.max_contour.value())
        self.zrange = r
        self.redraw()
        self.sb_actor.SetNumberOfLabels(self.ui.num_contour.value())
        self.ui.vtkWidget.update()

    
    def update_z_aspect(self):
        '''
        Updates z_aspect and redraws points displayed based on new z_aspect
        '''
        z_aspect = self.ui.z_aspect_sb.value()
        #scale point cloud
        self.point_actor.SetScale((1,1,z_aspect))
        self.point_actor.Modified()
        #now do axis_actor, can't scale in the same way as polydata . . .
        try:
            self.ren.RemoveActor(self.axis_actor) #it will need to be replaced
        except: pass
        #local limits for z axis - don't scale x & y limits, scale the z axis according to z aspect and the current limits
        nl=np.append(self.limits[0:4],([self.limits[-2]*z_aspect,self.limits[-1]*z_aspect]))
        self.axis_actor = get_axis(self.ren, nl, z_aspect)
        self.ren.AddActor(self.axis_actor)
        self.ui.vtkWidget.update()

    def trans_mean_z(self):
        '''
        Translate points to mean of z value of all points
        '''
        if hasattr(self,'outline'):
            a = np.vstack((self.raw_pts[self.active_pnt,:],self.outline))
            #update transformation matrix to move back with local matrix
            self.trans[2,-1] = - np.mean(a[:,2])
            self.apply_transformation()

    def trans_centroid(self):
        '''
        Translates active points to centroid of point cloud and outline
        '''
        if hasattr(self,'outline'):
            a = np.vstack((self.raw_pts[self.active_pnt,:],self.outline))
            self.trans[0:3,-1] = - np.mean(a, axis=0)
            self.apply_transformation()

    def svd(self):
        '''
        Get svd of active points, draw relevant vectors
        '''
        #only use active points for calculating SVD, but for underlying calc, the size will be length of active points by length of active points
        max_svd_pts = 10000
        svd_pts = self.raw_pts[self.active_pnt,:]
        if len(svd_pts) > max_svd_pts:
            ind = reduce_pnts(svd_pts,max_svd_pts)
            svd_pts = svd_pts[ind]

        #get centroid
        c = np.mean(svd_pts, axis = 0)

        #get plane of best fit normal and rotation matrix to get normal to align to z axis
        R = get_svd_orientation(svd_pts)

        self.trans[0:3,0:3] = R
        self.trans[0:3,-1] = -c

        #translate svd_pts to the centroid, rotate, and then translate back
        svd_pts = svd_pts - c
        svd_pts = np.dot(svd_pts, R)
        svd_pts = svd_pts + c

        #update transformation matrix to move back with local matrix
        m = np.eye(4)
        m[0:3,-1] = c
        self.trans = m @ self.trans

        self.apply_transformation()

    
    def rotate(self):
        '''
        Manually rotates data about z axis based on input from ui
        '''
        a = np.deg2rad(self.ui.rotate_z.value()) #negative for counterclockwise
        self.trans[0:2,0:2]=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        self.apply_transformation()
    
    def auto_rotate(self):
        '''
        Tries to find corners of the outline, and aligns longest side to x axis
        '''
        if not hasattr(self,'outline'):
            return
        #reorder outline and get corner index
        corner_ind, self.outline = get_corner_ind(self.outline)

        #get coordinates of corners
        corners = self.outline[corner_ind,:]

        #calculate side lengths - follow standard 2D element face numbering
        s1 = corners[1,:] - corners[0,:]
        s2 = corners[2,:] - corners[1,:]
        s3 = corners[3,:] - corners[2,:]
        s4 = corners[0,:] - corners[3,:]
        s = np.vstack((s1,s2,s3,s4))
        mag = np.sqrt((s*s).sum(axis=1))
        #find u (x axis, longest) 
        u = s[mag == np.amax(mag),:][0]
        u = u/np.linalg.norm(u) #normalize
        #v vector will be the cross product of u and z axis
        v = np.cross(u,[0,0,1])
        v = v/np.linalg.norm(v)#normalize

        #make rotation matrix & apply
        self.trans[0:2,0:2] = np.array([[u[0],v[0]],[u[1],v[1]]] )
        self.apply_transformation()

    def apply_transformation(self):
        '''
        Apply/update transformation
        '''
        T = self.trans.copy()
        self.raw_pts = do_transform(self.raw_pts,T)
        if hasattr(self,'outline'):
            self.outline = do_transform(self.outline,T)
        self.c_trans = T @ self.c_trans #pre-multiply, same as np.dot(self.c_trans,T.T)

        self.trans = np.eye(4)#reset the local transformation matrix

        self.reset_display()
        self.update_decimate()
        self.reset_dir()

    def preview_outline(self):
        '''
        Generates an outline on the basis of the semiperimeter set in the ui
        '''
        self.ui.outline_status_label.setText('Processing outline . . .')
        QtWidgets.QApplication.processEvents()

        #if there isn't an active triangulation, get one.
        if not hasattr(self,'tri'):
            self.actuate_triangulate()

        #calculate the hull with the shapely alpha shape routine
        self.outline = alpha_shape(\
        self.raw_pts[self.active_pnt,:], \
        self.tri,\
        self.ui.alpha_cutoff.value())

        if self.outline is not None:
            #draw it with a different color and point size than a 'registered' outline
            if hasattr(self,'outline_actor'):
                    self.ren.RemoveActor(self.outline_actor)
            self.outline_actor, _ = gen_outline(self.outline,(255,127,80),self.point_size*2)
            self.ren.AddActor(self.outline_actor)
            self.ui.vtkWidget.update()
            self.ui.outline_status_label.setText('Generated preview.')
        else:
            self.ui.outline_status_label.setText('Outline processing failed.')

    def apply_outline(self):
        '''
        Accepts the outline by calling redraw and changing the registered state
        '''
        if not hasattr(self,'outline'):
            return
        self.registered = True
        self.ui.rotate_z_box.setEnabled(True)
        self.ui.outline_offset.setEnabled(True)
        self.reset_display() #

    def reset_outline(self):
        '''
        Removes the outline actor if it exists and changes the registered state
        '''
        if not hasattr(self,'outline'):
            return
        else:
            del self.outline
        self.registered = False
        self.actuate_triangulate()
        self.ui.rotate_z_box.setEnabled(False)
        self.ui.outline_offset.setEnabled(False)
        if hasattr(self,'outline_actor'):
            self.ren.RemoveActor(self.outline_actor)
            
        self.ui.vtkWidget.update()

    def choose_path(self):
        '''
        Permits selection of the cutting path
        '''
        self.draw_selection_directions()
        self.current_arrow_color = vtk.vtkNamedColors().GetColor3d("LightSalmon")
        self.iren.AddObserver("EndPickEvent",self.check_pick)
        self.ui.choose_cut_path_button.setEnabled(False)

    def chose_dir(self):
        '''
        Permits the selection of cutting direction
        '''
        self.draw_selection_directions()
        self.current_arrow_color = vtk.vtkNamedColors().GetColor3d("Tomato")
        self.iren.AddObserver("EndPickEvent",self.check_pick)
        self.ui.choose_cut_dir_button.setEnabled(False)

    def check_pick(self,object,event):
        """
        Activates a pick event for arrows
        """
        self.picked = None
        self.ui.vtkWidget.setFocus()
        picker=vtk.vtkPropPicker()

        click_pos=self.iren.GetEventPosition()

        picker.Pick(click_pos[0],click_pos[1],0,self.ren)
        picked_actor = picker.GetActor()
        count=0
        if picked_actor:
            self.picked = identify_actor(picked_actor,self.a)
            self.a[self.picked].GetProperty().SetColor(self.current_arrow_color)
            
            #remove remaining actors
            mask = np.full(len(self.a), True, dtype=bool)
            mask[self.picked] = False
            ind = np.arange(len(self.a),dtype=int)[mask]
            for i in ind:
                self.ren.RemoveActor(self.a[i])
        
            if tuple(self.current_arrow_color) == tuple(vtk.vtkNamedColors().GetColor3d("LightSalmon")):
                self.cut_path = self.directions[self.picked,:]
                self.cut_path_arrow = self.a[self.picked]
            
            if tuple(self.current_arrow_color) == tuple(vtk.vtkNamedColors().GetColor3d("Tomato")):
                self.cut_dir = self.directions[self.picked,:]
                self.cut_dir_arrow = self.a[self.picked]
        
        self.iren.RemoveObservers("EndPickEvent")
        self.ui.vtkWidget.update()
    
    def draw_active_directions(self):
        '''
        Draws direction actors based on the presence of self.cut_dir and self.cut_path. Called on load/display reset
        '''
        #if there's a cutting direction
        if not hasattr(self,'cut_dir'):
            return
        
        if hasattr(self,'cut_orient_actor'):
            self.ren.RemoveActor(self.cut_orient_actor)
        
        if hasattr(self,'cut_path'):
            self.cut_orient_actor = gen_cutting_orientation_actor(\
            self.limits,self.cut_dir,self.cut_path)
        else:
            self.cut_orient_actor = gen_cutting_orientation_actor(\
            self.limits,self.cut_dir)

        self.ren.AddActor(self.cut_orient_actor)

        #remove selection arrows
        if hasattr(self,'cut_path_arrow'):
            self.ren.RemoveActor(self.cut_path_arrow)

        if hasattr(self,'cut_dir_arrow'):
            self.ren.RemoveActor(self.cut_dir_arrow)
            
        self.ui.vtkWidget.update()

    def reset_dir(self):
        '''
        Resets all directions
        '''
        if hasattr(self,'a'):
            #remove all selectable actors that might be drawn
            for a in self.a:
                self.ren.RemoveActor(a)
        
        if hasattr(self,'picked'):
            del self.picked
            del self.a
        if hasattr(self,'cut_path_arrow'):
            self.ren.RemoveActor(self.cut_path_arrow)
        
        if hasattr(self,'cut_dir_arrow'):
            self.ren.RemoveActor(self.cut_dir_arrow)
            
        if hasattr(self,'cut_orient_actor'):
            self.ren.RemoveActor(self.cut_orient_actor)
        
        #delete existing directions
        try: 
            del self.cut_dir
        except: pass
        
        try: del self.cut_path
        except: pass
        
        self.ui.choose_cut_dir_button.setEnabled(True)
        self.ui.choose_cut_path_button.setEnabled(True)
        self.ui.vtkWidget.update()

    def draw_selection_directions(self):
        '''
        Draws a series of selectable arrows for specifying directions. Specific directions are dependant on two cases:
        1) No directions have been selected
        2) Either the cutting path or direction have been specified, so show directions that apply to those that have not been selected
        '''
        #get size of arrows
        asize=np.maximum(self.limits[1]-self.limits[0],self.limits[3]-self.limits[2])*0.10
                
        #up, down, left, right
        self.directions = np.array([
        [0,1,0], 
        [0,-1,0],
        [-1,0,0],
        [1,0,0]
        ])
        start = np.array([
        [0,asize,0],
        [0,-asize,0],
        [-asize,0,0],
        [asize,0,0]
        ])
        
        mask = np.full(len(self.directions), True, dtype=bool)

        if hasattr(self,'picked'):
            if self.picked is not None:
                if self.picked <=1:
                    #only show the remaining directions
                    mask[:2] = False
                else:
                    mask[-2:] = False
        
        self.directions = self.directions[mask,:]
        start = start[mask,:]
        
        #translate starts to centroid of bounding box
        cent = np.array([
        (self.limits[1]+self.limits[0])/2,
        (self.limits[3]+self.limits[2])/2,
        0
        ])
        
        start = start + cent
        
        #create a list of arrow actors
        self.a = []
        for i in range(len(self.directions)):
             self.a.append(draw_arrow(start[i,:],\
             asize,\
             self.directions[i,:],\
             self.ren,\
             False)[0])
        
        #update
        self.ui.vtkWidget.update()

    def check_save_state(self, id, entry):
        '''
        1) Checks to make sure there's an active outline and point cloud, warn user if there isn't.
        2) In the instance that there is existing data according to 'id' in the current save file, warn the user that it will be over-written.
        Returns false if checks fail
        '''
        if not self.registered:
            info_msg('Saving the current data requires both an outline and an outlined surface. Either import an outline for the current point data or process one from an unregistered point cloud.')
            return False
        
        #check contents of output_filename:
        if self.output_filename is None:
            return False
        
        with h5py.File(self.output_filename, 'r') as f:
            existing_keys = list(f.keys())
            if not existing_keys or id not in existing_keys:
                return True
        
        with h5py.File(self.output_filename, 'r') as f:
            #check the id's entry keys
            g = f['%s'%id]
            if entry in list(g.keys()):
                overwrite = warning_msg(self, \
                'There is existing data in the specified file for the target field, overwrite?'
                )
                if not overwrite: #get a new file to save to
                    new_save_filename, _ = get_save_file('*.pyCM')
                    if new_save_filename is None:
                        return False
                    if new_save_filename == self.output_filename:
                        info_msg('Currently active save file selected; specify another.')
                    else:
                        self.output_filename = initialize_HDF5(new_save_filename)
                        return True
                else:
                    #overwrite
                    return True
            else:
                #there isn't an entry in the id
                return True

    def write_output(self):
        '''
        Writes output to a file if it exists, gets a save location and initializes a save file, checking first with self.check_save_state if any data will be overwritten.
        '''
        if self.ui.save_reference.isChecked():
            id = 'ref'
            name = 'reference'
        else:
            id = 'float' #button group is exclusive
            name = 'floating'
        #get entry & cast as string as it's going to be a key
        entry = str(self.ui.entry_spec.currentIndex())
        
        #initialized as None, and initialize_HDF5 returns None if no file is selected
        if self.output_filename is None:
            self.output_filename = initialize_HDF5()
        
        passed = self.check_save_state(id, entry)
        
        if not passed:
            self.ui.save_label.setText('Ready')
            return
        else:
            #write data - if it doesn't exist, create entries, otherwise assign values to existing
            self.ui.save_label.setText('Saving . . .')
            QtWidgets.QApplication.processEvents()
            entry_dict = {'points':self.raw_pts, 'outline':self.outline, 'active':self.active_pnt}

            with h5py.File(self.output_filename, 'r+') as f:
                if entry == '0':
                    if id in list(f.keys()):
                        g = f['%s'%id]
                        del g['transform']
                    else:
                        g = f.create_group("%s"%id)
                    g.create_dataset('transform', data = self.c_trans)
                    if hasattr(self,'cut_dir'):
                            f['%s'%id].attrs['cut_dir'] = self.cut_dir 
                    if hasattr(self,'cut_path'):
                            f['%s'%id].attrs['cut_path'] = self.cut_path

                g = f['%s'%id]
                if entry in list(g.keys()):
                    gg = g[entry]
                    for k, v in entry_dict.items():
                        del gg[k]
                        gg.create_dataset(k, data = v)
                        # data = gg[k]
                        # data[()] = v
                else:
                    gg = g.create_group(entry)
                    for k, v in entry_dict.items():
                        gg.create_dataset(k, data = v)

                f.attrs['date_modified'] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                #write paths as attributes

        self.ui.save_label.setText('Saved %s entry %s to %s'%(name,entry,os.path.basename(self.output_filename)))

    def keypress(self,obj,event):
        '''
        VTK interactor specific keypress binding from the editor tab
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
        elif key == "a":
            flip_visible(self.axis_actor)
        elif key == "r":
            self.actuate_pick()
        self.ui.vtkWidget.update()

    def update_preview(self, value):
        '''
        Reads the output file and (re)draws contents in the preview tab
        If the prewiew tab has a selected dataset, call load_from_preview with it
        '''

        if self.output_filename is None:
            return
        self.ui.preview_widget.ref_pnts, \
        self.ui.preview_widget.ref_outlines, \
        self.ui.preview_widget.float_pnts, \
        self.ui.preview_widget.float_outlines, \
        self.ui.preview_widget.cut_attr = read_file(self.output_filename)
        
        self.ui.preview_widget.draw()

        #check if something has been picked from the preview_widget, then 'picked' will have entries
        if value == 0:
            if self.ui.preview_widget.picked:
                self.load_from_preview(self.ui.preview_widget.picked)
                self.ui.preview_widget.picked = {}

    def load_from_preview(self,sel):
        
        #check if there's anything currently being displayed
        if hasattr(self,'outline_actor') or hasattr(self,'point_actor'):
            clear_display = warning_msg(self, \
                'Clear data currently loaded in the editor to load the selected dataset?'
                )
            if not clear_display:
                return
        
        self.ren.RemoveAllViewProps()
        
        #check if sel is either 'Ref' or 'Float'
        if list(sel.keys())[0] == 'Ref':
            entry = sel['Ref']
            ref_pnts, ref_outlines, _, _, \
            cut_attr,ref_active_list, _, ref_trans, _ \
            = read_file(self.output_filename, return_active = True)
            self.c_trans = ref_trans
            self.outline = ref_outlines[entry]
            self.raw_pts = ref_pnts[entry]
            active_pnt = ref_active_list[entry]

            self.ui.save_reference.setChecked(True)
        else:
            entry = sel['Float']
            _, _, float_pnts, float_outlines, \
            cut_attr, _, float_active_list, _, float_trans \
            = read_file(self.output_filename, return_active = True)
            self.c_trans = float_trans
            self.outline = float_outlines[entry]
            self.raw_pts = float_pnts[entry]
            active_pnt = float_active_list[entry]
            
            self.ui.save_floating.setChecked(True)
            
        self.bool_pnt = np.zeros(len(self.raw_pts), dtype=bool)
        self.bool_pnt[active_pnt] = True
        self.active_pnt = np.arange(0, len(self.bool_pnt), 1, dtype=int)
        self.registered = True
        self.apply_decimate()
        
        #manage entries
        self.insert_entries(entry)
        self.ui.entry_spec.setCurrentIndex(entry)
        self.change_entries() #to enable/disable cut orientation
        
        #turn on decimation, reorientation, outline processing & saving
        self.ui.point_processing_box.setEnabled(True)
        self.ui.save_box.setEnabled(True)
        self.ui.reorient_box.setEnabled(True)
        self.ui.outliner_box.setEnabled(True)
        

def get_svd_orientation(points):
    '''
    Returns a 3x3 rotation matrix required to take z component of the orthonormal matrix of points to either 0,0,1 or 0,0,-1 at the centroid of points, depending on concavity.
    '''
    #get singular vectors, v, with an origin of the centroid of points
    _,_,v = np.linalg.svd(points - np.mean(points, axis=0))
    #normal of plane of best fit is the right singular vector
    normal = v[2]
    
    #handles the case if the dataset is net convex vs. concave relative to +z
    if normal[2] < 0:
        target = np.array([0,0,-1])
    else: 
        target = np.array([0,0,1])
        
    #solve for angle and axis between normal and target
    angle = np.arccos(np.dot(target,normal))
    axis = np.cross(target,normal)
    axis = axis / np.linalg.norm(axis)#normalize

    #convenience variables for the rotation matrix
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    ux = axis[0]
    uy = axis[1]
    uz = axis[2]

    R = np.array([[c+ux**2*t, ux*uy*t-uz*s, ux*uz*t+uy*s],
    [uy*ux*t+uz*s, c+uy**2*t, uy*uz*t-ux*s],
    [uz*ux*t-uy*s, uz*uy*t+ux*s, c+uz**2*t]])
    
    #return the inverted rotation for pre-multiplication
    return np.linalg.inv(R)

def tri_normal_z(points,tri):
    '''
    Returns an array of normalized z components of each triangle in a 2D triangulation, tri - an index of points
    '''
    triangles = points[tri.vertices]
    V = triangles[:,2,:] - triangles[:,0,:]
    U = triangles[:,1,:] - triangles[:,0,:]
    n = np.cross(U,V)
    
    mag = np.sqrt((n ** 2).sum(-1)) #inner1d can make faster
    z_norm = np.abs(n[:,-1])/mag
    
    return z_norm, np.mean(mag)

def alpha_shape(points, tri, alpha=1):
    """
    Compute the alpha shape (concave hull) of a set
    of points using shapely routine.
    @param points: Iterable container of points.
    @param alpha: alpha multiplier. The higher the
    multiplier, the more points remain. Setting too low will fail and will return None
    """
    from shapely.ops import unary_union, polygonize
    import shapely.geometry as geometry
    
    coords = points.copy()
    
    triangles = coords[tri.vertices]
    a = ((triangles[:,0,0] - triangles[:,1,0]) ** 2 + (triangles[:,0,1] - triangles[:,1,1]) ** 2) ** 0.5
    b = ((triangles[:,1,0] - triangles[:,2,0]) ** 2 + (triangles[:,1,1] - triangles[:,2,1]) ** 2) ** 0.5
    c = ((triangles[:,2,0] - triangles[:,0,0]) ** 2 + (triangles[:,2,1] - triangles[:,0,1]) ** 2) ** 0.5
    s = ( a + b + c ) / 2.0 #semiperimeter
    areas = (s*(s-a)*(s-b)*(s-c)) ** 0.5 #Heron's formula
    triangles = coords[tri.vertices]
    
    #make sure the areas are valid
    valid = np.isreal(areas) & ~np.isnan(areas) & ~np.isinf(areas)

    np.seterr(divide='ignore', invalid='ignore')
    co = a[valid] * b[valid] * c[valid] / (4.0 * areas[valid]) #circumradius
    np.seterr(divide='warn', invalid='warn')
    
    cutoff = np.mean(co[np.isreal(co) & ~np.isnan(co) & ~np.isinf(co)])
    filtered = triangles[valid][co < (alpha*cutoff)] #as opposed to the published 1 / alpha

    edge1 = filtered[:,(0,1)]
    edge2 = filtered[:,(1,2)]
    edge3 = filtered[:,(2,0)]
    edge_points = np.unique(np.concatenate((edge1,edge2,edge3)), axis = 0).tolist()
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    chull = unary_union(triangles)
    try:
        x,y = chull.exterior.coords.xy
        outline = np.column_stack((x,y,np.zeros(len(x)))) #return outline appearing at z=0
        return outline
    except:
        #if chull doesn't have the exterior attribute (alpha was too small)
        return None

def read_data(file):
    _, ext = os.path.splitext(file)
    if ext.lower() == '.csv':
        data = np.genfromtxt(file,delimiter=',')
        print('NAMRC data type recognised')
    elif ext.lower() == '.dat':
        data = np.genfromtxt(file,skip_header=1) / 1e3 #convert from micron to mm
    elif ext.lower() == '.txt':
        data = np.genfromtxt(file)
    return data

def read_file(file, return_active = False):
    '''
    Reads output file generated by the editor from this module. Returns a reference and floating list of outlines and active points.
    Active set to true will return active points instead of applying them to pnts
    TO DO - turn into a dictionary for re-editing
    '''
    ref_pnts = []
    ref_outlines = []
    float_pnts = []
    float_outlines = []
    ref_active_list = []
    float_active_list = []
    ref_trans = None
    float_trans = None
    fields = ['ref', 'float']
    cut_attr = dict.fromkeys(fields)
    for k in cut_attr.keys():
        cut_attr[k] = {}
    with h5py.File(file, 'r') as f:
        try: 
            gr = f[fields[0]]
            ref_trans = gr['transform'][()]
            for k in gr.keys():
                if k.isdigit():
                    pnts = gr['%s/points'%k][()]
                    active = gr['%s/active'%k][()]
                    if not return_active:
                        ref_pnts.insert(int(k),pnts[active])
                    else:
                        ref_pnts.insert(int(k),pnts)
                        ref_active_list.insert(int(k),active)
                    ref_outlines.insert(int(k), gr['%s/outline'%k][()])

        except: pass
        try: 
            gf = f[fields[1]]
            float_trans = gf['transform'][()]
            for k in gf.keys():
                if k.isdigit():
                    pnts = gf['%s/points'%k][()]
                    active = gf['%s/active'%k][()]
                    if not return_active:
                        float_pnts.insert(int(k),pnts[active])
                    else:
                        float_pnts.insert(int(k),pnts)
                        float_active_list.insert(int(k),active)
                    float_outlines.insert(int(k), gf['%s/outline'%k][()])
        except: pass
        try:
            for k in fields:
                cut_attr[k].update(f[k].attrs)
        except: pass
    
    if not return_active:
        return ref_pnts, ref_outlines, float_pnts, float_outlines, cut_attr
    else: 
        return ref_pnts, ref_outlines, float_pnts, float_outlines, cut_attr,ref_active_list, float_active_list, ref_trans, float_trans

if __name__ == "__main__":
    if len(sys.argv)>1:
        launch(sys.argv[1])
    else:
        launch()