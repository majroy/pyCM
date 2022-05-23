#!/usr/bin/env python
'''
Uses VTK Python to allow for alignment and averaging point clouds associated with the contour method. Full interaction requires a 3-button mouse and keyboard. See documentation for keyboard/mouse bindings.
1.6 - Updated for overall version 2.
'''

__author__ = "M.J. Roy"
__version__ = "1.6"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014--"

import sys, os
from pkg_resources import Requirement, resource_filename
import numpy as np
from scipy.interpolate import griddata
import vtk
import vtk.util.numpy_support as v2n
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
from pyCMcommon import *
from icp import *
from registration import read_file as reg_read_file


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

class aa_main_window(object):
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
            MainWindow.setWindowTitle("pyCM - Alignment and averaging tool v%s" %__version__)
        
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
        
        #make display box
        self.display_box = QtWidgets.QGroupBox('Display')
        display_layout = QtWidgets.QGridLayout()
        self.display_box.setLayout(display_layout)
        self.z_aspect_sb = QtWidgets.QSpinBox()
        self.z_aspect_sb.setPrefix("Scaling: ")
        self.z_aspect_sb.setToolTip("Scaling factor applied to Z axis")
        self.z_aspect_sb.setValue(1)
        self.z_aspect_sb.setMinimum(1)
        self.z_aspect_sb.setMaximum(1000)
        self.draw_directions_button = QtWidgets.QPushButton('Show cut')
        self.draw_directions_button.setToolTip('Show cutting direction(s)')
        self.draw_directions_button.setCheckable(True)
        self.caption_rb = QtWidgets.QCheckBox('Caption entries')
        self.caption_rb.setToolTip('Label entries in viewport')
        self.caption_rb.setChecked(False)
        
        ref_op_slider_label = QtWidgets.QLabel("Reference opacity:")
        self.ref_op_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.ref_op_slider.setStyleSheet("QSlider::handle:horizontal {background-color: rgb(242, 101, 34);}")
        self.ref_op_slider.setRange(0,100)
        self.ref_op_slider.setSliderPosition(100)
        float_op_slider_label = QtWidgets.QLabel("Floating opacity:")
        self.float_op_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.float_op_slider.setStyleSheet("QSlider::handle:horizontal {background-color: rgb(255, 205, 52);}")
        self.float_op_slider.setRange(0,100)
        self.float_op_slider.setSliderPosition(100)
        avg_op_slider_label = QtWidgets.QLabel("Averaged opacity:")
        self.avg_op_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.avg_op_slider.setStyleSheet("QSlider::handle:horizontal {background-color: rgb(25, 25, 112);}")
        self.avg_op_slider.setRange(0,100)
        self.avg_op_slider.setSliderPosition(100)
        self.avg_op_slider.setEnabled(False)
        
        local_display_layout = QtWidgets.QHBoxLayout()
        local_display_layout.addWidget(self.z_aspect_sb)
        local_display_layout.addWidget(self.caption_rb)
        local_display_layout.addWidget(self.draw_directions_button)
        display_layout.addLayout(local_display_layout,0,0,1,2)
        display_layout.addWidget(ref_op_slider_label,1,0,1,1)
        display_layout.addWidget(self.ref_op_slider,1,1,1,1)
        display_layout.addWidget(float_op_slider_label,2,0,1,1)
        display_layout.addWidget(self.float_op_slider,2,1,1,1)
        display_layout.addWidget(avg_op_slider_label,3,0,1,1)
        display_layout.addWidget(self.avg_op_slider,3,1,1,1)
        
        self.display_box.setEnabled(False)
        
        #make mirror box
        self.mirror_box = QtWidgets.QGroupBox('Mirroring')
        mirror_layout = QtWidgets.QHBoxLayout()
        self.mirror_box.setLayout(mirror_layout)
        self.mirror_x_button = QtWidgets.QPushButton('ZX')
        self.mirror_x_button.setToolTip('Mirror floating data on ZY (x = 0) plane')
        self.mirror_y_button = QtWidgets.QPushButton('ZY')
        self.mirror_y_button.setToolTip('Mirror floating data on ZY (y = 0) plane')
        self.align_reset_button = QtWidgets.QPushButton("Reset")
        self.align_reset_button.setToolTip('Undo all transformations to the floating point dataset')
        mirror_layout.addWidget(self.mirror_x_button)
        mirror_layout.addWidget(self.mirror_y_button)
        mirror_layout.addWidget(self.align_reset_button)
        
        self.mirror_box.setEnabled(False)
        
        #make alignment box
        self.align_box = QtWidgets.QGroupBox('Alignment')
        align_layout = QtWidgets.QGridLayout()
        self.align_box.setLayout(align_layout)
        cent_layout = QtWidgets.QHBoxLayout()
        cent_label=QtWidgets.QLabel("Move to centroids:")
        self.centroid_reference = QtWidgets.QRadioButton('Reference')
        self.centroid_reference.setChecked(True)
        self.centroid_floating = QtWidgets.QRadioButton('Floating')
        centroid_button_group = QtWidgets.QButtonGroup()
        centroid_button_group.addButton(self.centroid_reference)
        centroid_button_group.addButton(self.centroid_floating)
        centroid_button_group.setExclusive(True)
        self.cent_move_button = QtWidgets.QPushButton("Apply")
        cent_layout.addWidget(cent_label)
        cent_layout.addWidget(self.centroid_reference)
        cent_layout.addWidget(self.centroid_floating)
        
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
        self.move_button.setSizePolicy(sizePolicy)
        
        align_algo_layout = QtWidgets.QHBoxLayout()
        self.k_neighbour_button = QtWidgets.QPushButton("K-neighbour ICP")
        self.k_neighbour_button.setToolTip('Generate and apply a 2D transformation based on outlines')
        self.corner_best_fit_button = QtWidgets.QPushButton("Corner SVD")
        self.corner_best_fit_button.setToolTip('Solve and apply the best-fit 2D transform to outline corners with single value decomposition')
        self.vtk_icp_button = QtWidgets.QPushButton("VTK ICP")
        self.vtk_icp_button.setToolTip('Generate and apply a VTK 3D transformation based on outlines')
        align_algo_layout.addWidget(self.k_neighbour_button)
        align_algo_layout.addWidget(self.corner_best_fit_button)
        align_algo_layout.addWidget(self.vtk_icp_button)
        
        align_layout.addLayout(cent_layout,0,0,1,2)
        align_layout.addWidget(self.cent_move_button,0,2,1,1)
        align_layout.addWidget(trans_x_label,1,0,1,1)
        align_layout.addWidget(self.trans_x,1,1,1,1)
        align_layout.addWidget(trans_y_label,2,0,1,1)
        align_layout.addWidget(self.trans_y,2,1,1,1)
        align_layout.addWidget(rotate_label,3,0,1,1)
        align_layout.addWidget(self.rotate_z,3,1,1,1)
        align_layout.addWidget(self.move_button,1,2,3,1)
        align_layout.addLayout(align_algo_layout,4,0,1,3)
        align_layout.setRowStretch(align_layout.rowCount(), 1)
        
        self.align_box.setEnabled(False)
        
        #make average layout
        self.average_box = QtWidgets.QGroupBox('Averaging')
        average_layout = QtWidgets.QGridLayout()
        self.average_box.setLayout(average_layout)
        # grid_label=QtWidgets.QLabel("Grid spacing")
        self.grid_size = QtWidgets.QDoubleSpinBox()
        self.grid_size.setPrefix('Grid size: ')
        self.grid_size.setToolTip("Spacing of grid which will be used to average the datasets")
        self.grid_size.setValue(0)
        self.grid_size.setMaximum(300)
        self.grid_size.setMinimum(-300)
        self.grid_size.setDecimals(4)
        self.mask_average_rb = QtWidgets.QCheckBox('Outline crop')
        self.mask_average_rb.setToolTip("Mask/crop display of averaged data to the reference outline")
        self.nearest_neigh_rb = QtWidgets.QCheckBox('N-neighbour')
        self.nearest_neigh_rb.setToolTip("Set undefined values to nearest neighbour as opposed to lowest defined value - requires re-averaging to effect")
        self.mask_average_rb.setChecked(True)
        self.nearest_neigh_rb.setChecked(True)
        self.average_button = QtWidgets.QPushButton('Average')
        self.average_button.setToolTip('Start averaging')
        self.reset_average_button = QtWidgets.QPushButton('Reset')
        self.reset_average_button.setToolTip('Invalidate averaging and return to alignment')
        self.averaging_status_label = QtWidgets.QLabel("Ready")
        
        # average_layout.addWidget(grid_label,0,0,1,1)
        average_layout.addWidget(self.grid_size,0,0,1,1)
        average_layout.addWidget(self.mask_average_rb,0,1,1,1)
        average_layout.addWidget(self.nearest_neigh_rb,0,2,1,1)
        average_layout.addWidget(self.average_button,1,1,1,1)
        average_layout.addWidget(self.reset_average_button,1,2,1,1)
        average_layout.addWidget(self.averaging_status_label,2,0,1,3)
        
        self.average_box.setEnabled(False)
        
        #make save box
        self.save_box = QtWidgets.QGroupBox('Write current')
        save_layout = QtWidgets.QHBoxLayout()
        self.save_box.setLayout(save_layout)
        #make save buttons
        self.save_button = QtWidgets.QPushButton("Save")
        self.save_button.setToolTip('Save data to hdf5-formatted results file')
        self.save_label = QtWidgets.QLabel('Ready')
        self.save_label.setWordWrap(True)
        save_layout.addWidget(self.save_label)
        verticalSpacer = QtWidgets.QSpacerItem(50, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        save_layout.addItem(verticalSpacer)
        save_layout.addWidget(self.save_button)
        self.save_box.setEnabled(False)
        
        
        lvlayout=QtWidgets.QVBoxLayout()
        lvlayout.minimumSize()
        lvlayout.addWidget(self.display_box)
        lvlayout.addWidget(self.mirror_box)
        lvlayout.addWidget(self.align_box)
        lvlayout.addWidget(self.average_box)
        lvlayout.addWidget(self.save_box)
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
        self.ui = aa_main_window()
        self.ui.setup(self)
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(vtk.vtkNamedColors().GetColor3d("slategray"))
        self.ren.GradientBackgroundOn()

        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style=vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        # self.iren.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ui.vtkWidget.Initialize()
        
        self.file = None
        self.active_dir = os.getcwd()
        self.point_size = 2
        self.limits = np.empty(6)
        self.gsize = 0
        self.trans = np.eye(4) #default
        self.c_trans = np.eye(4) #initialise cumulative transformation
        
        
        #connect functions/methods to ui here
        self.ui.z_aspect_sb.valueChanged.connect(self.update_z_aspect)
        self.ui.caption_rb.toggled.connect(self.hide_caption)
        self.ui.draw_directions_button.clicked.connect(self.draw_cut_directions)
        self.ui.ref_op_slider.valueChanged[int].connect(self.change_ref_opacity)
        self.ui.float_op_slider.valueChanged[int].connect(self.change_float_opacity)
        self.ui.avg_op_slider.valueChanged[int].connect(self.change_avg_opacity)
        self.ui.mirror_x_button.clicked.connect(self.trans_mirror_x)
        self.ui.mirror_y_button.clicked.connect(self.trans_mirror_y)
        self.ui.align_reset_button.clicked.connect(self.reset_trans)
        self.ui.trans_x.editingFinished.connect(self.clear_var_yz)
        self.ui.trans_y.editingFinished.connect(self.clear_var_xz)
        self.ui.rotate_z.editingFinished.connect(self.clear_var_xy)
        self.ui.cent_move_button.clicked.connect(self.centroid_move)
        self.ui.move_button.clicked.connect(self.trans_from_ui)
        self.ui.vtk_icp_button.clicked.connect(self.trans_from_vtk_icp)
        self.ui.k_neighbour_button.clicked.connect(self.trans_from_k_neighbour)
        self.ui.corner_best_fit_button.clicked.connect(self.trans_from_corners)
        self.ui.average_button.clicked.connect(self.average)
        self.ui.reset_average_button.clicked.connect(self.reset_average)
        self.ui.mask_average_rb.toggled.connect(self.draw_average)
        self.ui.save_button.clicked.connect(self.write_output)

    def get_data(self):
        '''
        Calls load_h5 method to load data. Calls reset_display, reset_align, reset_average
        '''
        if self.file is None:
            self.file, self.active_dir = get_file("*.pyCM",self.active_dir)
        
        #make sure valid file was selected
        if self.file is None or not(os.path.isfile(self.file)):
            return
        
        #delete any existing objects
        if hasattr(self,'active_pnt'):
            del self.active_pnt
        
        #clear renderer
        self.ren.RemoveAllViewProps()
        
        #otherwise start reading in registration data
        self.rp, \
        self.ro, \
        self.fp, \
        self.fo, \
        self.cut_attr = reg_read_file(self.file)
        
        #move outlines and data to the mean z value of point clouds
        rm = np.mean(np.vstack(self.rp)[:,-1])
        fm = np.mean(np.vstack(self.fp)[:,-1])
        for i in range(len(self.rp)):
            self.rp[i][:,-1] = self.rp[i][:,-1] - rm
            self.ro[i][:,-1] = rm
            self.fp[i][:,-1] = self.fp[i][:,-1] - fm
            self.fo[i][:,-1] = fm
        
        self.trans = np.eye(4) #default
        self.c_trans = np.eye(4) #initialise
        
        #try reading in any alignment/averaged data
        avg, gsize, T = read_file(self.file)
        if avg:
            self.avg = avg
            self.gsize = gsize
            self.ui.grid_size.setValue(gsize)
            self.draw(0)
            self.trans = T
            self.apply_transformation()
            self.draw_average()
        else:
            self.draw(0)

        self.ui.display_box.setEnabled(True)
        if not self.cut_attr['ref']:
            self.ui.draw_directions_button.setEnabled(False)
        else:
            self.ui.draw_directions_button.setEnabled(True)
        self.ui.mirror_box.setEnabled(True)
        self.ui.align_box.setEnabled(True)
        self.ui.average_box.setEnabled(True)
        
            
    def draw(self,option=1):
        '''
        Draws reference and floating point clouds depending on option:
        0 - from load, draws both
        1 - draws floating only
        '''
        if option == 0:
        #remove all actors
            self.ren.RemoveAllViewProps()
            self.ref_actor_list = []
            self.ref_outline_polydata_list = []
            self.ref_caption_actor_list = []
            
            color=(242, 101, 34)
            for i in range(len(self.rp)):
                rp_actor, _, _, _, = gen_point_cloud(self.rp[i],color,self.point_size)
                rp_actor.GetProperty().SetOpacity(\
                self.ui.ref_op_slider.value() / 100)
                self.ren.AddActor(rp_actor)
                self.ref_actor_list.append(rp_actor)
                ro_actor, ro_pc = gen_outline(self.ro[i],color,self.point_size)
                self.ren.AddActor(ro_actor)
                self.ref_outline_polydata_list.append(ro_pc)
                self.ref_actor_list.append(ro_actor)
                cap_actor = gen_caption_actor('%s'%i, ro_actor)
                self.ren.AddActor(cap_actor)
                self.ref_caption_actor_list.append(cap_actor)
        
        elif option == 1:
            for actor in self.float_actor_list:
                self.ren.RemoveActor(actor)
            for actor in self.float_caption_actor_list:
                self.ren.RemoveActor(actor)

        self.float_actor_list = []
        self.float_outline_polydata_list = []
        self.float_caption_actor_list = []
        
        color = (255, 205, 52)
        for i in range(len(self.fp)):
            fp_actor, _, _, _, = gen_point_cloud(self.fp[i],color,self.point_size)
            fp_actor.GetProperty().SetOpacity(\
            self.ui.float_op_slider.value() / 100)
            self.ren.AddActor(fp_actor)
            self.float_actor_list.append(fp_actor)
            fo_actor, fo_pc = gen_outline(self.fo[i],color,self.point_size)
            self.ren.AddActor(fo_actor)
            self.float_outline_polydata_list.append(fo_pc)
            self.float_actor_list.append(fo_actor)
            cap_actor = gen_caption_actor('%s'%i, fo_actor)
            self.ren.AddActor(cap_actor)
            self.float_caption_actor_list.append(cap_actor)
        
        self.limits = get_limits(np.vstack(self.rp + self.fp))
        self.update_z_aspect()
        
        if self.ui.draw_directions_button.isChecked():
            self.draw_cut_directions()
        
        if not self.ui.caption_rb.isChecked():
            self.hide_caption()
        
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
    
    def draw_average(self):
        '''
        -(Re)draws averaged data set
        -flips the visibility of the floating outline
        -turns the reference outline white
        '''
        if hasattr(self,'avg_actor_list'):
            for actor in self.avg_actor_list:
                self.ren.RemoveActor(actor)
            self.ren.RemoveActor(self.sb_actor)
        self.avg_actor_list = []
        
        #change opacity of floating and reference data
        self.ui.ref_op_slider.setValue(0) 
        self.ui.float_op_slider.setValue(0)
        self.ui.avg_op_slider.setEnabled(True)
        
        #handle outlines
        float_outline_actors = self.float_actor_list[1::2]
        for actor in float_outline_actors:
            if actor.GetVisibility():
                flip_visible(actor)
        ref_outline_actors = self.ref_actor_list[1::2]
        for actor in ref_outline_actors:
            actor.GetProperty().SetColor((255,255,255))
        
        for i in range(len(self.avg)):
            if hasattr(self,'active_pnt'):
                if self.ui.mask_average_rb.isChecked():
                    pnts = self.avg[i][self.active_pnt[i]]
            else:
                pnts = self.avg[i]
            
            if i == 0: 
                avg_actor, \
                _, \
                _, lut = \
                gen_point_cloud(pnts,'blues')
                self.zrange = (self.limits[-2],self.limits[-1])
            else:
                avg_actor, \
                _, \
                _, _ = \
                gen_point_cloud(pnts,'blues')
                self.zrange = (self.limits[-2],self.limits[-1])

            avg_actor.GetProperty().SetOpacity(\
            self.ui.avg_op_slider.value() / 100)
            self.ren.AddActor(avg_actor)
            self.avg_actor_list.append(avg_actor)
            
        #handle scalebar
        sb_widget = gen_scalar_bar()
        sb_widget.SetInteractor(self.iren)
        sb_widget.On()
        self.sb_actor = sb_widget.GetScalarBarActor()
        self.sb_actor.SetLookupTable(lut)
        self.ren.AddActor(self.sb_actor)
        
        self.update_z_aspect()
        
        
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
    
    def hide_caption(self):
        state = self.ui.caption_rb.isChecked()
        for actor in self.ref_caption_actor_list:
            if not state:
                actor.VisibilityOff()
            else:
                actor.VisibilityOn()
        for actor in self.float_caption_actor_list:
            if not state:
                actor.VisibilityOff()
            else:
                actor.VisibilityOn()
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
            
        if hasattr(self,'float_cut_orient_actor'):
            self.ren.RemoveActor(self.float_cut_orient_actor)
        
        if self.ui.draw_directions_button.isChecked():
            pass
        else:
            self.ui.vtkWidget.update()
            return
            
        if cut_attr['ref']:
            ld = cut_attr['ref']
            if ld['cut_path'].any():
                self.ref_cut_orient_actor = gen_cutting_orientation_actor(\
                get_limits(np.vstack(self.ro)),\
                ld['cut_dir'],\
                ld['cut_path'])
            else:
                self.ref_cut_orient_actor = gen_cutting_orientation_actor(\
                get_limits(np.vstack(self.ro)),\
                ld['cut_dir'])
            
            self.ren.AddActor(self.ref_cut_orient_actor)
            self.ref_cut_orient_actor.GetProperty().SetColor(\
            (242/255, 101/255, 34/255)\
            )
            
        if cut_attr['float']:
            ld = cut_attr['float']
            if ld['cut_path'].any():
                self.float_cut_orient_actor = gen_cutting_orientation_actor(\
                get_limits(np.vstack(self.fo)),\
                ld['cut_dir'],\
                ld['cut_path'])
            else:
                self.float_cut_orient_actor = gen_cutting_orientation_actor(\
                get_limits(np.vstack(self.fo)),\
                ld['cut_dir'])
            self.ren.AddActor(self.float_cut_orient_actor)
            self.float_cut_orient_actor.GetProperty().SetColor(\
            (255/255, 205/255, 52/255)\
            )
        
        self.ui.vtkWidget.update()
    
    def reset_average(self):
        #change opacity of floating and reference data
        self.ui.ref_op_slider.setValue(100) 
        self.ui.float_op_slider.setValue(100)
        self.ui.avg_op_slider.setEnabled(False)
        self.ui.save_box.setEnabled(False)
        
        self.draw(0)
        
        self.ui.mirror_box.setEnabled(True)
        self.ui.align_box.setEnabled(True)
        
    
    def reset_trans(self):
        '''
        inverts the cumulative transformation matrix and applies it to the floating data. Resets self.trans and c_trans.
        '''
        self.trans = np.linalg.inv(self.c_trans)
        self.apply_transformation()
        
        
    def centroid_move(self):
        '''
        Translates either reference data directly or floating data by apply_transformation, depending on what radio button in the ui is selected.
        TO DO: UPDATE REF entry and transformation matrix in data record
        '''
        
        if self.ui.centroid_reference.isChecked():
            T = np.eye(4)
            T[0:3,-1] = - np.mean(np.vstack((self.rp + self.ro)), axis=0)
            for i in range(len(self.rp)):
                self.rp[i] = do_transform(self.rp[i],T)
                self.ro[i] = do_transform(self.ro[i],T)
            self.draw(0)
        elif self.ui.centroid_floating.isChecked():
            self.trans[0:3,-1] = - np.mean(np.vstack((self.fp + self.fo)), axis=0)
            self.apply_transformation()
    
    def trans_mirror_x(self):
        '''
        Mirrors floating on ZX (X=0) plane
        '''
        
        self.trans[1,1] = -1
        self.apply_transformation()
        
    def trans_mirror_y(self):
        '''
        Mirrors floating on ZY (Y=0) plane
        '''
        self.trans[0,0] = -1
        self.apply_transformation()

    def trans_from_ui(self):
        '''
        Builds and applies a transformation from entries on the ui
        '''
        self.trans[0,-1] = self.ui.trans_x.value()
        self.trans[1,-1] = self.ui.trans_y.value()
        a = np.deg2rad(self.ui.rotate_z.value()) #negative for counterclockwise
        self.trans[0:2,0:2]=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        self.apply_transformation()

    def trans_from_vtk_icp(self):
        '''
        Get icp transformation matrix from outline polydata
        '''
        fo_pd = self.float_outline_polydata_list
        ro_pd = self.ref_outline_polydata_list
        
        T = [vtk_icp(\
        fo_pd[i],\
        ro_pd[i]) for i in \
        range(len(fo_pd))\
        ]
        self.trans = sum(T)/len(T)
        self.apply_transformation()
        
    
    def trans_from_k_neighbour(self):
        '''
        Get 2D k-neighbour transformation matrix from icp algo, build and apply 3D according to this.
        '''
        #build array of all outlines with each having the same number of points
        respaced_float_outline = []
        for i in range(len(self.fo)):
            local_outline, _, _ = respace_equally(self.fo[i][:,:2],len(self.ro[i]))
            respaced_float_outline.append(local_outline)
        #A (fo) and B (ro) need to have the same shape. Respace fo to ro.
        A = np.vstack(respaced_float_outline)
        reference_outline = [x[:,:2] for x in self.ro]
        B = np.vstack(reference_outline)
        two_d_trans, _, _ = icp(A,B)
        #pack relevant values into self.trans from two_d_trans
        self.trans[:2,-1] = two_d_trans[:2,-1]
        self.trans[0:2,0:2] = two_d_trans[0:2,0:2]
        self.apply_transformation()
    
    def trans_from_corners(self):
        '''
        Uses get_corner_ind to find corners, and aligns according to best_fit_transform imported from icp
        '''
        #get fo corners
        A = []
        for outline in self.fo:
            i, new_fo = get_corner_ind(outline)
            A.append(new_fo[i,:2])
        B = []
        for outline in self.ro:
            j, new_ro = get_corner_ind(outline)
            B.append(new_ro[j,:2])
        
        _, R, t = best_fit_transform(np.vstack(A), np.vstack(B))
        #pack relevant values into self.trans
        self.trans[:2,-1] = t
        self.trans[0:2,0:2] = R
        self.apply_transformation()
    
    def apply_transformation(self):
        '''
        Apply/update transformation
        - NB translations are not captured correctly in c_trans
        '''
        T = self.trans.copy()
        #assumes that the number of entries in points and outline are the same
        for i in range(len(self.fp)):
            self.fp[i] = do_transform(self.fp[i],T)
            self.fo[i] = do_transform(self.fo[i],T)
        self.c_trans = T @ self.c_trans #pre-multiply, same as np.dot(self.c_trans,T.T)
        
        #update cutting directions of floating
        if self.cut_attr['float']:
            self.cut_attr['float']['cut_dir'] = T[:3,:3] @ self.cut_attr['float'].get('cut_dir')
            self.cut_attr['float']['cut_path'] = T[:3,:3] @ self.cut_attr['float'].get('cut_path')
                
        self.draw(1) #calls draw updating floating only
        self.trans = np.eye(4)#reset the local transformation matrix

    def average(self):
        '''
        Averages point cloud and calls draw_average with the result
        '''
        
        self.ui.averaging_status_label.setText("Averaging, setting grid . . .")
        QtWidgets.QApplication.processEvents()
        
        if self.ui.grid_size.value() == 0:
            #create unique combinations of reference points and calculate their mean separation distance based on the first entry
            target = self.rp[0]
            combination = np.array(np.meshgrid(target[:,:2])).reshape(-1,2)
            distance = np.sqrt(np.sum(np.diff(combination,axis=0)**2,axis=1))
            self.ui.grid_size.setValue(np.mean(distance))
            QtWidgets.QApplication.processEvents()
            
        num_entries = len(self.ro)
        self.avg = [None] * num_entries
        self.active_pnt = [None] * num_entries
        self.gsize = self.ui.grid_size.value()
        
        #start loop based on entries
        for i in range(num_entries):
            #create meshgrid based on limits of reference outline (x,y = meshgrid = y,x = mgrid)
            ref_limits = get_limits(self.ro[i], 0.01)
            
            x = np.linspace(ref_limits[0],ref_limits[1],\
            int((ref_limits[1]-ref_limits[0]) / self.gsize))
            y = np.linspace(ref_limits[2],ref_limits[3],\
            int((ref_limits[3]-ref_limits[2]) / self.gsize)) 
            
            grid_x, grid_y = np.meshgrid(\
            x,y,
            indexing='xy',\
            sparse = False)
            
            #apply grid to reference data
            rp = self.rp[i]#.copy()
            ref_grid = griddata(\
            rp[:,:2],\
            rp[:,-1],\
            (grid_x,grid_y),\
            method='linear')
            
            fp = self.fp[i]#.copy()
            #apply grid to the (hopefully) aligned floating data
            aligned_grid = griddata(\
            fp[:,:2],\
            fp[:,-1],\
            (grid_x,grid_y),\
            method='linear')
            
            if self.ui.nearest_neigh_rb.isChecked():
                self.ui.averaging_status_label.setText("Applying nearest neighbours to outliers . . .")
                QtWidgets.QApplication.processEvents()        
                ref_grid_nearest = griddata(\
                rp[:,:2],\
                rp[:,-1],\
                (grid_x,grid_y),\
                method='nearest')
                
                ref_mask = np.isnan(ref_grid).any(axis=1)
                ref_grid[ref_mask] = ref_grid_nearest[ref_mask]
            
                aligned_grid_nearest = griddata(\
                fp[:,:2],\
                fp[:,-1],\
                (grid_x,grid_y),\
                method='nearest')
                
                aligned_mask = np.isnan(aligned_grid).any(axis=1)
                aligned_grid[aligned_mask] = aligned_grid_nearest[aligned_mask]
            
            self.ui.averaging_status_label.setText("Averaging entry %d with grid . . ."%i)
            QtWidgets.QApplication.processEvents()
            
            avg_grid = (ref_grid + aligned_grid)/2
            
            #flatten the grids (if not using sparse grids)
            avg_xy_flattened = np.hstack(\
            (np.ravel(grid_x.T)[np.newaxis].T,\
            np.ravel(grid_y.T)[np.newaxis].T))

            #re-tile x & y to match average if gridx &y are sparse
            # avg_xy_flattened = np.hstack((
            # np.tile(grid_x,(grid_y.shape[1],1)),\
            # np.repeat(grid_y,grid_x.shape[0])[np.newaxis].T))
            
            raw_avg = np.hstack(\
            (avg_xy_flattened, \
            np.ravel(avg_grid.T)[np.newaxis].T))
        
            #generate mask with in_poly
            self.active_pnt[i] = in_poly(self.ro[i],raw_avg)

            #plug potential nan values
            plug = np.nanmin(raw_avg[:,2])
            raw_avg[:,2] = np.nan_to_num(raw_avg[:,2], nan=plug)
            
            self.avg[i] = raw_avg

        self.ui.averaging_status_label.setText("Rendering . . .")
        QtWidgets.QApplication.processEvents()
        #call draw_average
        self.draw_average()
        
        self.ui.mirror_box.setEnabled(False)
        self.ui.align_box.setEnabled(False)
        self.ui.save_box.setEnabled(True)
        
        self.ui.averaging_status_label.setText("Ready")

    def update_z_aspect(self):
        '''
        Updates z_aspect and redraws points displayed based on new z_aspect
        '''
        z_aspect = self.ui.z_aspect_sb.value()
        #scale points, not outlines, so skipping every other entry in the actor list
        for actor in (self.float_actor_list + self.ref_actor_list)[::2]:
            actor.SetScale((1,1,z_aspect))
            actor.Modified()

        if hasattr(self,'avg_actor_list'):
            for actor in self.avg_actor_list:
                actor.SetScale((1,1,z_aspect))
                actor.Modified()
        
        #now do axis_actor, can't scale in the same way as polydata
        try:
            self.ren.RemoveActor(self.axis_actor) #it will need to be replaced
        except: pass
        #local limits for z axis - don't scale x & y limits, scale the z axis according to z aspect and the current limits
        nl=np.append(self.limits[0:4],([self.limits[-2]*z_aspect,self.limits[-1]*z_aspect]))
        self.axis_actor = get_axis(self.ren, nl, z_aspect)
        self.ren.AddActor(self.axis_actor)
        self.ui.vtkWidget.update()

    def clear_var_yz(self):
        '''
        following clear_var_* functions clear other inputs to ensure serial manual manipulation.
        '''
        self.ui.trans_y.setValue(0)
        self.ui.rotate_z.setValue(0)
        
    def clear_var_xz(self):
        self.ui.trans_x.setValue(0)
        self.ui.rotate_z.setValue(0)
        
    def clear_var_xy(self):
        self.ui.trans_x.setValue(0)
        self.ui.trans_y.setValue(0)
        
        
    def change_ref_opacity(self,value):
        if hasattr(self,'ref_actor_list'):
            for actor in self.ref_actor_list[::2]:
                actor.GetProperty().SetOpacity(value/100)
        self.ui.vtkWidget.update()
    
    def change_float_opacity(self,value):
        if hasattr(self,'float_actor_list'):
            for actor in self.float_actor_list[::2]:
                actor.GetProperty().SetOpacity(value/100)
        self.ui.vtkWidget.update()
        
    def change_avg_opacity(self,value):
        if hasattr(self,'avg_actor_list'):
            for actor in self.avg_actor_list:
                actor.GetProperty().SetOpacity(value/100)
        self.ui.vtkWidget.update()

    def check_save_state(self, id):
        '''
        Checks to make sure that there are data objects pertaining to an averaged surface within the interactor against those that might be present in the specified file
        '''
        if not hasattr(self,'avg'):
            info_msg('Saving the current step requires an averaged dataset.')
            return False
        
        if self.file is None:
            return False

        with h5py.File(self.file, 'r') as f:
            existing_keys = list(f.keys())
            if not existing_keys or id not in existing_keys:
                return True
        
        with h5py.File(self.file, 'r') as f:
            #check the id's entry keys
            g = f['%s'%id]
            if 'aa' in list(g.keys()):
                overwrite = warning_msg(self, \
                'There is existing data in the specified file for the target field, overwrite?'
                )
                if not overwrite: #fail check
                    return False
                else:
                    #overwrite
                    return True
            else:
                #there isn't an entry in the id
                return True

    def write_output(self):
        
        passed = self.check_save_state('aa')
        
        if not passed:
            self.ui.save_label.setText('Ready')
            return
        else:
            self.ui.save_label.setText('Saving . . .')
            QtWidgets.QApplication.processEvents()
            
            with h5py.File(self.file, 'r+') as f:
                if 'aa' in list(f.keys()):
                    del f['aa']
                g = f.create_group('aa')
                g.attrs['grid_size'] = self.gsize
                g.create_dataset('transform', data = self.c_trans)
                for i in range(len(self.avg)):
                    gg = g.create_group(str(i))
                    gg.create_dataset('points',data = self.avg[i][self.active_pnt[i]])
                f.attrs['date_modified'] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                
        self.ui.save_label.setText('Saved to %s'%(os.path.basename(self.file)))
        
        
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
        elif key == "Up":
            self.ren.GetActiveCamera().Roll(30)
        elif key == "Down":
            self.ren.GetActiveCamera().Roll(-30)
        elif key == "1":
            xyview(self.ren)
        elif key == "2":
            yzview(self.ren)
        elif key == "3":
            xzview(self.ren)
        elif key == "a":
            flip_visible(self.axis_actor)
        if key == "l":
            self.file = None
            self.get_data()
            
        self.ui.vtkWidget.update()

def read_file(file):
    '''
    Reads output file generated by the editor from this module. Returns a list of avg and active entries ready for plotting, along with transform and grid size
    '''
    avg = []
    gsize = None
    T = np.eye(4)
    with h5py.File(file, 'r') as f:
        try:
            g = f['aa']
            gsize = g.attrs['grid_size']
            for k in g.keys():
                if k.isdigit():
                    local_avg = g['%s/points'%k][()]
                    #insert should generate same result as append as keys are read sequentially
                    avg.insert(int(k), local_avg)
            if k == 'transform':
                T = g['%s'%k][()]
        except: pass
    
    return avg, gsize, T

if __name__ == "__main__":
    if len(sys.argv)>1:
        launch(sys.argv[1])
    else:
        launch()