#!/usr/bin/env python
'''
Uses VTK Python to allow for fitting an averaged dataset associated with the
 contour method. Full interaction requires a 3-button mouse and keyboard.
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
from scipy.interpolate import griddata, bisplrep, bisplev
from scipy.spatial import Delaunay
import vtk
import vtk.util.numpy_support as v2n
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import Qt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib import rc
from pyCM.pyCMcommon import *
from pyCM.align_average import read_file as aa_read_file
from pyCM.registration import read_file as reg_read_file
from pyCM.decimate import decimate_ui

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

class sf_main_window(object):
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
            MainWindow.setWindowTitle("pyCM - Surface fitting tool v%s" %__version__)
        
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
        self.mask_average_rb = QtWidgets.QCheckBox('Outline crop')
        self.mask_average_rb.setToolTip("Mask/crop display of fitted surface to outline")
        self.mask_average_rb.setChecked(False)
        self.caption_rb = QtWidgets.QCheckBox('Caption entries')
        self.caption_rb.setToolTip('Label entries in viewport')
        self.caption_rb.setChecked(False)
        self.show_outlines_rb = QtWidgets.QCheckBox('Hide outlines')
        self.show_outlines_rb.setToolTip("Hides outlines if selected, shows if deselected")
        self.show_outlines_rb.setChecked(True)
        
        avg_op_slider_label = QtWidgets.QLabel("Averaged:")
        self.avg_op_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.avg_op_slider.setStyleSheet("QSlider::handle:horizontal {background-color: rgb(25, 25, 112);}")
        self.avg_op_slider.setRange(0,100)
        self.avg_op_slider.setSliderPosition(100)
        self.avg_op_slider.setEnabled(False)
        self.hide_avg_sb = QtWidgets.QCheckBox('Hide legend')
        self.hide_avg_sb.setToolTip("Hides colour map/contour legend for averaged data")
        self.hide_avg_sb.setEnabled(False)
        
        fit_op_slider_label = QtWidgets.QLabel("Fit:")
        self.fit_op_slider = QtWidgets.QSlider(Qt.Horizontal)
        self.fit_op_slider.setStyleSheet("QSlider::handle:horizontal {background-color: rgb(145,33,158);}")
        self.fit_op_slider.setRange(0,100)
        self.fit_op_slider.setSliderPosition(100)
        self.fit_op_slider.setEnabled(False)
        self.hide_fit_sb = QtWidgets.QCheckBox('Hide legend')
        self.hide_fit_sb.setToolTip("Hides colour map/contour legend for fitted surface")
        self.hide_fit_sb.setEnabled(False)

        local_display_layout = QtWidgets.QHBoxLayout()
        local_display_layout.addWidget(self.z_aspect_sb)
        local_display_layout.addWidget(self.mask_average_rb)
        local_display_layout.addWidget(self.caption_rb)
        local_display_layout.addWidget(self.show_outlines_rb)
        display_layout.addLayout(local_display_layout,0,0,1,3)
        display_layout.addWidget(avg_op_slider_label,1,0,1,1)
        display_layout.addWidget(self.avg_op_slider,1,1,1,1)
        display_layout.addWidget(self.hide_avg_sb,1,2,1,1)
        display_layout.addWidget(fit_op_slider_label,2,0,1,1)
        display_layout.addWidget(self.fit_op_slider,2,1,1,1)
        display_layout.addWidget(self.hide_fit_sb,2,2,1,1)

        self.display_box.setEnabled(False)
        
        #make active profile display
        entry_spec_layout = QtWidgets.QHBoxLayout()
        self.entry_spec = QtWidgets.QComboBox()
        self.entry_spec.setToolTip('Set focus on this entry')
        self.show_all_entries_button = QtWidgets.QPushButton("Show all")
        self.show_all_entries_button.setToolTip('Show all entries')
        entry_spec_layout.addWidget(self.entry_spec)
        entry_spec_layout.addWidget(self.show_all_entries_button)
        
        #make decimation box
        self.decimate_box = QtWidgets.QGroupBox('Decimation')
        decimate_layout = QtWidgets.QGridLayout()
        self.decimate_box.setLayout(decimate_layout)
        self.active_picking_indicator = QtWidgets.QLabel('Active')
        self.active_picking_indicator.setStyleSheet("background-color : gray; color : darkGray;")
        self.active_picking_indicator.setAlignment(QtCore.Qt.AlignCenter)
        self.active_picking_indicator.setFixedSize(50, 20)
        self.active_picking_indicator.setToolTip('Press R with interactor in focus to activate/deactivate manual point selection')
        self.active_picking_indicator.setEnabled(False)
        self.undo_last_pick_button=QtWidgets.QPushButton('Undo last')
        self.undo_last_pick_button.setToolTip('Undo last manual selection')
        self.apply_pick_button=QtWidgets.QPushButton('Apply')
        self.apply_pick_button.setToolTip('Remove selected points and update')
        self.reset_pick_button=QtWidgets.QPushButton('Reset')
        self.reset_pick_button.setToolTip('Reset all selected points and update')
        
        decimate_layout.addWidget(self.active_picking_indicator,0,0,1,1)
        decimate_layout.addWidget(self.undo_last_pick_button,0,1,1,1)
        decimate_layout.addWidget(self.apply_pick_button,0,2,1,1)
        decimate_layout.addWidget(self.reset_pick_button,0,3,1,1)
        
        self.decimate_box.setEnabled(False)
        
        #make bv box
        self.bv_box = QtWidgets.QGroupBox('Bivariate spline fitting')
        bv_layout = QtWidgets.QGridLayout()
        self.bv_box.setLayout(bv_layout)
        spline_space_label = QtWidgets.QLabel("Spacing")
        self.ksx = QtWidgets.QDoubleSpinBox()
        self.ksx.setPrefix("x = ")
        self.ksx.setToolTip("Knot spacing in x direction")
        self.ksx.setValue(0)
        self.ksx.setMinimum(0.0000001)
        self.ksx.setMaximum(1000)
        self.ksy = QtWidgets.QDoubleSpinBox()
        self.ksy.setPrefix("y = ")
        self.ksy.setToolTip("Knot spacing in y direction")
        self.ksy.setValue(0)
        self.ksy.setMinimum(0.0000001)
        self.ksy.setMaximum(10000)
        spline_order_label = QtWidgets.QLabel("Order")
        self.kox = QtWidgets.QSpinBox()
        self.kox.setPrefix("x: ")
        self.kox.setToolTip("Spline order in x direction")
        self.kox.setValue(3)
        self.kox.setMinimum(1)
        self.kox.setMaximum(5)
        self.koy = QtWidgets.QSpinBox()
        self.koy.setPrefix("y: ")
        self.koy.setToolTip("Spline order in y direction")
        self.koy.setValue(3)
        self.koy.setMinimum(1)
        self.koy.setMaximum(5)
        self.fit_all_rb = QtWidgets.QCheckBox('Fit all')
        self.fit_all_rb.setToolTip("Fit all entries with these parameters if selected, otherwise only the focussed entry.")
        self.fit_all_rb.setChecked(False)
        self.fit_spline_button = QtWidgets.QPushButton('Fit')
        self.fit_spline_button.setToolTip('Fit data to bivariate spline and update display')
        self.spline_status_label = QtWidgets.QLabel("Ready")
        
        bv_layout.addWidget(spline_space_label,0,0,1,1)
        bv_layout.addWidget(self.ksx,0,1,1,1)
        bv_layout.addWidget(self.ksy,0,2,1,1)
        bv_layout.addWidget(spline_order_label,1,0,1,1)
        bv_layout.addWidget(self.kox,1,1,1,1)
        bv_layout.addWidget(self.koy,1,2,1,1)
        bv_layout.addWidget(self.fit_all_rb,2,1,1,1)
        bv_layout.addWidget(self.fit_spline_button,2,2,1,1)
        bv_layout.addWidget(self.spline_status_label,3,0,1,3)
        
        
        # line extraction from surface
        self.extract_box = QtWidgets.QGroupBox('Evaluate fit on plane') #layout assigned further on
        #main layout for buttons and canvas
        extract_layout = QtWidgets.QGridLayout()

        # labels for axes
        x_label = QtWidgets.QLabel("x")
        y_label = QtWidgets.QLabel("y")

        # x, y, z of first point
        start_label = QtWidgets.QLabel("Centre")
        start_label.setToolTip('Centre of the target plane')
        self.point1_x_coord = QtWidgets.QDoubleSpinBox()
        self.point1_x_coord.setPrefix('x = ')
        self.point1_x_coord.setMinimum(-100000)
        self.point1_x_coord.setMaximum(100000)
        self.point1_y_coord = QtWidgets.QDoubleSpinBox()
        self.point1_y_coord.setPrefix('y = ')
        self.point1_y_coord.setMinimum(-100000)
        self.point1_y_coord.setMaximum(100000)

        # x, y, z of second point
        end_label = QtWidgets.QLabel("Normal")
        end_label.setToolTip('Normal of the target plane')
        self.point2_x_coord = QtWidgets.QDoubleSpinBox()
        self.point2_x_coord.setPrefix('x = ')
        self.point2_x_coord.setMinimum(-100000)
        self.point2_x_coord.setMaximum(100000)
        self.point2_y_coord = QtWidgets.QDoubleSpinBox()
        self.point2_y_coord.setPrefix('y = ')
        self.point2_y_coord.setMinimum(-100000)
        self.point2_y_coord.setMaximum(100000)
        
        self.extract_button = QtWidgets.QPushButton('Update')
        self.extract_button.setToolTip('Show/update line trace')
        self.export_line_status_label = QtWidgets.QLabel("Ready")
        self.export_line_status_label.setWordWrap(True)
        self.export_line_button = QtWidgets.QPushButton('Export')
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
        self.extract_box.setEnabled(False)
        
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
        save_layout.addItem(verticalSpacer)
        save_layout.addWidget(self.save_button)
        self.save_box.setEnabled(False)
        
        lvlayout=QtWidgets.QVBoxLayout()
        lvlayout.minimumSize()
        lvlayout.addWidget(self.display_box)
        lvlayout.addLayout(entry_spec_layout)
        lvlayout.addWidget(self.decimate_box)
        lvlayout.addWidget(self.bv_box)
        lvlayout.addWidget(self.extract_box)
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
        self.ui = sf_main_window()
        self.ui.setup(self)
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(vtk.vtkNamedColors().GetColor3d("slategray"))
        self.ren.GradientBackgroundOn()

        self.file = None
        self.active_dir = os.getcwd()
        self.point_size = 2
        self.picking = False
        self.picked = None

        
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style=vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        # self.iren.AddObserver("MouseMoveEvent", self.on_mouse_move)
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.ui.vtkWidget.Initialize()
        
        #display functions
        self.ui.z_aspect_sb.valueChanged.connect(self.update_z_aspect)
        self.ui.show_outlines_rb.toggled.connect(self.hide_outlines)
        self.ui.caption_rb.toggled.connect(self.hide_caption)
        self.ui.avg_op_slider.valueChanged[int].connect(self.change_avg_opacity)
        self.ui.fit_op_slider.valueChanged[int].connect(self.change_fit_opacity)
        self.ui.hide_avg_sb.toggled.connect(self.hide_avg_legend)
        self.ui.hide_fit_sb.toggled.connect(self.hide_fit_legend)
        
        
        self.ui.entry_spec.activated.connect(self.focus_entry)
        self.ui.show_all_entries_button.clicked.connect(self.draw)
        #decimate functions
        self.ui.undo_last_pick_button.clicked.connect(self.undo_pick)
        self.ui.apply_pick_button.clicked.connect(self.apply_decimate)
        self.ui.reset_pick_button.clicked.connect(self.reset_decimate)
        
        self.ui.fit_spline_button.clicked.connect(self.fit_bv)
        self.ui.extract_button.clicked.connect(self.extract_line)
        self.ui.export_line_button.clicked.connect(self.export_line)
        
        self.ui.save_button.clicked.connect(self.write_output)
        
    def get_data(self):
        if self.file is None:
            self.file, self.active_dir = get_file("*.pyCM",self.active_dir)
        
        
        #make sure valid file was selected
        if self.file is None or not(os.path.isfile(self.file)):
            return
        
        #clear the renderer completely
        self.ren.RemoveAllViewProps()
        
        #try reading data from this step
        self.active_pnt, self.eval_points, self.tri, self.rse, self.bv_tck_list = read_file(self.file)
        self.avg, gsize, _ = aa_read_file(self.file)
        if gsize is None:
            return
        if not self.active_pnt:
            #populate ui entries & initialize remainder
            self.ui.ksx.setValue(6 * gsize)
            self.ui.ksy.setValue(6 * gsize)
            for i in range(len(self.avg)):
                self.active_pnt.append(None)
                self.bv_tck_list.append(None)
                self.eval_points.append(None)
                self.tri.append(None)
                self.rse.append(None)
        else:
            #populate ui entries from tck list of the first entry
            self.ui.ksx.setValue(self.bv_tck_list[0][-1])
            self.ui.ksy.setValue(self.bv_tck_list[0][-2])
        #get outlines
        _, self.outlines, _, _, _, = reg_read_file(self.file)
        

        #move points as needed and initialize actor/polydata lists
        self.fit_actor_list = []
        self.fit_polydata_list = []
        self.avg_surf_polydata_list = []
        self.bool_pnt = []
        self.ui.entry_spec.clear()
        for i in range(len(self.avg)):
            #update average pnts to be only in outline
            self.avg[i][:,-1] = self.avg[i][:,-1] - np.mean(self.avg[i][:,-1])
            self.outlines[i][:,-1] = np.mean(self.avg[i][:,-1])
            self.ui.entry_spec.insertItem(i,'Entry %d'%i)
            #update active to be only in outline - active pnt changes from boolean to an index referencing the total points in self.avg
            self.bool_pnt.append(np.ones(len(self.avg[i]), dtype=bool))
            self.active_pnt[i] = np.arange(0, len(self.bool_pnt[i]), 1, dtype=int)
            #create entries for spline entities as the user may elect to fit entry non-sequentially
            self.fit_actor_list.append(vtk.vtkActor())
            self.fit_polydata_list.append(vtk.vtkPolyData())
            self.avg_surf_polydata_list.append(vtk.vtkPolyData())
        
        self.limits = get_limits(np.vstack(self.avg + self.outlines))
        self.draw_splines()
        self.draw() #resets the interactor to show all
        
        
        self.ui.display_box.setEnabled(True)
        self.ui.decimate_box.setEnabled(True)
    
    def hide_outlines(self):
        
        for actor in self.outline_actor_list:
            flip_visible(actor)
        self.ui.vtkWidget.update()
    
    def hide_caption(self):
        for actor in self.caption_actor_list:
            flip_visible(actor)
        self.ui.vtkWidget.update()
    
    def hide_avg_legend(self):
        if hasattr(self,'sb_actor'):
            flip_visible(self.sb_actor)
            self.ui.vtkWidget.update()
    
    def hide_fit_legend(self):
        if hasattr(self,'sb_fit_actor'):
            flip_visible(self.sb_fit_actor)
            self.ui.vtkWidget.update()
            
    def draw(self):
        
        #remove all avg actors
        if hasattr(self,'avg_actor_list'):
            for actor in self.avg_actor_list:
                self.ren.RemoveActor(actor)
            for actor in self.outline_actor_list:
                self.ren.RemoveActor(actor)
            self.ren.RemoveActor(self.sb_actor)
        #make all potential spline objects visible
        for actor in self.fit_actor_list:
            try: actor.VisibilityOn()
            except: pass
        self.avg_actor_list = []
        self.avg_polydata_list = []
        self.color_list = []
        self.outline_actor_list = []
        self.caption_actor_list = []

        for i in range(len(self.avg)):
            avg = self.avg[i]
            active = self.active_pnt[i]
            pnts = avg[active,:]

            if i == 0: 
                avg_actor, \
                avg_polydata, \
                color, lut = \
                gen_point_cloud(pnts,'blues')
                self.zrange = (self.limits[-2],self.limits[-1])
            else:
                avg_actor, \
                avg_polydata, \
                color, _ = \
                gen_point_cloud(pnts,'blues')
                self.zrange = (self.limits[-2],self.limits[-1])
            
            self.ren.AddActor(avg_actor)
            self.avg_polydata_list.append(avg_polydata)
            self.color_list.append(color)
            self.avg_actor_list.append(avg_actor)
            
            o_actor, _ = gen_outline(self.outlines[i],(255, 255, 255),self.point_size)
            o_actor.SetPickable(0)
            self.outline_actor_list.append(o_actor)
            self.ren.AddActor(o_actor)
            cap_actor = gen_caption_actor('%s'%i, o_actor)
            self.ren.AddActor(cap_actor)
            self.caption_actor_list.append(cap_actor)
            
            
        #handle scalebar
        sb_widget = gen_scalar_bar()
        sb_widget.SetInteractor(self.iren)
        sb_widget.On()
        self.sb_actor = sb_widget.GetScalarBarActor()
        self.sb_actor.SetLookupTable(lut)
        self.ren.AddActor(self.sb_actor)
        if self.ui.hide_avg_sb.isChecked():
            self.sb_actor.VisibilityOff()
        
        self.ui.avg_op_slider.setEnabled(True)
        self.ui.hide_avg_sb.setEnabled(True)
        self.update_z_aspect()

        if self.ui.show_outlines_rb.isChecked():
            self.hide_outlines()
        if not self.ui.caption_rb.isChecked():
            self.hide_caption()

        # self.ui.avg_op_slider.setValue(100)

        self.ren.ResetCamera()
        self.ui.vtkWidget.update()

    def update_z_aspect(self):
        '''
        Updates z_aspect and redraws points displayed based on new z_aspect
        '''
        z_aspect = self.ui.z_aspect_sb.value()

        if hasattr(self,'avg_actor_list'):
            for actor in self.avg_actor_list:
                actor.SetScale((1,1,z_aspect))
                actor.Modified()
        
        if hasattr(self,'fit_actor_list'):
            for actor in self.fit_actor_list:
                actor.SetScale((1,1,z_aspect))
                actor.Modified()
                
        #now do axis_actor, can't scale in the same way as polydata
        if hasattr(self,'axis_actor'):
            if not self.axis_actor.GetVisibility():
                return
        try:
            self.ren.RemoveActor(self.axis_actor) #it will need to be replaced
        except: pass
        #local limits for z axis - don't scale x & y limits, scale the z axis according to z aspect and the current limits
        nl=np.append(self.limits[0:4],([self.limits[-2]*z_aspect,self.limits[-1]*z_aspect]))
        self.axis_actor = get_axis(self.ren, nl, z_aspect)
        self.ren.AddActor(self.axis_actor)
        
        
        self.ui.vtkWidget.update()

    def change_avg_opacity(self,value):
        if hasattr(self,'avg_actor_list'):
            for actor in self.avg_actor_list:
                actor.GetProperty().SetOpacity(value/100)
        self.ui.vtkWidget.update()

    def change_fit_opacity(self,value):
        if hasattr(self,'fit_actor_list'):
            for actor in self.fit_actor_list:
                actor.GetProperty().SetOpacity(value/100)
        self.ui.vtkWidget.update()

    def choose_actor(self):
        '''
        starts pick observer for picking an actor, gets the picked index from the observer
        '''
        
        picker=vtk.vtkPropPicker()
        self.iren.SetPicker(picker)
        self.iren.AddObserver("EndPickEvent",self.check_actor_pick)
        
        if self.picked is not None:
            self.ui.entry_spec.setCurrentIndex(self.picked)

    def check_actor_pick(self,object,event):
        """
        Activates a pick event for actors, returns the index of the picked actor and hides the rest
        """

        click_pos=self.iren.GetEventPosition()
        
        picker = self.iren.GetPicker()
        picker.Pick(click_pos[0],click_pos[1],0,self.ren)
        picked_actor = picker.GetActor()
        highlight_color = tuple(vtk.vtkNamedColors().GetColor3d("Tomato"))
        if picked_actor:
            self.picked = identify_actor(picked_actor,self.avg_actor_list)

        self.iren.RemoveObservers("EndPickEvent")
        
    def focus_entry(self):
        '''
        Hides all actors and resets limits according to the entry in the combo box. Calls draw to reverse any entries
        '''
        # self.draw()
        target = self.ui.entry_spec.currentIndex()
        
        #remove remaining actors
        mask = np.full(len(self.avg_actor_list), True, dtype=bool)
        mask[target] = False
        ind = np.arange(len(self.avg_actor_list),dtype=int)[mask]
        for i in ind:
            
            self.avg_actor_list[i].VisibilityOff()
            self.fit_actor_list[i].VisibilityOff()
            if not self.ui.show_outlines_rb.isChecked():
                self.outline_actor_list[i].VisibilityOff()
        
        self.avg_actor_list[target].VisibilityOn()
        self.fit_actor_list[target].VisibilityOn()
        if not self.ui.show_outlines_rb.isChecked():
            self.outline_actor_list[target].VisibilityOn()

        self.limits = get_limits(np.vstack((self.avg[target],self.outlines[target])))
        self.update_z_aspect()
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()

    def actuate_decimation_pick(self):
        '''
        Starts picking and handles ui button display for decimation
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

    def picker_callback(self,obj,event):
        '''
        Manual picking callback function
        '''
        target = self.ui.entry_spec.currentIndex()
        polydata = self.avg_polydata_list[target]
        
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
        extract.SetInputData(polydata)
        extract.Update()
        extracted = extract.GetOutput()

        ids = vtk.vtkIdTypeArray()
        ids = extracted.GetPointData().GetArray("vtkOriginalPointIds")

        if ids:
            self.last_selected_colors = []
            for i in range(ids.GetNumberOfTuples()):
                self.last_selected_colors.append(self.color_list[target].GetTuple(i))
                self.bool_pnt[target][ids.GetValue(i)] = False
            #store them in an array for an undo operation, and update self.bool_pnt
            self.last_selected_ids=ids
            self.update_mask()

    def update_mask(self):
        '''
        Updates the coloration of poly data on the basis of what's been masked.
        '''
        target = self.ui.entry_spec.currentIndex()
        polydata = self.avg_polydata_list[target]

        localind = np.arange(0, len(self.bool_pnt[target]), 1, dtype=int)
        localind = localind[np.where(np.logical_not(self.bool_pnt[target]))]

        cl = self.color_list[target]
        for i in localind:
            cl.SetTuple(i,(255,99,71))

        polydata.GetPointData().SetScalars(cl)
        polydata.Modified()
        self.avg_actor_list[target].Modified()
        self.ui.vtkWidget.update()

    def undo_pick(self):
        '''
        Undoes the last pick registered
        '''
        target = self.ui.entry_spec.currentIndex()
        if hasattr(self,"last_selected_ids"):
            cl = self.color_list[target]
            polydata = self.avg_polydata_list[target]
            #update colours
            for i in range(self.last_selected_ids.GetNumberOfTuples()):
                cl.SetTuple(self.last_selected_ids.GetValue(i),self.last_selected_colors[i])
                self.bool_pnt[target][self.last_selected_ids.GetValue(i)] = True
            
            polydata.GetPointData().SetScalars(cl)
            polydata.Modified()
            self.avg_actor_list[target].Modified()
            self.ui.vtkWidget.update()
            del self.last_selected_ids, self.last_selected_colors
        else:
            return
    
    def apply_decimate(self):
        '''
        Updates avg pnts based on active_pnt and draws
        '''
        target = self.ui.entry_spec.currentIndex()
        #update active_pnt for redraw
        self.active_pnt[target] = self.active_pnt[target][self.bool_pnt[target]]
        self.bool_pnt[target] = np.ones(len(self.active_pnt[target]), dtype=bool)
        self.draw()
        self.focus_entry()
        
    def reset_decimate(self):
        '''
        Resets points of active entry
        '''
        target = self.ui.entry_spec.currentIndex()
        self.bool_pnt[target]=np.ones(len(self.avg[target]), dtype=bool)
        self.active_pnt[target] = np.arange(0, len(self.bool_pnt[target]), 1, dtype=int)
        self.draw()
        self.focus_entry()
    
    def fit_bv(self):
        '''
        Fits data to bivariate spline according to settings in ui
        FIX - entryspec is added to each time on load causing this to crash.
        '''
        #get input from ui for local processing
        gx = self.ui.ksx.value()
        gy = self.ui.ksy.value()
        kx = self.ui.kox.value()
        ky = self.ui.koy.value()
        
        #get flag for fit entity in focus or all
        if self.ui.fit_all_rb.isChecked():
            ind = [i for i in range(len(self.avg))]
        else:
            ind = [self.ui.entry_spec.currentIndex()]
        
        for i in ind:
            limits = get_limits(self.outlines[i],0)

            tx = np.linspace(\
            limits[0], limits[1],\
            int((limits[1] - limits[0]) / gx))
            ty = np.linspace(\
            limits[2], limits[3],\
            int((limits[3] - limits[2]) / gy))
            if len(tx) < 3:
                self.ui.spline_status_label.setText('Knot spacing too sparse in x for entry %d'%i)
                return
            if len(ty) < 3:
                self.ui.spline_status_label.setText('Knot spacing too sparse in y for entry %d'%i)
                return
            #insert repeated knots at periphery
            tx = np.insert(tx, 0, [limits[0]] * kx)
            tx = np.insert(tx, -1, [limits[1]] * kx)
            ty = np.insert(ty, 0, [limits[2]] * ky)
            ty = np.insert(ty, -1, [limits[3]] * ky)
            x = self.avg[i][self.active_pnt[i],0]
            y = self.avg[i][self.active_pnt[i],1]
            z = self.avg[i][self.active_pnt[i],2]
            try:
                tck = bisplrep(\
                x,y,z,\
                kx=kx, ky=ky, tx=tx, ty=ty, task=-1)
            except ValueError as ve:
                self.ui.spline_status_label.setText('Fit failed on entry %d'%i)
                return
            #append target spacing to tck
            tck.append(gx)
            tck.append(gy)
            self.bv_tck_list[i] = tck
            self.ui.spline_status_label.setText('Calculated spline for entry %d . . .'%i)
            QtWidgets.QApplication.processEvents()
        
            self.evaluate_spline(i)

        self.draw_splines()
        

    def evaluate_spline(self, ind):
        
        self.ui.spline_status_label.setText('Evaluating spline for entry %d . . .'%ind)
        QtWidgets.QApplication.processEvents()
        pts = self.avg[ind]
        tri = Delaunay(pts[:,:2])
        self.tri[ind] = tri.vertices
        tck = self.bv_tck_list[ind]
        znew = np.empty(len(pts))
        for i in range(len(pts)):
            znew[i] = bisplev(pts[i,0],pts[i,1], tck[:-2])
        
        self.eval_points[ind] = np.hstack((self.avg[ind][:,:2],znew[np.newaxis].T))
        self.rse[ind] = np.sqrt((pts[:,-1] - znew)**2)
        
    def draw_splines(self):
        
        for actor in self.fit_actor_list:
            try:
                self.ren.RemoveActor(actor)
            except:
                pass
        lut = None
        for i in range(len(self.eval_points)):
            if self.eval_points[i] is not None:
                if self.ui.mask_average_rb.isChecked():
                    actor, pd, lut = gen_surface(\
                    self.eval_points[i],\
                    self.tri[i],\
                    self.outlines[i],\
                    self.rse[i])
                    #get avg surface polydata with same triangulation
                    _, avg_pd = gen_surface(\
                    self.avg[i],\
                    self.tri[i],\
                    self.outlines[i], vals = None)
                else:
                    actor, pd, lut, = gen_surface(self.eval_points[i],\
                    self.tri[i],\
                    vals = self.rse[i])
                    _, avg_pd = gen_surface(self.avg[i],self.tri[i])
                actor.SetPickable(0)
            else:
                return
            self.fit_actor_list[i] = actor
            self.fit_polydata_list[i] = pd
            self.avg_surf_polydata_list[i] = avg_pd
            self.ren.AddActor(actor)
        #scalebar
        if hasattr(self,'sb_fit_actor'):
            self.ren.RemoveActor(self.sb_fit_actor)
        if lut is not None:
            sb_widget = gen_scalar_bar(side = 'right')
            sb_widget.SetInteractor(self.iren)
            sb_widget.On()
            self.sb_fit_actor = sb_widget.GetScalarBarActor()
            self.sb_fit_actor.SetLookupTable(lut)
            self.ren.AddActor(self.sb_fit_actor)
        
        self.ui.avg_op_slider.setValue(0)
        self.change_avg_opacity(0)
        self.ui.hide_avg_sb.setChecked(True)
        self.hide_avg_legend()
        
        self.update_z_aspect()
        self.ui.vtkWidget.update()
        self.ui.spline_status_label.setText('Ready')
        QtWidgets.QApplication.processEvents()
        
        self.ui.fit_op_slider.setEnabled(True)
        self.ui.hide_fit_sb.setEnabled(True)
        self.ui.extract_box.setEnabled(True)
        self.ui.save_box.setEnabled(True)

    def extract_line(self):
        '''
        Method to extract a line from both avg data and fit, update canvas and generate a line actor
        '''
        
        try: self.ren.RemoveActor(self.line_actor)
        except: pass
        
        target = self.ui.entry_spec.currentIndex()
        
        #get points from ui
        c = np.asarray([self.ui.point1_x_coord.value(), self.ui.point1_y_coord.value(), 0])
        n = np.asarray([self.ui.point2_x_coord.value(), self.ui.point2_y_coord.value(), 0])

        #define a plane on the basis of the cross product of line segment p1-p2 and the z axis
        plane = vtk.vtkPlane()
        plane.SetOrigin(c)
        plane.SetNormal(n)
        
        self.extracted_avg_pnts, _, = \
        extract_contour(plane, self.avg_surf_polydata_list[target])
        self.extracted_fit_pnts, self.line_actor = \
        extract_contour(plane, self.fit_polydata_list[target])
        
        #project points on plane
        v = np.cross(n,[0,0,1])
        avg_proj = np.dot(v[:2],self.extracted_avg_pnts[:,:2].T).T /np.linalg.norm(v[:2])
        fit_proj = np.dot(v[:2],self.extracted_fit_pnts[:,:2].T).T /np.linalg.norm(v[:2])
        
        self.ui.figure.clear()
        ax = self.ui.figure.add_subplot(111)
        ax.plot(avg_proj,self.extracted_avg_pnts[:,2],'.',\
        color=(25/255,25/255,112/255), alpha = 0.5, label="Averaged")
        ax.plot(fit_proj,self.extracted_fit_pnts[:,2],'.',\
        color=(145/255,33/255,158/255), alpha = 0.5, label="Fitted")
        ax.legend(loc="best")
        ax.set_ylabel("z (mm)")
        ax.set_xlabel("Distance along trace (mm)")
        ax.grid(visible=True, which='major', color='#666666', linestyle='-')
        ax.minorticks_on()
        ax.grid(visible=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        self.ui.figure.tight_layout()
        self.ui.canvas.draw()
        
        self.ren.AddActor(self.line_actor)
        self.ui.export_line_button.setEnabled(True)
        self.ui.export_line_status_label.setText("Ready")
        
        self.update_z_aspect()
        self.ui.vtkWidget.update()

    def export_line(self):
        """
        Collects data from ui, writes to a valid file
        """
        
        fileo, _ = get_save_file('*.csv')
        if fileo is None:
            return
        
        data = np.column_stack((self.extracted_avg_pnts, self.extracted_fit_pnts[:,-1]))
        
        np.savetxt(fileo,
        data,
        delimiter=',',
        header = "x, y, z_avg, z_fit",
        fmt='%.4f, %.4f, %.4f, %.4f')
        
        self.ui.export_line_status_label.setText('Exported to %s'%fileo)


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
        elif key == "a":
            flip_visible(self.axis_actor)
        elif key == "p":
            #needs to be replaced with a button
            self.choose_actor()
        elif key == "r":
            self.actuate_decimation_pick()
        if key == "l":
            self.file = None
            self.get_data()
            
        self.ui.vtkWidget.update()

    def check_save_state(self, id):
        '''
        Checks to make sure that there are data objects pertaining to an averaged surface within the interactor against those that might be present in the specified file
        '''
        if not hasattr(self,'eval_points'):
            info_msg('Saving the current step requires a fitted dataset.')
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
            if 'eval_points' in list(g.keys()):
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
        
        passed = self.check_save_state('fit')
        
        tck_entries = ['tx', 'ty', 'c', 'kx', 'ky', 'ksx', 'ksy']
        
        if not passed:
            self.ui.save_label.setText('Ready')
            return
        else:
            self.ui.save_label.setText('Saving . . .')
            QtWidgets.QApplication.processEvents()
            
            with h5py.File(self.file, 'r+') as f:
                if 'fit' in list(f.keys()):
                    del f['fit']
                g = f.create_group('fit')
                for i in range(len(self.avg)):
                    gg = g.create_group(str(i))
                    gg.create_dataset('active',data = self.active_pnt[i])
                    if self.eval_points[i] is not None:
                        #then this fit hasn't been conducted
                        gg.create_dataset('eval_points',data = self.eval_points[i])
                        gg.create_dataset('tri',data = self.tri[i])
                        gg.create_dataset('rse',data = self.rse[i])
                        ggg = gg.create_group('bv_tck')
                        #create dictionary of spline entities for each entry
                        d = {key:value for key,value in zip(tck_entries,self.bv_tck_list[i])}
                        for key, value in d.items():
                            ggg.create_dataset(key, data=value)
                    
                f.attrs['date_modified'] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
                
        self.ui.save_label.setText('Saved to %s'%(os.path.basename(self.file)))

def extract_contour(plane, polydata):
    '''
    Based on vtkPlane and a surface vtkPolyData passed, return points which are on the intersection of the polydata and plane. Also returns a trace actor corresponding to subject polydata
    '''
    
    pd = vtk.vtkPolyData()
    #operate on a deep copy as scalars will be changed
    pd.DeepCopy(polydata)
    
    data = vtk.vtkDoubleArray()
    data.SetNumberOfTuples(pd.GetNumberOfPoints())
    pts = pd.GetPoints()
    for i in range(pd.GetNumberOfPoints()):
        point = pts.GetPoint(i)
        data.SetTuple1(i, plane.EvaluateFunction(point))
    pd.GetPointData().SetScalars(data)
    
    cf = vtk.vtkContourFilter()
    cf.SetInputData(pd)
    cf.ComputeScalarsOff()
    cf.ComputeNormalsOff()
    cf.GenerateValues(1,0,0)
    cf.Update()
    
    extracted_pnts = v2n.vtk_to_numpy(cf.GetOutput().GetPoints().GetData())
    
    a_limits = pd.GetBounds()
    rad = np.maximum(a_limits[1]-a_limits[0],a_limits[3]-a_limits[2])*0.005
    
    tf = vtk.vtkTubeFilter()
    tf.SetInputConnection(cf.GetOutputPort())
    tf.SetRadius(rad) #1% of the maximum dimension
    tf.SetNumberOfSides(25)
    tf.Update()
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tf.GetOutputPort())
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor((255,255,255))
    
    return v2n.vtk_to_numpy(cf.GetOutput().GetPoints().GetData()), actor
    

def gen_surface(pts,tri,poly = None,vals = None):
    '''
    Returns a polydata actor based on points pts and triangulation tri
    '''

    vtk_pnts = vtk.vtkPoints()
    vtk_triangles = vtk.vtkCellArray()
    
    #load up points
    vtk_pnts.SetData(v2n.numpy_to_vtk(pts))
    if poly is not None:
        #calculate centroids of triangles
        centroids = np.array([np.mean(i, axis = 0) for i in pts[tri]])
        inside = in_poly(poly,centroids)
        tri = tri[inside]
    
    for i in tri:
        triangle=vtk.vtkTriangle()
        for j in range(0,3):
            triangle.GetPointIds().SetId(j,i[j])
        vtk_triangles.InsertNextCell(triangle)
        
    triangle_pd = vtk.vtkPolyData()
    triangle_pd.SetPoints(vtk_pnts)
    triangle_pd.SetPolys(vtk_triangles)
    
    if vals is not None:
        vtk_z_array = v2n.numpy_to_vtk(vals)#v2n.numpy_to_vtk(pts[:,-1])
        lut = get_diverging_lut('reds') #add options here for other baseline color series
        lut.SetTableRange(np.amin(vals), np.amax(vals))
        lut.Build()
        colors = lut.MapScalars(vtk_z_array,vtk.VTK_COLOR_MODE_DEFAULT,-1,vtk.VTK_RGB)
        triangle_pd.GetPointData().SetScalars(colors)

    #filter so that edges are shared - not a valid polydata object to pass . . .
    clean_pd = vtk.vtkCleanPolyData()
    clean_pd.SetInputData(triangle_pd)
    
    # Create a mapper and actor for smoothed dataset
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(clean_pd.GetOutputPort())

    surf_actor = vtk.vtkActor()
    surf_actor.SetMapper(mapper)
    surf_actor.GetProperty().SetInterpolationToFlat()
    surf_actor.GetProperty().SetRepresentationToSurface()
    # fit_actor.GetProperty().SetLighting(False)
    if vals is None:
        surf_actor.GetProperty().SetColor(1,0.804,0.204)
        return surf_actor, triangle_pd
    else:
        return surf_actor, triangle_pd, lut

def read_file(file):
    '''
    Reads output file generated by the editor from this module. Returns a list of avg and active entries ready for plotting, along with transform and grid size
    '''
    #create empty lists which are returned if reading fails
    active = []
    eval_points = []
    tri = []
    rse = []
    bv_tck = []

    tck_entries = ['tx', 'ty', 'c', 'kx', 'ky', 'ksx', 'ksy']
    
    with h5py.File(file, 'r') as f:
        try:
            g = f['fit']
            active = [None] * len(g.keys())
            eval_points = [None] * len(g.keys())
            tri = [None] * len(g.keys())
            rse = [None] * len(g.keys())
            bv_tck = [None] * len(g.keys())
            
            for k in g.keys():
                if k.isdigit():
                    active[int(k)] = g['%s/active'%k][()]
                    if '%s/eval_points'%k in g:
                        eval_points[int(k)] = g['%s/eval_points'%k][()]
                        tri[int(k)] = g['%s/tri'%k][()]
                        rse[int(k)] = g['%s/rse'%k][()]
                        
                        #build spline list entry from fields
                        read_tck = []
                        local_bv_tck = g['%s/bv_tck'%k]
                        for entry in tck_entries:
                            read_tck.append(local_bv_tck[entry][()])
                        bv_tck[int(k)] = read_tck
        except: pass
    
    return active, eval_points, tri, rse, bv_tck

if __name__ == "__main__":
    if len(sys.argv)>1:
        launch(sys.argv[1])
    else:
        launch()