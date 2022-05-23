#!/usr/bin/env python
'''
Decimate widget for pyCM
'''

import sys
from PyQt5 import QtCore, QtGui, QtWidgets

class decimate_ui(QtWidgets.QWidget):
    def __init__(self, parent = None):
        super(decimate_ui, self).__init__(parent)
        
        hl = QtWidgets.QHBoxLayout()
        
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(100)
        sizePolicy.setVerticalStretch(100)
        
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
        
        hl.addWidget(self.point_processing_box)
        self.setLayout(hl)
        

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    widget = decimate_ui()
    widget.show()
    sys.exit(app.exec_())