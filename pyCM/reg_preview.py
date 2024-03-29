#!/usr/bin/env python
'''
Preview widget for pyCM point/outline registration
'''

import sys
import numpy as np
import vtk
from PyQt5 import QtCore, QtGui, QtWidgets
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from pyCM.pyCMcommon import *

class registration_viewer(QtWidgets.QWidget):
    def __init__(self, parent = None):
        super(registration_viewer, self).__init__(parent)
        
        vl = QtWidgets.QVBoxLayout()
        
        self.vtkWidget = QVTKRenderWindowInteractor(parent)
        
        left_viewport = [0, 0, 0.5, 1] #[xmin max ymin ymax]
        right_viewport = [0.5, 0, 1, 1]
        self.preview1_ren = vtk.vtkRenderer()
        self.preview2_ren = vtk.vtkRenderer()
        self.preview1_ren.SetBackground(vtk.vtkNamedColors().GetColor3d("aliceblue"))
        self.preview1_ren.GradientBackgroundOn()
        self.preview2_ren.SetBackground(vtk.vtkNamedColors().GetColor3d("aliceblue"))
        self.preview2_ren.GradientBackgroundOn()
        self.preview1_ren.SetViewport(left_viewport)
        self.preview2_ren.SetViewport(right_viewport)
        self.vtkWidget.GetRenderWindow().AddRenderer(self.preview1_ren)
        self.vtkWidget.GetRenderWindow().AddRenderer(self.preview2_ren)
        
        self.picking = False
        self.picked = {}
        
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.iren.SetInteractorStyle(style)
        self.preview1_ren.GetActiveCamera().ParallelProjectionOn()
        self.preview2_ren.SetActiveCamera(self.preview1_ren.GetActiveCamera())
        self.setLayout(vl)
        vl.addWidget(self.vtkWidget)
        self.initialize_vtk()
        
    def initialize_vtk(self):
        '''
        Run-once setup/definition of interactor specifics
        '''
        
        self.ref_pnts = self.ref_outlines = self.float_pnts = self.float_outlines = []
        
        axes_actor = vtk.vtkAxesActor()
        self.axes = vtk.vtkOrientationMarkerWidget()
        self.axes.SetOrientationMarker(axes_actor)
        
        self.iren.AddObserver("KeyPressEvent",self.keypress)
        self.axes.SetInteractor(self.iren)
        
        self.axes.EnabledOn()
        self.axes.InteractiveOn()

        self.preview1_ren.ResetCamera()
        self.iren.Initialize()
        
    def draw(self):
        '''
        Draws whatever has been passed as ref_pnts etc
        '''
        if not self.ref_pnts:
            info_msg('Need to have at least reference points registered for a preview.')
            return
        
        self.preview1_ren.RemoveAllViewProps()
        self.preview2_ren.RemoveAllViewProps()
        
        #get limits of ref points:
        limits = get_limits(np.vstack(self.ref_pnts+self.float_pnts))
        
        self.ref_point_actor_list = [] 
        self.float_point_actor_list = []
        for e in range(len(self.ref_pnts)):
            point_actor, \
            _, \
            _, lut = \
            gen_point_cloud(self.ref_pnts[e], None, (limits[-2],limits[-1]))
            outline_actor, _ = gen_outline(self.ref_outlines[e])
            outline_actor.SetPickable(0)
            cap_actor = gen_caption_actor('%s'%e, outline_actor)
            cap_actor.SetPickable(0)
            self.ref_point_actor_list.append(point_actor)
            self.ref_point_actor_list.append(outline_actor)
            self.ref_point_actor_list.append(cap_actor)
        
        for e in range(len(self.float_pnts)):
            point_actor, \
            _, \
            _, _ = \
            gen_point_cloud(self.float_pnts[e], None, (limits[-2],limits[-1]))
            outline_actor, _ = gen_outline(self.float_outlines[e])
            outline_actor.SetPickable(0)
            cap_actor = gen_caption_actor('%s'%e, outline_actor)
            cap_actor.SetPickable(0)
            self.float_point_actor_list.append(point_actor)
            self.float_point_actor_list.append(outline_actor)
            self.float_point_actor_list.append(cap_actor)
        
        for actor in self.ref_point_actor_list:
            self.preview1_ren.AddActor(actor)
        for actor in self.float_point_actor_list:
            self.preview2_ren.AddActor(actor)

        ref_label_actor = gen_info_actor('Reference', self.preview1_ren)
        float_label_actor = gen_info_actor('Floating', self.preview2_ren)
        self.preview1_ren.AddActor(ref_label_actor)
        self.preview2_ren.AddActor(float_label_actor)
        
        #cutting orientations
        cut_attr = self.cut_attr.copy()
        if cut_attr['ref']:
            ld = cut_attr['ref']
            if ld['cut_path'].any():
                ref_cut_orient_actor = gen_cutting_orientation_actor(\
                limits,\
                ld['cut_dir'],\
                ld['cut_path'])
            else:
                ref_cut_orient_actor = gen_cutting_orientation_actor(\
                limits,\
                ld['cut_dir'])
            self.preview1_ren.AddActor(ref_cut_orient_actor)
            
        if cut_attr['float']:
            ld = cut_attr['float']
            if ld['cut_path'].any():
                float_cut_orient_actor = gen_cutting_orientation_actor(\
                limits,\
                ld['cut_dir'],\
                ld['cut_path'])
            else:
                float_cut_orient_actor = gen_cutting_orientation_actor(\
                limits,\
                ld['cut_dir'])
            self.preview2_ren.AddActor(float_cut_orient_actor)
        
        self.preview1_ren.ResetCamera()
        
        self.vtkWidget.update()

    def actuate_entry_pick(self):
        
        if self.picking:
            self.iren.RemoveObservers('LeftButtonPressEvent')
            self.iren.AddObserver('LeftButtonPressEvent',self.default_left_button)
            self.picking = False
            QtWidgets.QApplication.processEvents()
            
        else:
            self.iren.AddObserver('LeftButtonPressEvent', self.picker_callback)
            self.picking = True
    
    def picker_callback(self,obj,event):
        
        #reset opacity of all actors pending a new pick
        for actor in self.ref_point_actor_list:
            actor.GetProperty().SetOpacity(1)
        for actor in self.float_point_actor_list:
            actor.GetProperty().SetOpacity(1)
        
        self.picked = {}
        float_picker=vtk.vtkPropPicker()
        ref_picker = vtk.vtkPropPicker()
        
        click_pos=self.iren.GetEventPosition()
        
        ref_picker.Pick(click_pos[0],click_pos[1],0,self.preview1_ren)
        float_picker.Pick(click_pos[0],click_pos[1],0,self.preview2_ren)
        float_picked_actor = float_picker.GetActor()
        ref_picked_actor = ref_picker.GetActor()
        
        if ref_picked_actor:
           # picked = ref_picked_actor
           picked = identify_actor(ref_picked_actor,self.ref_point_actor_list)
           self.ref_point_actor_list[picked].GetProperty().SetOpacity(0.2)
           self.picked = {'Ref': picked//3}
           self.actuate_entry_pick()

        if float_picked_actor:
           # picked = float_picked_actor
           picked = identify_actor(float_picked_actor,self.float_point_actor_list)
           self.float_point_actor_list[picked].GetProperty().SetOpacity(0.2)
           self.picked = {'Float': picked//3}
           self.actuate_entry_pick()
        
        self.iren.RemoveObservers("EndPickEvent")
        self.vtkWidget.update()

    def default_left_button(self, obj, event):
        #forward standard events according to the default style`
        self.iren.GetInteractorStyle().OnLeftButtonDown()

    def keypress(self,obj,event):
        key = obj.GetKeyCode()
        if key == "1":
            xyview(self.preview1_ren)
        elif key == "p":
            self.actuate_entry_pick()
        self.vtkWidget.update()

    def closeEvent(self,event):
        '''
        Finalize all VTK widgets to negate OpenGL messages/errors
        '''
        self.vtkWidget.close()

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    widget = registration_viewer()
    widget.show()
    sys.exit(app.exec_())