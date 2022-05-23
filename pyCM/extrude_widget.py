#!/usr/bin/env python
'''
pyCM mesh_widget - runs external thread to extrude and mesh an outline, either with GMSH or Abaqus external executables.
'''

__author__ = "M.J. Roy"
__version__ = "0.1"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014--"

import os, io
import subprocess as sp
import numpy as np
import vtk
from PyQt5 import QtGui, QtWidgets, QtCore, QtSvg
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from pkg_resources import Requirement, resource_filename
import yaml
from pyCMcommon import get_file, get_save_file, write_dxf, gen_filtered_ugrid

class run_gmsh(QThread):
    '''
    Sets up and runs external thread, emits 100 when done.
    '''
    _signal = pyqtSignal(int)
    def __init__(self,input_file,exe,quads):
        super(run_gmsh, self).__init__()
        #variables passed here
        self.input_file = input_file
        self.exe = exe
        self.quads = quads

    def run(self):
        current_dir = os.getcwd()
        output_dir = os.path.dirname(self.input_file)
        base = os.path.basename(self.input_file)
        self.output_file = os.path.splitext(base)[0]+'.vtk'
        os.chdir(output_dir)
        try:
            if self.quads:
                out=sp.check_output([self.exe,"-3",self.input_file,"-o",self.output_file], shell=True)
                main_cell_type = 12 #linear quads
            else:
                # make sure second order tets are generated
                out=sp.check_output([self.exe,"-3","-order","2",self.input_file,"-o",self.output_file], shell=True)
                main_cell_type = 24 #quadratic tets
            print(out.decode("utf-8"))
            print('pyCM: filtering non-volumetric elements . . .')
            #filter resulting vtk file with threshold to get rid of non-volumetric parts
            mesh_source = vtk.vtkUnstructuredGridReader()
            mesh_source.SetFileName(self.output_file)
            mesh_source.Update()
            mesh_as_read = mesh_source.GetOutput()

            cell_type_array = vtk.vtkIntArray()
            cell_type_array.SetName("Type")
            cell_type_array.SetNumberOfComponents(1)
            cell_type_array.SetNumberOfTuples(mesh_as_read.GetNumberOfCells())
            mesh_as_read.GetCellData().AddArray(cell_type_array)

            threshold = vtk.vtkThreshold()
            threshold.SetInputData(mesh_as_read)
            threshold.SetUpperThreshold(main_cell_type)
            threshold.SetThresholdFunction(2) #upper
            threshold.SetInputArrayToProcess(0,0,0,1,"Type")
            threshold.Update()
            mesh = threshold.GetOutput()
            mesh.GetCellData().RemoveArray("Type")
            
            clean_mesh = gen_filtered_ugrid(mesh)
            
            writer = vtk.vtkUnstructuredGridWriter()
            writer.SetInputData(clean_mesh)
            writer.SetFileName(self.output_file)
            writer.Update()
            print('pyCM: Done.')
        except sp.CalledProcessError as e:
            print("Gmsh command failed for some reason.")
            print(e)
        self._signal.emit(100)
        os.chdir(current_dir)
    
class extrude_widget(QtWidgets.QDialog):

    def __init__(self, parent, outline, dist, nn=2, length=10):
        super(extrude_widget, self).__init__(parent)
        
        self.outline = outline
        self.dist = dist
        self.nn = nn
        self.length = length
        
        self.setWindowTitle("pyCM - Extruded mesh generator" )
        self.setWindowFlag(Qt.WindowContextHelpButtonHint, False)
        self.setMinimumSize(QtCore.QSize(450, 400))
        
        svg = QtSvg.QSvgWidget(resource_filename("pyCM","meta/extrude_def.svg"))
        svg.renderer().setAspectRatioMode(Qt.KeepAspectRatio)
        param_layout = QtWidgets.QGridLayout()
        num_nodes_label = QtWidgets.QLabel('Node count:')
        self.num_nodes_sb = QtWidgets.QSpinBox()
        self.num_nodes_sb.setPrefix('N = ')
        self.num_nodes_sb.setToolTip('Number of nodes along extrusion direction')
        self.num_nodes_sb.setValue(nn)
        self.extrude_depth_sb = QtWidgets.QDoubleSpinBox()
        extrude_depth_label = QtWidgets.QLabel('Extrusion length:')
        self.extrude_depth_sb.setPrefix('L = ')
        self.extrude_depth_sb.setSuffix(' mm')
        self.extrude_depth_sb.setToolTip('Length of extrusion')
        self.extrude_depth_sb.setValue(length)
        self.min_dist_sb = QtWidgets.QDoubleSpinBox()
        elem_length_label = QtWidgets.QLabel('Minimum element length:')
        self.min_dist_sb.setPrefix('l = ')
        self.min_dist_sb.setSuffix(' mm')
        self.min_dist_sb.setMinimum(0.01)
        self.min_dist_sb.setDecimals(3)
        self.min_dist_sb.setValue(dist)
        self.min_dist_sb.setToolTip('Minimum element length')
        
        self.quad_rb=QtWidgets.QRadioButton("Quads")
        self.tet_rb=QtWidgets.QRadioButton("Tets")
        self.quad_rb.setChecked(True)

        mtype_button_group = QtWidgets.QButtonGroup()
        mtype_button_group.addButton(self.tet_rb)
        mtype_button_group.addButton(self.quad_rb)
        mtype_button_group.setExclusive(True)
        
        param_layout.addWidget(num_nodes_label,0,0,1,1)
        param_layout.addWidget(self.num_nodes_sb,0,1,1,1)
        param_layout.addWidget(extrude_depth_label,1,0,1,1)
        param_layout.addWidget(self.extrude_depth_sb,1,1,1,1)
        param_layout.addWidget(elem_length_label,2,0,1,1)
        param_layout.addWidget(self.min_dist_sb,2,1,1,1)
        param_layout.addWidget(self.quad_rb,3,0,1,1)
        param_layout.addWidget(self.tet_rb,3,1,1,1)
        
        
        self.pbar = QtWidgets.QProgressBar(self, textVisible=True)
        self.pbar.setAlignment(Qt.AlignCenter)
        self.pbar.setFormat("Idle")
        self.pbar.setFont(QtGui.QFont("Helvetica",italic=True))
        self.pbar.setValue(0)

        self.run_button = QtWidgets.QPushButton('Run')
        gmsh_exec_path_label = QtWidgets.QLabel('Path to Gmsh executable:')
        self.gmsh_exec_path = QtWidgets.QLineEdit()
        gmsh_choose_path = QtWidgets.QPushButton('...')
        gmsh_choose_path.setMaximumWidth(20)
        work_dir_path_label = QtWidgets.QLabel('Working directory:')
        self.work_dir_path = QtWidgets.QLineEdit()
        wd_choose_path = QtWidgets.QPushButton('...')
        wd_choose_path.setMaximumWidth(20)

        run_layout = QtWidgets.QGridLayout()
        run_layout.addWidget(self.run_button,2,0,1,1)
        run_layout.addWidget(gmsh_exec_path_label,0,0,1,1)
        run_layout.addWidget(self.gmsh_exec_path,0,1,1,2)
        run_layout.addWidget(gmsh_choose_path,0,3,1,1)
        run_layout.addWidget(work_dir_path_label,1,0,1,1)
        run_layout.addWidget(self.work_dir_path,1,1,1,2)
        run_layout.addWidget(wd_choose_path,1,3,1,1)
        run_layout.addWidget(self.pbar,2,1,1,3)

        self.run_button.clicked.connect(self.run_mesh)
        gmsh_choose_path.clicked.connect(self.set_gmsh)
        wd_choose_path.clicked.connect(self.set_wd)
        
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addWidget(svg)
        vertical_spacer1 = QtWidgets.QSpacerItem(20, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.layout.addItem(vertical_spacer1)
        self.layout.addLayout(param_layout)
        vertical_spacer2 = QtWidgets.QSpacerItem(10, 10, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.layout.addItem(vertical_spacer2)
        self.layout.addLayout(run_layout)

        self.setLayout(self.layout)
        self.read_config()
        self.show()

    def run_mesh(self):
        self.make_config_change() #save anything that the user might have put into config boxes
        
        self.input_file, _ = get_save_file("*.geo", self.work_dir_path.text())
        if self.input_file is None:
            return
        
        if self.tet_rb.isChecked():
            quads = False
        else:
            quads = True
            
        write_geo(self.outline,\
        self.num_nodes_sb.value(),\
        self.extrude_depth_sb.value(),\
        self.min_dist_sb.value(), self.input_file, quads)
        
        self.thread = run_gmsh(self.input_file,self.gmsh_exec_path.text(),quads)
        self.thread._signal.connect(self.signal_accept)
        self.thread.start()
        self.pbar.setTextVisible(True)
        self.pbar.setStyleSheet("")
        self.pbar.setRange(0,0)
        
    def signal_accept(self, msg):
        if int(msg) == 100:
            self.pbar.setRange(0,100)
            self.pbar.setValue(0)
            self.pbar.setFormat("Complete")
            self.pbar.setStyleSheet("QProgressBar"
              "{"
              "background-color: lightgreen;"
              "border : 1px"
              "}")

    def set_gmsh(self):
        f,_ = get_file("*.*")
        if f is None or not(os.path.isfile(f)):
            return
        self.gmsh_exec_path.setText(f)
        self.make_config_change()

    def set_wd(self):
        dir = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory"))
        if dir != '':
            self.work_dir_path.setText(dir)
            self.make_config_change()
        else:
            return

    def read_config(self):
        fname=resource_filename("pyCM","meta/pyCMconfig.yml")
        with open(fname, 'r') as f:
            cur_yaml = yaml.safe_load(f)
        
        self.gmsh_exec_path.setText(cur_yaml['FEA']['gmsh_exec'])
        self.work_dir_path.setText(cur_yaml['FEA']['work_dir'])


    def make_config_change(self):
        new_entries = dict(
        gmsh_exec = str(self.gmsh_exec_path.text()),
        work_dir = str(self.work_dir_path.text())
        )
        fname=resource_filename("pyCM","meta/pyCMconfig.yml")
        with open(fname, 'r') as yamlfile:
            cur_yaml = yaml.safe_load(yamlfile)
            cur_yaml['FEA'].update(new_entries)
        if cur_yaml:
            with open(fname, 'w') as yamlfile:
                yaml.safe_dump(cur_yaml, yamlfile)

    def closeEvent(self, event):
        '''
        Method from Qdialog
        '''
        try:
            self.vtk_file = os.path.join(self.work_dir_path.text(),self.thread.output_file)
        except:
            self.vtk_file = None

def write_geo(outline, nn, extrude_depth, dist, fileo, quads = True):
    '''
    Writes a gmsh geo file giving instructions to extrude the given outline
    '''
    
    fid=io.StringIO()
    
    n = len(outline)
    #calculate bias of node spacing
    bias_coeff = (extrude_depth / float(dist))**(1/float(nn - 1)) #upper bound
    bias_range = np.linspace(bias_coeff/2,bias_coeff,1000)
    intersection = dist*(1 - np.power(bias_range, nn))/(1 - bias_range)
    bias_index = np.where(intersection > extrude_depth)
    bias = bias_range[bias_index[0][0]]
    
    l = np.array([])
    for i in range(1, nn+1):
        l = np.append(l, dist*(1-bias**i)/(1-float(bias)))
    l[-1] = extrude_depth
    l = l/float(extrude_depth)
    # l = np.flip(l)
    
    #write points and lines corresponding to outline
    pc = 0
    for i in range(n):
        pc += 1
        output_line  = "Point(%i) = {%8.8f,%8.8f,%8.8f};\n" \
                % (pc, outline[i,0],outline[i,1],0)
        fid.write(output_line)
    lc = 0
    for i in range(n-1):
        lc+=1
        output_line = "Line(%i) = {%i,%i};\n"%(lc,lc,lc+1)
        fid.write(output_line)
    lc+=1
    fid.write("Line(%i) = {%i,%i};\n"%(lc,lc,lc-n+1)) #last line to enclose
    
    #indicate that the outline is looped
    fid.write("Line Loop(%i) = {%i:%i};\n" %(n+1,1,lc))
    
    sec=n+1
    #generate plane
    fid.write("Plane Surface(%i) = {%i};\n" %(n+2,n+1))
    
    #recombine triangles if quads is true
    if quads:
        fid.write("Recombine Surface {%i};\n\n" %(n+2)) #for quads, otherwise tets
    sec+=1
    
    #extrusion keywords
    fid.write("OutOfPlane[]= Extrude {0, 0, %8.8f} {\n Surface{%i};\n Layers{ {"%(extrude_depth, sec)) 

    for i in range(len(l)-1):
        fid.write("1,")
    fid.write("1}, {")
    for i in range(len(l)-1):
        fid.write("%2.4f,"%l[i])
    if quads: #quads vs. tets
        fid.write("%2.4f} };\n Recombine;};\n \n//EOF"%l[-1])
    else:
        fid.write("%2.4f} };};\n \n//EOF"%l[-1])

    with open(fileo, 'w+') as f: f.write(fid.getvalue())

def write_abaqus_py(outline, nn, extrude_depth, dist, fileo, quads = True):

    if quads:
        elem_type = "C3D8"
    else:
        elem_type = "C3D10"
    

    base = os.path.basename(fileo)
    dxf_file = os.path.splitext(base)[0]+'.dxf'
    write_dxf(dxf_file,outline)

    fid = io.StringIO()

    s1="""
# RawInpWriter.py
# Abaqus python script to automatically generate an extruded mesh 
# Geometry and mesh is based on *.dxf file
# Produced by pyCM
####################################################################
import os
from abaqus import *
from abaqusConstants import *
from abaqus import backwardCompatibility
backwardCompatibility.setValues(reportDeprecated=False)
from caeModules import *
from driverUtils import executeOnCaeStartup
from dxf2abq import importdxf
# Parameters employed:\n"""

    s2="""
DXF_file=os.path.normpath(DXF_file)
executeOnCaeStartup()
Mdb()
importdxf(fileName=DXF_file)
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches[os.path.basename(OutputFname)])
p = mdb.models['Model-1'].Part(name='Part-1', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=Depth)
s.unsetPrimaryObject()
f = p.faces
faces = f.findAt((CentPoint, ))
p.Set(faces=faces, name='SURFACE')
del mdb.models['Model-1'].sketches['__profile__']
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=OFF)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
a = mdb.models['Model-1'].rootAssembly
e1 = a.instances['Part-1-1'].edges
pickedEdges2 = e1.findAt((EdgePoint, ))
a.seedEdgeByBias(biasMethod=SINGLE, end1Edges=pickedEdges2, minSize=MinLength, 
    maxSize=MaxLength, constraint=FINER)
elemType1 = mesh.ElemType(elemCode=ElemType, elemLibrary=STANDARD)
a = mdb.models['Model-1'].rootAssembly
c1 = a.instances['Part-1-1'].cells
cells1 = c1.getSequenceFromMask(mask=('[#1 ]', ), )
pickedRegions =(cells1, )
"""

    s3="""
a.setElementType(regions=pickedRegions, elemTypes=(elemType1,))
a = mdb.models['Model-1'].rootAssembly
partInstances =(a.instances['Part-1-1'], )
a.generateMesh(regions=partInstances)
mdb.models['Model-1'].setValues(noPartsInputFile=ON)
mdb.Job(name=OutputFname, model='Model-1', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', parallelizationMethodExplicit=DOMAIN, numDomains=1, 
    activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=1)
mdb.jobs[OutputFname].writeInput(consistencyChecking=OFF)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])
session.viewports['Viewport: 1'].view.fitView()\n"""


    fid.write("%s"%s1)
    fid.write("DXF_file=r'%s'\n" %dxf_file)
    fid.write("Depth=%8.8f\n" %extrude_depth);
    fid.write("NN=%i\n" %nn);
    fid.write("MinLength=%8.8f\n" %(dist));
    fid.write("MaxLength=%8.8f\n" %(dist*extrude_depth));
    fid.write("EdgePoint=(%8.8f,%8.8f,%8.8f)\n" %(outline[0,0],outline[0,1],dist*extrude_depth));
    cent = np.mean(outline,axis=0)
    fid.write("CentPoint=(%8.8f,%8.8f,%8.8f)\n" %(cent[0],cent[1],0));
    fid.write("OutputFname='%s'\n"
    %(os.path.splitext(os.path.basename(fileo))[0]))
    fid.write("ElemType=%s\n"%elem_type)
    fid.write("%s"%s2)
    
    if quads:
        fid.write("a.setMeshControls(regions=cells1, algorithm=ADVANCING_FRONT)\n")
    else:
        fid.write("a.setMeshControls(regions=cells1, elemShape=TET, technique=FREE)\n")
    fid.write("%s"%s3)
    
    with open(fileo, 'w+') as f: f.write(fid.getvalue())

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    w = extrude_widget(None, None, 0, 0, 0)
    sys.exit(app.exec_())