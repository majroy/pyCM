#!/usr/bin/env python
'''
pyCM fea_widget - runs external thread to run an FEA from a pyCM output file.
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
import vtk.util.numpy_support as v2n
from PyQt5 import QtGui, QtWidgets, QtCore
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from pkg_resources import Requirement, resource_filename
import yaml
from pyCMcommon import get_file, get_save_file, gen_filtered_ugrid, read_file_for_fea, vtkug_writer
from frd_access import postprocess as frd_post

class execute_fea(QThread):
    '''
    Sets up and runs external thread, emits 100 when done.
    '''
    _signal = pyqtSignal(int)
    def __init__(self,input_file,exe,ccx, post, output):
        super(execute_fea, self).__init__()
        #variables passed here
        self.input_file = input_file[:-4] #strip of suffix
        self.exe = exe #executable path
        self.which_exe = ccx #t/f = ccx or abaqus
        self.run_post = post
        self.results_file = output

    def run(self):
        current_dir = os.getcwd()
        output_dir = os.path.dirname(self.input_file)
        base = os.path.basename(self.input_file)
        os.chdir(output_dir)
        if self.which_exe:
            try:
                print('pyCM exec: %s -i %s'%(self.exe,base))
                out=sp.check_output([self.exe,"-i",base], shell=True)
                print("Calculix output log:")
                print("----------------")
                print(out.decode("utf-8"))
                print("----------------")
                print("pyCM: Calculix run completed . . . Idle")
            except sp.CalledProcessError as e:
                print("Calculix command failed for some reason.")
                print(e)
            if self.run_post:
                frd_process(base)
                replace_grid(base+'.vtu',self.results_file)
        else:
            try:
                print('pyCM exec: %s job=%s int ask_delete=OFF'%(self.exe,base))
                out=sp.check_output([self.exe,"job="+base,"int", "ask_delete=OFF"],shell=True)

                print("Abaqus output log:")
                print("----------------")
                print(out.decode("utf-8"))
                print("----------------")
                print("pyCM: Abaqus run completed . . . Idle")
            except sp.CalledProcessError as e:
                print("Abaqus command failed for some reason.")
                print(e)
            if self.run_post:
                print('pyCM: Running post-processing . . .')
                script_name = resource_filename("pyCM","odb_access.py")
                odb_process(self.exe,script_name,base)
                replace_grid(base+'.vtu',self.results_file)
        
        self._signal.emit(100)
        os.chdir(current_dir)
        
    
class fea_widget(QtWidgets.QDialog):

    def __init__(self, parent, file):
        super(fea_widget, self).__init__(parent)
        self.file = file

        self.setWindowTitle("pyCM - FEA execution: %s"%os.path.basename(self.file))
        self.setWindowFlag(Qt.WindowContextHelpButtonHint, False)
        self.setMinimumSize(QtCore.QSize(450, 200))

        param_layout = QtWidgets.QGridLayout()
        
        self.abq_rb=QtWidgets.QRadioButton("Abaqus")
        self.ccx_rb=QtWidgets.QRadioButton("CalculiX")
        self.ccx_rb.setChecked(True)

        mtype_button_group = QtWidgets.QButtonGroup()
        mtype_button_group.addButton(self.ccx_rb)
        mtype_button_group.addButton(self.abq_rb)
        mtype_button_group.setExclusive(True)
        
        self.merge_results_rb = QtWidgets.QCheckBox('Extract results after run')
        self.merge_results_rb.setToolTip('Extract results and update pyCM file after concluding run.')
        self.merge_results_rb.setChecked(True)
        
        param_layout.addWidget(self.ccx_rb,0,0,1,1)
        param_layout.addWidget(self.abq_rb,0,1,1,1)
        param_layout.addWidget(self.merge_results_rb,0,2,1,1)
        
        self.pbar = QtWidgets.QProgressBar(self, textVisible=True)
        self.pbar.setAlignment(Qt.AlignCenter)
        self.pbar.setFormat("Idle")
        self.pbar.setFont(QtGui.QFont("Helvetica",italic=True))
        self.pbar.setValue(0)

        self.run_button = QtWidgets.QPushButton('Run')
        ccx_exec_path_label = QtWidgets.QLabel('CalculiX executable:')
        self.ccx_exec_path = QtWidgets.QLineEdit()
        ccx_choose_path = QtWidgets.QPushButton('...')
        ccx_choose_path.setMaximumWidth(30)
        ccx_choose_path.setAutoDefault(False)
        abq_exec_path_label = QtWidgets.QLabel('Abaqus executable:')
        self.abq_exec_path = QtWidgets.QLineEdit()
        abq_set_path = QtWidgets.QPushButton('Set')
        abq_set_path.setMaximumWidth(30)
        abq_set_path.setAutoDefault(False)
        work_dir_path_label = QtWidgets.QLabel('Working directory:')
        self.work_dir_path = QtWidgets.QLineEdit()
        wd_choose_path = QtWidgets.QPushButton('...')
        wd_choose_path.setMaximumWidth(30)
        wd_choose_path.setAutoDefault(False)

        run_layout = QtWidgets.QGridLayout()

        run_layout.addWidget(abq_exec_path_label,0,0,1,1)
        run_layout.addWidget(self.abq_exec_path,0,1,1,2)
        run_layout.addWidget(abq_set_path,0,3,1,1)
        
        run_layout.addWidget(ccx_exec_path_label,1,0,1,1)
        run_layout.addWidget(self.ccx_exec_path,1,1,1,2)
        run_layout.addWidget(ccx_choose_path,1,3,1,1)
        
        run_layout.addWidget(work_dir_path_label,2,0,1,1)
        run_layout.addWidget(self.work_dir_path,2,1,1,2)
        run_layout.addWidget(wd_choose_path,2,3,1,1)
        run_layout.addWidget(self.run_button,3,0,1,1)
        run_layout.addWidget(self.pbar,3,1,1,3)

        self.run_button.clicked.connect(self.run_fea)
        ccx_choose_path.clicked.connect(self.set_ccx)
        abq_set_path.clicked.connect(self.set_abq)
        wd_choose_path.clicked.connect(self.set_wd)
        
        self.layout = QtWidgets.QVBoxLayout()
        self.layout.addLayout(param_layout)
        self.layout.addLayout(run_layout)

        self.setLayout(self.layout)
        self.read_config()
        self.show()

    
    
    def run_fea(self):
        self.make_config_change() #save anything that the user might have put into config boxes
        
        self.input_file, _ = get_save_file("*.inp", self.work_dir_path.text())
        if self.input_file is None:
            return
        
        run_ccx = False
        if self.ccx_rb.isChecked():
            run_ccx = True

        run_post = False
        if self.merge_results_rb.isChecked():
            run_post = True
        
        #get the mesh from self.file
        _, bc_disp_nodes, bc_disp, bc_nodes, _, mod, pr, mesh = read_file_for_fea(self.file)
        
        write_inp(bc_disp_nodes,\
        bc_disp,\
        bc_nodes,\
        mod, pr, mesh, run_ccx, self.input_file)
        
        if run_ccx:
            self.thread = execute_fea(self.input_file,self.ccx_exec_path.text(),run_ccx, run_post, self.file)
        else:
            self.thread = execute_fea(self.input_file,self.abq_exec_path.text(),run_ccx, run_post, self.file)
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
    
    def set_abq(self):
        self.make_config_change()
    
    def set_ccx(self):
        f,_ = get_file("*.*")
        if f is None or not(os.path.isfile(f)):
            return
        self.ccx_exec_path.setText(f)
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
            read = yaml.safe_load(f)
        
        self.ccx_exec_path.setText(read['FEA']['ccx_exec'])
        self.abq_exec_path.setText(read['FEA']['abq_exec'])
        self.work_dir_path.setText(read['FEA']['work_dir'])


    def make_config_change(self):
        new_entries = dict(
        ccx_exec = str(self.ccx_exec_path.text()),
        abq_exec = str(self.abq_exec_path.text()),
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
        Not implemented
        '''
        pass

def odb_process(executable,script_name,base):
    '''
    Calls odb_access using sp as it is an Abaqus Python script and needs to be called by the Abaqus Python executable.
    '''
    
    try:
        print('pyCM exec: %s python %s %s %s'%(executable, script_name, base+'.odb', base+'.vtu'))
        out=sp.check_output([executable,"python",script_name, base+'.odb', base+'.vtu'],shell=True)
        print("----------------")
        print(out.decode("utf-8"))
        print("----------------")
    except sp.CalledProcessError as e:
        print("Abaqus post-processing failed for some reason.")
        print(e)
        
def frd_process(base):
    '''
    Calls frd_access directly
    '''
    try:
        frd_post(base+'.frd',base+'.vtu')
    except Exception as e:
        print("CalculiX post-processing failed for some reason.")
        print(e)
    
def replace_grid(vtu_file, pyCM_file):
    '''
    Reads a vtu file and replaces the unstructured grid object in a pyCM file with it
    '''

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_file)
    reader.Update()
    
    w = vtkug_writer()
    w.SetInputConnection(reader.GetOutputPort())
    w.SetFileName(pyCM_file)
    w.Update()
    print('pyCM: Merged FEA results to file.')

def write_inp(bc_disp_nodes, bc_disp, bc_nodes, mod, pr, mesh, ccx, file):
    '''
    Writes a contour method input deck either for ccx or abaqus based on parameters
    '''
    
    mesh = gen_filtered_ugrid(mesh)
    
    nodes = v2n.vtk_to_numpy(mesh.GetPoints().GetData())
    nodes=np.column_stack((np.arange(1,len(nodes)+1),nodes))

    cells = v2n.vtk_to_numpy(mesh.GetCells().GetData())
    #determine element type based on first entry in cells, 8-C3D8, 10-C3D10
    el_type = cells[0]
    
    cells=np.resize(cells+1,(len(cells)//(el_type+1),el_type+1))
    cells=np.column_stack((np.arange(1,len(cells)+1),cells[:,1::]))
    
    fid = io.BytesIO()
    fid.write(str.encode('*HEADING\n'))
    fid.write(str.encode('**pyCM input deck\n'))
    fid.write(str.encode('**%s\n'%file))
    
    fid.write(str.encode('*NODE\n'))
    np.savetxt(fid,nodes,fmt='%i,%.6f,%.6f,%.6f',delimiter=',')
    fid.write(str.encode('*ELEMENT, TYPE=C3D%i\n'%el_type))
    np.savetxt(fid,cells,fmt='%i',delimiter=',')
    fid.write(str.encode('*ELSET, ELSET=DOMAIN, GENERATE\n'))
    fid.write(str.encode('%i,%i,%i\n'%(1,len(cells),1)))
    fid.write(str.encode('*SOLID SECTION, ELSET=DOMAIN, MATERIAL=USERSPEC\n'))
    fid.write(str.encode('*MATERIAL, NAME=USERSPEC\n'))
    fid.write(str.encode('*ELASTIC, TYPE=ISO\n'))
    fid.write(str.encode('%7.0f,%.3f\n'%(mod, pr)))
    if ccx:
        fid.write(str.encode('*STEP\n'))
    else:
        fid.write(str.encode('*STEP, NAME=CONFORM\n'))
    fid.write(str.encode('*STATIC\n'))
    fid.write(str.encode('*BOUNDARY\n'))
    fid.write(str.encode('%i, 1,2, 0\n'%(bc_nodes[0]+1)))
    fid.write(str.encode('%i, 2, 0\n'%(bc_nodes[1]+1)))
    fid.write(str.encode('*BOUNDARY\n'))
    for entry in range(len(bc_disp_nodes)):
        if bc_disp_nodes[entry] is None:
            continue
        node_list = bc_disp_nodes[entry]
        node_disp = bc_disp[entry]
        for ind in range(len(node_list)):
            fid.write(str.encode('%i, 3,, %6.6f\n'%(node_list[ind]+1,node_disp[ind])))
    if ccx:
        fid.write(str.encode('*NODE FILE\n'))
        fid.write(str.encode('U\n'))#Coords by default
        fid.write(str.encode('*EL FILE\n'))
        fid.write(str.encode('S\n'))#Coords by default
    else:
        fid.write(str.encode('*EL PRINT\n'))
        fid.write(str.encode('COORD,S\n'))#have to specify coords
    fid.write(str.encode('*ENDSTEP'))

    with open(file, 'w+') as f: 
        f.write(fid.getvalue().decode("utf-8"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    if len(sys.argv)>1:
        fw = fea_widget(None, sys.argv[1])
    else:
        fw = fea_widget(None, 'No file specified')
    sys.exit(app.exec_())