#!/usr/bin/env python
'''
VTK and QT to allow for preprocessing FEAs associated with the contour method. 
See documentation for the following main methods in this package:
-point_cloud
-align_average
-fit_surface
-preprocess
-------------------------------------------------------------------------------
1.4 - loading between tabs improved
'''
__author__ = "M.J. Roy"
__version__ = "1.4"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2019"

import sys, os.path, shutil
import subprocess as sp
from pkg_resources import Requirement, resource_filename
import numpy as np
import scipy.io as sio
import vtk
import vtk.util.numpy_support as vtk_to_numpy
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
#Change following to local import for dev
from pyCM.pyCMcommon import *

import pyCM.point_cloud as pc
import pyCM.align_average as aa
import pyCM.fit_surface as sf
import pyCM.preprocess as pre
import pyCM.postprocess as post




class Ui_MainWindow(object):
    def setupUi(self, MainWindow, width, height):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setWindowTitle("pyCM - main v%s" %__version__)
        MainWindow.setWindowIcon(QtGui.QIcon(resource_filename("pyCM","meta/pyCM_icon.png")))
        MainWindow.setEnabled(True)
        MainWindow.resize(width, height)
        MainWindow.showMaximized()

        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.mainLayout = QtWidgets.QHBoxLayout(self.centralwidget)
        self.mainLayout.addWidget(self.tabWidget)
        MainWindow.setCentralWidget(self.centralwidget)
        self.tabLayout = QtWidgets.QHBoxLayout(self.tabWidget)

        #add tabs
        self.pctab = QtWidgets.QWidget(self.tabWidget)
        self.tabWidget.addTab(self.pctab, "Point cloud editor")

        self.aatab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.aatab, "Alignment/averaging")
        self.sftab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.sftab, "Surface fitting")
        self.pretab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.pretab, "Preprocessing")
        self.posttab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.posttab, "Postprocessing")
        
        MainWindow.setCentralWidget(self.centralwidget)
        #detect changes between tabs
        self.tabWidget.currentChanged.connect(self.on_change)
        #add menubar
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        #do file menu first
        fileMenu = self.menubar.addMenu('&File')
        
        loadButton = QtWidgets.QAction('Load', MainWindow)
        loadButton.setShortcut('Ctrl+L')
        loadButton.setStatusTip('Load pyCM data file')
        loadButton.triggered.connect(lambda: self.populate(None))
        
        copyButton = QtWidgets.QAction('Copy', MainWindow)
        copyButton.setShortcut('Ctrl+C')
        copyButton.setStatusTip('Copy current results file.')
        copyButton.triggered.connect(self.copy)
        
        # clearallButton = QtWidgets.QAction('Restart', MainWindow)
        # clearallButton.setShortcut('Ctrl+Shift+R')
        # clearallButton.setStatusTip('Restart interactor without affecting current analysis.')
        # clearallButton.triggered.connect(self.restart)
        
        exitButton = QtWidgets.QAction('Exit', MainWindow)
        exitButton.setShortcut('Ctrl+Q')
        exitButton.setStatusTip('Exit application')
        exitButton.triggered.connect(MainWindow.close)
        
        #add buttons to menus
        fileMenu.addAction(loadButton)
        fileMenu.addAction(copyButton)
        # fileMenu.addAction(clearallButton) #debug
        fileMenu.addAction(exitButton)

        
        #options menu
        optMenu = self.menubar.addMenu('&Options')
        setFEA = QtWidgets.QAction('Set FEA configuration', MainWindow)
        setFEA.setShortcut('Ctrl+E')
        setFEA.setStatusTip('Identify/set FEA executable locations')
        setFEA.triggered.connect(self.getFEAconfig)
        
        setWorkDir = QtWidgets.QAction('Set FEA working directory', MainWindow)
        setWorkDir.setShortcut('Ctrl+D')
        setWorkDir.setStatusTip('Select working directory to perform FEA')
        setWorkDir.triggered.connect(self.setFEAworkdir)
        
        
        optMenu.addAction(setWorkDir)
        optMenu.addAction(setFEA)
        
        ##To be implemented
        #reporting menu
        # optMenu = self.menubar.addMenu('&Reporting')

        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        #set up tabs
        self.initialize_all()
        
        
    def setup_pc(self):
        lhLayout=QtWidgets.QHBoxLayout(self.pctab)
        self.pctab.setLayout(lhLayout)
        self.pcui=pc.pnt_interactor(self.centralwidget)
        lhLayout.addWidget(self.pcui)
        self.pcui.iren.Initialize()
        self.pcui.unsaved_changes=False

    def setup_aa(self):
        lhLayout=QtWidgets.QHBoxLayout(self.aatab)
        self.aatab.setLayout(lhLayout)
        self.aaui=aa.aa_interactor(self.centralwidget)
        lhLayout.addWidget(self.aaui)
        self.aaui.iren.Initialize()
        self.aaui.averaged=False
        self.aaui.unsaved_changes=False
        
    def setup_sf(self):
        lhLayout=QtWidgets.QHBoxLayout(self.sftab)
        self.sftab.setLayout(lhLayout)
        self.sfui=sf.surf_int(self.centralwidget)
        lhLayout.addWidget(self.sfui)
        self.sfui.iren.Initialize()
        self.sfui.fitted=False
        self.sfui.unsaved_changes=False

    def setup_pre(self):
        lhLayout=QtWidgets.QHBoxLayout(self.pretab)
        self.pretab.setLayout(lhLayout)
        self.preui=pre.msh_interactor(self.centralwidget)
        lhLayout.addWidget(self.preui)
        self.preui.iren.Initialize()
        self.preui.preprocessed = False
        self.preui.unsaved_changes=False
        
    def setup_post(self):
        lhLayout=QtWidgets.QHBoxLayout(self.posttab)
        self.posttab.setLayout(lhLayout)
        self.postui=post.pp_interactor(self.centralwidget)
        lhLayout.addWidget(self.postui)
        self.postui.iren.Initialize()
        self.postui.unsaved_changes=False

    def populate(self, file):
        """
        Populates all ui's with contents of mat file
        """
        if file == None:
            self.activeFile, _, = get_file('*.mat')
        else: self.activeFile=file
        
        if self.activeFile != None:
            MainWindow.setWindowTitle("%s  -  pyCM v%s" %(self.activeFile,__version__))
            self.pcui.fileo=self.activeFile
            self.pcui.load_mat()
            self.aaui.get_input_data(self.activeFile)
            self.sfui.get_input_data(self.activeFile)
            self.preui.get_input_data(self.activeFile)
            self.postui.get_input_data(self.activeFile)
        else: return
        
    def initialize_all(self):
        #run set up on all tabs again with null argument
        self.setup_pc()
        self.setup_aa()
        self.setup_sf()
        self.setup_pre()
        self.setup_post()
    
    def on_change(self):
        
        #check if there's an activeFile & change title bar for new analyses
        
        #perform checks from each step and load as necessary
        
        QtWidgets.QApplication.processEvents()
        
        #check the status of alignment/averaging against the surface input file
        if not self.pcui.unsaved_changes and self.tabWidget.currentIndex()==1:
            try:
                
                if not hasattr(self,'activeFile'): #otherwise will reload data un-necessarily
                    self.activeFile = self.pcui.fileo
                    self.aaui.get_input_data(self.activeFile)
                    MainWindow.setWindowTitle("%s  -  pyCM v%s" %(self.activeFile,__version__))

            except:
                ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Either no point cloud data, or loading alignment and averaging data failed.", \
                QtWidgets.QMessageBox.Ok)
                self.tabWidget.setCurrentIndex(0)
        elif not (self.pcui.refWritten and self.pcui.floatWritten) and self.tabWidget.currentIndex()>0:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Insufficient point cloud data available for this step.", \
                QtWidgets.QMessageBox.Ok)
            self.tabWidget.setCurrentIndex(0)
            
            
        #check if there are unsaved changes from the point editor
        if self.pcui.unsaved_changes and self.tabWidget.currentIndex()!=0:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "There are unsaved changes pending to a point cloud. Ignore?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if ret == QtWidgets.QMessageBox.Yes: #don't incorporate
                #reload relevant parts of results file & clear unsaved changes.
                self.pcui.load_mat()
                return
            else: 
                #change unsaved flag temporarily to suppress on_change dialog
                self.tabWidget.setCurrentIndex(0)


        #check if there are unsaved changes from alignment/averaging
        if self.pcui.unsaved_changes and self.tabWidget.currentIndex()!=1:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "There are unsaved changes pending on alignment/averaging. Ignore?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if ret == QtWidgets.QMessageBox.Yes: #don't incorporate
                #reload relevant parts of results file & clear unsaved changes.
                self.aaui.get_input_data(self.activeFile)
                return
            else: 
                self.tabWidget.setCurrentIndex(1)


        #check if there are unsaved changes from surface fitting
        if self.sfui.unsaved_changes and self.tabWidget.currentIndex()!=2:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "There are unsaved changes pending to fitting. Ignore?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if ret == QtWidgets.QMessageBox.Yes: #don't incorporate
                #reload relevant parts of results file & clear unsaved changes.
                self.sfui.get_input_data(self.activeFile)
                return
            else: 
                self.tabWidget.setCurrentIndex(2)
                
        #check if there are unsaved changes from preprocessing
        if self.preui.unsaved_changes and self.tabWidget.currentIndex()!=3:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "There are unsaved changes pending to preprocessing. Ignore?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if ret == QtWidgets.QMessageBox.Yes: #don't incorporate
                #reload relevant parts of results file & clear unsaved changes.
                self.preui.get_input_data(self.activeFile)
                return
            else: 
                self.tabWidget.setCurrentIndex(3)
        
        #check the status of alignment/averaging against the surface input file
        if self.aaui.averaged and self.tabWidget.currentIndex()==2:
            try:
                if not self.sfui.unsaved_changes:
                    self.sfui.get_input_data(self.activeFile)
                    # print('reloaded fitted surface') #debug

            except:
                ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Loading results from alignment/averaging failed.", \
                QtWidgets.QMessageBox.Ok)
        elif not (self.aaui.averaged) and self.tabWidget.currentIndex()==2:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Need to have an averaged surface saved for this step.", \
                QtWidgets.QMessageBox.Ok)
            self.tabWidget.setCurrentIndex(1)
            
        #check the status of surface fitting 
        if self.sfui.fitted and self.tabWidget.currentIndex()==3:
            try:
                if not self.preui.unsaved_changes:
                    self.preui.get_input_data(self.activeFile)
            except:
                ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Loading results from fitting failed.", \
                QtWidgets.QMessageBox.Ok)
        elif not (self.sfui.fitted) and self.tabWidget.currentIndex()==3:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Need to have a fitted surface saved for this step.", \
                QtWidgets.QMessageBox.Ok)
            self.tabWidget.setCurrentIndex(2)
            
        #check the status of the FEA
        if self.preui.preprocessed and self.tabWidget.currentIndex()==4:
            try:
                if not self.preui.unsaved_changes:
                    self.postui.get_input_data(self.activeFile)
            except:
                ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Loading results from fitting failed.", \
                QtWidgets.QMessageBox.Ok)
        elif not (self.preui.preprocessed) and self.tabWidget.currentIndex()==4:
            ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Need to have a completed FEA for this step.", \
                QtWidgets.QMessageBox.Ok)
            self.tabWidget.setCurrentIndex(3)

    def getFEAconfig(self):
        #check config file
        try:
            self.filec = resource_filename("pyCM","pyCMconfig.yml")#needs to be pyCM
        except: #resource_filename will inform if the directory doesn't exist
            print("Did not find config file in the pyCM installation directory.")
        try:
            with open(self.preui.filec,'r') as ymlfile:
                readcfg = yaml.load(ymlfile, Loader=yaml.FullLoader)    
                l=[readcfg['FEA']['abaqusExec'],readcfg['FEA']['gmshExec'],readcfg['FEA']['ccxExec']]
                self.preui.cfg=pre.GetFEAconfig(l,self.filec)
        except:
            try:
                pre.GetFEAconfig(['','',''],self.filec)
            except:
                print("Failed to set config file.")
                return
    
    def setFEAworkdir(self):
        #sets directory to run FEA from - writes outputd variable to preprocessor interactor
        self.preui.outputd = str(QtWidgets.QFileDialog.getExistingDirectory(None, "Select directory", "",QtWidgets.QFileDialog.ShowDirsOnly))

    def copy(self):
    
        if hasattr(self,'activeFile'):
            if self.activeFile != None:
                launchlocation, _ = os.path.split(self.activeFile)
                newFile, _, = get_open_file('*.mat',launchlocation)
                if newFile == self.activeFile:
                    ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                    "Overwriting in this manner isn't supported. Modify current results file.", \
                    QtWidgets.QMessageBox.Ok)
                    return
                #copy it
                try:
                    shutil.copyfile(self.activeFile, newFile)
                except:
                    return
                
                ret=QtWidgets.QMessageBox.warning(MainWindow, "pyCM Warning", \
                "Results file copied. Open this copy?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
                if ret == QtWidgets.QMessageBox.Yes: 
                    #populate with copy
                    self.populate(newFile)
                
    # def restart(self):
        # self.centralwidget.close
        # os.execl(sys.executable, sys.executable, *sys.argv) #debug
        # sp.call(sys.executable + ' "' + os.path.realpath(__file__) + '"')

if __name__ == "__main__":

    app = QtWidgets.QApplication(sys.argv)
    screen_res = app.desktop().screenGeometry()
    
    
    spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
    splash_pix = QtGui.QPixmap(spl_fname,'PNG')
    splash = QtWidgets.QSplashScreen(splash_pix)
    splash.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint | QtCore.Qt.FramelessWindowHint)
    splash.setMask(splash_pix.mask())
    

    splash.show()
    app.processEvents()
    
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow,0.75*screen_res.width(),0.75*screen_res.height())
    
    MainWindow.show()
    splash.finish(MainWindow)
    
    ow = vtk.vtkFileOutputWindow();
    ow.SetFileName("vtk_errors.txt");
    ow.SetInstance(ow);
    #send to console instead
    # ow = vtk.vtkOutputWindow();
    # ow.SendToStdErrOn()
    
    ret = app.exec_()
    
    if sys.stdin.isatty() and not hasattr(sys,'ps1'):
        del MainWindow
        sys.exit(ret)


