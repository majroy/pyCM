#!/usr/bin/env python
'''
VTK and QT to allow for preprocessing FEAs associated with the contour method. 
See documentation for the following main methods in this package:
-point_cloud
-align_average
-fit_surface
-preprocess
-------------------------------------------------------------------------------
1.5 - associated with overall pyCM 2.0 release
'''
__author__ = "M.J. Roy"
__version__ = "1.5"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014--"

import sys,os,ctypes, time
from PyQt5 import QtCore, QtGui, QtWidgets
import vtk
from pkg_resources import Requirement, resource_filename
from importlib.metadata import version
from pyCMcommon import make_splash, get_file
import registration as reg
import align_average as aa
import fit_surface as fs
import preprocess as pre
import postprocess as post

class main_window(QtWidgets.QMainWindow):
    '''
    Inherits most attributes from QMainWindow
    '''
    def __init__(self, parent=None):
        super(main_window, self).__init__(parent)
        self.setWindowTitle("pyCM v%s" %version('pyCM'))
        self.setWindowIcon(QtGui.QIcon(resource_filename("pyCM","meta/pyCM_icon.png")))

        if os.name == 'nt':
            myappid = 'pyCM.main.%s'%version('pyCM') # arbitrary string
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(myappid) #windows taskbar icon
        
        self.setMinimumSize(QtCore.QSize(1000, 1000))
        
        #create menubar
        menubar = QtWidgets.QMenuBar(self)
        self.setMenuBar(menubar)
        #file menu
        file_menu = menubar.addMenu('&File')
        load_button = QtWidgets.QAction('Load', self)
        load_button.setShortcut('Ctrl+L')
        load_button.setStatusTip('Load pyCM data file')
        
        file_menu.addAction(load_button)
        
        #create tabwidget
        self.tabWidget = QtWidgets.QTabWidget()
        
        
        #add tabs
        self.regtab = QtWidgets.QWidget(self.tabWidget)
        self.tabWidget.addTab(self.regtab, "Registration")
        self.aatab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.aatab, "Alignment/averaging")
        self.fstab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.fstab, "Surface fitting")
        self.pretab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.pretab, "Preprocessing")
        self.posttab = QtWidgets.QWidget()
        self.tabWidget.addTab(self.posttab, "Postprocessing")
        self.setCentralWidget(self.tabWidget)

        #add a status bar
        self.statusbar = QtWidgets.QStatusBar(self)
        # self.statusbar.setObjectName("statusbar")
        self.setStatusBar(self.statusbar)

        #connect functions
        load_button.triggered.connect(self.populate)

        self.tabWidget.setCurrentIndex(0)
        self.initialize_all()
        self.file = None #active datafile
        self.active_dir = os.getcwd()
    
    def center(self):
        frame = self.frameGeometry()
        center = QtWidgets.QDesktopWidget().availableGeometry().center()
        frame.moveCenter(center)
        self.move(frame.topLeft())
    
    def closeEvent(self,event):
        '''
        Finalize all VTK widgets to negate OpenGL messages/errors
        '''
        self.regui.ui.vtkWidget.close()
        self.regui.ui.preview_widget.vtkWidget.close()
        self.aaui.ui.vtkWidget.close()
        self.fsui.ui.vtkWidget.close()
        self.preui.ui.vtkWidget.close()
        self.postui.ui.vtkWidget.close()
    
    def setup_reg(self):
        lhLayout = QtWidgets.QHBoxLayout(self.regtab)
        self.regui=reg.interactor(self.tabWidget)
        lhLayout.addWidget(self.regui)
    
    def setup_aa(self):
        lhLayout = QtWidgets.QHBoxLayout(self.aatab)
        self.aaui=aa.interactor(self.tabWidget)
        lhLayout.addWidget(self.aaui)
    
    def setup_fs(self):
        lhLayout = QtWidgets.QHBoxLayout(self.fstab)
        self.fsui = fs.interactor(self.tabWidget)
        lhLayout.addWidget(self.fsui)
    
    def setup_pre(self):
        lhLayout = QtWidgets.QHBoxLayout(self.pretab)
        self.preui = pre.interactor(self.tabWidget)
        lhLayout.addWidget(self.preui)
        
    def setup_post(self):
        lhLayout = QtWidgets.QHBoxLayout(self.posttab)
        self.postui = post.interactor(self.tabWidget)
        lhLayout.addWidget(self.postui)
    
    def initialize_all(self):
        self.setup_reg()
        self.setup_aa()
        self.setup_fs()
        self.setup_pre()
        self.setup_post()

    def populate(self):
        """
        Populates all ui with contents of pyCM file
        """
        
        self.file, self.active_dir = get_file("*.pyCM",self.active_dir)
        
        #make sure valid file was selected
        if self.file is None or not(os.path.isfile(self.file)):
            return
        
        self.setWindowTitle("%s  -  pyCM v%s" %(self.file,version('pyCM')))
        self.regui.output_filename = self.file
        
        self.regui.ui.tab_widget.setCurrentIndex(1)
        self.aaui.file = self.file
        self.aaui.get_data()
        self.fsui.file = self.file
        self.fsui.get_data()
        self.preui.file = self.file
        self.preui.get_data()
        self.postui.file = self.file
        self.postui.get_data()


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    app.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps)

    splash = make_splash()
    splash.show()

    app_main_window = main_window()
    app_main_window.center()
    app_main_window.show()

    splash.finish(app_main_window)
    
    sys.exit(app.exec_())