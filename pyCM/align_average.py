#!/usr/bin/env python
'''
Uses VTK python to allow for editing point clouds associated with the contour 
method. Full interaction requires a 3-button mouse and keyboard.
-------------------------------------------------------------------------------
Current mapping is as follows:
LMB - rotate about point cloud centroid.
MMB - pan
RMB - zoom/refresh window extents
1 - view 1, default, looks down z axis onto xy plane
2 - view 2, looks down x axis onto zy plane
3 - view 3, looks down y axis onto zx plane
z - increase z-aspect ratio
x - decrease z-aspect ratio
c - return to default z-aspect
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
a - remove compass/axes
o - hide/restore outlines
l - load a results file
-------------------------------------------------------------------------------
ver 19-01-08
1.1 - Initial release
1.2 - Refactored to use PyQt interface and eliminated global variables
1.3 - Refactored to use PyQt5, Python 3
1.4 - added option to remove start at centroid
1.5 - added additional tools for alignment
'''
__author__ = "M.J. Roy"
__version__ = "1.5"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"
__copyright__ = "(c) M. J. Roy, 2014-2019"

import sys
import os.path
import vtk
import vtk.util.numpy_support as VN
import numpy as np
import numpy.matlib
import scipy.io as sio
from scipy.interpolate import griddata
from scipy.spatial.distance import pdist, squareform
from matplotlib import path
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
from pyCM.pyCMcommon import *
from pkg_resources import Requirement, resource_filename


def aa_def(*args,**kwargs):
    """
    Main function, builds qt interaction
    """    
    app = QtWidgets.QApplication.instance()
    if app is None:
        app = QtWidgets.QApplication(sys.argv)
    
    spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
    splash_pix = QtGui.QPixmap(spl_fname,'PNG')
    splash = QtWidgets.QSplashScreen(splash_pix)
    splash.setMask(splash_pix.mask())

    splash.show()
    app.processEvents()
    
    window = aa_interactor(None)

    if len(args)==1: 
        aa_interactor.get_input_data(window,args[0])
    else: 
        aa_interactor.get_input_data(window,None)

    window.show()
    splash.finish(window)
    window.iren.Initialize() # Need this line to actually show the render inside Qt

    ret = app.exec_()
    
    if sys.stdin.isatty() and not hasattr(sys,'ps1'):
        sys.exit(ret)
    else:
        return window

class ali_avg(object):
    """
    Class to build qt interaction, including VTK widget
    setupUi builds, initialize starts VTK widget
    """
    
    def setupUi(self, MainWindow):
        MainWindow.setWindowTitle("pyCM - Alignment and averaging tool v%s" %__version__)
        MainWindow.setWindowIcon(QtGui.QIcon(resource_filename("pyCM","meta/pyCM_icon.png")))
        self.centralWidget = QtWidgets.QWidget(MainWindow)
        if hasattr(MainWindow,'setCentralWidget'):
            MainWindow.setCentralWidget(self.centralWidget)
        else:
            self.centralWidget=MainWindow
        self.mainlayout=QtWidgets.QGridLayout(self.centralWidget)
        
        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        
        mainUiBox = QtWidgets.QGridLayout()
        
        self.vtkWidget.setMinimumSize(QtCore.QSize(1050, 600))
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(10)
        sizePolicy.setVerticalStretch(10)
        sizePolicy.setHeightForWidth(self.vtkWidget.sizePolicy().hasHeightForWidth())
        self.vtkWidget.setSizePolicy(sizePolicy)
        
        self.statLabel=QtWidgets.QLabel("Idle")
        self.statLabel.setWordWrap(True)
        self.statLabel.setFont(QtGui.QFont("Helvetica",italic=True))
        self.statLabel.setMinimumWidth(100)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.statLabel.sizePolicy().hasHeightForWidth())
        self.statLabel.setSizePolicy(sizePolicy)

        headFont=QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold)

        
        # #define buttons/widgets
        # self.reloadButton = QtWidgets.QPushButton('Load')
        scalingLabel=QtWidgets.QLabel("Active axis for scaling")
        scalingLabel.setFont(headFont)
        self.xsButton=QtWidgets.QRadioButton("x")
        self.ysButton=QtWidgets.QRadioButton("y")
        self.zsButton=QtWidgets.QRadioButton("z")
        self.zsButton.setChecked(True)
        self.scalingButtonGroup = QtWidgets.QButtonGroup()
        self.scalingButtonGroup.addButton(self.xsButton)
        self.scalingButtonGroup.addButton(self.ysButton)
        self.scalingButtonGroup.addButton(self.zsButton)
        self.scalingButtonGroup.setExclusive(True)
        scaleBoxlayout = QtWidgets.QGridLayout()
        scaleBoxlayout.addWidget(self.xsButton,1,1)
        scaleBoxlayout.addWidget(self.ysButton,1,2)
        scaleBoxlayout.addWidget(self.zsButton,1,3)
        
        horizLine1=QtWidgets.QFrame()
        horizLine1.setFrameStyle(QtWidgets.QFrame.HLine)
        mirrorLabel=QtWidgets.QLabel("Mirroring")
        mirrorLabel.setFont(headFont)
        self.mirrorXbutton = QtWidgets.QPushButton('ZY')
        self.mirrorYbutton = QtWidgets.QPushButton('ZX')

        horizLine2=QtWidgets.QFrame()
        horizLine2.setFrameStyle(QtWidgets.QFrame.HLine)
        alignLabel=QtWidgets.QLabel("Alignment")
        alignLabel.setFont(headFont)
        self.centRefButton=QtWidgets.QPushButton("Move reference centroid to origin")
        self.centFloatButton=QtWidgets.QPushButton("Move float centroid to origin")

        self.transXlabel=QtWidgets.QLabel("Translate x:")
        self.transX = QtWidgets.QDoubleSpinBox()
        self.transX.setValue(0)
        self.transX.setMaximum(300)
        self.transX.setMinimum(-300)
        self.transYlabel=QtWidgets.QLabel("Translate y:")
        self.transY = QtWidgets.QDoubleSpinBox()
        self.transY.setValue(0)
        self.transY.setMaximum(300)
        self.transY.setMinimum(-300)
        self.rotateZlabel=QtWidgets.QLabel("Rotate about z (deg):")
        self.rotateZ= QtWidgets.QDoubleSpinBox()
        self.rotateZ.setValue(0)
        self.rotateZ.setMaximum(180)
        self.rotateZ.setMinimum(-180)
        self.transButton=QtWidgets.QPushButton('Transform floating')
        self.numPntsOutline = QtWidgets.QSpinBox()
        self.numPntsOutline.setMaximum(10000)
        self.reduceOutlineButton = QtWidgets.QPushButton('Decimate outlines')

        
        alignAlgoButtonGroup = QtWidgets.QButtonGroup()
        self.useVTKalignButton=QtWidgets.QRadioButton("VTK ICP")
        self.useICPalignButton=QtWidgets.QRadioButton("K-neighbour ICP")
        self.useVTKalignButton.setChecked(True)
        alignAlgoButtonGroup.addButton(self.useVTKalignButton)
        alignAlgoButtonGroup.addButton(self.useICPalignButton)
        alignAlgoButtonGroup.setExclusive(True)
        self.X180Button = QtWidgets.QPushButton("Flip X")
        self.Y180Button = QtWidgets.QPushButton("Flip Y")
        self.alignButton = QtWidgets.QPushButton("Align")
        self.acceptAlignButton = QtWidgets.QPushButton("Accept")
        
        self.alignButton.setStyleSheet("background-color : None ")
        
        horizLine3=QtWidgets.QFrame()
        horizLine3.setFrameStyle(QtWidgets.QFrame.HLine)
        averageLabel=QtWidgets.QLabel("Averaging")
        averageLabel.setFont(headFont)
        
        #widgets for setting grid
        gridLabel=QtWidgets.QLabel("Grid spacing:")
        self.gridInd = QtWidgets.QDoubleSpinBox()
        self.gridInd.setValue(0)
        self.gridInd.setMaximum(5)
        self.gridInd.setMinimum(0.001)
        
        self.averageButton = QtWidgets.QPushButton('Average')
        self.averageButton.setStyleSheet("background-color : None ")
        
        horizLine4=QtWidgets.QFrame()
        horizLine4.setFrameStyle(QtWidgets.QFrame.HLine)
        self.writeButton=QtWidgets.QPushButton('Write')
        
        horizLine5=QtWidgets.QFrame()
        horizLine5.setFrameStyle(QtWidgets.QFrame.HLine)


        #add widgets to ui
        mainUiBox.addWidget(scalingLabel,0,0,1,2)
        mainUiBox.addLayout(scaleBoxlayout,1,0,1,2)
        mainUiBox.addWidget(horizLine1,2,0,1,2)
        mainUiBox.addWidget(mirrorLabel,3,0,1,2)
        mainUiBox.addWidget(self.mirrorYbutton,4,0,1,1)
        mainUiBox.addWidget(self.mirrorXbutton,4,1,1,1)
        mainUiBox.addWidget(horizLine2,5,0,1,2)
        mainUiBox.addWidget(alignLabel,6,0,1,2)
        mainUiBox.addWidget(self.centRefButton,7,0,1,2)
        mainUiBox.addWidget(self.centFloatButton,8,0,1,2)
        mainUiBox.addWidget(self.transXlabel,9,0,1,1)
        mainUiBox.addWidget(self.transX,9,1,1,1)
        mainUiBox.addWidget(self.transYlabel,10,0,1,1)
        mainUiBox.addWidget(self.transY,10,1,1,1)
        mainUiBox.addWidget(self.rotateZlabel,11,0,1,1)
        mainUiBox.addWidget(self.rotateZ,11,1,1,1)
        mainUiBox.addWidget(self.transButton,12,0,1,2)
        mainUiBox.addWidget(self.X180Button,13,0,1,1)
        mainUiBox.addWidget(self.Y180Button,13,1,1,1)
        mainUiBox.addWidget(self.numPntsOutline,14,0,1,1)
        mainUiBox.addWidget(self.reduceOutlineButton,14,1,1,1)
        
        mainUiBox.addWidget(self.useVTKalignButton,15,0,1,1)
        mainUiBox.addWidget(self.useICPalignButton,15,1,1,1)

        mainUiBox.addWidget(self.alignButton,16,0,1,1)
        mainUiBox.addWidget(self.acceptAlignButton,16,1,1,1)
        mainUiBox.addWidget(horizLine3,17,0,1,2)
        mainUiBox.addWidget(averageLabel,18,0,1,2)
        mainUiBox.addWidget(gridLabel,19,0,1,1)
        mainUiBox.addWidget(self.gridInd,19,1,1,1)
        mainUiBox.addWidget(self.averageButton,20,0,1,2)
        mainUiBox.addWidget(horizLine4,21,0,1,2)
        mainUiBox.addWidget(self.writeButton,22,0,1,2)
        mainUiBox.addWidget(horizLine5,23,0,1,2)
        # mainUiBox.addWidget(self.statusLabel,18,0,1,2)

        mainUiBox.setColumnMinimumWidth(0,mainUiBox.columnMinimumWidth(0))
        mainUiBox.setColumnMinimumWidth(1,mainUiBox.columnMinimumWidth(0))
        
        lvLayout=QtWidgets.QVBoxLayout()
        lhLayout=QtWidgets.QHBoxLayout()
        lvLayout.addLayout(mainUiBox)
        lvLayout.addStretch(1)
        lhLayout.addLayout(lvLayout)
        lhLayout.addStretch(2)



        self.mainlayout.addWidget(self.vtkWidget,0,0,1,1)
        self.mainlayout.addLayout(lhLayout,0,1,1,1)
        self.mainlayout.addWidget(self.statLabel,1,0,1,2)
        
    def initialize(self):
        self.vtkWidget.start()

class aa_interactor(QtWidgets.QWidget):
    """
    Sets up the main VTK window, reads file and sets connections between UI and interactor
    """
    def __init__(self, parent):
        super(aa_interactor,self).__init__(parent)
        self.ui = ali_avg()
        self.ui.setupUi(self)
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(0.1, 0.2, 0.4)
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        style=vtk.vtkInteractorStyleTrackballCamera()
        style.AutoAdjustCameraClippingRangeOn()
        self.iren.SetInteractorStyle(style)
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.cp=self.ren.GetActiveCamera().GetPosition()
        self.fp=self.ren.GetActiveCamera().GetFocalPoint()
        self.iren.AddObserver("KeyPressEvent", self.keypress)
        
        self.PointSize=2
        self.LineWidth=1
        self.Zaspect=1.0
        self.limits=np.empty(6)
        self.picking=False
        self.Offset=0
        self.mirrored=False
        self.aligned=False
        
        # self.ui.reloadButton.clicked.connect(lambda: self.get_input_data(None))
        self.ui.mirrorXbutton.clicked.connect(lambda: self.flipside('x'))
        self.ui.mirrorYbutton.clicked.connect(lambda: self.flipside('y'))
        self.ui.centRefButton.clicked.connect(lambda: self.zero_pos('ref'))
        self.ui.centFloatButton.clicked.connect(lambda: self.zero_pos('float'))
        self.ui.transButton.clicked.connect(lambda: self.shift())
        self.ui.X180Button.clicked.connect(lambda: self.flip('x'))
        self.ui.Y180Button.clicked.connect(lambda: self.flip('y'))
        self.ui.alignButton.clicked.connect(lambda: self.align())
        self.ui.acceptAlignButton.clicked.connect(lambda: self.accept_align())
        self.ui.averageButton.clicked.connect(lambda: self.average())
        self.ui.writeButton.clicked.connect(lambda: self.write())
        self.ui.reduceOutlineButton.clicked.connect(lambda: self.reduce_outline())
    
    def update_float(self):
    
        if hasattr(self,'fActor'):
            self.ren.RemoveActor(self.fActor)
            self.ren.RemoveActor(self.fOutlineActor)
        
        color=(255, 205, 52)
        self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
        self.ren.AddActor(self.fActor)
        
        self.fOutlineActor, self.fOPC = gen_outline(self.fO_local,color,self.PointSize)
        self.ren.AddActor(self.fOutlineActor)
    
    def update_limits(self):
        # recalculate extents of interactor
        RefMin=np.amin(np.vstack((self.flp,self.rp)),axis=0)
        RefMax=np.amax(np.vstack((self.flp,self.rp)),axis=0)

        extents=RefMax-RefMin #extents
        rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
        self.limits=[RefMin[0]-rl, \
        RefMax[0]+rl, \
        RefMin[1]-rl, \
        RefMax[1]+rl, \
        RefMin[2],RefMax[2]]

        #add axes
        self.add_axis(self.limits,[1,1,1])
        
        s,nl,axs=self.get_scale()

        self.fActor.SetScale(s)
        self.fActor.Modified()

        self.add_axis(nl,axs)

        #update
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()    
    
    def zero_pos(self,p):
        '''
        Moves the outline and point cloud from the centroidof the outline of p to 0,0,0 
        '''
        self.unsaved_changes=True
        
        if self.averaged == True: #then set it false and change the button
                self.averaged = False
                self.ui.averageButton.setStyleSheet("background-color : None ")
                
        if self.aligned == True: #then set it false and change the button
            self.aligned = False
            self.ui.alignButton.setStyleSheet("background-color : None ")
            
        local_trans=np.identity(4)
        if p == "ref":
            self.ren.RemoveActor(self.rActor)
            #perform move on datasets
            centroid = np.mean(self.rO, axis = 0)
            #update the homogeneous transformation matrix
            local_trans[0:3,3]=-centroid
            self.refTrans.append(local_trans)
            self.rp = self.rp - centroid
            self.rO = self.rO - centroid
            self.rO_local = self.rO_local - centroid

            color=(242, 101, 34)
            self.ren.RemoveActor(self.rOutlineActor)
            self.rPC, self.rActor, _, = gen_point_cloud(self.rp,color,self.PointSize)
            self.ren.AddActor(self.rActor)
            self.rOutlineActor, self.rOPC = gen_outline(self.rO_local,color,self.PointSize)
            self.ren.AddActor(self.rOutlineActor)
            
        if p == "float":
            
            # perform move on datasets
            centroid = np.mean(self.fO, axis = 0)
            #update homogeneous transformation matrix
            local_trans[0:3,3]=-centroid
            self.floatTrans.append(local_trans)
            self.flp = self.flp - centroid
            self.fO = self.fO - centroid
            self.fO_local = self.fO_local - centroid
            self.update_float()

        self.update_limits()
        self.ren.ResetCamera()
    
    def reduce_outline(self):
        '''
        Decimates alias outlines used for alignment using reduce_equally.
        '''

        self.ren.RemoveActor(self.rOutlineActor)
        
        if self.ui.numPntsOutline.value() > len(self.rO_local):
            self.rO_local=self.rO
            self.fO_local=self.fO
        
        #Do reference first
        color=(242, 101, 34)

        X = respace_equally(self.rO_local,self.ui.numPntsOutline.value())[0]
        self.rO_local=np.zeros((self.ui.numPntsOutline.value(),3))
        self.rO_local[:,:-1]=X
        self.rOutlineActor, self.rOPC = gen_outline(self.rO_local,color,self.PointSize)
        self.ren.AddActor(self.rOutlineActor)
        
        color=(255, 205, 52)

        X = respace_equally(self.fO_local,self.ui.numPntsOutline.value())[0]
        self.fO_local=np.zeros((self.ui.numPntsOutline.value(),3))
        self.fO_local[:,:-1]=X
        
        self.update_float()
        self.update_limits()
        
    def shift(self):
        '''
        Applies rigid body transformations to the floating dataset
        '''
        self.unsaved_changes=True
        
        if self.averaged == True: #then set it false and change the button
            self.averaged = False
            self.ui.averageButton.setStyleSheet("background-color : None ")
            
        if self.aligned == True: #then set it false and change the button
            self.aligned = False
            self.ui.alignButton.setStyleSheet("background-color : None ")

        #get x and y translations and z rotation and update float transformation matrix
        local_trans = np.identity(4)
        a=np.deg2rad(float(self.ui.rotateZ.value()))
        local_trans[0:2,0:2]=np.array([[np.cos(a),-np.sin(a)],[np.sin(a),np.cos(a)]])
        
        local_trans[0,-1]=float(self.ui.transX.value())
        local_trans[1,-1]=float(self.ui.transY.value())
        self.floatTrans.append(local_trans)
        
        #apply operation
        self.flp=apply_trans(self.flp,local_trans)
        self.fO=apply_trans(self.fO,local_trans)
        self.fO_local=apply_trans(self.fO_local,local_trans)

        self.update_float()
        self.update_limits()
    
    def flip(self,axis):
        '''
        Applies a rotation of 180 degrees about 'axis'
        '''
        self.unsaved_changes=True
        
        if self.averaged == True: #then set it false and change the button
            self.averaged = False
            self.ui.averageButton.setStyleSheet("background-color : None ")
        
        if self.aligned == True: #then set it false and change the button
            self.aligned = False
            self.ui.alignButton.setStyleSheet("background-color : None ")

        
        local_trans = np.identity(4)
        if axis == 'x':
            local_trans[1,1] = -1
            local_trans[2,2] = -1
        if axis == 'y':
            local_trans[0,0] = -1
            local_trans[2,2] = -1
        #update overall homogeneous transformation matrix
        self.floatTrans.append(local_trans)
        
        #apply operation
        self.flp=np.dot(self.flp,local_trans[0:3,0:3])
        self.fO=np.dot(self.fO,local_trans[0:3,0:3])
        self.fO_local=np.dot(self.fO_local,local_trans[0:3,0:3])
        
        
        self.update_float()
        self.update_limits()
        
        
    def write(self):
        
        mat_vars=sio.whosmat(self.fileo)
        if not set(['aa', 'trans']).isdisjoint([item for sublist in mat_vars for item in sublist]): #tell the user that they might overwrite their data
            ret=QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
                "There is already data associated with this analysis step saved. Overwrite and invalidate subsequent steps?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if ret == QtWidgets.QMessageBox.No: #don't overwrite
                return
            else:
                #delete fitting parameters with pyCMcommon helper function, which negates FEA pre-processing as well.
                clear_mat(self.fileo,['x_out','aa_mask','spline_x']) 

        mat_contents=sio.loadmat(self.fileo)
        
        new={'trans': {'ref':self.refTrans, 'float':self.floatTrans},'aa': {'pnts': self.ap, 'gsize': self.gsize}}
        
        mat_contents.update(new) #update the dictionary
            
        sio.savemat(self.fileo,mat_contents)    
        self.ui.statLabel.setText("Wrote data.")
        self.unsaved_changes=False
        
    def average(self):
        
        self.unsaved_changes=True
        
        if hasattr(self,'aActor'):
            self.ren.RemoveActor(self.aActor)
        
        self.ui.statLabel.setText("Averaging, applying grid . . .")
        QtWidgets.QApplication.processEvents()
        
        #temporarily shift all data such that it appears in the first cartesian quadrant
        tT=np.amin(self.rO,axis=0)
        self.rO, self.fO, self.rp, self.flp=self.rO-tT, self.fO-tT, self.rp-tT, self.flp-tT
        
        #use max to get a 'window' for assessing grid spacing
        RefMax=np.amax(self.rO,axis=0)
        RefMin=np.amin(self.rO,axis=0)
        windowVerts=np.matrix([[0.25*RefMin[0], 0.25*RefMin[1]],
        [0.25*RefMin[0], 0.25*(RefMax[1])],
        [0.25*(RefMax[1]), 0.25*(RefMax[1])],
        [0.25*(RefMax[0]), 0.25*(RefMin[1])]]);
        
        p=path.Path(windowVerts)
        inWindow=p.contains_points(self.rp[:,:2]) #first 2 columns of RefPoints is x and y
        
        windowed=self.rp[inWindow,:2]
        
        #populate grid size if attribute doesn't exist
        if not hasattr(self,'gsize'):
            gs=squareform(pdist(windowed,'euclidean')) 
            self.gsize = np.mean(np.sort(gs)[:,1])
            self.ui.gridInd.setValue(self.gsize)
        else:
            self.gsize=self.ui.gridInd.value()
        
        #grid the reference based on gsize, bumping out the grid by 10% in either direction
        grid_x, grid_y = np.meshgrid(
        np.linspace(1.1*RefMin[0],1.1*RefMax[0],int((1.1*RefMax[0]-1.1*RefMin[0])/self.gsize)),
        np.linspace(1.1*RefMin[1],1.1*RefMax[1],int((1.1*RefMax[1]-1.1*RefMin[1])/self.gsize)), 
        indexing='xy')
        
        #apply the grid to the reference data
        grid_Ref=griddata(self.rp[:,:2],self.rp[:,-1],(grid_x,grid_y),method='linear')
        
        #apply the grid to the aligned data
        grid_Align=griddata(self.flp[:,:2],self.flp[:,-1],(grid_x,grid_y),method='linear')
        
        self.ui.statLabel.setText("Averaging using grid . . .")
        QtWidgets.QApplication.processEvents()
        
        #average z values
        grid_Avg=(grid_Ref+grid_Align)/2
        
        #make sure that there isn't anything averaged outside the floating outline
        p=path.Path(self.rO[:,:2])
        inTest=np.hstack((np.ravel(grid_x.T)[np.newaxis].T,np.ravel(grid_y.T)[np.newaxis].T))
        inOutline=p.contains_points(inTest)
        
        #averaged points
        self.ap = np.hstack((inTest[inOutline,:], \
                    np.ravel(grid_Avg.T)[np.newaxis].T[inOutline]))
                    
        #move everything back to original location
        self.rO, self.fO, self.rp, self.flp, self.ap = \
        self.rO+tT, self.fO+tT, self.rp+tT, self.flp+tT, self.ap+tT
        
        self.ui.statLabel.setText("Rendering . . .")
        QtWidgets.QApplication.processEvents()
        
        #show it
        color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
        _, self.aActor, _, = gen_point_cloud(self.ap,color,self.PointSize)
        self.ren.AddActor(self.aActor)
        
        s,nl,axs=self.get_scale()

        self.aActor.SetScale(s)
        self.aActor.Modified()
        
        #update
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()
        self.ui.statLabel.setText("Averaging complete.")
        self.averaged=True
        self.ui.averageButton.setStyleSheet("background-color :rgb(77, 209, 97);")
        
    
    def flipside(self,flipDirection):
        self.ui.statLabel.setText("Starting mirroring . . .")
        self.ui.vtkWidget.update()
        local_trans=np.identity(4)

        self.unsaved_changes=True
        
        if self.averaged == True: #then set it false and change the button
            self.averaged = False
            self.ui.averageButton.setStyleSheet("background-color : None ")
        
        if flipDirection == "x":
            local_trans[0,0]=-1
            self.flp[:,0]=-self.flp[:,0]
            # self.fO[:,0]=-self.fO[:,0]
            self.fO_local[:,0]=-self.fO_local[:,0]
            
            if not self.mirrored:
                self.ui.mirrorXbutton.setStyleSheet("background-color :rgb(77, 209, 97);")
                self.mirrored=True
                self.mirror_plane="x"
            elif self.mirror_plane == "x":
                self.ui.mirrorXbutton.setStyleSheet("background-color: None")
                self.mirrored=False
            elif self.mirror_plane == "y":
                self.flp[:,1]=-self.flp[:,1]
                # self.fO[:,1]=-self.fO[:,1]
                self.fO_local[:,1]=-self.fO_local[:,1]
                self.mirrored=True
                self.mirror_plane="x"
                self.ui.mirrorYbutton.setStyleSheet("background-color: None")
                self.ui.mirrorXbutton.setStyleSheet("background-color :rgb(77, 209, 97);")
                
        elif flipDirection == "y":
            local_trans[1,1]=-1
            self.flp[:,1]=-self.flp[:,1]
            # self.fO[:,1]=-self.fO[:,1]
            self.fO_local[:,1]=-self.fO_local[:,1]
            if not self.mirrored:
                self.ui.mirrorYbutton.setStyleSheet("background-color :rgb(77, 209, 97);")
                self.mirrored=True
                self.mirror_plane="y"
            elif self.mirror_plane == "y":
                self.ui.mirrorYbutton.setStyleSheet("background-color: None")
                self.mirrored=False
            elif self.mirror_plane == "x":
                self.flp[:,0]=-self.flp[:,0]
                # self.fO[:,0]=-self.fO[:,0]
                self.fO_local[:,0]=-self.fO_local[:,0]
                self.mirrored=True
                self.mirror_plane="y"
                self.ui.mirrorXbutton.setStyleSheet("background-color: None")
                self.ui.mirrorYbutton.setStyleSheet("background-color :rgb(77, 209, 97);")

        self.floatTrans.append(local_trans)
        self.update_float()
        self.update_limits()
        self.ren.ResetCamera()
        self.ui.statLabel.setText("Mirror operation complete.")
    
    def accept_align(self):
        '''
        Accepts the current alignment and allows analysis to proceed if the profile has not been algorithmically aligned with the align button being pressed.
        '''
        self.aligned = True
        self.ui.alignButton.setStyleSheet("background-color :rgb(77, 209, 97);")
    
    def align(self):
        '''
        Uses the built-in icp landmark transformation provided by vtk to outline actors, and then updates the renderer and the stored homogeneous transformation matrix transM
        '''
        self.unsaved_changes=True
        
        if self.averaged == True: #then set it false and change the button
            self.averaged = False
            self.ui.averageButton.setStyleSheet("background-color : None ")
            
        if self.aligned == True: #then set it false and change the button
            self.aligned = False
            self.ui.alignButton.setStyleSheet("background-color : None ")

            
        self.ui.statLabel.setText("Starting alignment . . .")
        QtWidgets.QApplication.processEvents()
        
        
        
        if self.ui.useVTKalignButton.isChecked():
            icp_trans=vtk.vtkIterativeClosestPointTransform()
            icp_trans.SetSource(self.fOPC)
            icp_trans.SetTarget(self.rOPC)


            icp_trans.StartByMatchingCentroidsOn()
            icp_trans.GetLandmarkTransform().SetModeToRigidBody()
            icp_trans.SetMeanDistanceModeToRMS()
            icp_trans.CheckMeanDistanceOn()
            icp_trans.SetMeanDistanceModeToAbsoluteValue()
            # icp_trans.SetMaximumNumberOfLandmarks(200)
            icp_trans.DebugOn()
            icp_trans.Modified()
            icp_trans.Update()
            icp_trans.Inverse()
            
            T=np.ones(shape=(4,4))
            for i in range(4):
                for j in range(4):
                    T[i,j]=icp_trans.GetMatrix().GetElement(i, j)
            T=np.linalg.inv(T)

        if self.ui.useICPalignButton.isChecked():
            self.reduce_outline()
            T,_,_ = icp(self.fO_local,self.rO_local)
            
            
            
            
        #apply operation
        self.flp=apply_trans(self.flp,T)
        self.fO_local=apply_trans(self.fO_local,T)
        self.floatTrans.append(T)

        
        self.update_float()
        self.update_limits()
        
        
        if self.mirrored==False:
            self.ui.statLabel.setText("WARNING alignment proceeded without a mirror operation. Alignment complete.")
        else:
            self.ui.statLabel.setText("Alignment complete.")
        
        
    def get_input_data(self,filem):
        """
        Loads the content of a *.mat file pertaining to this particular step
        """
        
        if hasattr(self,'rActor'): #then remove everything
            self.ren.RemoveActor(self.rActor)
            self.ren.RemoveActor(self.fActor)
            self.ren.RemoveActor(self.rOutlineActor)
            self.ren.RemoveActor(self.fOutlineActor)
            
        if hasattr(self,'aActor'):
            self.ren.RemoveActor(self.aActor)
        
        if filem == None:
            filem, _, =get_file('*.mat')
        
        if filem: #check variables
            mat_contents = sio.loadmat(filem)
            self.fileo=filem
            if 'aa' in mat_contents:

                
                #draw floating and reference datasets
                
                self.rp=mat_contents['ref']['rawPnts'][0][0]
                ind=mat_contents['ref']['mask'][0][0][0]
                self.rO=mat_contents['ref']['x_out'][0][0]
                self.rO_local=self.rO
                
                self.refTrans=mat_contents['trans']['ref'][0][0]
                
                
                self.rp=self.rp[np.where(ind)]
                
                #apply the transform with post multiplication
                
                for transformation in self.refTrans:
                    self.rp=apply_trans(self.rp,transformation)
                    self.rO=apply_trans(self.rO,transformation)
                    
                self.rO_local=self.rO

                
                color=(242, 101, 34)
                self.rPC, self.rActor, _, = gen_point_cloud(self.rp,color,self.PointSize)
                self.ren.AddActor(self.rActor)
                self.rOutlineActor, self.rOPC = gen_outline(self.rO_local,color,self.PointSize)
                self.ren.AddActor(self.rOutlineActor)
                
                s,nl,axs=self.get_scale()
                
                self.rActor.SetScale(s)
                self.rActor.Modified()
                
                #do other one, but with transformed floating points
                self.flp=mat_contents['float']['rawPnts'][0][0]
                ind=mat_contents['float']['mask'][0][0][0]
                self.fO=mat_contents['float']['x_out'][0][0]
                self.fO_local = self.fO
                
                self.flp=self.flp[np.where(ind)]


                self.floatTrans=mat_contents['trans']['float'][0][0]
                #read in as np array
                
                for transformation in self.floatTrans:
                    self.flp=apply_trans(self.flp,transformation)
                    self.fO=apply_trans(self.fO,transformation)
                
                
                self.fO_local=self.fO
                
                self.update_float()
                
                #after applied, convert transformation arrays to lists
                self.refTrans = self.refTrans.tolist()
                self.floatTrans = self.floatTrans.tolist()
                
                #show aligned and averaged data
                self.ap=mat_contents['aa']['pnts'][0][0]
                self.gsize=mat_contents['aa']['gsize'][0][0]

                #do grid
                self.ui.gridInd.setValue(self.gsize)

                
            
                color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
                _, self.aActor, _, = gen_point_cloud(self.ap,color,self.PointSize)
                self.ren.AddActor(self.aActor)

                self.aActor.SetScale(s)
                self.aActor.Modified()
                
                self.update_limits()
                self.ren.ResetCamera()

                self.ui.numPntsOutline.setValue(np.min([len(self.rO),len(self.fO)]))
                
                self.ui.statLabel.setText("This dataset has already been aligned and averaged.")
                self.aligned = True
                self.ui.alignButton.setStyleSheet("background-color :rgb(77, 209, 97);")
                self.mirrored=True
                self.averaged=True
                self.ui.averageButton.setStyleSheet("background-color :rgb(77, 209, 97);")
            else:
                self.ui.statLabel.setText("This dataset has not been previously aligned.")
                self.averaged = False

                try:
                    self.rp=mat_contents['ref']['rawPnts'][0][0]
                    ind=mat_contents['ref']['mask'][0][0][0]
                    self.rO=mat_contents['ref']['x_out'][0][0]
                    self.rO_local=self.rO
                    
                    
                    self.rp=self.rp[np.where(ind)]
                    
                    color=(242, 101, 34)
                    self.rPC, self.rActor, _, = gen_point_cloud(self.rp,color,self.PointSize)
                    self.ren.AddActor(self.rActor)
                    self.rOutlineActor, self.rOPC = gen_outline(self.rO_local,color,self.PointSize)
                    self.ren.AddActor(self.rOutlineActor)
                    
                    #do other one
                    self.flp=mat_contents['float']['rawPnts'][0][0]
                    ind=mat_contents['float']['mask'][0][0][0]
                    self.fO=mat_contents['float']['x_out'][0][0]
                    self.fO_local=self.fO
                    
                    
                    self.flp=self.flp[np.where(ind)]
                    
                    #populate outline
                    self.ui.numPntsOutline.setValue(np.min([len(self.rO),len(self.fO)]))
                    

                    
                    color=(255, 205, 52)
                    self.fPC, self.fActor, _, = gen_point_cloud(self.flp,color,self.PointSize)
                    self.ren.AddActor(self.fActor)
                    self.fOutlineActor, self.fOPC = gen_outline(self.fO_local,color,self.PointSize)
                    self.ren.AddActor(self.fOutlineActor)
                    
                    RefMin=np.amin(np.vstack((self.flp,self.rp)),axis=0)
                    RefMax=np.amax(np.vstack((self.flp,self.rp)),axis=0)

                    extents=RefMax-RefMin #extents
                    rl=0.1*(np.amin(extents)) #linear 'scale' to set up interactor
                    self.limits=[RefMin[0]-rl, \
                    RefMax[0]+rl, \
                    RefMin[1]-rl, \
                    RefMax[1]+rl, \
                    RefMin[2],RefMax[2]]

                    #add axes
                    self.add_axis(self.limits,[1,1,1])
                    
                    #initialize both transformation matrices
                    self.refTrans=[]
                    self.floatTrans=[]

                except Exception as e:
                    print("Couldn't read in both sets of data.")
                    print(e)
                
        else:
            print("Invalid *.mat file")
            return
            
        self.unsaved_changes=False
        #update
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()
    
    def keypress(self,obj,event):
        key = obj.GetKeyCode()

        if key =="1":
            xyview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
        elif key =="2":
            yzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
        elif key =="3":
            xzview(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
        elif key=="z":
            self.Zaspect=self.Zaspect*2
            s,nl,axs=self.get_scale()
            if hasattr(self,'pointActor'):
                self.pointActor.SetScale(s)
                self.pointActor.Modified()
            if hasattr(self,'rActor'):
                self.rActor.SetScale(s)
                self.rActor.Modified()
            if hasattr(self,'fActor'):
                self.fActor.SetScale(s)
                self.fActor.Modified()
            if hasattr(self,'aActor'):
                self.aActor.SetScale(s)
                self.aActor.Modified()
            
            self.add_axis(nl,axs)

        elif key=="x":
            self.Zaspect=self.Zaspect*0.5
            s,nl,axs=self.get_scale()
            if hasattr(self,'pointActor'):
                self.pointActor.SetScale(s)
            if hasattr(self,'rActor'):
                self.rActor.SetScale(s)
                self.rActor.Modified()
            if hasattr(self,'fActor'):
                self.fActor.SetScale(s)
                self.fActor.Modified()
            if hasattr(self,'aActor'):
                self.aActor.SetScale(s)
                self.aActor.Modified()

            self.add_axis(nl,axs)


        elif key=="c":
            self.Zaspect=1.0
            s,_,_,=self.get_scale()
            if hasattr(self,'pointActor'):
                self.pointActor.SetScale(s)
            if hasattr(self,'rActor'):
                self.rActor.SetScale(s)
                self.rActor.Modified()
            if hasattr(self,'fActor'):
                self.fActor.SetScale(s)
                self.fActor.Modified()
            if hasattr(self,'aActor'):
                # self.fActor.SetScale(1,1,self.Zaspect)
                self.aActor.SetScale(s)
                self.aActor.Modified()
            self.add_axis(self.limits,[1,1,1])
            self.ren.ResetCamera()

        elif key=="i":
            im = vtk.vtkWindowToImageFilter()
            writer = vtk.vtkPNGWriter()
            im.SetInput(self.ui.vtkWidget._RenderWindow)
            im.Update()
            writer.SetInputConnection(im.GetOutputPort())
            writer.SetFileName("Avg_aligned.png")
            writer.Write()
            print("Screen output saved to %s" %os.path.join(os.getcwd(),'Avg_aligned.png'))

        elif key=="a":
            if hasattr(self,'ax3D'):
                flip_visible(self.ax3D)
            
        elif key == "o":
            if hasattr(self,'outlineActor'):
                flip_visible(self.outlineActor)
        
        elif key == "f":
            if hasattr(self,'ax3D'):
                flip_colors(self.ren,self.ax3D)
                
                
        elif key=="l":
            self.get_input_data(None)
        
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()

    def get_scale(self):
        '''
        Returns array for the keypress function based on what radio button is selected.
        '''
        if self.ui.xsButton.isChecked():
            s=np.array([self.Zaspect,1,1])
            nl=np.append([self.limits[0]*self.Zaspect,self.limits[1]*self.Zaspect],self.limits[2:])
            axs=np.array([1/self.Zaspect,1,1])
            
        elif self.ui.ysButton.isChecked():
            s=np.array([1,self.Zaspect,1])
            nl=np.append(self.limits[0:2],([self.limits[2]*self.Zaspect,self.limits[3]*self.Zaspect],self.limits[4:]))
            axs=np.array([1,1/self.Zaspect,1])
        else:
            s=np.array([1,1,self.Zaspect])
            nl=np.append(self.limits[0:4],([self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect]))
            axs=np.array([1,1,1/self.Zaspect])
        return s,nl,axs

    def add_axis(self,limits,scale):
        if hasattr(self,"ax3D"):
            self.ren.RemoveActor(self.ax3D)
        self.ax3D = vtk.vtkCubeAxesActor()
        self.ax3D.ZAxisTickVisibilityOn()
        self.ax3D.SetXTitle('X')
        self.ax3D.SetXUnits('mm')
        self.ax3D.SetYTitle('Y')
        self.ax3D.SetYUnits('mm')
        self.ax3D.SetZTitle('Z')
        self.ax3D.SetZUnits('mm')
        self.ax3D.SetBounds(limits)
        self.ax3D.SetZAxisRange(limits[-2]*scale[2],limits[-1]*scale[2])
        self.ax3D.SetXAxisRange(limits[0]*scale[0],limits[1]*scale[0])
        self.ax3D.SetYAxisRange(limits[2]*scale[1],limits[3]*scale[1])
        self.ax3D.SetCamera(self.ren.GetActiveCamera())
        self.ren.AddActor(self.ax3D)
        if not(self.ren.GetBackground()==(0.1, 0.2, 0.4)):
            flip_colors(self.ren,self.ax3D)

def apply_trans(P,T):
    '''
    Apply rotation/reflection and translation in homogeneous matrix T on a discrete basis to a Nx3 point cloud P
    '''

    return np.dot(P,T[0:3,0:3])+T[0:3,-1]


if __name__ == '__main__':
    if len(sys.argv)>1:
        aa_def(sys.argv[1])
    else:
        aa_def()