#!/usr/bin/env python
'''
Uses VTK python to allow for fitting an averaged dataset associated with the
 contour method. Full interaction requires a 3-button mouse and keyboard.
-------------------------------------------------------------------------------
Current mapping is as follows:
LMB - rotate about point cloud centroid.
MMB - pan
RMB - zoom
1 - view 1, default, looks down z axis onto xy plane
2 - view 2, looks down x axis onto zy plane
3 - view 3, looks down y axis onto zx plane
z - increase z-aspect ratio
x - decrease z-aspect ratio
c - return to default z-aspect
f - flip colors from white on dark to dark on white
i - save output to .png in current working directory
r - remove/reinstate compass/axes
o - remove/reinstate outline
-------------------------------------------------------------------------------
ver 1.1 17-17-03
1.1 - Initial release
1.2 - Refactored for PyQt5 & Python 3.x
1.3 - Refactored to handle self-restraint features
1.4 - Fixed deprecated scipy/numpy 'list-like' issues when saving the spline
1.5 - Fixed wrong index assigned to 'tck_c' when saving the spline
'''
__author__ = "M.J. Roy"
__version__ = "1.5"
__email__ = "matthew.roy@manchester.ac.uk"
__status__ = "Experimental"

import os,sys,time
import vtk
from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
import scipy.io as sio
from scipy.interpolate import griddata,bisplrep,bisplev
from scipy.spatial.distance import pdist, squareform
from scipy.spatial import Delaunay
from matplotlib import path
from pkg_resources import Requirement, resource_filename
from pyCM.pyCMcommon import *


def sf_def(*args, **kwargs):
    app = QtWidgets.QApplication(sys.argv)
    
    spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
    splash_pix = QtGui.QPixmap(spl_fname,'PNG')
    splash = QtWidgets.QSplashScreen(splash_pix)
    splash.setMask(splash_pix.mask())

    splash.show()
    app.processEvents()
    
    window = surf_int(None)

    if len(args)==1:
        surf_int.get_input_data(window,args[0])
    else:
        surf_int.get_input_data(window,None)
    
    window.show()
    splash.finish(window)
    window.iren.Initialize() # Need this line to actually show the render inside Qt
    
    ret = app.exec_()
    
    if sys.stdin.isatty() and not hasattr(sys,'ps1'):
        sys.exit(ret)
    else:
        return window

class sf_main_window(object):

    def setupUi(self, MainWindow):
        MainWindow.setWindowTitle("pyCM - surface fitting v%s" %__version__)
        MainWindow.setWindowIcon(QtGui.QIcon(resource_filename("pyCM","meta/pyCM_icon.png")))
        if hasattr(MainWindow,'setCentralWidget'):
            MainWindow.setCentralWidget(self.centralWidget)
        else:
            self.centralWidget=MainWindow
        self.mainlayout=QtWidgets.QGridLayout(self.centralWidget)

        self.vtkWidget = QVTKRenderWindowInteractor(self.centralWidget)
        mainUiBox = QtWidgets.QGridLayout()
        sectionBox = QtWidgets.QGridLayout()

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

        self.pickLabel=QtWidgets.QLabel("Remove points")
        self.pickLabel.setFont(headFont)
        self.pickHelpLabel=QtWidgets.QLabel("Press R to activate")
        self.pickActiveLabel=QtWidgets.QLabel("Pick active")
        self.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }");
        self.pickActiveLabel.setFont(QtGui.QFont("Helvetica",italic=True))
        self.undoLastPickButton=QtWidgets.QPushButton('Undo last pick')
        self.reloadButton = QtWidgets.QPushButton('Undo all/reload')

        splineLabel=QtWidgets.QLabel("Bivariate spline fitting")
        splineLabel.setFont(QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold))
        self.numLabel1=QtWidgets.QLabel("Knot spacing (x)")
        self.numEdit1 = QtWidgets.QDoubleSpinBox()
        self.numEdit1.setMaximum(10000)
        self.numEdit1.setMinimum(0.0000001)
        self.numEdit1.setValue(4)
        self.numEdit1.setDecimals(3)
        
        
        self.numLabel2=QtWidgets.QLabel("Knot spacing (y)")
        self.numEdit2 = QtWidgets.QDoubleSpinBox()
        self.numEdit2.setMaximum(10000)
        self.numEdit2.setMinimum(0.0000001)
        self.numEdit2.setValue(4)
        self.numEdit2.setDecimals(3)
        
        self.numLabel3=QtWidgets.QLabel("Spline order (x)")
        self.numEdit3 = QtWidgets.QSpinBox()
        self.numEdit3.setValue(3)
        self.numEdit3.setMinimum(1)
        self.numEdit3.setMaximum(5)
        # self.numEdit3.setMinimumWidth(50)
        
        self.numLabel4=QtWidgets.QLabel("Spline order (y)")
        self.numEdit4 = QtWidgets.QSpinBox()
        self.numEdit4.setValue(3)
        self.numEdit4.setMinimum(1)
        self.numEdit4.setMaximum(5)
        # self.numEdit4.setMinimumWidth(50)
        

        numLabel6=QtWidgets.QLabel("Display resolution:")
        self.drx=QtWidgets.QDoubleSpinBox()
        self.drx.setMaximum(10000)
        self.drx.setMinimum(0.0000001)
        self.drx.setDecimals(3)
        drx_label=QtWidgets.QLabel("x")

        self.yMin=QtWidgets.QLineEdit()
        self.dry=QtWidgets.QDoubleSpinBox()
        self.dry.setMaximum(10000)
        self.dry.setMinimum(0.0000001)
        self.dry.setDecimals(3)
        dry_label=QtWidgets.QLabel("y")
        
        self.updateButton = QtWidgets.QPushButton('Fit')

        #statLabel is what the ui is doing, not to be confused with statLabel
        self.statLabel=QtWidgets.QLabel("Idle")
        self.statLabel.setWordWrap(True)
        self.statLabel.setFont(QtGui.QFont("Helvetica",italic=True))

        horizLine1=QtWidgets.QFrame()
        horizLine1.setFrameStyle(QtWidgets.QFrame.HLine)
        #set size of line
        horizLine2=QtWidgets.QFrame()
        horizLine2.setFrameStyle(QtWidgets.QFrame.HLine)
        horizLine3=QtWidgets.QFrame()
        horizLine3.setFrameStyle(QtWidgets.QFrame.HLine)
        horizLine4=QtWidgets.QFrame()
        horizLine4.setFrameStyle(QtWidgets.QFrame.HLine)        

        sectionLabel=QtWidgets.QLabel("Data sectioning")
        sectionLabel.setFont(QtGui.QFont("Helvetica [Cronyx]",weight=QtGui.QFont.Bold))
        sectionIntList=[]
        self.xMin=QtWidgets.QLineEdit()
        sectionIntList.append(self.xMin)
        x0_label=QtWidgets.QLabel("x0")

        self.xMax=QtWidgets.QLineEdit()
        sectionIntList.append(self.xMax)
        x1_label=QtWidgets.QLabel("x1")

        self.yMin=QtWidgets.QLineEdit()
        y0_label=QtWidgets.QLabel("y0")
        sectionIntList.append(self.yMin)
        
        self.yMax=QtWidgets.QLineEdit()
        y1_label=QtWidgets.QLabel("y1")
        sectionIntList.append(self.yMax)
        self.sectionButton = QtWidgets.QPushButton('Section')
        self.revertButton = QtWidgets.QPushButton('Revert')
        self.writeButton = QtWidgets.QPushButton('Write')
        # self.writeButton.setMinimumWidth(50)

        #mainUiBox is the main container, the sectionBox is nested within
        mainUiBox.addWidget(horizLine1,0,0,1,2)
        mainUiBox.addWidget(self.pickLabel,1,0,1,2)
        mainUiBox.addWidget(self.pickHelpLabel,2,0,1,1)
        mainUiBox.addWidget(self.pickActiveLabel,2,1,1,1)
        mainUiBox.addWidget(self.undoLastPickButton,3,0,1,2)
        mainUiBox.addWidget(self.reloadButton,4,0,1,2)
        mainUiBox.addWidget(horizLine2,5,0,1,2)

        mainUiBox.addWidget(splineLabel,6,0,1,2)
        mainUiBox.addWidget(self.numLabel1,7,0,1,1)
        mainUiBox.addWidget(self.numEdit1,7,1,1,1)
        mainUiBox.addWidget(self.numLabel2,8,0,1,1)
        mainUiBox.addWidget(self.numEdit2,8,1,1,1)
        mainUiBox.addWidget(self.numLabel3,9,0,1,1)
        mainUiBox.addWidget(self.numEdit3,9,1,1,1)
        mainUiBox.addWidget(self.numLabel4,10,0,1,1)
        mainUiBox.addWidget(self.numEdit4,10,1,1,1)
        mainUiBox.addWidget(self.updateButton,12,0,1,2)
        mainUiBox.addWidget(numLabel6,13,0,1,2)
        
        displayLayout = QtWidgets.QGridLayout()
        displayLayout.addWidget(drx_label,0,0)
        displayLayout.addWidget(self.drx,0,1)
        displayLayout.addWidget(dry_label,0,2)
        displayLayout.addWidget(self.dry,0,3)
        
        mainUiBox.addLayout(displayLayout,14,0,1,2)
        mainUiBox.addWidget(horizLine2,15,0,1,2)
        mainUiBox.addLayout(sectionBox,16,0,5,2)
        mainUiBox.addWidget(horizLine3,21,0,1,2)
        mainUiBox.addWidget(self.writeButton,22,0,1,2)
        mainUiBox.addWidget(horizLine4,23,0,1,2)

        
        sectionBox.addWidget(sectionLabel,0,0,1,4)
        sectionBox.addWidget(x0_label,1,0,1,1)
        sectionBox.addWidget(sectionIntList[0],1,1,1,1)
        sectionBox.addWidget(x1_label,1,2,1,1)
        sectionBox.addWidget(sectionIntList[1],1,3,1,1)
        sectionBox.addWidget(y0_label,2,0,1,1)
        sectionBox.addWidget(sectionIntList[2],2,1,1,1)
        sectionBox.addWidget(y1_label,2,2,1,1)
        sectionBox.addWidget(sectionIntList[3],2,3,1,1)
        sectionBox.addWidget(self.sectionButton,3,0,1,2)
        sectionBox.addWidget(self.revertButton,3,2,1,2)

        lvLayout=QtWidgets.QVBoxLayout()
        lvLayout.addLayout(mainUiBox)
        lvLayout.addStretch(1)
    
        self.mainlayout.addWidget(self.vtkWidget,0,0,1,1)
        self.mainlayout.addLayout(lvLayout,0,1,1,1)
        self.mainlayout.addWidget(self.statLabel,1,0,1,2)

    
    def initialize(self):
        self.vtkWidget.start()


class surf_int(QtWidgets.QWidget):

    def __init__(self, parent):
        super(surf_int,self).__init__(parent)
        self.ui = sf_main_window()
        self.ui.setupUi(self)
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(0.1, 0.2, 0.4)
        self.ui.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.ui.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.ren.GetActiveCamera().ParallelProjectionOn()
        self.cp=self.ren.GetActiveCamera().GetPosition()
        self.fp=self.ren.GetActiveCamera().GetFocalPoint()
        self.iren.AddObserver("KeyPressEvent", self.keypress)

        self.PointSize=2
        self.LineWidth=1
        self.Zaspect=1.0
        self.limits=np.empty(6)
        self.picking=False
        self.fitted=False

        self.ui.updateButton.clicked.connect(lambda: self.onUpdateSpline())
        self.ui.sectionButton.clicked.connect(lambda: self.Cut())
        self.ui.revertButton.clicked.connect(lambda: self.RemoveCut())
        self.ui.writeButton.clicked.connect(lambda: self.write())
        self.ui.reloadButton.clicked.connect(lambda: self.get_input_data())
        self.ui.undoLastPickButton.clicked.connect(lambda: self.undo_pick())
        self.ui.numEdit1.valueChanged.connect(self.changeUpdateBackground)
        self.ui.numEdit2.valueChanged.connect(self.changeUpdateBackground)
        self.ui.numEdit3.valueChanged.connect(self.changeUpdateBackground)
        self.ui.numEdit4.valueChanged.connect(self.changeUpdateBackground)
        self.ui.drx.valueChanged.connect(self.changeUpdateBackground)
        self.ui.dry.valueChanged.connect(self.changeUpdateBackground)
        
    def changeUpdateBackground(self):
        self.ui.updateButton.setStyleSheet("background-color : None")

    def get_input_data(self,filem):
        """
        Loads the content of a *.mat file pertaining to this particular step
        """
        
        if hasattr(self,'pointActor'): #then remove everything
            self.ren.RemoveActor(self.pointActor)
            self.ren.RemoveActor(self.outlineActor)
        
        if hasattr(self,'splineActor'):
            self.ren.RemoveActor(self.splineActor)

        
        if filem == None:
            filem, _, =get_file('*.mat')
            
        if filem:
            mat_contents = sio.loadmat(filem)
            self.fileo=filem
            
            try:
                pts=mat_contents['aa']['pnts'][0][0]
                refTrans=mat_contents['trans']['ref'][0][0]
                

                self.pts=pts[~np.isnan(pts).any(axis=1)] #remove all nans
                self.RefOutline=np.concatenate(mat_contents['ref']['x_out'],axis=0)[0]
                for transformation in refTrans:
                    self.RefOutline = np.dot(self.RefOutline,transformation[0:3,0:3])+transformation[0:3,-1]
                
                self.RefMin=np.amin(self.RefOutline,axis=0)
                self.RefMax=np.amax(self.RefOutline,axis=0)
                self.limits=get_limits(np.vstack((self.pts,self.RefOutline)))
                self.ui.xMin.setText('%.3f'%self.limits[0])
                self.ui.xMax.setText('%.3f'%self.limits[1])
                self.ui.yMin.setText('%.3f'%self.limits[2])
                self.ui.yMax.setText('%.3f'%self.limits[3])
                
                self.bool_pnt=np.ones(len(self.pts), dtype=bool) #initialize mask
                
                #Generate actors
                color=(int(0.2784*255),int(0.6745*255),int(0.6941*255))
                self.vtkPntsPolyData, self.pointActor, self.colors, = gen_point_cloud(self.pts,color,self.PointSize)
                

                
                self.ren.AddActor(self.pointActor)
                self.outlineActor, _, = gen_outline(self.RefOutline,color,self.PointSize)
                
                self.ren.AddActor(self.outlineActor)

                #add axes
                self.axisActor = add_axis(self.ren,self.limits,[1,1,1])
                
                if 'spline_x' in mat_contents: #then it can be displayed & settings displayed
                    order=mat_contents['spline_x']['order'][0][0][0]
                    #build tck list from mat contents
                    self.tck = [mat_contents['spline_x']['tck_x'][0][0][0], \
                    mat_contents['spline_x']['tck_y'][0][0][0], \
                    mat_contents['spline_x']['tck_c'][0][0][0], \
                    order[0], order[1]]
                    spacing=mat_contents['spline_x']['kspacing'][0][0][0]
                    # smoothing=mat_contents['spline_x']['smooth'][0][0]
                    self.gx,self.gy=spacing[0],spacing[1]
                    self.ui.numEdit1.setValue(self.gx)
                    self.ui.numEdit2.setValue(self.gy)
                    self.ui.numEdit3.setValue(int(order[0]))
                    self.ui.numEdit4.setValue(int(order[1]))
                    self.bool_pnt=mat_contents['aa_mask'][0]

                    #find points to be painted red
                    localind=np.asarray(range(len(self.bool_pnt)))
                    localind=localind[np.where(np.logical_not(self.bool_pnt))]
                    
                    for i in localind:
                        #turn them red
                        self.colors.SetTuple(i,(255,0,0))

                    self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
                    self.vtkPntsPolyData.Modified()
                    
                    
                    self.DisplayFit()
                    self.fitted=True
                else: #there is no fitted surface
                    self.fitted=False
                    self.ui.updateButton.setStyleSheet("background-color : None")
                    
            except Exception as e:
                print('Load from pyCM fit_surface failed; error message returned:')
                print(str(e))
                
                
            
        else:
            self.ui.statLabel.setText("Invalid file.")
            return
            
        self.unsaved_changes=False
        #update
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()
            
    def onUpdateSpline(self):
        p,ro,rmin,rmax=self.pts[np.where(self.bool_pnt)],self.RefOutline,self.RefMin,self.RefMax
        self.ui.statLabel.setText("Fitting . . .")
        QtWidgets.qApp.processEvents()
        #read input parameters
        self.gx=self.ui.numEdit1.value()
        self.gy=self.ui.numEdit2.value()
        kx=self.ui.numEdit3.value()
        ky=self.ui.numEdit4.value()

        tx=np.linspace(rmin[0],rmax[0],int((rmax[0]-rmin[0])/self.gx))
        ty=np.linspace(rmin[1],rmax[1],int((rmax[1]-rmin[1])/self.gy))

        #make sure both x & y have enough values in either direction
        if len(tx)<3 or len(ty)<3:
            self.ui.statLabel.setText("Grid too large . . .")
            self.ui.updateButton.setEnabled(True)
            return
        tx=np.insert(tx,0,[rmin[0]] * kx) #to make sure knots are repeated at edges
        tx=np.insert(tx,-1,[rmax[0]] * kx)
        ty=np.insert(ty,0,[rmin[1]] * ky) #to make sure knots are repeated at edges
        ty=np.insert(ty,-1,[rmax[1]] * ky)
        

        try:
            self.tck = bisplrep(p[:,0], p[:,1], p[:,2], kx=kx, ky=ky, tx=tx, ty=ty, task=-1) #get spline representation
        except ValueError as ve:
            self.ui.statLabel.setText("Last fit failed.")
            return
        
        self.DisplayFit()
            
    def DisplayFit(self):
        
        try:
            #now evaluate for show
            if not hasattr(self,'dryval'): #then it won't have drxval
                rx=np.linspace(self.RefMin[0],self.RefMax[0],int((self.RefMax[0]-self.RefMin[0])/(self.gx/4)))
                ry=np.linspace(self.RefMin[1],self.RefMax[1],int((self.RefMax[1]-self.RefMin[1])/(self.gy/4)))
                self.drxval,self.dryval=rx[1]-rx[0],ry[1]-ry[0]
                self.ui.drx.setValue(self.drxval)
                self.ui.dry.setValue(self.dryval)
            else: #read the res directly from UI, apply and update
                self.drxval,self.dryval=self.ui.drx.value(),self.ui.dry.value()
                rx=np.linspace(self.RefMin[0],self.RefMax[0],int((self.RefMax[0]-self.RefMin[0])/(self.drxval)))
                ry=np.linspace(self.RefMin[1],self.RefMax[1],int((self.RefMax[1]-self.RefMin[1])/(self.dryval)))
                self.drxval,self.dryval=rx[1]-rx[0],ry[1]-ry[0]
                self.ui.drx.setValue(self.drxval)
                self.ui.dry.setValue(self.dryval)
                
            grid_x, grid_y = np.meshgrid(rx,ry,indexing='xy')
            grid_x=np.transpose(grid_x)
            grid_y=np.transpose(grid_y)

            points2D=np.vstack([grid_x.flatten(),grid_y.flatten()]).T
            tri=Delaunay(points2D)
            tri= tri.simplices
            
            znew = bisplev(grid_x[:,0], grid_y[0,:], self.tck)

            points3D=np.vstack([points2D[:,0],points2D[:,1],znew.flatten()]).T
            
            #crop according to what's in the outline
            #temporarily move the outline to the first quadrant
            t_offset=np.array([0.1*self.RefMin[0],0.1*self.RefMin[1]])
            
            for j in range(0,2):
                if t_offset[j]<0:
                    points3D[:,j]=points3D[:,j]-t_offset[j]
                    self.RefOutline[:,j]=self.RefOutline[:,j]-t_offset[j]


            pth=path.Path(self.RefOutline[:,:2])
            inOutline=pth.contains_points(points3D[:,:2])

            ind=np.array(range(0,len(inOutline)))
            ind_out=ind[np.invert(inOutline)]
            
            #remove all triangles that contain points outside the perimeter
            tri = tri[(np.isin(tri,ind_out)==(0,0,0)).all(axis=1),:]

            #move things back
            for j in range(0,2):
                if t_offset[j]<0:
                    points3D[:,j]=points3D[:,j]+t_offset[j]
                    self.RefOutline[:,j]=self.RefOutline[:,j]+t_offset[j]

            self.ui.statLabel.setText("Rendering . . .")
            self.DisplaySplineFit(points3D,tri)

            zeval = np.empty(np.size(self.pts[:,2]))
            for i in range(len(zeval)):
                zeval[i] = bisplev(self.pts[i,0], self.pts[i,1], self.tck)

            a=(zeval-self.pts[:,2])
            RSME=(np.sum(a*a)/len(a))**0.5

            self.ui.statLabel.setText("RSME: %2.2f micron."%(RSME*1000))
            
        except Exception as e:
            print(str(e))
            self.ui.statLabel.setText("Failed to show fit.")
        
        self.ui.updateButton.setEnabled(True)
        self.ui.updateButton.setStyleSheet("background-color :rgb(77, 209, 97);")

        
    def DisplaySplineFit(self,p,t):

        if hasattr(self,'splineActor'):
            self.ren.RemoveActor(self.splineActor)
        self.SplinePoints = vtk.vtkPoints()
        triangles = vtk.vtkCellArray()
        
        #load up points
        for i in p:
            self.SplinePoints.InsertNextPoint(i)

        for i in t:
            triangle=vtk.vtkTriangle()
            for j in range(0,3):
                triangle.GetPointIds().SetId(j,i[j])
            triangles.InsertNextCell(triangle)

        trianglePolyData = vtk.vtkPolyData()
        trianglePolyData.SetPoints(self.SplinePoints)
        trianglePolyData.SetPolys(triangles)
        #filter so that edges are shared
        self.cSplinePolyData = vtk.vtkCleanPolyData()
        self.cSplinePolyData.SetInputData(trianglePolyData)

        # Create a mapper and actor for smoothed dataset
        self.Smapper = vtk.vtkPolyDataMapper()
        self.Smapper.SetInputConnection(self.cSplinePolyData.GetOutputPort())

        self.splineActor = vtk.vtkActor()
        self.splineActor.SetMapper(self.Smapper)
        self.splineActor.GetProperty().SetInterpolationToFlat()
        self.splineActor.GetProperty().SetRepresentationToSurface()
        self.splineActor.GetProperty().SetColor(1,0.804,0.204)
        self.splineActor.GetProperty().SetOpacity(0.75)
        self.splineActor.GetProperty().SetLighting(False)
        self.splineActor.SetScale(1,1,self.Zaspect)
        self.ren.AddActor(self.splineActor)
        self.ui.vtkWidget.update()


    def Cut(self):
        pts=np.array([float(self.ui.xMin.text()),float(self.ui.xMax.text()),float(self.ui.yMin.text()),float(self.ui.yMax.text())])

        planex1 = vtk.vtkPlane()
        planex1.SetOrigin(pts[0],0,0)
        planex1.SetNormal(1,0,0)
        
        planex2 = vtk.vtkPlane()
        planex2.SetOrigin(pts[1],0,0)
        planex2.SetNormal(-1,0,0)
        
        planey1 = vtk.vtkPlane()
        planey1.SetOrigin(0,pts[2],0)
        planey1.SetNormal(0,1,0)
        
        planey2 = vtk.vtkPlane()
        planey2.SetOrigin(0,pts[3],0)
        planey2.SetNormal(0,-1,0)
        
        planeCollection = vtk.vtkPlaneCollection()
        planeCollection.AddItem(planex1)
        planeCollection.AddItem(planex2)
        planeCollection.AddItem(planey1)
        planeCollection.AddItem(planey2)
        
        Omapper=self.outlineActor.GetMapper()
        Pmapper=self.pointActor.GetMapper()
        Omapper.SetClippingPlanes(planeCollection)
        Pmapper.SetClippingPlanes(planeCollection)

        if hasattr(self,'Smapper'):
            self.Smapper.SetClippingPlanes(planeCollection)
        nl=np.array(self.axisActor.GetBounds())
        self.ren.RemoveActor(self.axisActor)
        #add axes
        self.axisActor = add_axis(self.ren,self.limits,[1,1,1])

        
        #update
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()
    
    def RemoveCut(self):
        self.ui.xMin.setText('%.3f'%self.limits[0])
        self.ui.xMax.setText('%.3f'%self.limits[1])
        self.ui.yMin.setText('%.3f'%self.limits[2])
        self.ui.yMax.setText('%.3f'%self.limits[3])
        Omapper=self.outlineActor.GetMapper()
        Pmapper=self.pointActor.GetMapper()
        Omapper.RemoveAllClippingPlanes()
        Pmapper.RemoveAllClippingPlanes()
        if hasattr(self,'Smapper'):
            self.Smapper.RemoveAllClippingPlanes()
        self.ren.RemoveActor(self.axisActor)
        #add axes
        self.axisActor = add_axis(self.ren,self.limits,[1,1,1])

        
        #update
        self.ren.ResetCamera()
        self.ui.vtkWidget.update()
        self.ui.vtkWidget.setFocus()


    def keypress(self,obj,event):
        key = obj.GetKeyCode()

        if key =="1":
            XYView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
        elif key =="2":
            YZView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)
        elif key =="3":
            XZView(self.ren, self.ren.GetActiveCamera(),self.cp,self.fp)

        elif key=="z":
            self.Zaspect=self.Zaspect*2
            s,nl,axs=self.get_scale()
            if hasattr(self,'splineActor'):
                self.splineActor.SetScale(s)
                self.splineActor.Modified()
            if hasattr(self,'pointActor'):
                self.pointActor.SetScale(s)
                self.pointActor.Modified()
            self.ren.RemoveActor(self.axisActor)
            self.axisActor = add_axis(self.ren,nl,axs)

        elif key=="x":
            self.Zaspect=self.Zaspect*0.5
            s,nl,axs=self.get_scale()
            if hasattr(self,'splineActor'):
                self.splineActor.SetScale(s)
                self.splineActor.Modified()
            if hasattr(self,'pointActor'):
                self.pointActor.SetScale(s)
                self.pointActor.Modified()
            self.ren.RemoveActor(self.axisActor)
            self.axisActor = add_axis(self.ren,nl,axs)

        elif key=="c":
            self.Zaspect=1.0
            s,_,_,=self.get_scale()
            if hasattr(self,'splineActor'):
                self.splineActor.SetScale(s)
                self.splineActor.Modified()
            if hasattr(self,'pointActor'):
                self.pointActor.SetScale(s)
                self.pointActor.Modified()

            self.ren.RemoveActor(self.axisActor)
            self.axisActor = add_axis(self.ren,self.limits,[1,1,1])
            self.ren.ResetCamera()

        elif key=="i":
            im = vtk.vtkWindowToImageFilter()
            writer = vtk.vtkPNGWriter()
            im.SetInput(self.ui.vtkWidget._RenderWindow)
            im.Update()
            writer.SetInputConnection(im.GetOutputPort())
            writer.SetFileName("spline_fit.png")
            writer.Write()
            print("Screen output saved to %s" %os.path.join(os.getcwd(),'spline_fit.png'))
        
        elif key=="a":
            FlipVisible(self.axisActor)
            
        elif key =="f": #flip color scheme for printing
            FlipColors(self.ren,self.axisActor,None)
            FlipColors(self.ren,self.pointActor,1)
            if hasattr(self,'splineActor'):
                FlipColors(self.ren,self.splineActor,0)
                
        elif key == "o":
            FlipVisible(self.outlineActor)
            
        elif key == "Z":
            self.PointSize=updatePointSize(self.pointActor,self.PointSize*2)
            if hasattr(self,'splineActor'):
                self.LineWidth=updateLineWidth(self.splineActor,self.LineWidth*2)
            
        elif key == "X":
            self.PointSize=updatePointSize(self.pointActor,self.PointSize*0.5)
            if hasattr(self,'splineActor'):
                self.LineWidth=updateLineWidth(self.splineActor,self.LineWidth*0.5)
            
        elif key == "C":
            self.PointSize=updatePointSize(self.pointActor,1)
            if hasattr(self,'splineActor'):
                self.LineWidth=updateLineWidth(self.splineActor,1)

        elif key == "e":
            self.write()
        
        elif key == "q":
            if sys.stdin.isatty():
                sys.exit("Surface fitting complete.")
            else:
                print("Surface fitting completed.")
                return
                
        elif key=="l":
            self.get_input_data(None)
            
            
        elif key == "r":
            if self.picking == True:
                self.picking =False
                self.show_picking()
            else:
                self.picking =True
                self.show_picking()
                self.start_pick()

        self.ui.vtkWidget.update()

    def show_picking(self):
        #Updates when the 'r' button is pressed to provide a link between VTK & Qt hooks
        if self.picking == True:
            self.ui.pickActiveLabel.setStyleSheet("QLabel { background-color : red; color : white; }");
        else:
            self.ui.pickActiveLabel.setStyleSheet("QLabel { background-color : gray; color : darkGray; }");
    
    def start_pick(self):
        #Required to change interactor
        style=vtk.vtkInteractorStyleRubberBandPick()
        self.iren.SetInteractorStyle(style)
        picker = vtk.vtkAreaPicker()
        self.iren.SetPicker(picker)
        picker.AddObserver("EndPickEvent", self.picker_callback)
    
    def picker_callback(self,obj,event):
        
        extract = vtk.vtkExtractSelectedFrustum()
    
        fPlanes=obj.GetFrustum() #collection of planes based on unscaled display
    
        #scale frustum to account for the zaspect
        scaledPlanes=vtk.vtkPlanes()
        scaledNormals=vtk.vtkDoubleArray()
        scaledNormals.SetNumberOfComponents(3)
        scaledNormals.SetNumberOfTuples(6)
        scaledOrigins=vtk.vtkPoints()
        for j in range(6):
            i=fPlanes.GetPlane(j)
            k=i.GetOrigin()
            q=i.GetNormal()
            scaledOrigins.InsertNextPoint(k[0],k[1],k[2]/float(self.Zaspect))
            scaledNormals.SetTuple(j,(q[0],q[1],q[2]*float(self.Zaspect)))
        scaledPlanes.SetNormals(scaledNormals)
        scaledPlanes.SetPoints(scaledOrigins)
            
        
        extract.SetFrustum(scaledPlanes)
        extract.SetInputData(self.vtkPntsPolyData)
        extract.Update()
        extracted = extract.GetOutput()
        
        ids = vtk.vtkIdTypeArray()
        ids = extracted.GetPointData().GetArray("vtkOriginalPointIds")

        
        if ids:
            #store them in an array for an undo operation
            self.lastSelectedIds=ids
            for i in range(ids.GetNumberOfTuples()):
                #turn them red
                self.colors.SetTuple(ids.GetValue(i),(255,0,0))
                self.bool_pnt[ids.GetValue(i)]=False
        
            self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
            self.vtkPntsPolyData.Modified()
        
        
        self.ui.vtkWidget.update()
        #set flag on ui to show that data has been modified
        self.unsaved_changes=True

    def undo_pick(self):
        if hasattr(self,"lastSelectedIds"):
            for i in range(self.lastSelectedIds.GetNumberOfTuples()):
                #turn them from red to starting color
                self.colors.SetTuple(self.lastSelectedIds.GetValue(i),(70, 171, 176))
                self.bool_pnt[self.lastSelectedIds.GetValue(i)]=True
            self.vtkPntsPolyData.GetPointData().SetScalars(self.colors)
            self.vtkPntsPolyData.Modified()
            self.ui.vtkWidget.update()
        else:
            self.ui.statLabel.setText("No picked selection to revert.")
            
    def write(self):
        
        mat_vars=sio.whosmat(self.fileo)
        if not set(['spline_x']).isdisjoint([item for sublist in mat_vars for item in sublist]): #tell the user that they might overwrite their data
            ret=QtWidgets.QMessageBox.warning(self, "pyCM Warning", \
                "There is already data associated with this analysis step saved. Overwrite and invalidate subsequent steps?", \
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, QtWidgets.QMessageBox.No)
            if ret == QtWidgets.QMessageBox.No: #don't overwrite
                return
            else:
                #delete fitting parameters with pyCMcommon helper function, which negates key FEA parameters.
                clear_mat(self.fileo,['vtk','pickedCornerInd','FEA']) 
        
        if hasattr(self,'tck'): #then spline fitting has been done
            mat_contents=sio.loadmat(self.fileo)
            coefs=[np.reshape(self.tck[2],(len(self.tck[0])-self.tck[3]-1,-1))]
            number=np.array([len(self.tck[0]),len(self.tck[1])])

            #16/07/2020 - scipy.io.savemat does not support 'list-like' arrays. Therefore the spline object now contains the elements of the tck list.
            new={'spline_x': {'form': 'B-', 'kspacing': [self.gx, self.gy], 'coefs': coefs, 'number': number, 'tck_x': self.tck[0],'tck_y': self.tck[1], 'tck_c': self.tck[2], 'order':np.array([self.tck[3],self.tck[4]])},  'x_out':self.RefOutline, 'aa_mask':self.bool_pnt}
            
            
            
            
            mat_contents.update(new)
            sio.savemat(self.fileo,mat_contents)
            
            self.ui.statLabel.setText("Output written.")
            self.fitted=True
            self.unsaved_changes=False
        else:
            self.ui.statLabel.setText("Nothing to write.")

    def get_scale(self):
        '''
        Returns array for the keypress function
        '''
        s=np.array([1,1,self.Zaspect])
        nl=np.append(self.limits[0:4],([self.limits[-2]*self.Zaspect,self.limits[-1]*self.Zaspect]))
        axs=np.array([1,1,1/self.Zaspect])
        return s,nl,axs        
            



def XYView(renderer, camera,cp,fp):
    camera.SetPosition(0,0,cp[2]+0)
    camera.SetFocalPoint(fp)
    camera.SetViewUp(0,1,0)
    camera.OrthogonalizeViewUp()
    camera.ParallelProjectionOn()
    renderer.ResetCamera()

def YZView(renderer, camera,cp,fp):
    camera.SetPosition(cp[2]+0,0,0)
    camera.SetFocalPoint(fp)
    camera.SetViewUp(0,0,1)
    camera.OrthogonalizeViewUp()
    camera.ParallelProjectionOn()
    renderer.ResetCamera()


def XZView(renderer,camera,cp,fp):
    vtk.vtkObject.GlobalWarningDisplayOff() #otherwise there's crystal eyes error . . .
    camera.SetPosition(0,cp[2]+0,0)
    camera.SetFocalPoint(fp)
    camera.SetViewUp(0,0,1)
    camera.OrthogonalizeViewUp()
    camera.ParallelProjectionOn()
    renderer.ResetCamera()
    


def FlipVisible(actor):
    if actor.GetVisibility():
        actor.VisibilityOff()
    else:
        actor.VisibilityOn()

def updatePointSize(actor,NewPointSize):
    actor.GetProperty().SetPointSize(NewPointSize)
    actor.Modified()
    return NewPointSize

def updateLineWidth(actor,NewLineWidth):
    actor.GetProperty().SetLineWidth(NewLineWidth)
    actor.Modified()
    return NewLineWidth

def FlipColors(ren,actor,contrast):
    if ren.GetBackground()==(0.1, 0.2, 0.4):
        if hasattr(actor,'GetXAxesLinesProperty'):
            actor.GetTitleTextProperty(0).SetColor(0,0,0)
            actor.GetLabelTextProperty(0).SetColor(0,0,0)
            actor.GetXAxesLinesProperty().SetColor(0,0,0)
            actor.SetXTitle('x') #there's a vtk bug here . . .
            
            actor.GetTitleTextProperty(1).SetColor(0,0,0)
            actor.GetLabelTextProperty(1).SetColor(0,0,0)
            actor.GetYAxesLinesProperty().SetColor(0,0,0)

            actor.GetTitleTextProperty(2).SetColor(0,0,0)
            actor.GetLabelTextProperty(2).SetColor(0,0,0)
            actor.GetZAxesLinesProperty().SetColor(0,0,0)
            ren.SetBackground(1, 1, 1)
        else:
            if contrast == 1:
                actor.GetProperty().SetColor(0.2784,0.6745,0.6941)
            else:
                actor.GetProperty().SetColor(1,0.804,0.204)

    else:
        if hasattr(actor,'GetXAxesLinesProperty'):
            actor.GetTitleTextProperty(0).SetColor(1,1,1)
            actor.GetLabelTextProperty(0).SetColor(1,1,1)
            actor.GetXAxesLinesProperty().SetColor(1,1,1)
            actor.SetXTitle('X')
            
            actor.GetTitleTextProperty(1).SetColor(1,1,1)
            actor.GetLabelTextProperty(1).SetColor(1,1,1)
            actor.GetYAxesLinesProperty().SetColor(1,1,1)
            actor.SetYTitle('Y')
            
            actor.GetTitleTextProperty(2).SetColor(1,1,1)
            actor.GetLabelTextProperty(2).SetColor(1,1,1)
            actor.GetZAxesLinesProperty().SetColor(1,1,1)
            ren.SetBackground(0.1, 0.2, 0.4)
        else:
            if contrast == 1:
                actor.GetProperty().SetColor(0.0353, 0.1922, 0.2706)
            else:
                actor.GetProperty().SetColor(0.8039, 0.3490, 0.2902)

if __name__ == "__main__":
    if len(sys.argv)>1:
        sf_def(sys.argv[1])
    else:
        sf_def()
