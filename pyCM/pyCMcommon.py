import os,io,sys,yaml
import vtk
from vtkmodules.vtkCommonColor import vtkNamedColors, vtkColorSeries
import numpy as np
import vtk.util.numpy_support as v2n
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.numpy_interface import dataset_adapter as dsa
from PyQt5 import QtCore, QtGui, QtWidgets
from pkg_resources import Requirement, resource_filename
from importlib.metadata import version
from scipy.interpolate import interp1d
from pyclipper import PyclipperOffset, scale_to_clipper, scale_from_clipper, JT_SQUARE, ET_CLOSEDPOLYGON
from matplotlib import path
import h5py

from datetime import datetime

def make_splash():
    '''
    Makes and returns the pyCM Qt splash window object
    '''
    spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
    splash_pix = QtGui.QPixmap(spl_fname,'PNG')
    splash = QtWidgets.QSplashScreen(splash_pix, QtCore.Qt.SplashScreen)
    splash.setWindowFlags(QtCore.Qt.WindowStaysOnTopHint | QtCore.Qt.FramelessWindowHint)
    
    font = splash.font()
    font.setPixelSize(20)
    font.setWeight(QtGui.QFont.Bold)
    splash.setFont(font)
    
    splash.showMessage(version('pyCM'),QtCore.Qt.AlignRight | QtCore.Qt.AlignBottom, QtCore.Qt.darkGray)
    return splash

def make_logo(ren):
    spl_fname=resource_filename("pyCM","meta/pyCM_logo.png")
    img_reader = vtk.vtkPNGReader()
    img_reader.SetFileName(spl_fname)
    img_reader.Update()
    logo = vtk.vtkLogoRepresentation()
    logo.SetImage(img_reader.GetOutput())
    logo.ProportionalResizeOn()
    logo.SetPosition( 0.1, 0.1 ) #lower left
    logo.SetPosition2( 0.8, 0.8 ) #upper right
    logo.GetImageProperty().SetDisplayLocationToBackground()
    ren.AddViewProp(logo)
    logo.SetRenderer(ren)
    return logo

def get_file(*args):
    '''
    Returns absolute path to filename and the directory it is located in from a PyQt5 filedialog. First value is file extension, second is a string which overwrites the window message.
    '''
    ext = args[0]
    if len(args)>1:
        launchdir = args[1]
    else: launchdir = os.getcwd()
    ftypeName={}
    ftypeName['*.txt']=["Select point cloud data file:", "*.txt", "TXT File"]
    ftypeName['*.csv']=["Select point cloud data file:", "*.csv", "PC-DMIS point measurement file"]
    ftypeName['*.dat']=["Select unregistered point cloud data file:", "*.dat", "NanoFocus Origin format"]
    ftypeName['*.pyCM'] = ["pyCM HDF5-format data file", "*.pyCM", "pyCM file"]
    ftypeName['*.vtk']=["Select the legacy VTK file:", "*.vtk", "VTK File"]
    ftypeName['*.inp']=["Select inp file:", "*.inp", "INP File"]
    ftypeName['*.*'] = ["pyCM external executable", "*.*", "..."]
    
    if ext=='*.txt':
        filer = QtWidgets.QFileDialog.getOpenFileName(None, ftypeName[ext][0], 
         launchdir, \
         (ftypeName[ext][2]+' ('+ftypeName[ext][1]+');;' \
         +ftypeName['*.csv'][2]+' ('+ftypeName['*.csv'][1]+');;'\
         +ftypeName['*.dat'][2]+' ('+ftypeName['*.dat'][1]+');;All Files (*.*)'))
    elif ext=='mesh':
        filer = QtWidgets.QFileDialog.getOpenFileName(None, "Select either legacy VTK file or FEA input file:", 
         launchdir, \
         (ftypeName['*.vtk'][2]+' ('+ftypeName['*.vtk'][1]+');;'\
         +ftypeName['*.inp'][2]+' ('+ftypeName['*.inp'][1]+');;All Files (*.*)'))
    else:
        filer = QtWidgets.QFileDialog.getOpenFileName(None, ftypeName[ext][0], 
         os.getcwd(),(ftypeName[ext][2]+' ('+ftypeName[ext][1]+');;All Files (*.*)'))

    if filer[0] == '':
            return None, None
    else:
        return filer[0], os.path.dirname(filer[0])

def get_save_file(*args):
    '''
    Returns absolute path to filename and the directory it is located in from a PyQt5 filedialog. First value is file extension, second is a string which overwrites the window message.
    '''
    ext = args[0]
    if len(args)>1:
        dir = args[1]
    else: dir = os.getcwd()
    
    file_type_names = {}
    file_type_names['*.pyCM'] = 'HDF5-formatted pyCM output file'
    file_type_names['*.csv'] = 'pyCM comma delimited output file'
    file_type_names['*.dxf'] = 'pyCM Drawing eXchange Format file'
    file_type_names['*.geo'] = 'GMSH instructions file'
    file_type_names['*.inp'] = 'pyCM FEA input file'

    filer = QtWidgets.QFileDialog.getSaveFileName(None, "Save as:", \
    dir, \
    str(file_type_names[ext]+' ('+ext+')') \
    )
    
    if filer[0] == '':
        return None, None
    else:
        return filer[0], os.path.dirname(filer[0])

def gen_outline(pts, color = (1,1,1), size = 2):
    '''
    Returns an outline actor with specified pts, color and size. Incoming pnts should be ordered.
    '''
    if color[0]<=1 and color != None:
        color=(int(color[0]*255),int(color[1]*255),int(color[2]*255))
    if color[0]>=1 and color != None:
        color=(color[0]/float(255),color[1]/float(255),color[2]/float(255))
    points=vtk.vtkPoints()

    points.SetData(v2n.numpy_to_vtk(pts))

    lineseg=vtk.vtkPolygon()
    lineseg.GetPointIds().SetNumberOfIds(len(pts))
    for i in range(len(pts)):
        lineseg.GetPointIds().SetId(i,i)
    linesegcells=vtk.vtkCellArray()
    linesegcells.InsertNextCell(lineseg)
    outline=vtk.vtkPolyData()
    outline.SetPoints(points)
    outline.SetVerts(linesegcells)
    outline.SetLines(linesegcells)
    Omapper=vtk.vtkPolyDataMapper()
    Omapper.SetInputData(outline)
    outlineActor=vtk.vtkActor()
    outlineActor.SetMapper(Omapper)
    outlineActor.GetProperty().SetColor(color)
    outlineActor.GetProperty().SetPointSize(size)
    return outlineActor, outline

def gen_point_cloud(pts,color=None,r=None,size=2):
    '''
    Returns vtk objects and actor for a point cloud having size points, returns color array associated with the actor/polydata object as well as a lookuptable for rendering a scalebar if colouration is applied based on height. color needs to be specified as an RGB 0-255 tuple.
    '''
    
    lut=None
    vtkPnts = vtk.vtkPoints()
    vtkVerts = vtk.vtkCellArray()
    
    #load up points
    vtkPnts.SetData(v2n.numpy_to_vtk(pts))

    for i in np.arange(len(pts)):
        vtkVerts.InsertNextCell(1)
        vtkVerts.InsertCellPoint(i)

    pC = vtk.vtkPolyData()
    pC.SetPoints(vtkPnts)
    pC.SetVerts(vtkVerts)

    mapper = vtk.vtkDataSetMapper()

    if color is None:
        vtk_z_array = v2n.numpy_to_vtk(pts[:,-1])
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)
        if r is None:
            lut.SetTableRange(np.amin(pts[:,-1]), np.amax(pts[:,-1]))
        else:
            lut.SetTableRange(r[0],r[1])
        lut.Build()
        colors = lut.MapScalars(vtk_z_array,vtk.VTK_COLOR_MODE_DEFAULT,-1,vtk.VTK_RGB)

        
    elif isinstance(color,tuple):
        #if assigning single color every point, lut will return as None
        colors=vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        for i in np.arange(len(pts)):
            colors.InsertNextTuple(color)
        pC.GetPointData().SetScalars(colors)
        mapper.SetInputData(pC)
        
    elif isinstance(color, str):
        vtk_z_array = v2n.numpy_to_vtk(pts[:,-1])
        lut = get_diverging_lut(color) #add options here for other baseline color series
        if r is None:
            lut.SetTableRange(np.amin(pts[:,-1]), np.amax(pts[:,-1]))
        else:
            lut.SetTableRange(r[0],r[1])
        lut.Build()
        colors = lut.MapScalars(vtk_z_array,vtk.VTK_COLOR_MODE_DEFAULT,-1,vtk.VTK_RGB)

    
    pC.GetPointData().SetScalars(colors)
    mapper.SetInputData(pC)

    actor=vtk.vtkActor()
    actor.SetMapper(mapper)

    actor.GetProperty().SetPointSize(size)
    return actor, pC, colors, lut

    
def get_limits(pts, factor = 0.1):
    '''
    Returns a bounding box with x,y values bumped out by factor
    '''
    RefMin = np.amin(pts,axis=0)
    RefMax = np.amax(pts,axis=0)

    extents=RefMax-RefMin #extents
    rl=factor*(np.amin(extents[0:2])) #linear 'scale' to set up interactor
    return [RefMin[0]-rl, \
      RefMax[0]+rl, \
      RefMin[1]-rl, \
      RefMax[1]+rl, \
      RefMin[2],RefMax[2]]

def get_axis(renderer, limits, scale, z = False):#(ren,limits,scale):
    '''
    Adds returns a cubeaxesactor based on limits passed with blank Z axis
    '''
    ax3D = vtk.vtkCubeAxesActor()
    if z:
        ax3D.ZAxisTickVisibilityOn()
        ax3D.SetZTitle('Z')
        ax3D.ZAxisMinorTickVisibilityOn()
        ax3D.ZAxisLabelVisibilityOn()
        ax3D.ZAxisTickVisibilityOn()
    else:
        ax3D.ZAxisMinorTickVisibilityOff()
        ax3D.ZAxisLabelVisibilityOff()
        ax3D.ZAxisTickVisibilityOff()
    ax3D.SetXTitle('X')
    ax3D.SetYTitle('Y')
    ax3D.SetBounds(limits)
    ax3D.SetZAxisRange(limits[-2]*scale,limits[-1]*scale)
    ax3D.SetCamera(renderer.GetActiveCamera())

    return ax3D

def gen_cutting_orientation_actor(limits,cut_dir,cut_path=None):
    '''
    Returns a vtk actor corresponding to cutting directions and limits
    '''
    #generate an arrow in the bottom right of the bounding box
    arrow_source=vtk.vtkArrowSource()
    arrow_source.SetShaftRadius(0.12)
    arrow_source.SetTipResolution(25)
    arrow_source.SetShaftResolution(25)
    arrow_source.SetTipRadius(0.35)
    arrow_source.SetTipLength(0.7/2)
    length=np.maximum(limits[1]-limits[0],limits[3]-limits[2])*0.05
    start = np.array([limits[0], limits[2], 0])
    
    end = start + (length*cut_dir)
    norm_x =(end - start)/length

    arbitrary=np.array([1,1,1]) #can be replaced with a random vector
    norm_z=np.cross(norm_x,arbitrary/np.linalg.norm(arbitrary))
    norm_y=np.cross(norm_z,norm_x)

    # Create the direction cosine matrix by writing values directly to an identity matrix
    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    for i in range(3):
        matrix.SetElement(i, 0, norm_x[i])
        matrix.SetElement(i, 1, norm_y[i])
        matrix.SetElement(i, 2, norm_z[i])

    #Apply transforms
    transform = vtk.vtkTransform()
    transform.Translate(start)
    transform.Concatenate(matrix)
    transform.Scale(length*2, length, length)
    
    # Transform the polydata
    transform_pd = vtk.vtkTransformPolyDataFilter()
    transform_pd.SetTransform(transform)
    transform_pd.SetInputConnection(arrow_source.GetOutputPort())
    transform_pd.Update()
    
    #Create mapper
    mapper = vtk.vtkDataSetMapper()
    
    #generate linesource if there's a cut_path
    if cut_path is not None:
        line_source = vtk.vtkLineSource()
        line_source.SetPoint1(start)
        #point 2 will be a projection along cut_path, length away
        line_source.SetPoint2(start + (length * (-cut_path)))
        line_source.Update()
        tube_filter = vtk.vtkTubeFilter()
        tube_filter.SetInputConnection(line_source.GetOutputPort())
        tube_filter.SetRadius(0.12*length)
        tube_filter.SetNumberOfSides(25)
        tube_filter.Update()
        
        #create append filter to merge them together
        append_filter = vtk.vtkAppendFilter()
        append_filter.AddInputData(tube_filter.GetOutput())
        append_filter.AddInputData(transform_pd.GetOutput())#GetOutputPort())
        append_filter.Update()
        mapper.SetInputConnection(append_filter.GetOutputPort())
    else:
        mapper.SetInputConnection(transform_pd.GetOutputPort())
    
    cut_orient_actor = vtk.vtkActor()
    cut_orient_actor.SetMapper(mapper)
    colors = vtk.vtkNamedColors()
    cut_orient_actor.GetProperty().SetColor(colors.GetColor3d("Violet"))
    return cut_orient_actor

def flip_visible(actor):
    '''
    Convenience function for changing the visibility of actors
    '''
    if actor.GetVisibility():
        actor.VisibilityOff()
    else:
        actor.VisibilityOn()

def get_diverging_lut(c = 'blues'):
    '''
    builds lut from named colors: https://htmlpreview.github.io/?https://github.com/Kitware/vtk-examples/blob/gh-pages/VTKNamedColorPatches.html
    '''
    colors = vtkNamedColors()
    # Colour transfer function.
    ctf = vtk.vtkColorTransferFunction()
    ctf.SetColorSpaceToDiverging()
    if c == 'blues':
        p1 = [0.0] + list(colors.GetColor3d('MidnightBlue'))
        p2 = [0.5] + list(colors.GetColor3d('Gainsboro'))
        p3 = [1.0] + list(colors.GetColor3d('cobalt_green'))
    elif c == 'reds':
        p1 = [0.0] + list(colors.GetColor3d('cobalt_violet_deep'))
        p2 = [0.5] + list(colors.GetColor3d('wheat'))
        p3 = [1.0] + list(colors.GetColor3d('cadmium_yellow_light'))
    ctf.AddRGBPoint(*p1)
    ctf.AddRGBPoint(*p2)
    ctf.AddRGBPoint(*p3)

    table_size = 256
    lut = vtk.vtkLookupTable()
    lut.SetNumberOfTableValues(table_size)
    lut.Build()

    for i in range(0, table_size):
        rgba = list(ctf.GetColor(float(i) / table_size))
        rgba.append(1)
        lut.SetTableValue(i, rgba)

    return lut
    
def xyview(ren):
    camera = ren.GetActiveCamera()
    camera.SetPosition(0,0,1)
    camera.SetFocalPoint(0,0,0)
    camera.SetViewUp(0,1,0)
    ren.ResetCamera()

def yzview(ren):
    camera = ren.GetActiveCamera()
    camera.SetPosition(1,0,0)
    camera.SetFocalPoint(0,0,0)
    camera.SetViewUp(0,0,1)
    ren.ResetCamera()

def xzview(ren):
    vtk.vtkObject.GlobalWarningDisplayOff() #mapping from '3' triggers an underlying stereoview that most displays do not support for trackball interactors
    camera = ren.GetActiveCamera()
    camera.SetPosition(0,1,0)
    camera.SetFocalPoint(0,0,0)
    camera.SetViewUp(0,0,1)
    ren.ResetCamera()
    
def reduce_pnts(pnts, val, z_cutoff = None, mode=0):
    '''
    returns an index of reduced points based on params mode and val:
    mode 0 - val is the number of points to retain (default)
    mode 1 - val is a percentage of points to retain
    mode 2 - val is a z value that the third column of val has to be greater than in order to keep
    '''
    localind = np.arange(0, len(pnts), 1, dtype=int)
    #re-index based on z_cutoff
    if z_cutoff is not None:
        localind = localind[pnts[:,-1] > z_cutoff]
    
    #re-index based on either percentage or number of remaining points above the cutoff
    if mode == 0:
        ind = localind[np.arange(0,len(localind),int(len(localind)/val),dtype=int)]
    elif mode == 1:
        red = float(val)/100
        ind = localind[np.arange(0,len(localind),int(1/red),dtype=int)]
    
    return ind

def draw_arrow(start,length,direction,renderer,invert = True):
    """
    Draws and scales an arrow with a defined starting point, direction and length, adds to the renderer, returns the actor
    """
    arrow_source=vtk.vtkArrowSource()
    arrow_source.SetShaftRadius(0.12)
    arrow_source.SetTipResolution(12)
    arrow_source.SetTipRadius(0.35)
    arrow_source.SetTipLength(0.7/2)
    arrow_source.SetTipResolution(25)
    arrow_source.SetShaftResolution(25)
    
    if invert:
        arrow_source.InvertOn()
    else:
        arrow_source.InvertOff()
    end = start + (length*direction)
    norm_x =(end - start)/length

    arbitrary=np.array([1,1,1]) #can be replaced with a random vector
    norm_z=np.cross(norm_x,arbitrary/np.linalg.norm(arbitrary))
    norm_y=np.cross(norm_z,norm_x)
    
    # Create the direction cosine matrix by writing values directly to an identity matrix
    matrix = vtk.vtkMatrix4x4()
    matrix.Identity()
    for i in range(3):
        matrix.SetElement(i, 0, norm_x[i])
        matrix.SetElement(i, 1, norm_y[i])
        matrix.SetElement(i, 2, norm_z[i])
        
    #Apply transforms
    transform = vtk.vtkTransform()
    transform.Translate(start)
    transform.Concatenate(matrix)
    transform.Scale(length*2, length, length)
 
    # Transform the polydata
    transform_pd = vtk.vtkTransformPolyDataFilter()
    transform_pd.SetTransform(transform)
    transform_pd.SetInputConnection(arrow_source.GetOutputPort())
    transform_pd.Update()
    arrow_polydata = transform_pd.GetOutput()
    
    #Create mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(arrow_polydata)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    if renderer is not None:
        renderer.AddActor(actor)
    return actor, arrow_polydata
    
def do_transform(points, T):
    '''
    Applies 4x4 transformation matrix to points and returns the result
    @Param - points, Nx3 matrix of points; T - 4x4 homologous matrix
    '''
    X = points.copy()
    X = X.transpose()
    X = np.append(X, np.ones((1, X.shape[1])), axis=0) #pad with 1's
    X = T @ X #apply by matrix multiplication
    return X[0:3].transpose() #return an Nx3

def respace_equally(X,input):
    '''
    Takes X an array of points, respaces them on the basis of input, either a floating point value of what the target interval between points is, or an integer which is the total number of points. Returns the new array of points, the perimeter and the number of points.
    '''
    distance=np.sqrt(np.sum(np.diff(X,axis=0)**2,axis=1))
    s=np.insert(np.cumsum(distance),0,0)
    Perimeter=np.sum(distance)

    if not isinstance(input,(int)):
        nPts=round(Perimeter/input)
    else:
        nPts=input
    
    sNew=np.linspace(0,s[-1],nPts)
    fx = interp1d(s,X[:,0])
    fy = interp1d(s,X[:,1])
    
    Xnew=fx(sNew)
    Ynew=fy(sNew)
    
    X_new=np.stack((Xnew,Ynew),axis=-1)
    return X_new,Perimeter,nPts

def get_corner_ind(X):
    '''
    Returns the index of corner points of 3D point array X
    Performs a check and re-orders such that the returned index is ccw, with the first point occuring in the bottom left corner of the outline.
    '''

    #calculate limits
    limits = get_limits(X, 0.01)
    
    #Calculate 2D corners, starting with finding the bottom left corner by minimizing euclidean distance between bounding box corner and outline
    d = np.sqrt((X[:,0]-limits[0])**2+(X[:,1]-limits[2])**2)
    ind=np.where(d==np.amin(d))[0][0] #index of value that corresponds

    #reorder the points so that ind is first
    X=np.vstack((X[ind::,:],X[0:ind+1,:]))

    c_target=np.array([
    [limits[0],limits[3]], #xmin,ymax
    [limits[1],limits[3]], #xmax,ymax
    [limits[1],limits[2]] #xmax,ymin
    ])
    ind=np.array([])
    for i in c_target:
        d=np.array([])
        for j in range(len(X[:,0])):
            d=np.append(d,
            np.sqrt((i[0]-X[j,0])**2+(i[1]-X[j,1])**2)
                )
        ind=np.append(ind,np.where(d==np.amin(d))[0][0])
    corner_ind = np.sort(np.append(ind,0)).astype(int)
    corners = X[corner_ind,:]
    
    #ensure ccw ordering with orientation matrix of triangle 0,1,2
    O = np.ones((3,3))
    O[:,-2] = corners[:3,0]
    O[:,-1] = corners[:3,1]
    det = np.linalg.det(O)
    
    if det > 0:
        # c3 = corners[3,:].copy()
        # c1 = corners[1,:].copy()
        i1 = corner_ind[1]
        i3 = corner_ind[3]
        # corners[1,:] = c3
        # corners[3,:] = c1
        corner_ind[1] = i3
        corner_ind[3] = i1
        
    return corner_ind, X

def offset_poly(poly, offset):
    '''
    Uses pyclipper to offset param poly and return an offset polygon according to offset param. poly can be 3D (XYZ by N), but will return a 2D (XY by N)
    '''
    
    X = tuple(map(tuple, poly[:,:2]))
    scaled_poly = scale_to_clipper(X)
    pco = PyclipperOffset()
    pco.AddPath(scaled_poly, JT_SQUARE, ET_CLOSEDPOLYGON)
    scaled_result = pco.Execute(scale_to_clipper(offset))
    
    result = scale_from_clipper(scaled_result)
    
    return np.asarray(result[0])

def in_poly(poly,pnts):
    '''
    Convenience function for matplotlib's contains_points, returning a boolean array corresponding to pnts being inside poly as True.
    '''
    p = path.Path(poly[:,:2])
    return p.contains_points(pnts[:,:2])

def gen_scalar_bar(title = None, num_contours = 13, side = 'left'):
    '''
    Returns a formatted scalebar widget based on the incoming lookup table, title and number of labels to the left or right of the interacator depending on 'side'
    '''
    bar_widget = vtk.vtkScalarBarWidget()
    scalarBarRep = bar_widget.GetRepresentation()
    if side == 'left':
        scalarBarRep.GetPositionCoordinate().SetValue(0.005,0.01) #bottom left
        scalarBarRep.GetPosition2Coordinate().SetValue(0.095,0.98) #top right
    elif side == 'right':
        scalarBarRep.GetPositionCoordinate().SetValue(0.903,0.01) #bottom left
        scalarBarRep.GetPosition2Coordinate().SetValue(0.095,0.98) #top right
    sb_actor=bar_widget.GetScalarBarActor()

    sb_actor.SetTitle(title)
    sb_actor.SetNumberOfLabels(num_contours)

    #attempt to change scalebar properties
    sb_actor.GetLabelTextProperty().SetColor(1,1,1)
    sb_actor.GetTitleTextProperty().SetColor(1,1,1)
    sb_actor.GetLabelTextProperty().SetFontSize(1)
    sb_actor.GetTitleTextProperty().SetFontSize(1)
    sb_actor.SetLabelFormat("%.4f")

    return bar_widget

def gen_line_actor(p1, p2, res = 50, color = None, radius = 0.5):
    '''
    Returns a line running from p1 to p2, glyphed for visibility
    '''

    line_source = vtk.vtkLineSource()
    line_source.SetResolution(res)
    line_source.SetPoint1(p1)
    line_source.SetPoint2(p2)
    line_source.Update()
    
    tube_filter = vtk.vtkTubeFilter()
    tube_filter.SetInputConnection(line_source.GetOutputPort())
    tube_filter.SetRadius(radius)
    tube_filter.SetNumberOfSides(res)
    tube_filter.Update()
    
    ph_res = res
    th_res = res
    
    sphere1 = vtk.vtkSphereSource()
    sphere1.SetCenter(p1)
    sphere1.SetPhiResolution(ph_res)
    sphere1.SetThetaResolution(th_res)
    sphere1.SetRadius(radius * 2)
    sphere1.Update()
    
    sphere2 = vtk.vtkSphereSource()
    sphere2.SetCenter(p2)
    sphere2.SetPhiResolution(ph_res)
    sphere2.SetThetaResolution(th_res)
    sphere2.SetRadius(radius * 2)
    sphere2.Update()
    
    appendFilter = vtk.vtkAppendPolyData()
    appendFilter.AddInputData(sphere1.GetOutput())
    appendFilter.AddInputData(tube_filter.GetOutput())
    appendFilter.AddInputData(sphere2.GetOutput())
    appendFilter.Update()
    
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(appendFilter.GetOutput())
    line_actor = vtk.vtkActor()
    line_actor.SetMapper(mapper)

    if color == None:
        line_actor.GetProperty().SetColor(vtk.vtkNamedColors().GetColor3d("Violet"))
    else:
        line_actor.GetProperty().SetColor(color)
    return line_actor



def gen_info_actor(message,ren, color=(0,0,0)):
    '''
    Returns an information actor comprised of the incoming message string positioned correctly according to the incoming renderer
    '''
    
    textmapper = vtk.vtkTextMapper()
    textmapper.SetInput(message)
    textProperty = vtk.vtkTextProperty()
    textProperty.SetFontSize(16)
    textProperty.SetJustificationToCentered()
    textProperty.SetColor(color)
    # textProperty.SetColor(vtk.vtkNamedColors().GetColor3d('tomato'))
    textmapper.SetTextProperty(textProperty)
    info_actor = vtk.vtkActor2D()
    info_actor.SetMapper(textmapper)
    #get size of renderwindow
    size = ren.GetSize() #(width,height)
    info_actor.SetPosition(int(0.5*size[0]), int(0.001*size[1]))
    return info_actor

def gen_caption_actor(message, actor = None, color = (0,0,0)):
    '''
    Captions an actor
    '''
    caption_actor = vtk.vtkCaptionActor2D()
    b = actor.GetBounds()
    caption_actor.SetAttachmentPoint((b[0],b[2],b[4]))
    caption_actor.SetCaption(message)
    caption_actor.SetThreeDimensionalLeader(False)
    caption_actor.BorderOff()
    caption_actor.LeaderOff()
    caption_actor.SetWidth(0.25 / 3.0)
    caption_actor.SetHeight(0.10 / 3.0)
    
    p = caption_actor.GetCaptionTextProperty()
    p.SetColor(color)
    p.BoldOn()
    p.ItalicOff()
    p.SetFontSize(16)
    p.ShadowOn()
    return caption_actor
    
def identify_actor(actor, list_of_actors):
    """
    Returns the index of actor in list_of_actors by casting address as a 16-bit integer
    """
    
    ai = int(actor.GetAddressAsString('vtkPolyData')[5:], 16)
    for a in range(len(list_of_actors)):
        if ai == int(list_of_actors[a].GetAddressAsString('vtkPolyData')[5:], 16):
            break
    return a

def initialize_HDF5(file=None):
    '''
    Create an HDF5 file with relevant empty structures.
    '''
    if file is None:
        #launch get_save_file
        file, _ = get_save_file('*.pyCM')
    
    if file is None:
            return

    with h5py.File(file, 'w') as f:
        f.attrs['date_created'] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        
        f.attrs['date_modified'] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    
    return file


class vtkug_writer(VTKPythonAlgorithmBase):
    '''
    Adapted Berk Geveci's class to write VTK unstructured grids to HDF5: from 'Developing HDF5 readers using vtkPythonAlgorithm', 2014 
    https://blog.kitware.com/hdf5-reader-and-writer-for-unstructured-grids/
    '''
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, \
        nInputPorts=1, \
        inputType='vtkUnstructuredGrid', \
        nOutputPorts=0)
 
        self.__FileName = ""
        self.__NumberOfPieces = 1
        self.__CurrentPiece = 0


    def RequestData(self, request, inInfo, outInfo):
        info = inInfo[0].GetInformationObject(0)
        inp = dsa.WrapDataObject(vtk.vtkDataSet.GetData(info))
 
        if self.__CurrentPiece == 0:
            self.__File = h5py.File(self.__FileName, 'r+')
            if "mesh" in self.__File:
                del self.__File["mesh"]
        
        model = self.__File.create_group("mesh")
        
        model.attrs['bounds'] = inp.GetBounds()
 
        model.create_dataset("cells", data=inp.Cells)
        model.create_dataset("cell_types", data=inp.CellTypes)
        model.create_dataset("cell_locations", data=inp.CellLocations)
        model.create_dataset("points", data=inp.Points)
 
        pdata = model.create_group("point_data")
        for name in inp.PointData.keys():
            pdata.create_dataset(name, data=inp.PointData[name])
        
        self.__File.attrs['date_modified'] = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
        
        if self.__CurrentPiece < self.__NumberOfPieces - 1:
            # If we are not done, ask the pipeline to re-execute us.
            self.__CurrentPiece += 1
            request.Set(
                vtk.vtkStreamingDemandDrivenPipeline.CONTINUE_EXECUTING(),
                1)
        else:
            # Stop execution
            request.Remove(
                vtk.vtkStreamingDemandDrivenPipeline.CONTINUE_EXECUTING())
            self.__File.close()
            del self.__File
        return 1
 
    def RequestInformation(self, request, inInfo, outInfo):
        # Reset values.
        self.__CurrentPiece = 0
        return 1
 
    def RequestUpdateExtent(self, request, inInfo, outInfo):
        info = inInfo[0].GetInformationObject(0)
        info.Set(
            vtk.vtkStreamingDemandDrivenPipeline.UPDATE_NUMBER_OF_PIECES(),
            self.__NumberOfPieces)
        info.Set(
            vtk.vtkStreamingDemandDrivenPipeline.UPDATE_PIECE_NUMBER(),
            self.__CurrentPiece)
        return 1
 
    def SetFileName(self, fname):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname
 
    def GetFileName(self):
        return self.__FileName
 
    def SetNumberOfPieces(self, npieces):
        if npieces != self.__NumberOfPieces:
            self.Modified()
            self.__NumberOfPieces = npieces
 
    def GetNumberOfPieces(self):
        return self.__NumberOfPieces


class vtkug_reader(VTKPythonAlgorithmBase):
    '''
    Companion reader to vtkug_writer. Returns an unstructured grid on block 0 of a multiblock dataset, i.e. vtkug_reader.GetOutputDataObject(0).GetBlock(0). Adapted from Berk Geveci's class 'HDF5 Reader and Writer for Unstructured Grids', 2015.
    '''
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
            nInputPorts=0,
            nOutputPorts=1, outputType='vtkMultiBlockDataSet')
 
        self.__FileName = ""
 
    def RequestData(self, request, inInfo, outInfo):
        output = dsa.WrapDataObject(vtk.vtkMultiBlockDataSet.GetData(outInfo))
        f = h5py.File(self.__FileName, 'r')
        idx = 0
        grp = f["mesh"]
        ug = vtk.vtkUnstructuredGrid()
        output.SetBlock(idx, ug)
        idx += 1
        ug = dsa.WrapDataObject(ug)
        cells = grp['cells'][:]
        locations = grp['cell_locations'][:]
        types = grp['cell_types'][:]
        ug.SetCells(types, locations, cells)
        pts = grp['points'][:]
        ug.Points = pts
        pt_arrays = grp['point_data']
        for pt_array in pt_arrays:
            array = pt_arrays[pt_array][:]
            ug.PointData.append(array, pt_array)
 
        return 1
 
    def SetFileName(self, fname):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname
 
    def GetFileName(self):
        return self.__FileName

def warning_msg(parent, message):
    '''
    Generates a warning box to the parent qt widget with the message, returns true/false for yes/no
    '''
    ret = QtWidgets.QMessageBox.warning(parent, "pyCM Warning", \
    message, \
    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No, \
    QtWidgets.QMessageBox.No)
    if ret == QtWidgets.QMessageBox.No: #don't overwrite
        return False
    else: return True

def info_msg(message):
    '''
    Generates an info box
    '''
    msg = QtWidgets.QMessageBox()
    msg.setIcon(QtWidgets.QMessageBox.Information)
    msg.setText(message)
    msg.setWindowTitle("pyCM Error")
    msg.exec_()
    
def write_dxf(file, outline):
    '''
    writes a closed dxf polygon based on points from outline to file
    '''
    n = len(outline)
    fid = io.StringIO()

    fid.write("0\nSECTION\n2\nENTITIES\n0\n")
    for i in range(n-1):
        outputline = "LINE\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n0\n" \
                %(outline[i,0],outline[i,1],0,outline[i+1,0],outline[i+1,1],0)
        fid.write(outputline)
    #last line
    fid.write("LINE\n10\n%f\n20\n%f\n30\n%f\n11\n%f\n21\n%f\n31\n%f\n0\n" \
                %(outline[-1,0],outline[-1,1],0,outline[0,0],outline[0,1],0))
    fid.write("ENDSEC\n0\nEOF\n")
    
    with open(file, 'w+') as f: f.write(fid.getvalue())
    
    fid.close()

def gen_filtered_ugrid(ugrid):
    '''
    Eliminates low-order elements from an input unstructured mesh. Returns an unstructured grid with only volumetric cells (elements). Retains original no=de count and numbering.
    '''

    cell_type = vtk.vtkCellTypes()
    ugrid.GetCellTypes(cell_type)
    
    if cell_type.IsType(24) == 1: #if it contains 2nd order tets
        vtk_type = 24
        n_nodes_per_element = 10
    elif cell_type.IsType(12) == 1: #if it contains linear quads
        vtk_type = 12
        n_nodes_per_element = 8

    cells  = ugrid.GetCells()
    raw_point_ids = v2n.vtk_to_numpy(cells.GetData()) #1d list of cell_type & connectivity

    cell_offsets = v2n.vtk_to_numpy(cells.GetOffsetsArray()) #1d list of offsets of ids
    cell_connectivity = []

    ind = 0 #index of raw_points
    for i in range(1,len(cell_offsets)):
        local = raw_point_ids[ind:cell_offsets[i]+i]
        if len(local) == n_nodes_per_element+1: #only get high order (volumetric) elements
            cell_connectivity.append(local)
        ind = cell_offsets[i]+i

    #cell_points now contains <node number, cell_node_1 . . . cell_node_n>
    cell_connectivity = np.array(cell_connectivity)
    #remove node number from cell_points
    cell_connectivity = cell_connectivity[:,1::]
    
    nodes = ugrid.GetPoints()
    new_cells = vtk.vtkCellArray()

    new_ugrid = vtk.vtkUnstructuredGrid()
    for element in cell_connectivity:
        new_ugrid.InsertNextCell(vtk_type, n_nodes_per_element, element)
    new_ugrid.SetPoints(nodes)
    
    return new_ugrid

def convert_inp_to_vtk(infile,outfile):
    """
    Converts abaqus inp file into a legacy ASCII vtk file. First order quads (C3D8) and third order tets (C3D10) are supported.
    """
    fid = open(infile)
    
    #flags for identifying sections of the inp file
    inpKeywords=["*Node", "*Element", "*Nset", "*Elset", "*NODE", "*ELEMENT", "*NSET", "*ELSET"]
    
    #map abaqus mesh types to vtk objects
    vtkType={'C3D8': 12, 'C3D10': 24}

    #create counter for all lines in the inp file, and array to store their location
    i=0
    lineFlag=[];
    #read file and find both where keywords occur as well as the element type used
    while 1:
        lines = fid.readlines(100000)
        if not lines:
            break
        for line in lines:
            i+=1
            for keyword in inpKeywords:
                if line[0:len(keyword)]==keyword:
                    lineFlag.append(i)
                    if keyword=="*Element" or keyword=="*ELEMENT":
                        line = line.replace("\n", "")
                        CellNum=vtkType[line.split("=")[-1]]

    fid.close()
    #use genfromtxt to read between lines id'ed by lineFlag to pull in nodes and elements
    Nodes=np.genfromtxt(infile,skip_header=lineFlag[0],skip_footer=i-lineFlag[1]+1,delimiter=",")
    Elements=np.genfromtxt(infile,skip_header=lineFlag[1],skip_footer=i-lineFlag[2]+1,delimiter=",")
    #Now write it in VTK format to a new file starting with header

    fid=open(outfile,'wb+')
    fid.write(str.encode('# vtk DataFile Version 2.0\n'))
    fid.write(str.encode('%s,created by pyCM\n'%outfile[:-4]))
    fid.write(str.encode('ASCII\n'))
    fid.write(str.encode('DATASET UNSTRUCTURED_GRID\n'))
    fid.write(str.encode('POINTS %i double\n'%len(Nodes)))
    
    #dump nodes
    np.savetxt(fid,Nodes[:,1::],fmt='%.6f')
    fid.write(str.encode('\n'))
    fid.write(str.encode('CELLS %i %i\n'%(len(Elements),len(Elements)*len(Elements[0,:]))))
    #Now elements, stack the number of nodes in the element instead of the element number
    Cells=np.hstack((np.ones([len(Elements[:,0]),1])*len(Elements[0,1::]),Elements[:,1::]-1))
    np.savetxt(fid,Cells,fmt='%i')
    fid.write(str.encode('\n'))

    #Write cell types
    fid.write(str.encode('CELL_TYPES %i\n'%len(Elements)))
    CellType=np.ones([len(Elements[:,0]),1])*CellNum
    np.savetxt(fid,CellType,fmt='%i')

    fid.close()

def read_file_for_fea(file):
    '''
    Reads a pyCM h5 file for data pertinent for input file and viewing purposes
    '''
    outlines = []
    bc_disp_nodes = []
    bc_disp = []
    bc_nodes = []
    trans = None
    mod = None
    pr = None
    mesh = None
    read_mesh = False
    
    with h5py.File(file, 'r') as f:
        try:
            g = f['bc_prop']
            #get number of entries
            n = sum(k.isdigit() for k in g.keys())
            outlines = [None] * n
            bc_disp = [None] * n
            bc_disp_nodes = [None] * n
            for k in g.keys():
                if k.isdigit():
                    outlines[int(k)] = g['%s/outline'%k][()]
                    if '%s/bc_disp_val'%k in g:
                        bc_disp[int(k)] = g['%s/bc_disp_val'%k][()]
                        bc_disp_nodes[int(k)] = g['%s/surface_nodes'%k][()]
            if 'rigid_body_nodes' in list(g.keys()):
                bc_nodes = g['rigid_body_nodes'][()]
            if 'transform' in list(g.keys()):
                trans = g['transform'][()]
            mod = g.attrs['modulus']
            pr = g.attrs['poisson_ratio']
        except: pass
        
        if 'mesh' in f.keys():
            read_mesh = True
    
    if read_mesh:
        try:
            r = vtkug_reader()
            r.SetFileName(file)
            r.Update()
            mesh = r.GetOutputDataObject(0).GetBlock(0)
        except:
            pass
    
    return outlines, bc_disp_nodes, bc_disp, bc_nodes, trans, mod, pr, mesh