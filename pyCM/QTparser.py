""" Parses the .dat file and extrapolates to the nodes in the .inp file """
import os
import numpy as np

# locations for the ABAQUS .dat and .inp files
DAT_FILE = "../data/abaqus_output.dat"
INP_FILE = "../data/Nik_abq.inp"

# lookup keywords
DAT_FILE_LOOKUP_STR = "E L E M E N T   O U T P U T"
INP_FILE_NODE_LOOKUP_STR = "*NODE"
INP_FILE_ELEM_LOOKUP_STR = "*ELEMENT, TYPE=C3D8"
INP_FILE_ELEM_END_LOOKUP_STR = "*ELSET, ELSET=DOMAIN, GENERATE"

# numpy data types for the quadrature point extraction
# the columns are: node id, quadrature id, x coord, y coord, z coord, S33
QP_VAR_TYPE = "int32, float64, float64, float64, float64"

# numpy data types for nodal point extraction
NP_VAR_TYPE = "int32, float64, float64, float64"

# numpy data types for elements extraction
EL_VAR_TYPE = "int32, int32, int32, int32, int32, int32, int32, int32, int32"

def get_quadrature_data(file_name):
    """ Reads the quadrature point coordinates and stress values.
        Returns a numpy array. """
    # initialize
    curr_line = 0
    read_flag = False
    feed_flag = False
    num_steps = 0
    row_start = 0
    row_end = 1

    # locate the start and end line of the quadrature block
    with open(file_name) as dat_file:
        p_lines = dat_file.readlines()
        for line in p_lines:
            num_steps = num_steps + 1
            if line.rstrip() and not feed_flag:
                line = line.strip()
                # find the element output section
                if line.find(DAT_FILE_LOOKUP_STR) >= 0:
                    read_flag = True
                # we can read the quadrature point data
                # but we need to skip 3 lines
                if read_flag:
                    curr_line = curr_line + 1
                    # coordinate and stress data reached
                    if curr_line == 5:
                        feed_flag = True
                        row_start = num_steps - 1
                        curr_line = 0
            else:
                if not line.rstrip():
                    if read_flag is True and curr_line is 0:
                        # count the number of lines at the end of the .dat file
                        row_end = row_end + 1
    dat_file.close()

    # extract quadrature data for
    # node id, quadrature id, x coord, y coord, z coord, S33
    element_data = np.genfromtxt(file_name, skip_header=row_start, skip_footer=11, \
                                usecols=(0, 2, 3, 4, 7), autostrip=True,             \
                                dtype=QP_VAR_TYPE)
    return element_data

def get_node_data(file_name):
    """ Reads the nodal point coordinates.
        Returns a numpy array. """
    # initialize
    curr_line = 0
    node_start = 0
    node_end = 0
    elem_start = 0
    elem_end = 0

    with open(file_name) as inp_file:
        p_lines = inp_file.readlines()
        for line in p_lines:
            curr_line = curr_line + 1
            if line.find(INP_FILE_NODE_LOOKUP_STR) >= 0:
                node_start = curr_line
            if line.find(INP_FILE_ELEM_LOOKUP_STR) >= 0:
                node_end = curr_line - 1
                elem_start = curr_line
            if line.find(INP_FILE_ELEM_END_LOOKUP_STR) >= 0:
                elem_end = curr_line - 1

    node_end = curr_line - node_end
    elem_end = curr_line - elem_end

    inp_file.close()

    # extract nodal point data for
    # node id, x coord, y coord, z coord, S33
    node_data = np.genfromtxt(file_name, skip_header=node_start, skip_footer=node_end, \
                                delimiter=',', dtype=NP_VAR_TYPE)

    element_data = np.genfromtxt(file_name, skip_header=elem_start, skip_footer=elem_end, \
                                delimiter=',', dtype=EL_VAR_TYPE)

    # numpy provides a structured array which is not useful for our purposes
    # we need a 2d array
    node_data = node_data.view().reshape(len(node_data), -1)
    return node_data, element_data

def CreateVTK_file(vtk_file):
    """ Creates a legacy VTK file from an ABAQUS .inp file """
    fid = open(file_name)

    #flags for identifying sections of the inp file
    inp_keywords = ["*Node", "*Element", "*Nset", "*Elset"]

    #map abaqus mesh types to vtk objects
    vtkType = {}
    vtkType['C3D8'] = 12
    vtkType['C3D10'] = 24

    #create counter for all lines in the inp file, and array to store their location
    i=0
    lineFlag=[]

    #read file and find both where keywords occur as well as the element type used
    while 1:
        lines = fid.readlines(100000)
        if not lines:
            break
        for line in lines:
            i+=1
            for keyword in inp_keywords:
                if line[0:len(keyword)] == keyword:
                    lineFlag.append(i)
                    if keyword=="*Element":
                        line = line.replace("\n", "")
                        CellNum=vtkType[line.split("=")[-1]]
    fid.close()
    #use genfromtxt to read between lines id'ed by lineFlag to pull in nodes and elements
    Nodes = np.genfromtxt(file_name, skip_header=lineFlag[0], skip_footer=i-lineFlag[1]+1, delimiter=",")
    Elements = np.genfromtxt(file_name, skip_header=lineFlag[1], skip_footer=i-lineFlag[2]+1, delimiter=",")

    #Now write it in VTK format to a new file starting with header
    fid = open(vtk_file, 'w+')
    fid.write("# vtk DataFile Version 2.0\n")
    fid.write("%s,created by pyCM\n"%vtk_file[:-4])
    fid.write("ASCII\n")
    fid.write("DATASET UNSTRUCTURED_GRID\n")
    fid.write("POINTS %i double\n"%len(Nodes))

    #dump nodes
    np.savetxt(fid, Nodes[:,1::], fmt='%.6f')
    fid.write("\n")
    fid.write("CELLS %i %i\n"%(len(Elements),len(Elements)*len(Elements[0,:])))
    #Now elements, stack the number of nodes in the element instead of the element number
    Cells = np.hstack((np.ones([len(Elements[:,0]), 1]) * len(Elements[0, 1::]), Elements[:,1::] - 1))
    np.savetxt(fid,Cells,fmt='%i')
    fid.write("\n")

    #Write cell types
    fid.write("CELL_TYPES %i\n"%len(Elements))
    CellType = np.ones([len(Elements[:,0]), 1]) * CellNum
    np.savetxt(fid, CellType, fmt='%i')
    fid.close()

def file_path(rel_file_path):
    """ converts from relative to absolute file path """
    script_dir = os.path.dirname(__file__)
    abs_file_path = os.path.join(script_dir, rel_file_path)
    return abs_file_path

def main():
    """ Main function """

    # obtain quadrature data - this seems to hoard memory, however
    # we need the maximum number of elements to match nodes
    quadrature_data = get_quadrature_data(file_path(DAT_FILE))
    #print("Quadrature data")
    #print(quadrature_data)
    a = quadrature_data.shape
    #MaxElem = quadrature_data[:, 0]

    # obtain nodal and element data from .dat
    node_data, element_data = get_node_data(file_path(INP_FILE))
    b = node_data.shape
    c = element_data.shape
    #print("Nodal data")
    #print(node_data)

    #print("Element data")
    #print(element_data)

    # match each element with its corresponding nodes
    match_nodes(node_data, element_data)

    for element in element_data:
        element_id = element[0]
        #node_ids = element[1:]
        #print(element)
        print(element_id)
        #print(node_ids)

    # print vtk file without the table
    #CreateVTK_file(file_path(INP_FILE), "vtk_output.vtk")

def inspect_cubes():
    """draws cubes to show the nodal and quadrature points"""
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # plot quadrature points
    quad_points = np.array([[4.177, 1.443, 1.415],
                            [4.235, 2.666, 1.415],
                            [3.147, 2.718, 1.415],
                            [3.059, 1.443, 1.415],
                            [4.177, 1.443, 2.549],
                            [4.235, 2.666, 2.549],
                            [3.147, 2.718, 2.549],
                            [3.059, 1.443, 2.549]])

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot vertices
    ax.scatter3D(quad_points[:, 0], quad_points[:, 1], quad_points[:, 2])

    # list the sides of polygons of figure
    verts = [[quad_points[0],quad_points[1],quad_points[2],quad_points[3]],
    [quad_points[4],quad_points[5],quad_points[6],quad_points[7]],
    [quad_points[0],quad_points[1],quad_points[5],quad_points[4]],
    [quad_points[2],quad_points[3],quad_points[7],quad_points[6]],
    [quad_points[1],quad_points[2],quad_points[6],quad_points[5]],
    [quad_points[4],quad_points[7],quad_points[3],quad_points[0]],
    [quad_points[2],quad_points[3],quad_points[7],quad_points[6]]]

    # plot sides
    quad_coll = (Poly3DCollection(verts, linewidths=1, \
                                        edgecolors='r', alpha=0.5))
    face_color = mpl.colors.rgb2hex([0.5, 0.5, 1])
    quad_coll.set_facecolor(face_color)
    ax.add_collection3d(quad_coll)

    # plot nodal points
    nodal_points = np.array([[4.569, 1.003, 1.000],
                            [4.650, 3.087, 1.000],
                            [2.784, 3.211, 1.000],
                            [2.614, 0.969, 1.000],
                            [4.569, 1.003, 2.964],
                            [4.650, 3.087, 2.964],
                            [2.784, 3.211, 2.964],
                            [2.614, 0.969, 2.964]])

    # plot vertices
    ax.scatter3D(nodal_points[:, 0], nodal_points[:, 1], nodal_points[:, 2])

    # list the sides of polygons of figure
    verts = [[nodal_points[0],nodal_points[1],nodal_points[2],nodal_points[3]],
    [nodal_points[4],nodal_points[5],nodal_points[6],nodal_points[7]],
    [nodal_points[0],nodal_points[1],nodal_points[5],nodal_points[4]],
    [nodal_points[2],nodal_points[3],nodal_points[7],nodal_points[6]],
    [nodal_points[1],nodal_points[2],nodal_points[6],nodal_points[5]],
    [nodal_points[4],nodal_points[7],nodal_points[3],nodal_points[0]],
    [nodal_points[2],nodal_points[3],nodal_points[7],nodal_points[6]]]

    # plot sides

    node_coll = (Poly3DCollection(verts, linewidths=1, \
                                        edgecolors='b', alpha=0.25))
    face_color = mpl.colors.rgb2hex([0.25, 0.25, 0.5])
    node_coll.set_facecolor(face_color)
    ax.add_collection3d(node_coll)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

    input("Enter")

if __name__ == "__main__":
    #main()
    inspect_cubes()
