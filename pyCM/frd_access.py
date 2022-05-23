'''
Python 3 script for extracting stresses at nodes of a C3D8 and C3D10 formulation conducted by CalculiX. Output: ascii based *.vtu Consists of a deformed mesh and scalar fields at S11, S22, S33.
(c) M. J. Roy 2021
'''
import sys, io, shutil
import os.path
import numpy as np
from itertools import islice

def postprocess(frdname,outfile, displacement = False):


    if not os.path.isfile(frdname):
        sys.exit("Specified frd file to frd_access not valid.")

    #read file and get line numbers with key values identifying output
    i=0
    lineFlag=[]
    keyStrings=["    2C", "    3C", " -4  DISP", " -4  STRESS", " -3"]

    fid = open(frdname)
    while 1:
        lines = fid.readlines(100000)
        if not lines:
            break
        for line in lines:
            i+=1
            for keyString in keyStrings:
                if line[0:len(keyString)]==keyString:
                    lineFlag.append(i)
                
    fid.close()

    #read nodal data
    a = np.genfromtxt(frdname,
        delimiter = (3,10,12,12,12),
        skip_header=lineFlag[0], max_rows=lineFlag[1]-lineFlag[0]-1)
    node_array = a[:,1::] #node number, x,y,z coords

    n_nodes = np.size(node_array, 0)
    n_elements = int((lineFlag[3]-1 - lineFlag[2])/2) #Could also be read directly from '3C' flag
    
    #get element type on the basis of the length of the connectivity of the first element
    with open(frdname) as fid:
        first_element_line = list(islice(fid,lineFlag[2]+1,lineFlag[2]+2))[0]
        l = np.array([int(i) for i in first_element_line.split()])
    
    num_nodes_per_element = len(l)-1
    
    if num_nodes_per_element == 10:
        vtk_cell_type = 24
    elif num_nodes_per_element == 8:
        vtk_cell_type = 12

    #preallocate numpy array
    element_array = np.zeros([n_elements,num_nodes_per_element+1])
    #read element data using isslice, could be optimised with regex
    with open(frdname) as fid:
        lines = islice(fid,lineFlag[2],lineFlag[3]-1)
        e_num = 0
        for line in lines:
            l = np.array([int(i) for i in line.split()])
            #capture lines that once split have 8 values (C3D8)
            if l[0] == -1: #then element declaration
                e_num = l[1]
                element_array[e_num-1,0] = e_num-1 #cell numbering in vtk is 0 based
            if l[0] == -2 and e_num != 0:#then connectivity of element number declared the line above
                element_array[e_num-1,1::] = l[1::]-1 #node numbering in vtk is 0 based

    vtk_element_array = element_array
    vtk_element_array[:,0] = num_nodes_per_element #quads

    if displacement:
        #read displacement data and add to node_array
        a = np.genfromtxt(frdname,
            delimiter = (3,10,12,12,12),
            skip_header=lineFlag[4]+4, max_rows=lineFlag[5]-1-(lineFlag[4]+4))
            
        for row in a:
            node_array[int(row[1])-1,1::] = node_array[int(row[1])-1,1::] + row[2::]

    #read averaged stresses
    a = np.genfromtxt(frdname,
        delimiter = (3,10,12,12,12,12,12,12),
        skip_header=lineFlag[6]+6, max_rows=lineFlag[7]-1-(lineFlag[6]+6))
    stress_avg_array = a[:,2:5] #SXX/S11, SYY/S22, SZZ/S33

    #write output
    fid = io.StringIO()
    
    fid.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n<UnstructuredGrid><Piece NumberOfPoints="%i" NumberOfCells="%i">\n'%(n_nodes,n_elements))
    fid.write('<PointData>\n')
    #now point data (stresses at nodes)    
    for i in range(3):
        fid.write('<DataArray Name="S%s%s" type="Float64" format="ASCII">\n'%(i+1,i+1))
        np.savetxt(fid,stress_avg_array[:,i])
        fid.write('</DataArray>\n')
    fid.write('</PointData>\n')
    #write node coordinates
    fid.write('<Points>\n<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="%f" RangeMax="%f">\n'%(np.min(node_array[:,1::]),np.max(node_array[:,1::])))
    np.savetxt(fid,node_array[:,1::],fmt='%.6f')
    fid.write('</DataArray>\n</Points>')

    fid.write('<Cells>\n <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="%i" RangeMax="%i">\n'%(0,n_nodes+1))
    np.savetxt(fid,vtk_element_array[:,1::],fmt='%i')
    fid.write('</DataArray>\n<DataArray type="Int64" Name="offsets" format="ascii" RangeMin="%i" RangeMax="%i">\n'%(num_nodes_per_element,num_nodes_per_element*n_elements))
    offsets = np.asarray([i * num_nodes_per_element for i in range(n_elements+1)])
    np.savetxt(fid,offsets[1::],fmt='%i')
    fid.write('</DataArray>\n<DataArray type="UInt8" Name="types" format="ascii" RangeMin="%i" RangeMax="%i">\n'%(vtk_cell_type,vtk_cell_type))
    types = np.asarray([vtk_cell_type for i in range(n_elements+1)])
    np.savetxt(fid,types,fmt='%i')

    fid.write('</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>')
    
    with open(outfile, 'w+') as f:
        fid.seek(0)
        shutil.copyfileobj(fid, f)


if __name__ == "__main__":
    #check if incoming odb file is valid. Use try/catch for valid directory.
    frdname = sys.argv[1]
    outfile = sys.argv[2]
    postprocess(frdname,outfile)