'''
Abaqus Python script for extracting stresses at nodes. Output: ascii based *.vtu Consisting of deformed mesh and scalar fields of S11, S22, S33. Input odb and output file are arguments.
(c) M. J. Roy 2021
'''
import sys
import os.path
import numpy as np
from odbAccess import *


def postprocess(odbname, outfile, apply_displacement = False):
    if not os.path.isfile(odbname):
        sys.exit("Specified odb file to odb_access not valid.")


    odb = openOdb(odbname,readOnly=True)

    #access geometry and topology information ( odb->rootAssembly->instances->(nodes, elements) )
    rootassembly = odb.rootAssembly
    instance = rootassembly.instances

    #access attribute information
    step = odb.steps

    allinstancestr = str(instance)
    autoins = allinstancestr.split("'")
    instancename = autoins[1]
    node = instance[instancename].nodes
    element = instance[instancename].elements

    n_nodes = len(node)
    n_elements = len(element)

    allstepstr = str(step)
    autostep = allstepstr.split("'")
    stepname = autostep[1] #only 'named' step

    N_Frame = odb.steps[stepname].frames[-1] #last frame

    #create np array containing nodeLabel and baseline x, y, z
    node_array = np.zeros([n_nodes,4])
    #node numbering may not be continuous
    i = 0
    for n in node:
        node_array[i,:]=[n.label,n.coordinates[0],n.coordinates[1],n.coordinates[2]]
        i += 1

    #get displacements and add to base coordinates for deformed locations
    if apply_displacement:
        Displacement = N_Frame.fieldOutputs['U']
        i = 0
        for value in Displacement.values:
            node_array[i,1::] = node_array[value.nodeLabel-1,1::] + value.data
            i += 1

    #overwrite node labels but save actual node numbers
    actual_node_numbers = node_array[:,0]
    new_node_numbers = np.array(range(0,n_nodes))

    #get element data
    element_array = np.zeros([n_elements,len(element[0].connectivity)+1])
    j = 0
    for e in element:
        con = np.asarray([i for i in e.connectivity])
        #map new connectivity
        new_con = np.zeros(len(con))
        for jj in range(len(con)):
            new_con[jj] = new_node_numbers[np.where(actual_node_numbers == con[jj])]
        element_array[j,0] = e.label
        element_array[j,1::] = new_con #vtk node numbering starts at 0, not 1
        j += 1

    #get number of nodes per element based on the connectivity of the first element
    num_nodes_per_element = len(element[0].connectivity)

    vtk_element_array = element_array
    vtk_element_array[:,0] = num_nodes_per_element

    if num_nodes_per_element == 10:
        vtk_cell_type = 24
    elif num_nodes_per_element == 8:
        vtk_cell_type = 12

    Stress = N_Frame.fieldOutputs['S']
    node_Stress = Stress.getSubset(position=ELEMENT_NODAL)
    fieldValues = node_Stress.values

    #read stresses at nodes on elemental basis
    stress_array = np.zeros([len(node_array),8]) #node, Nxcount,culmulative S11,S22, S33
    stress_array[:,0] = actual_node_numbers
    i = 0
    for entry in fieldValues:
        if i % 10000 == 0:
            print 'working on stress entry ', entry.nodeLabel, ' : %d of %d'%(i,len(fieldValues))
        #find index
        ind = np.where(actual_node_numbers == entry.nodeLabel)
        stress_array[ind,1]+=1
        stress_array[ind,2]+=entry.data[0]
        stress_array[ind,3]+=entry.data[1]
        stress_array[ind,4]+=entry.data[2]
        # stress_array[ind,5]+=entry.data[3]
        # stress_array[ind,6]+=entry.data[4]
        # stress_array[ind,7]+=entry.data[5]
        i+=1

    #average stresses at nodes
    stress_avg_array = np.zeros([n_nodes,6])
    for i in range(6):
        stress_avg_array[:,i] = stress_array[:,i+2]/stress_array[:,1]

    labels = ['11', '22', '33', '12', '13', '23']

    fid=open(outfile,'w+')

    #write header
    fid.write('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">\n<UnstructuredGrid><Piece NumberOfPoints="%i" NumberOfCells="%i">\n'%(n_nodes,n_elements))

    #now point data (stresses at nodes)    
    fid.write('<PointData>\n')
    for i in range(6):
        fid.write('<DataArray Name="S%s" type="Float64" format="ASCII">\n'%(labels[i]))
        np.savetxt(fid,stress_avg_array[:,i])
        fid.write('</DataArray>\n')
    fid.write('</PointData>\n')

    #write node coordinates
    fid.write('<Points>\n<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii" RangeMin="%f" RangeMax="%f">\n'%(np.min(node_array[:,1::]),np.max(node_array[:,1::])))
    np.savetxt(fid,node_array[:,1::],fmt='%.6f')
    fid.write('</DataArray>\n</Points>')

    #write elements
    fid.write('<Cells>\n <DataArray type="Int64" Name="connectivity" format="ascii" RangeMin="%i" RangeMax="%i">\n'%(0,n_nodes+1))
    np.savetxt(fid,vtk_element_array[:,1::],fmt='%i')

    #write offsets
    fid.write('</DataArray>\n<DataArray type="Int64" Name="offsets" format="ascii" RangeMin="%i" RangeMax="%i">\n'%(num_nodes_per_element,num_nodes_per_element*n_elements))
    offsets = np.asarray([i * num_nodes_per_element for i in range(n_elements+1)])
    np.savetxt(fid,offsets[1::],fmt='%i')

    #write cell types
    fid.write('</DataArray>\n<DataArray type="UInt8" Name="types" format="ascii" RangeMin="%i" RangeMax="%i">\n'%(vtk_cell_type,vtk_cell_type))
    types = np.asarray([vtk_cell_type for i in range(n_elements+1)])
    np.savetxt(fid,types,fmt='%i')

    fid.write('</DataArray>\n</Cells>\n</Piece>\n</UnstructuredGrid>\n</VTKFile>')
    fid.close()

    #close odb file
    odb.close()

if __name__ == "__main__":
    #check if incoming odb file is valid. Use try/catch for valid directory.
    odbname = sys.argv[1]
    outfile = sys.argv[2]
    postprocess(odbname,outfile)