# postprocess

## Background
Tool for both post-processing the output from the linear elastic Finite Element Analysis of the contour method. Visualization is carried out with VTK and currently can support C3D8 element from an ABAQUS output (.dat) file.

## Initializing

**Input and output descriptors**

The input consists of a *.dat, *.inp and *.vtk files. All of them will be produced by the pre-processor. It is important to note that when selective the *.inp file the user must chose the intermediate input file containing *only the mesh* (not *.abq.inp!).

An analysis will be carried out at each input stage and it is normal to take a few seconds (especially on the *.inp step). Following the extrapolation from quadrature points to nodal points the data will be appended in the selected *.vtk file. This can then be visualized in the Qt window.

The tool is called from Python according to:
~~~
from pyCM import postprocess as pc
pc.post_process_tool()
~~~

##  Interaction functionality
The button "Extract S33" will start the process of asking the user for the relevant files. The button "Display S33" can then be used to visualize the *.vtk file.

<span>![<span>Main Window</span>](images/postprocess1.png)</span>
*<a name="fig1"></a> Figure 1: Postprocess*

## Known issues
- The vtk display seems to cut parts of the piece when rotating