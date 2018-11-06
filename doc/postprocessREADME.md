# postprocess

## Background
Tool for both post-processing the output from the linear elastic Finite Element Analysis of the contour method. Visualization is carried out with VTK and currently can support C3D8 element from an Abaqus output (.abq.dat) or Calculix output (.ccx.dat) file.

## Initializing

**Input and output descriptors**

This application takes information from the results file (.mat) and looks for the corresponding .dat file associated with the analysis. From the mesh already contained within the results file, the values at integration (gauss) points is translated to nodes using the relevant shape functions, with the complete output stored in a VTK unstructured mesh file in XML format, suitable for further postprocessing in [ParaView](https://www.paraview.org/).

The tool is called from Python according to:
~~~
from pyCM import postprocess as pp
pp.post_process_tool()
~~~

##  Interaction functionality
After loading a results file by pressing **l**, pressing the 'Extract' button will start reading the relevant .dat results file. Display, by default, is the longitudinal direction (S33) but can be changed to S11 and S22 (x and y directions). The number of contours on the scalebar can be changed, as well as the minimum and maximum stresses resolved.

<span>![<span>Main Window</span>](images/postprocess1.png)</span>
*<a name="fig1"></a> Figure 1: Postprocess*

**Keyboard and mouse mapping**

Key | Description
---  |---
Left mouse button 	|Rotate about the center of view
Middle mouse button 	|Pan
Right mouse button 	|Zoom/refresh window extents
1 	|View 1, default, looks down z axis onto xy plane
2 	|View 2, default, looks down x axis onto zy plane
3 	|View 3, default, looks down y axis onto zx plane
f | Flip colors from white on dark to dark on white
i | Save output to .png in current working directory
l | load/reload *.mat file to conduct/review/revise this analysis step


## Known issues
None at this time.