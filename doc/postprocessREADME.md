# postprocess

## Background
Tool for both post-processing the output from the linear elastic Finite Element Analysis of the contour method. Visualization is carried out with VTK and currently can support C3D8 element from an Abaqus output (.abq.dat) or Calculix output (.ccx.dat) file.

## Initializing

**Input and output descriptors**

This application takes information from the results file (.mat) and looks for the corresponding .dat file associated with the analysis. From the mesh already contained within the results file, the values at integration (gauss) points are translated to nodes using the relevant shape functions, with the complete output stored in a VTK unstructured mesh file in XML format, suitable for further postprocessing in [ParaView](https://www.paraview.org/).

The tool is called from Python according to:
~~~
from pyCM import postprocess as pp
pp.post_process_tool()
~~~

##  Interaction functionality
After loading a results file by pressing **l** or reloading the .mat file from `main`, the relevant .dat results file will be loaded. Display, by default, is the longitudinal direction (S33) but can be changed to S11 and S22 (x and y directions). The number of contours on the scalebar can be changed, as well as the minimum and maximum stresses resolved.

<span>![<span>Main Window</span>](images/postprocess1.png)</span>
*<a name="fig1"></a> Figure 1: Postprocess*

Stress profiles can also be extracted by inputting x and y locations of the start (x0, y0) and end positions (x1,y1) as well as the interval along the profile. The probe line will be displayed in the main viewport, along with a plot in a separate window. The values depicted are defined on the basis of the integration points of the underlying solution, interpolated on the basis of standard shape functions for the mesh type. These values are also written to the working directory as a `pyCM_line_probe.csv` file. The format of this file follows x,y,z,Sxx where Sxx is the selected stress component.

<span>![<span>Main Window</span>](images/postprocess2.png)</span>
*<a name="fig2"></a> Figure 2: Postprocess lineprobe extraction.*

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