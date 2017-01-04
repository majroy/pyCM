# fit_surface

## Background
Currently in beta - reads data from align_average and provides a GUI which is essentially a wrapper for scipy's FITPACK bivariate spline fitting function. Writes a .mat file containing the spline fit both from FITPACK as well as attempting to match Matlab's spline objects. Requires PyQt4.


## Hierarchy

* Module `fit_surface`
 - Function `interactor`

Import via `from pyCM import interactor`

## Initializing
GUI based.

**Input and output descriptors for the FEAtool function**
-Input from align_average
-Output in the form of a .mat file containing an outline and a spline object in both scipy and Matlab format.

##  Interaction functionality
Beyond the buttons located on the GUI, the following are available:

Key | Description
---  |---
Left mouse button 	|Rotate about the center of view
Middle mouse button 	|Pan
Right mouse button 	|Zoom/refresh window extents
1 	|View 1, default, looks down z axis onto xy plane
2 	|View 2, default, looks down x axis onto zy plane
3 	|View 3, default, looks down y axis onto zx plane
z | increase z-aspect ratio
x | decrease z-aspect ratio
c | return to default z-aspect
Shift-z | increase size of points
Shift-x | decrease size of points
Shift-c | return to default point size
f | flip colors from white on dark to dark on white
i | save output to .png in current working directory
r | remove/reinstate compass/axes
o | remove/reinstate outline
e | write output and exit

## Known issues
Largely untested.