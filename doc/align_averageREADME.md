# align_average

## Background
Performing a contour method generates two surfaces where one needs to be mirrored and then subsequently aligned and then averaged to eliminate/minimize cutting artefacts.

This module allows for viewing the two surfaces, with one being a reference and the other the surface to be aligned and averaged (floating). The floating surface can be then mirrored, and aligned with the reference surface. Upon alignment, then the two surfaces are averaged along a common grid dictated by the extents of the reference surface outline.

## Initializing

**Input and output descriptors for the `aa_def` function**

Input | Description
---  |---
Output file	| A *.mat file which at minimum, contains two data structures needed for subsequent processing, *ref* and *float* which contain the following:<ul><li>x,y,z: Nx1 arrays of the masked coordinate values. </li><li>rawPnts: Nx3 matrix of the points read in via the point cloud file.</li><li>mask: 1xN array of int8 values consisting of 0 and 1 where 0 indicates a masked point. Conversion to a boolean array will provide an index of rawPnts that were masked.</li><li>x_out: Nx3 matrix of the points that comprise the outline</li></ul> 

All output data is is added to the *.mat file which contained the originating data.

Output | Description
---  |---
Aligned and averaged points | A structure called `aa`, which contains the fields `pnts` - Nx3 matrix of points comprising the aligned and averaged data and `gsize` - the characteristic length of the grid that was used to average the data.
Transformation matrices | A structure called `trans` with fields `ref` and `float`, 4x4xN homogeneous transformation matrices of N operations carried out on each dataset, in the order in which they were performed.


The function can be called from interactive Python according to:
~~~
from pyCM import align_average as aa
aa.aa_def()
~~~
which will provide a GUI to locate the *.mat file with the floating and reference point cloud data, according to the description above. Alternatively, to specify the results file directly:
~~~
from pyCM import point_cloud as aa
aa.aa_def("PathToMatFile.mat")
~~~


##  Interaction functionality
The data will appear in a custom VTK interaction window after initializing. A set of axes is provided in the x and y directions to identify principle axes. The default view is looking down on the data in the z direction, after which the view can be rotated to show the data in perspective ([Fig. 1](#fig1)) by pressing the left mouse button. The middle mouse button provides a pan function while pressed, and the right mouse button zooms. There are three named views that are accessed via `1`, `2` and `3` looking down the z, x and y directions, respectively. 

<span>![<span>Main Window</span>](images/Avg_loaded.png)</span>  
*<a name="fig1"></a> Figure 1: Loaded data with the reference set in orange and the floating data set in yellow.*

For publication purposes, the ability to flip the default color scheme (dark on bright) has been provided. This is obtained by pressing f on the keyboard. Again, for publication purposes, a facility has been provided for printing the interaction window to file. Pressing `i` will print the interaction window to the current working directory as `Avg_aligned.png`.

As in other modules, a facility for increasing the aspect of the data has been provided, across all principle axes, according to the radio button selected. Pressing `z` increases the aspect ratio by 2x with each keypress, pressing `x` decreases by half, and `c` returns to the default aspect ratio.

Using the push buttons under the `Mirroring` pane will provide the ability to mirror the floating surface prior to alignment and averaging.

Automated alignment of the floating data to the reference is accomplished via an iterative closest point (ICP) technique, a [pre-existing VTK filter](https://www.vtk.org/doc/nightly/html/classvtkIterativeClosestPointTransform.html#details), or another K-neighbour ICP technique, depending on what was selected prior to pressing the `Align` button. Both of these algorithms seek to match the vertex in one surface with the closest surface point in the other, then applying the transformation which best matches in a least-square sense. Two ICP algorithms are provided as the VTK-sourced algorithm is fast, but closed to development. As there is often a threshold number of points deciding whether the K-neighbour ICP technique's success, one may alter the number of points on the outline to employ for alignment, by changing the number beside (and pressing) the `Decimate outlines` button.

For large point clouds which are nearly aligned, the iterative closest point has been noted to occasionally fail as the error between the starting transformation and the final is near floating point accuracy. Therefore, an intermediate translation is possible by entering values for x and y (or rotation about the z axis) translation and pressing the `Translate` button. Additionally, the ICP alignment can occasionally return an alignment that contains a 180 degree rotation about either the x and y axes. When this occurs, one may reverse it by pressing the relevant `Flip X` or `Flip Y` buttons and retrying alignment.

If automated alignment fails, then in some circumstances, the only route is to manually align the outlines. Far from ideal, best practice is to move both reference and floating centroids to the origin, and then perform incremental transformations until the outlines visually converge. Note that if alignment fails, this is likely because the 'hulls' of the profiles do not match; this is a keen indication that there is something wrong with the incoming data and the source will need to be reconsidered. Potential causes for this include i) the wrong files/data were selected or ii) the reference outline has a smaller included area than the floating. Some success has been found previously by reversing datasets *e.g.* which is the reference and which is floating. If one is satisfied with the alignment, then pressing the `Accept` button will move the analysis forward.

Once alignment has taken place, then averaging can take place by pressing the `Average` button. The data will initially be gridded and averaged according to the mean point spacing in the reference point cloud ([Fig. 2](#fig2)). Afterwards, the grid size can be modified - the GUI will accept values between 0.001 and 5 mm. However, if the underlying point cloud is too coarse or too dense, then other issues may arise *viz* such as memory issues or an artificial change in effective resolution of the final resolved stresses. 

<span>![<span>ZAspect</span>](images/Avg_averaged.png)</span>  
*<a name="fig2"></a> Figure 2: Reference, floating and aligned & averaged datasets shown in orange, yellow and light blue, respectively.*

Pressing the 'Write' button will write data to the originating *.mat file. 

A complete list of interaction keys is provided below. 

**Keyboard and mouse mapping**

Key | Description
---  |---
Left mouse button 	|Rotate about the center of view
Middle mouse button 	|Pan
Right mouse button 	|Zoom/refresh window extents
1 	|View 1, default, looks down z axis onto xy plane
2 	|View 2, default, looks down x axis onto zy plane
3 	|View 3, default, looks down y axis onto zx plane
z 	|Increase z-aspect ratio by a factor of 2
x 	|Decrease z-aspect ratio by a factor of 0.5
c 	|Return to default z-aspect ratio; x,y:z=1:1
f 	|Flip colour scheme from bright on dark to dark on bright.
i 	|Save visualization window as `Avg_aligned.png` to the current working directory
l	|Load a .mat file that contains a *ref* and *float* data structure specified above.

## Performance
The current version of this module has been tested with point clouds containing approximately 50,000 points each. Mirroring, alignment, averaging and output were all found to be near instantaneous. Extremely large datasets (on the order of 3 million points) were found to have Python attempt to demand too much memory.


## Known issues

Loading of extremely large datasets (on the order of a million points or more) has shown to create serious lag. Point clouds are better off sampled and reduced before using these tools. Not all hardware is supported; OpenGL errors have been noted when using 4k displays on older versions of VTK.