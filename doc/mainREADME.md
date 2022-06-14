# main

## Background
Tabbed interface which groups all steps of the contour method together to read and write a common datafile, along with collating all finite element analysis files.

## Initializing

Once pyCM has been installed, this GUI interface can be launched by typing
~~~
python -m pyCM.main
~~~
from a command line. This will launch the `main` script which in turn calls all other step-specific applications in discrete tabs. To read more about the functionality of each tab, see below where links to each respective step reside. 

## Menu
At this time, the drop down menu is sparsely populated, with the only functionality being the ability to:

* To load all pyCM results file with a *.pyCM extension (File)
* To load data from an active pyCM file with data pertaining to the current active step (File)

## Registration

The first tab is the point cloud editor which is required to start an analysis. To do so, all that is required is the same data types described in [registration](registrationREADME.md). 

## Alignment and averaging

The second tab uses the registered entries and allows teh user to mirror, align and average them. For more information, see  [align_average](align_averageREADME.md). 


## Surface fitting tab

The third tab permits the user to fit an analytical surface to the aligned and averaged surface generated in the previous step. See [fit_surface](fit_surfaceREADME.md). 

## FEA preprocessing tab

This tab currently permits a user to refactor (seed) an outline, generate a mesh, and subsequently assign boundary conditions as dictated by the analytical surface found in the previous step. See [preprocess](preprocessREADME.md).

## FEA postprocessing tab

This tab currently permits a user to generate a postprocessing file, and can display results from either Abaqus or Calculix-generated analyses across the entire mesh. It displays a contour plot of resolved longitudinal stresses by default, and stresses relieved by the cut in other principle directions. The number of contours as well as the maximum/minimum stresses can be set directly. See [postprocess](postprocessREADME.md).

## Known issues
Over-writing of data is still possible at this stage, it is recommended that the user works from a copy of a database file if any changes are to be experimented with, or if a sensitivity analysis is to be conducted. Changes in the analysis stream will not necessary negate results upstream, and therefore the user is advised to be cautious when using this tool. It has largely been written to a) perform an analysis in one direction along the analysis path, and b) playback analyses that have been conducted.
