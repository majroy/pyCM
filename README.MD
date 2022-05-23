[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0) 

- [Introduction](#introduction)
  * [Contour method](#contour-method)
  * [Analysis code](#analysis-code)
  * [Moving forward](#moving-forward)
- [Latest version](#latest-version)
- [Requirements](#requirements)
- [Installation](#installation)
- [Modules](#modules)
- [Contributions and extensions](#contributions-and-extensions)
- [Acknowledgements](#acknowledgements)

# Introduction
NB: Documentation is currently being updated - please see [Latest version](#latest-version).
## Contour method
The contour method for measuring residual stress is a cost effective, robust technique which has seen much interest both by academics and industry alike. First developed by Michael Prime at Los Alamos, it has seen a great deal of development both in the US and UK.

For more information on the contour method:
* Michael Prime's website containing an extensive background and literature overview: http://www.lanl.gov/contour/
* Hill Engineering http://hill-engineering.com/contour-method/
* StressMap, UK provider of contour method measurements based at the Open University http://www.stressmap.co.uk
* VEQTER, expertise in a wide range of residual stress measurement techniques http://www.veqter.co.uk/residual-stress-measurement/contour

## Analysis code
Recently it has been found that the visualization tools that are available in MATLAB are not the most robust for the purposes of the contour method. For example, large datasets developed with laser profilometers can often exceed the graphical capabilities of recent versions of MATLAB. The reasons for moving from MATLAB for contour method analysis are numerous:

* Proprietary nature of the built-in libraries, cost and licensing issues with external ones.
* Difficulties in dealing with path dependencies on managed desktops.
* Although 'compiled' MATLAB code can be run on computers without MATLAB installed, this is reliant on the target machine having an identical MATLAB component runtime installation, which is proprietary and quite large (250 Mb+).
* Diminished text file parsing built into MATLAB which makes working with FE codes difficult; this is particularly true for very large input files/meshes.

Each of these factors contribute to the overall portability of the code and present barriers to the residual stress characterisation community developing, validating and using analysis code.

Some of these have already been addressed via projects such as [GNU Octave](https://www.gnu.org/software/octave/), a free and open source language similar to MATLAB in terms of development environment and syntax. However, this language has been shown to be quite slow in some cases, as inter-operability with MATLAB was favoured in place of speed. Inter-operability between GNU Octave and MATLAB is still relatively incomplete as MATLAB is developed independently of Octave.

Many researchers have moved from using MATLAB to Python entirely for numerical analysis. Python is an extremely flexible language, with many packages available to perform the same work as MATLAB. Due to the open source nature of Python, there have been many contributions that allow for wrapping compiled science and engineering code written in other languages (Fortran, C, Java) to be called directly by Python. When coupled with [Numerical Python (NumPy)](http://www.numpy.org/)  and [SciPy](http://www.vtk.org/overview/), a significant amount of the same functionality is available as in MATLAB. Furthermore, there is a maintained library of Python bindings for VTK, the [Visualization Toolkit](http://www.scipy.org/), an extremely robust tool for visualising large datasets. For example, ParaView, an open source FEA post processor, makes extensive use of VTK - and the results of this application can be opened natively with this tool.

## Moving forward
It is for these reasons the overall aim is to future proof the analysis code by moving towards using Python as the language of choice. To this end, a Python package called pyCM (Python Contour Method) has been developed. The idea is to fully move over from using primarily MATLAB functions/scripts to Python equivalents in a step-by-step manner. The development strategy adopted is to build new tools where required, and allow for existing MATLAB tools to be used on the input and output from these Python tools. This will be accomplished using the SciPy packages which generate binary MATLAB files to store workspace variables, so processing can be partially performed in Python, and shifted back into MATLAB to use existing scripts.

# Latest version
Beta release 2.0.0 is the latest version of pyCM, which is currently being tested and documented. Note that this version is *not* backwards compatible with the last version, but the code for that version will remain here.

Alpha release 1.1.8post3 is last version which contains a complete toolchain for rectangular specimens with no self-restraint. Note that the post-processor is still in a beta form for testing purposes on *quadrilateral meshes only*; if tetrahedral meshes are employed, then the native FEA code's post processor should be employed (GraphiX or Abaqus Viewer).

**It has been confirmed that equivalent results are obtained with both CalculiX and Abaqus**.

# Requirements
Python 3.x with pip working correctly will install most functionality, as tested on Linux (Ubuntu) and Windows 10. If outline processing is required, then there is package that is needed that cannot be installed on Windows via pip - see below.

In order to optimize performance on Windows, it is recommended to get numpy, scipy and VTK binaries from Christoph Gohlke's [Unofficial Windows Binaries for Python Extension Packages](http://www.lfd.uci.edu/~gohlke/pythonlibs/). Install them via:

~~~
pip install <pkg_prefix>.whl
~~~

However, **the sources which arrive from PyPI are acceptable for most purposes and should automatically install**, with the exception of [shapely](https://pypi.org/project/Shapely/). Shapely is required for outline processing - see "Outline processing" over view for [point_cloud](doc/point_cloudREADME.md#Outline-processing).



# Installation

These instructions cover the steps needed to install all dependencies and then add the pyCM components to your Python PATH. This permits calling the pyCM components anywhere on a host machine.

1. Copy the appropriate pyCM-<version>.tar.gz file located in the *dist* folder.
2. From the command line, run `pip install pyCM-*` where `pyCM-*` file that was copied.

~~~
pip install pyCM-<version>.tar.gz
~~~

The remainder of the requisite packages will be installed from PyPI **with the exception of Shapely**. This latter package needs to be installed manually.

3. Check if it installed properly by trying to import the pyCM library from any working directory at the Python prompt/IDLE:
~~~
from pyCM import *
~~~
or a check of installed packages via the `pip freeze` command.

Uninstalling can be accomplished via `pip`:
~~~
pip uninstall pyCM
~~~

# Modules

* [main](doc/mainREADME.md) - a tabbed viewer which groups the following tools together to allow a more cohesive analysis route.
* [point_cloud](doc/point_cloudREADME.md) - utility to read and work with raw contour data both with and without perimeter files.
* [align_average](doc/align_averageREADME.md) - utility to manipulate, align and average two contour surface datasets with perimeter files.
* [fit_surface](doc/fit_surfaceREADME.md) - utility to fit averaged point clouds
* [preprocess](doc/preprocessREADME.md) - meshing/boundary condition building tool which works with Abaqus and OS FEA toolchains (Gmsh & Calculix).
* [postprocess](doc/postprocessREADME.md) - **BETA** currently only supports results from analyses with quadrilateral meshes, but for both Calculix and Abaqus.

# Contributions and extensions
This project is intended to best serve the overall residual stress community, so feel free to adjust as needed and add functionality. Please contact me directly by email prior to submitting pull requests so that we can discuss overall structure further.

Currently, this project has only been tested on Windows 10 x64 and limited Linux distributions. It is expected, but not tested to be compatible with all platforms.

# Acknowledgements
A number of these tools have been developed based on MATLAB code written by Greg Johnson, Chris Gill, Matt Fox, Philip Frankel and Yuan-Sheng Xiong (and probably others) at the University of Manchester between 2004-2013.

I'd also like to acknowledge Richard Moat at the Open University for outlining the architecture of Python and MATLAB interchangeability and highlighting the performance improvements in visualization garnered by VTK.