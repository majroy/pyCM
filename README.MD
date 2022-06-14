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
Invented by Michael Prime at Los Alamos, the contour method for measuring residual stress is a cost effective, robust technique which has seen much interest both by academics and industry alike.

For more information on the contour method:
* [Michael Prime's website](http://www.lanl.gov/contour/) containing an extensive background and literature overview. 
* [Hill Engineering](http://hill-engineering.com/contour-method/)
* [StressMap](http://www.stressmap.co.uk), UK provider of contour method measurements based at the Open University 
* [VEQTER](http://www.veqter.co.uk/residual-stress-measurement/contour), expertise in a wide range of residual stress measurement techniques 

To get a high level overview of the approach taken to contour method analysis taken here, please see our paper:
["pyCM: An open-source computational framework for residual stress analysis employing the Contour Method" by Roy et al. (2020).](https://www.sciencedirect.com/science/article/pii/S2352711019303668)

# Latest version
Beta release 2.0.0 is the latest version of pyCM. is *not* backwards compatible with previous versions. It offers substantial improvements over version 1.1.8, as well as rolling up many updates aimed at supporting the most recent packages. This release is an improvement over 1.1.8 in the following ways:
* The use of an HDF5 formatted storage file for speed and more transparent storage format
* The introduction of 'entries', meaning different domains can be processed independently. This means that it now supports items like pipes, which have two independent surfaces per cut face.
* Registration step includes more filters for decimation and height prviews
* Alignment and averaging has been improved, along with the ability to fill missing data with nearest neighbours.
* Surface fitting now shows a map of error between fitted surfaces and averaged data. A direct, exportable 2D comparison between fit and data on a plane has been implemented
* Finite element preprocessing allows the import and alignment of meshes.
* Post processing allows for an export of both line and all stresses at nodes. Both linear quadrilateral and second order tetrahedral meshes are supported.

Alpha release 1.1.8post3 is a legacy version which contains a complete toolchain for rectangular specimens with no self-restraint. Note that the post-processor is still in a beta form for testing purposes on *quadrilateral meshes only*; if tetrahedral meshes are employed, then the native FEA code's post processor should be employed (GraphiX or Abaqus Viewer).

**It has been confirmed that equivalent results are obtained with both CalculiX and Abaqus, CalculiX is significantly faster performing for smaller analyses**.

# Requirements
Python 3.9 for v2.0.0 has been used due to package dependencies at the time of development. Python 3.8 may work, 3.7 will not. The python package installer, pip, will install all package dependencies, as tested on Linux (Ubuntu) and Windows 10.

# Installation

These instructions cover the steps needed to install all dependencies and then add the pyCM components to your Python PATH. This permits calling the pyCM components anywhere on a host machine.

1. Copy the appropriate pyCM-<version>.tar.gz file located in the *dist* folder.
2. From the command line, run `pip install pyCM-*` where `pyCM-*` file that was copied, i.e. `pip install pyCM-<version>.tar.gz`. All required packages will be installed from PyPI.
3. Check if it installed properly by trying to import the pyCM library from any working directory at the Python prompt/IDLE with `from pyCM import *` or a check of installed packages via the `pip freeze` command. Uninstalling can be accomplished via pip and the command: `pip uninstall pyCM`

# Modules

* [main](doc/mainREADME.md) - a tabbed viewer which groups the following tools together to allow a more cohesive analysis route.
* [registration](doc/registrationREADME.md) - utility to read and work with raw contour data both with and without perimeter files.
* [align_average](doc/align_averageREADME.md) - utility to manipulate, align and average two contour surface datasets with perimeter files.
* [fit_surface](doc/fit_surfaceREADME.md) - utility to fit averaged point clouds
* [preprocess](doc/preprocessREADME.md) - meshing/boundary condition building tool. Also allows finite element analyses to run with either CalculiX or Abaqus.
* [postprocess](doc/postprocessREADME.md) - View and export results.

# Contributions and extensions
This project is intended to best serve the overall residual stress community, so feel free to adjust as needed and add functionality. Please contact me directly by email prior to submitting pull requests so that we can discuss overall structure further.

Currently, this project has only been tested on Windows 10 x64 and limited Linux distributions. It is expected, but not tested to be compatible with all platforms.

# Acknowledgements
A number of these tools have been developed based on MATLAB code written by Greg Johnson, Chris Gill, Matt Fox, Philip Frankel and Yuan-Sheng Xiong (and probably others) at the University of Manchester between 2004-2013.

I'd also like to acknowledge Richard Moat at the Open University for outlining the architecture of Python and MATLAB interchangeability and highlighting the performance improvements in visualization garnered by VTK.