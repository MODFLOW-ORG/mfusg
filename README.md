
# MODFLOW-USG

[![mfusg checks](https://github.com/MODFLOW-ORG/mfusg/actions/workflows/ci.yml/badge.svg)](https://github.com/MODFLOW-ORG/mfusg/actions/workflows/ci.yml)

[USGS Release Page](https://water.usgs.gov/ogw/mfusg)

## Introduction

A version of MODFLOW, called MODFLOW-USG (for UnStructured Grid), was developed to support a wide variety of structured and unstructured grid types, including nested grids and grids based on prismatic triangles, rectangles, hexagons, and other cell shapes. Flexibility in grid design can be used to focus resolution along rivers and around wells, for example, or to subdiscretize individual layers to better represent hydrostratigraphic units.

MODFLOW-USG is based on an underlying control volume finite difference (CVFD) formulation in which a cell can be connected to an arbitrary number of adjacent cells. To improve accuracy of the CVFD formulation for irregular grid-cell geometries or nested grids, a generalized Ghost Node Correction (GNC) Package was developed, which uses interpolated heads in the flow calculation between adjacent connected cells.

MODFLOW-USG includes a Groundwater Flow (GWF) Process, based on the GWF Process in MODFLOW-2005, as well as a new Connected Linear Network (CLN) Process to simulate the effects of multi-node wells, karst conduits, and tile drains, for example. The CLN Process is tightly coupled with the GWF Process in that the equations from both processes are formulated into one matrix equation and solved simultaneously. This robustness results from using an unstructured grid with unstructured matrix storage and solution schemes.

MODFLOW-USG also contains an optional Newton-Raphson formulation, based on the formulation in MODFLOW-NWT, for improving solution convergence and avoiding problems with the drying and rewetting of cells. Because the existing MODFLOW solvers were developed for structured and symmetric matrices, they were replaced with a new Sparse Matrix Solver (SMS) Package developed specifically for MODFLOW-USG. The SMS Package provides several methods for resolving nonlinearities and multiple symmetric and asymmetric linear solution schemes to solve the matrix arising from the flow equations and the Newton-Raphson formulation, respectively.


## How to Cite

### ***Citations for MODFLOW-USG:***

[Panday, Sorab, Langevin, C.D., Niswonger, R.G., Ibaraki, Motomu, and Hughes, J.D., 2013, MODFLOW-USG version 1: An unstructured grid version of MODFLOW for simulating groundwater flow and tightly coupled processes using a control volume finite-difference formulation: U.S. Geological Survey Techniques and Methods, book 6, chap. A45, 66 p.](https://pubs.usgs.gov/tm/06/a45)

#### ***Software/Code citation for MODFLOW 6:***

[Panday, Sorab, Langevin, C.D., Niswonger, R.G., Ibaraki, Motomu, and Hughes, J.D., 2017, MODFLOW-USG version 1.4.00, An unstructured grid version of MODFLOW for simulating groundwater flow and tightly coupled processes using a control volume finite-difference formulation: U.S. Geological Survey Software Release, 27 October 2017, https://dx.doi.org/10.5066/F7R20ZFJ](https://dx.doi.org/10.5066/F7R20ZFJ)


Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.

