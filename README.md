#BOUT++
<!---Build nice shields at shields.io-->
[![Build Status](https://travis-ci.org/boutproject/BOUT-dev.svg?branch=master)](https://travis-ci.org/boutproject/BOUT-dev)
[![License](https://img.shields.io/badge/license-GPL-blue.svg)](https://img.shields.io/badge/license-GPL-blue.svg)
[![py3comp](https://img.shields.io/badge/py3-compatible-brightgreen.svg)](https://img.shields.io/badge/py3-compatible-brightgreen.svg)
> ```
> .______     ______    __    __  .___________.
> |   _  \   /  __  \  |  |  |  | |           |  _     _
> |  |_)  | |  |  |  | |  |  |  | `---|  |----`_| |_ _| |_
> |   _  <  |  |  |  | |  |  |  |     |  |    |_   _|_   _|
> |  |_)  | |  `--'  | |  `--'  |     |  |      |_|   |_|
> |______/   \______/   \______/      |__|
> ```

A modular fluid simulation code written in C++.
Adapted from the original BOUT 3D 2-fluid simulation code.
Aims to be able to solve a wide variety of fluid models in
almost any curvilinear coordinate system.

Developed by B.Dudson (University of York) and M.Umansky (LLNL)
based on the original BOUT code by X.Xu (LLNL)
Further development by S.Farley (ANL) for coupling to PETSc and FACETS

Homepage found at [http://boutproject.github.io/](http://boutproject.github.io/)

## Table of Contents
* [Usage and installation](#usage-and-installation)
* [Overview of files](#overview-of-files)
* [Contributing](#contributing)
* [Terms of use](#terms-of-use)
* [License](#license)

## Usage and installation
Please see the users manual in the [manual](/manual/) section

## Overview of files

This directory contains

* **bin**                   Files for setting the BOUT++ configuration
* **examples**              Example models and test codes
* **externalpackages**      External packages needed for installing BOUT++
* **include**               Header files used in BOUT++
* **manual**                Manuals and documentation (also [doxygen](http://www.stack.nl/~dimitri/doxygen/) documentation)
* **src**                   The main code directory
* **CITATION**              Contains the paper citation for BOUT++
* **LICENSE**               GPL license
* **LICENSE.LESSER**        LGPL license
* **tools**                 Tools for helping with analysis, mesh generation, and data managment

  * **archiving**               Routines for managing input/output files e.g. compressing data, converting formats, and managing runs
  * **cyl_and_helimak_grids**   IDL codes for generating cylindrical and helimak grids
  * **eigensolver**             Matlab routines for solving eigenmodes
  * **idllib**                  Analysis codes in IDL. Add this to your IDL_PATH environment variable
  * **line_tracing**            IDL routines for line tracing of field lines
  * **line_tracing_v2**         Newer version of the IDL routines for line tracing of field lines
  * **mathematicalib**          Library for post processing using mathematica
  * **matlablib**               Library for post processing using mathlab
  * **numlib**                  Numerical IDL routines
  * **octave**                  Routines for post processing using octave
  * **plasmalib**               IDL routines for calculation of plasma parameters
  * **pdb2idl**                 Library to read Portable Data Binary (PDB) files into IDL
  * **pylib**                   Analysis codes in Python

    * **boutdata**        Routines to simplify accessing BOUT++ output
    * **boututils**       Some useful routines for accessing and plotting data
    * **bout_runners**    A python wrapper to submit several runs at once (either on a normal computer, or through a PBS system)
    * **post_bout**       Routines for post processing in BOUT++

  * **slab**              IDL routine for grid generation of a slab
  * **tokamak_grids**     Code to generate input grids for tokamak equilibria

    * **gridgen**         Grid generator in IDL. Hypnotoad GUI for converting G-EQDSK files into a flux-aligned orthogonal grid.
    * **elite**           Convert ELITE .eqin files into an intermediate binary file
    * **gato**            Convert DSKGATO files into intermediate binary format
    * **all**             Convert the intermediate binary file into BOUT++ input grid
    * **coils**           Routines for calculating the field due to external RMP coils and adding to existing equilibria
    * **cyclone**         Generate cyclone test cases (concentric circle "equilibrium" for local flux-surface calculations)
    * **py_gridgen**      Translation" into python of the corresponding IDL routines in the folder gridgen
    * **shifted_circle**  Produce shifted cirle equilibria input grids


## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Terms of use

BOUT++ is released under the LGPL, but since BOUT++ is a
scientific code we also ask that you show professional courtesy
when using this code:

1. Since you are benefiting from work on BOUT++, we ask that you
   submit any improvements you make to the code to us by emailing
   Ben Dudson at bd512@york.ac.uk
2. If you use BOUT++ results in a paper or professional publication,
   we ask that you send your results to one of the BOUT++ authors
   first so that we can check them. It is understood that in most cases
   if one or more of the BOUT++ team are involved in preparing results
   then they should appear as co-authors.
3. Publications or figures made with the BOUT++ code should acknowledge the
   BOUT++ code by citing B.Dudson et. al. Comp.Phys.Comm 2009 and/or
   other BOUT++ papers. See the file CITATION for details.



## License
Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu

BOUT++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BOUT++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

A copy of the LGPL license is in [LICENSE.LESSER](LICENSE.LESSER). Since this is based
on (and refers to) the GPL, this is included in [LICENSE](LICENSE).
