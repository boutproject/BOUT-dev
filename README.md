# BOUT++

<!---Build nice shields at shields.io-->
[![Build Status](https://github.com/boutproject/BOUT-dev/actions/workflows/tests.yml/badge.svg?branch=next)](https://github.com/boutproject/BOUT-dev/actions)
[![License](https://img.shields.io/badge/license-LGPL-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0.en.html)
[![py3comp](https://img.shields.io/badge/py3-compatible-brightgreen.svg)](https://img.shields.io/badge/py3-compatible-brightgreen.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13753882.svg)](https://doi.org/10.5281/zenodo.8369888)

```
.______     ______    __    __  .___________.
|   _  \   /  __  \  |  |  |  | |           |  _     _
|  |_)  | |  |  |  | |  |  |  | `---|  |----`_| |_ _| |_
|   _  <  |  |  |  | |  |  |  |     |  |    |_   _|_   _|
|  |_)  | |  `--'  | |  `--'  |     |  |      |_|   |_|
|______/   \______/   \______/      |__|
```

BOUT++ is a framework for writing fluid and plasma simulations in
curvilinear geometry. It is intended to be quite modular, with a
variety of numerical methods and time-integration solvers
available. BOUT++ is primarily designed and tested with reduced plasma
fluid models in mind, but it can evolve any number of equations, with
equations appearing in a readable form.

For example, the following set of equations for magnetohydrodynamics
(MHD):

![ddt_rho](http://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20%5Crho%7D%7B%5Cpartial%20t%7D%20%3D%20-%5Cmathbf%7Bv%7D%5Ccdot%5Cnabla%5Crho%20-%20%5Crho%5Cnabla%5Ccdot%5Cmathbf%7Bv%7D)

![ddt_p](http://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20p%7D%7B%5Cpartial%20t%7D%20%3D%20-%5Cmathbf%7Bv%7D%5Ccdot%5Cnabla%20p%20-%20%5Cgamma%20p%5Cnabla%5Ccdot%5Cmathbf%7Bv%7D)

![ddt_v](http://latex.codecogs.com/png.latex?%5Cfrac%7B%5Cpartial%20%5Cmathbf%7Bv%7D%7D%7B%5Cpartial%20t%7D%20%3D%20-%5Cmathbf%7Bv%7D%5Ccdot%5Cnabla%5Cmathbf%7Bv%7D%20&plus;%20%5Cfrac%7B1%7D%7B%5Crho%7D%28-%5Cnabla%20p%20&plus;%20%28%5Cnabla%5Ctimes%5Cmathbf%7BB%7D%29%5Ctimes%5Cmathbf%7BB%7D%29)

![ddt_B](http://latex.codecogs.com/png.latex?%7B%7B%5Cfrac%7B%5Cpartial%20%5Cmathbf%7BB%7D%7D%7B%5Cpartial%20t%7D%7D%7D%20%3D%20%5Cnabla%5Ctimes%28%5Cmathbf%7Bv%7D%5Ctimes%5Cmathbf%7BB%7D%29)

can be written simply as:

```cpp
ddt(rho) = -V_dot_Grad(v, rho) - rho*Div(v);
ddt(p)   = -V_dot_Grad(v, p) - g*p*Div(v);
ddt(v)   = -V_dot_Grad(v, v) + (cross(Curl(B),B) - Grad(p))/rho;
ddt(B)   = Curl(cross(v,B));
```

The full code for this example can be found in the [orszag-tang
example](examples/orszag-tang/mhd.cxx).

Jointly developed by University of York (UK), LLNL, CCFE, DCU, DTU,
and other international partners. See the Git logs for author details.


Homepage found at [http://boutproject.github.io/](http://boutproject.github.io/)

## Table of Contents
* [Requirements](#requirements)
* [Usage and installation](#usage-and-installation)
* [Terms of use](#terms-of-use)
* [Contributing](#contributing)
* [License](#license)

## Requirements

BOUT++ needs the following:

* A C++17 compiler
* MPI
* NetCDF

BOUT++ has the following optional dependencies:

* [FFTW3](https://www.fftw.org/) (strongly recommended!)
* [SUNDIALS](https://computing.llnl.gov/projects/sundials): CVODE, IDA, ARKODE
* [PETSc](https://petsc.org)
* [ADIOS2](https://adios2.readthedocs.io/)
* [SLEPc](https://slepc.upv.es/)
* LAPACK
* OpenMP
* Score-p (for performance diagnostics)

## Usage and installation
Please see the [users manual](http://bout-dev.readthedocs.io).

## Terms of use

BOUT++ is released under the LGPL, but since BOUT++ is a
scientific code we also ask that you show professional courtesy
when using this code:

1. Since you are benefiting from work on BOUT++, we ask that you
   submit any improvements you make to the code to us by submitting a
   pull request to this repository
2. If you use BOUT++ results in a paper or professional publication,
   we ask that you send your results to one of the BOUT++ authors
   first so that we can check them. It is understood that in most cases
   if one or more of the BOUT++ team are involved in preparing results
   then they should appear as co-authors.
3. If you use BOUT++ in your work, please help ensure that all the
   authors get the credit they deserve by citing BOUT++, preferably
   using the DOI of the version you used. See the file
   [CITATION.cff](CITATION.cff) for details. In addition, you may also
   cite either of the two main papers: B. Dudson et al,
   Comp. Phys. Comm. 2009, and B. Dudson et al, Phys. of Plasmas 2016

You can convert the CITATION.cff file into a Bibtex file as follows:

    pip3 install --user cffconvert
    cffconvert -if CITATION.cff -f bibtex -of CITATION.bib



## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) and the [manual page](https://bout-dev.readthedocs.io/en/stable/developer_docs/contributing.html)

## License
Copyright 2010-2024 BOUT++ contributors

BOUT++ is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BOUT++ is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

A copy of the LGPL license is in [LICENSE](LICENSE). Since this is based
on (and refers to) the GPL, this is included in [LICENSE.GPL](LICENSE.GPL).

BOUT++ links by default with some GPL licensed libraries. Thus if you
compile BOUT++ with any of them, BOUT++ will automatically be licensed
as GPL. Thus if you want to use BOUT++ with GPL non-compatible code,
make sure to compile without GPLed code.
