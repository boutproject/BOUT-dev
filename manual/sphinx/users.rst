Introduction
============

BOUT++ is a C++ framework for writing plasma fluid simulations with an
arbitrary number of equations in 3D curvilinear coordinates . It has
been developed from the original **BOU**\ ndary **T**\ urbulence 3D
2-fluid edge simulation code written by X.Xu and M.Umansky at LLNL.

Though designed to simulate tokamak edge plasmas, the methods used are
very general and almost any metric tensor can be specified, allowing the
code to be used to simulate (for example) plasmas in slab, sheared slab,
and cylindrical coordinates. The restrictions on the simulation domain
are that the equilibrium must be axisymmetric (in the z coordinate), and
that the parallelisation is done in the :math:`x` and :math:`y`
(parallel to :math:`\mathbf{B}`) directions.

The aim of BOUT++ is to automate the common tasks needed for simulation
codes, and to separate the complicated (and error-prone) details such as
differential geometry, parallel communication, and file input/output
from the user-specified equations to be solved. Thus the equations being
solved are made clear, and can be easily changed with only minimal
knowledge of the inner workings of the code. As far as possible, this
allows the user to concentrate on the physics, rather than worrying
about the numerics. This doesn’t mean that users don’t have to think
about numerical methods, and so selecting differencing schemes and
boundary conditions is discussed in this manual. The generality of the
BOUT++ of course also comes with a limitation: although there is a large
class of problems which can be tackled by this code, there are many more
problems which require a more specialised solver and which BOUT++ will
not be able to handle. Hopefully this manual will enable you to test
whether BOUT++ is suitable for your problem as quickly and painlessly as
possible.

This manual is written for the user who wants to run (or modify)
existing plasma models, or specify a new problem (grid and equations) to
be solved. In either case, it’s assumed that the user isn’t all that
interested in the details of the code. For a more detailed descriptions
of the code internals, see the developer and reference guides. After
describing how to install BOUT++ (section [sec:install]), run the test
suite (section [sec:runtestsuite]) and a few examples
(section [sec:running], more detail in section [sec:examples]),
increasingly sophisticated ways to modify the problem being solved are
introduced. The simplest way to modify a simulation case is by altering
the input options, described in section [sec:options]. Checking that the
options are doing what you think they should be by looking at the output
logs is described in section [sec:running], and an overview of the IDL
analysis routines for data post-processing and visualisation is given in
section [sec:output]. Generating new grid files, particularly for
tokamak equilibria, is described in section [sec:gridgen].

Up to this point, little programming experience has been assumed, but
performing more drastic alterations to the physics model requires
modifying C++ code. Section [sec:equations] describes how to write a new
physics model specifying the equations to be solved, using ideal MHD as
an example. The remaining sections describe in more detail aspects of
using BOUT++: section [sec:diffops] describes the differential operators
and methods available; section [sec:staggergrids] covers the
experimental staggered grid system.

Various sources of documentation are:

-  Most directories in the BOUT++ distribution contain a README file.
   This should describe briefly what the contents of the directory are
   and how to use them.

-  This user’s manual, which goes through BOUT++ from a user’s point of
   view

-  The developer’s manual, which gives details of the internal working
   of the code.

-  The reference guide, which summarises functions, settings etc.
   Intended more for quick reference rather than a guide.

- Most of the code contains Doxygen comment tags (which are slowly
  getting better). Running `doxygen <www.doxygen.org>`_ on these files
  should therefore generate an HTML reference. This is probably going
  to be the most up-to-date documentation.

License and terms of use
------------------------

::

    Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu

    BOUT++ is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    BOUT++ is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.

    A copy of the LGPL license is in COPYING.LESSER. Since this is based
    on (and refers to) the GPL, this is included in COPYING.

BOUT++ is free software, but since it is a scientific code we also ask
that you show professional courtesy when using this code:

#. Since you are benefiting from work on BOUT++, we ask that you submit
   any improvements you make to the code to us by emailing Ben Dudson at
   bd512@york.ac.uk

#. If you use BOUT++ results in a paper or professional publication, we
   ask that you send your results to one of the BOUT++ authors first so
   that we can check them. It is understood that in most cases if one or
   more of the BOUT++ team are involved in preparing results then they
   should appear as co-authors.

#. Publications or figures made with the BOUT++ code should
   acknowledge the BOUT++ code by citing `B.Dudson
   et. al. Comp.Phys.Comm 2009`_ and/or other BOUT++ papers. See the
   file CITATION for details.

   .. _B.Dudson et. al. Comp.Phys.Comm 2009: http://www.sciencedirect.com/science/article/B6TJ5-4VTCM95-3/2/ed200cd23916d02f86fda4ce6887d798

Getting started
===============

This section goes through the process of getting, installing, and
starting to run BOUT++. Only the basic functionality needed to use
BOUT++ is described here; the next section ([sec:advancedinstall]) goes
through more advanced options, and how to fix some common problems.

On large facilities (e.g NERSC or Archer), the compilers and libraries
needed should already be installed. It is common to organise libraries
using the ``modules`` system, so try typing “``module avail``” to get a
list of available modules. Some instructions for specific machines can
be found in appendix [apx:machineinstructions]. If you are installing on
your own machine, you may need to install an MPI compiler and the
libraries yourself.

This section will go through the following steps:

#. :ref:`Obtaining a copy of BOUT++ <sec-obtainbout>`

#. :ref:`Installing an MPI compiler <sec-installmpi>`

#. :ref:`Installing libraries <sec-libraries>`

#. :ref:`Configuring BOUT++ analysis codes <sec-configanalysis>`

#. :ref:`Compiling BOUT++ <sec-installbout>`

#. :ref:`Running the test suite <sec-runtestsuite>`

**Note**: In this manual commands to run in a BASH shell will begin with
’$’, and commands specific to CSH with a ’%’.

.. _sec-obtainbout:

Obtaining BOUT++
----------------

BOUT++ is hosted publicly on github at
http://github.com/boutproject/BOUT-dev. You can the latest stable
version from https://github.com/boutproject/BOUT-dev/releases. If you
want to develop BOUT++, you should use git to clone the repository. To
obtain a copy of the latest version, run

.. code-block:: bash

    $ git clone git://github.com/boutproject/BOUT-dev.git

which will create a directory ``BOUT-dev`` containing the code. To get
the latest changes later, go into the ``BOUT-dev`` directory and run

.. code-block:: bash

    $ git pull

Development is done on the “next” branch, which you can checkout with

.. code-block:: bash

    $ git checkout next

.. _sec-installmpi:

Installing an MPI compiler
--------------------------

To compile and run the examples BOUT++ needs an MPI compiler. If you are
installing on a cluster or supercomputer then the MPI C++ compilers will
already be installed, and on Cray or IBM machines will probably be
called ’CC’ and ’xlC’ respectively. If you’re installing on a smaller
server or your own machine then you need to check that you have an MPI
compiler by running

.. code-block:: bash

    $ mpicc

This should produce an error message along the lines of “no input
files”, but if you see something like “command not found” then you need
to install MPI first. There are several free MPI distributions
available, the main ones currently being MPICH2
(`www.mcs.anl.gov/mpi/mpich2 <www.mcs.anl.gov/mpi/mpich2>`__), OpenMPI
(`www.open-mpi.org/ <www.open-mpi.org/>`__), and LAM
(`www.lam-mpi.org/ <www.lam-mpi.org/>`__). On Ubuntu or Debian
distributions if you have administrator rights then you can install
MPICH2 by running

.. code-block:: bash

    $ sudo apt-get install mpich2 libmpich2-dev

If this works, and you now have a working ``mpicc`` command, skip to the
next section on installing libraries. If not, and particularly if you
don’t have administrator rights, you should install MPI in your home
directory by compiling it from source. In your home directory, create
two subdirectories: One called “install” where we’ll put the source
code, and one called “local” where we’ll install the MPI compiler:

.. code-block:: bash

    $ cd
    $ mkdir install
    $ mkdir local

Download the latest stable version of MPICH2 from
http://www.mcs.anl.gov/research/projects/mpich2/downloads/ and put the
file in the “install” subdirectory created above. At the time of writing
(June 2012), the file was called ``mpich2-1.4.1p1.tar.gz``. Untar the
file:

.. code-block:: bash

    $ tar -xzvf mpich2-1.4.1p1.tar.gz

which will create a directory containing the source code. ’cd’ into this
directory and run

.. code-block:: bash

    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

Each of which might take a while. This is the standard way of installing
software from source, and will also be used for installing libraries
later. The ``–prefix=`` option specifies where the software should be
installed. Since we don’t have permission to write in the system
directories (e.g. ``/usr/bin``), we just use a subdirectory of our home
directory. The ``configure`` command configures the install, finding the
libraries and commands it needs. ``make`` compiles everything using the
options found by ``configure``. The final ``make install`` step copies
the compiled code into the correct places under ``$HOME/local``.

To be able to use the MPI compiler, you need to modify the ``PATH``
environment variable. To do this, run

.. code-block:: bash

    $ export PATH=$PATH:$HOME/local/bin

and add this to the end of your startup file ``$HOME/.bashrc``. If
you’re using CSH rather than BASH, the command is

.. code-block:: bash

    % setenv PATH ${PATH}:${HOME}/local/bin

and the startup file is ``$HOME/.cshrc``. You should now be able to run
``mpicc`` and so have a working MPI compiler.

.. _sec-libraries:

Installing libraries
--------------------

After getting an MPI compiler, the next step is to make sure the
libraries BOUT++ needs are installed. At minimum BOUT++ needs the FFTW-3
library, and to run any of the examples you’ll also need NetCDF-4 or
HDF5 installed.

Most large machines (e.g. NERSC Hopper, HECToR, HPC-FF etc.) will have
these libraries and many more already installed, but you may need to
load a module to use them. To see a list of the available modules, try
running

.. code-block:: bash

    modules avail

which works on many systems, but not all. See your system’s
documentation on modules and which ones to load. If you don’t know, or
modules don’t work, you can still install libraries in your home
directory by following the instructions below.

If you’re installing on your own machine, then install the packages for
your distribution. On Ubuntu or Debian, the necessary packages can be
installed by running

.. code-block:: bash

    $ sudo apt-get install libfftw3-dev libnetcdf-dev

The easiest way to test if the libraries are installed correctly is try
configuring BOUT++. In the ``BOUT`` directory obtained previously, run

.. code-block:: bash

    $ ./configure

If this finishes by printing a summary, and paths for IDL, Python, and
Octave, then the libraries are set up and you can skip to the next
section. If you see a message
“``ERROR: FFTW not found. Required by BOUT++``” then you need to install
FFTW-3. If you haven’t already, create directories “install” and “local”
in your home directory:

.. code-block:: bash

    $ cd
    $ mkdir install
    $ mkdir local

Download the latest stable version from
http://www.fftw.org/download.html into the “install” directory. At the
time of writing, this was called ``fftw-3.3.2.tar.gz``. Untar this file,
and ’cd’ into the resulting directory. As with the MPI compiler,
configure and install the FFTW library into ``$HOME/local`` by running:

.. code-block:: bash

    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

Go back to the ``BOUT`` directory and re-run the configure script. If
you used ``$HOME/local`` as the prefix, BOUT++ configure should find the
FFTW library now. If you installed somewhere else, you can specify the
directory with the ``–with-fftw=`` option:

.. code-block:: bash

    $ ./configure --with-fftw=$HOME/local

Configure should now find FFTW, and search for the NetCDF library. If
configure finishes successfully, then skip to the next section, but if
you see a message ``NetCDF support disabled`` then configure couldn’t
find the NetCDF library. Unless you have PACT or pnetcdf installed, this
will be followed by a message
``ERROR: At least one file format must be supported``.

The latest versions of NetCDF have separated out the C++ API from the
main C library. As a result, you will need to download and install both.
Download the latest versions of the NetCDF-C and NetCDF-4 C++ libraries
from http://www.unidata.ucar.edu/downloads/netcdf. As of
January 2017, these are versions 4.4.1.1 and 4.3.0 respectively.

Untar the file and ’cd’ into the resulting directory:

.. code-block:: bash

    $ tar -xzvf netcdf-4.4.1.1.tar.gz
    $ cd netcdf-4.4.1.1

As with MPI compilers and FFTW, configure, then make and make install:

.. code-block:: bash

    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

Sometimes configure can fail, in which case try disabling Fortran:

.. code-block:: bash

    $ ./configure --prefix=$HOME/local --disable-fortran
    $ make
    $ make install

Similarly for the C++ API:

.. code-block:: bash

    $ tar -xzvf netcdf-cxx4-4.3.0.tar.gz
    $ cd netcdf-cxx4-4.3.0
    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

You may need to set a couple of environment variables as well:

.. code-block:: bash

    $ export PATH=$HOME/local/bin:$PATH
    $ export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

You should check where NetCDF actually installed its libraries. On some
systems this will be ``$HOME/local/lib``, but on others it may be, e.g.
``$HOME/local/lib64``. Check which it is, and set ``$LD_LIBRARY_PATH``
appropriately.

Go back to the BOUT directory and run the configure script again, this
time specifying both the location of FFTW (if you installed it from
source above), and the NetCDF library:

.. code-block:: bash

    $ ./configure --with-fftw=$HOME/local --with-netcdf=$HOME/local

which should now finish successfully, printing a summary of the
configuration:

.. code-block:: bash

    Configuration summary
      FACETS support: no
      PETSc support: no
      SLEPc support: no
      IDA support: yes
      CVODE support: yes
      ARKODE support: yes
      NetCDF support: yes
      Parallel-NetCDF support: no
      HDF5 support: yes (parallel: no)
      Hypre support: no
      MUMPS support: no

If not, see :ref:`sec-advancedinstall` for some things you can try to
resolve common problems.

.. _sec-configanalysis:

Configuring analysis routines
-----------------------------

The BOUT++ installation comes with a set of useful routines which can be
used to prepare inputs and analyse outputs. Most of this code is in IDL,
but an increasing amount is in Python. In particular all the test suite
scripts use Python, so to run these you’ll need this configured. If you
just want to compile BOUT++ then you can skip to the next section, but
make a note of what configure printed out.

When the configure script finishes, it prints out the paths you need to
get IDL, Python, and Octave analysis routines working. After running the
command which looks like

.. code-block:: bash

    $ export IDL_PATH=...

check that ``idl`` can find the analysis routines by running:

.. code-block:: bash

    $ idl
    IDL> .r collect
    IDL> help, /source

You should see the function ``COLLECT`` in the ``BOUT/tools/idllib``
directory. If not, something is wrong with your ``IDL_PATH`` variable.
On some machines, modifying ``IDL_PATH`` causes problems, in which case
you can try modifying the path inside IDL by running

.. code-block:: bash

    IDL> !path = !path + ":/path/to/BOUT/tools/idllib"

where you should use the full path. You can get this by going to the
``tools/idllib`` directory and typing ``pwd``. Once this is done
you should be able to use ``collect`` and other routines.

To use Python, you will need the NumPy and SciPy libraries. On Debian or
Ubuntu these can be installed with

.. code-block:: bash

    $ sudo apt-get install python-scipy

which should then add all the other dependencies like NumPy. To test if
everything is installed, run

.. code-block:: bash

    $ python -c "import scipy"

If not, see the SciPy website http://www.scipy.org for instructions on
installing.

To do this, the path to ``tools/pylib`` should be added to the
``PYTHONPATH`` environment variable. Instructions for doing this are
printed at the end of the configure script, for example:

.. code-block:: bash

    Make sure that the tools/pylib directory is in your PYTHONPATH
    e.g. by adding to your ~/.bashrc file

       export PYTHONPATH=/home/ben/BOUT/tools/pylib/:$PYTHONPATH

To test if this command has worked, try running

.. code-block:: bash

    $ python -c "import boutdata"

If this doesn’t produce any error messages then Python is configured
correctly.

.. _sec-installbout:

Compiling BOUT++
----------------

Once BOUT++ has been configured, you can compile the bulk of the code by
going to the ``BOUT`` directory (same as ``configure``) and running

.. code-block:: bash

    $ make

(on OS-X, FreeBSD, and AIX this should be ``gmake``). This should print
something like:

.. code-block:: bash

    ----- Compiling BOUT++ -----
    CXX      =  mpicxx
    CFLAGS   =  -O -DCHECK=2 -DSIGHANDLE \
     -DREVISION=13571f760cec446d907e1bbeb1d7a3b1c6e0212a \
     -DNCDF -DBOUT_HAS_PVODE
    CHECKSUM =  ff3fb702b13acc092613cfce3869b875
    INCLUDE  =  -I../include
      Compiling  field.cxx
      Compiling  field2d.cxx

At the end of this, you should see a file ``libbout++.a`` in the
``lib/`` subdirectory of the BOUT++ distribution. If you get an error,
please send an error report to a BOUT++ developer such as
benjamin.dudson@york.ac.uk containing

-  Which machine you’re compiling on

-  The output from make, including full error message

-  The ``make.config`` file in the BOUT++ root directory

.. _sec-runtestsuite:

Running the test suite
----------------------

In the ``examples/`` subdirectory there are a set of short test cases
which are intended to test portions of the BOUT++ code and catch any
bugs which could be introduced. To run the test cases, the Python
libraries must first be set up by following the instructions in
section [sec:configanalysis]. Go into the ``examples`` subdirectory and
run

.. code-block:: bash

    $ ./test_suite

This will go through a set of tests, each on a variety of different
processors. **Note:** currently this uses the ``mpirun`` command to
launch the runs, so won’t work on machines which use a job submission
system like PBS or SGE.

These tests should all pass, but if not please create an issue on
Github containing:

-  Which machine you’re running on

-  The ``make.config`` file in the BOUT++ root directory

-  The ``run.log.*`` files in the directory of the test which failed

If the tests pass, congratulations! You have now got a working
installation of BOUT++. Unless you want to use some experimental
features of BOUT++, skip to section [sec:running] to start running the
code.

.. _sec-advancedinstall:

Advanced installation options
=============================

This section describes some common issues encountered when configuring
and compiling BOUT++, and how to configure optional libraries like
SUNDIALS and PETSc.

File formats
------------

BOUT++ can currently use two different file formats: NetCDF-4_, and
HDF5_ and experimental support for parallel flavours of both. NetCDF
is a widely used format and so has many more tools for viewing and
manipulating files. In particular, the NetCDF-4 library can produce
files in either NetCDF3 “classic” format, which is backwards-compatible
with NetCDF libraries since 1994 (version 2.3), or in the newer NetCDF4
format, which is based on (and compatible with) HDF5. HDF5 is another
widely used format. If you have multiple libraries installed then BOUT++
can use them simultaneously, for example reading in grid files in NetCDF
format, but writing output data in HDF5 format.

.. _NetCDF-4: http://www.unidata.ucar.edu/software/netcdf/
.. _HDF5: https://www.hdfgroup.org/HDF5/

To enable NetCDF support, you will need to install NetCDF version 4.0.1
or later. Note that although the NetCDF-4 library is used for the C++
interface, by default BOUT++ writes the “classic” format. Because of
this, you don’t need to install zlib or HDF5 for BOUT++ NetCDF support
to work. If you want to output to HDF5 then you need to first install
the zlib and HDF5 libraries, and then compile NetCDF with HDF5 support.
When NetCDF is installed, a script ``nc-config`` should be put into
somewhere on the path. If this is found then configure should have all
the settings it needs. If this isn’t found then configure will search
for the NetCDF include and library files.

SUNDIALS
--------

The BOUT++ distribution includes a 1998 version of CVODE (then called
PVODE) by Scott D. Cohen and Alan C. Hindmarsh, which is the default
time integration solver. Whilst no serious bugs have been found in this
code (as far as the authors are aware of), several features such as
user-supplied preconditioners and constraints cannot be used with this
solver. Currently, BOUT++ also supports the SUNDIALS solvers CVODE, IDA
and ARKODE which are available from
https://computation.llnl.gov/casc/sundials/main.html.

| SUNDIALS is only downloadable from the home page, as submitting your
  name and e-mail is required for the download. As for the date of this
  typing, SUNDIALS version :math:`2.6.2` is the newest. In order for a
  smooth install it is recommended to install SUNDIALS from an install
  directory. The full installation guide is found in the downloaded
  ``.tar.gz``, but we will provide a step-by-step guide to install it
  and make it compatible with BOUT++ here.
|  

.. code-block:: bash

    $ cd ~
    $ mkdir local
    $ cd local
    $ mkdir examples
    $ cd ..
    $ mkdir install
    $ cd install
    $ mkdir sundials-install
    $ cd sundials-install
    $ # Move the downloaded sundials-2.6.2.tar.gz to sundials-install
    $ tar -xzvf sundials-2.6.2.tar.gz
    $ mkdir build
    $ cd build

    $ cmake \
      -DCMAKE_INSTALL_PREFIX=$HOME/local \
      -DEXAMPLES_INSTALL_PATH=$HOME/local/examples \
      -DCMAKE_LINKER=$HOME/local/lib \
      -DLAPACK_ENABLE=ON \
      -DOPENMP_ENABLE=ON \
      -DMPI_ENABLE=ON \
    $ ../sundials-2.6.2

    $ make
    $ make install

The SUNDIALS IDA solver is a Differential-Algebraic Equation (DAE)
solver, which evolves a system of the form
:math:`\mathbf{f}(\mathbf{u},\dot{\mathbf{u}},t) = 0`. This allows
algebraic constraints on variables to be specified.

To configure BOUT++ with SUNDIALS only (see section [sec:PETSc] on how
to build PETSc with SUNDIALS), go to the root directory of BOUT++ and
type

.. code-block:: bash

    $ ./configure --with-sundials

SUNDIALS will allow you to select at run-time which solver to use. See
:ref:`sec-timeoptions` for more details on how to do this.

PETSc
-----

BOUT++ can use PETSc http://www.mcs.anl.gov/petsc/ for time-integration
and for solving elliptic problems, such as inverting Poisson and
Helmholtz equations.

Currently, BOUT++ supports PETSc version :math:`3.1`, :math:`3.2`,
:math:`3.3` and :math:`3.4` (support for newer versions are planned for
the future). To install PETSc version :math:`3.4.5`, use the following
steps

.. code-block:: bash

    $ cd ~
    $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.5.tar.gz
    $ tar -xzvf petsc-3.4.5.tar.gz
    $ # Optional
    $ # rm petsc-3.4.5.tar.gz
    $ cd petsc-3.4.5

To build PETSc without SUNDIALS, configure with

.. code-block:: bash

    $ ./configure \
      --with-clanguage=cxx \
      --with-mpi=yes \
      --with-precision=double \
      --with-scalar-type=real \
      --with-shared-libraries=0

Add ``–with-debugging=yes`` to ``./configure`` in order to allow
debugging.

| To build PETSc with SUNDIALS, install SUNDIALS as explained in section
  :ref:`sec-sundials`, and append ``./configure`` with
  ``–with-sundials-dir=$HOME/local``
|  
| It is also possible to get PETSc to download and install MUMPS (see
  :ref:`sec-MUMPS`), by adding

.. code-block:: bash

    --download-mumps \
    --download-scalapack \
    --download-blacs \
    --download-f-blas-lapack=1 \
    --download-parmetis \
    --download-ptscotch \
    --download-metis

to ``./configure`` To make PETSc, type

.. code-block:: bash

    $ make PETSC_DIR=$HOME/petsc-3.4.5 PETSC_ARCH=arch-linux2-cxx-debug all

Should blas, lapack or any other packages be missing, you will get an
error, and a suggestion that you can append
``–download-name-of-package`` to the ``./configure`` line. You may want
to test that everything is configured properly. To do this, type

.. code-block:: bash

    $ make PETSC_DIR=$HOME/petsc-3.4.5 PETSC_ARCH=arch-linux2-cxx-debug test

To configure BOUT++ with PETSc, go to the BOUT++ root directory, and
type

.. code-block:: bash

    $ ./configure --with-petsc=$HOME/petsc-3.4.5

To configure BOUT++ with PETSc and sundials, type instead

.. code-block:: bash

    $ ./configure --with-petsc=$HOME/petsc-3.4.5 --with-sundials

Finally compile PETSc:

.. code-block:: bash

    $ make

To use PETSc, you have to define the variable ``PETSC_DIR`` to point to
the petsc directory, type

.. code-block:: bash

    $ export PETSC_DIR=$HOME/petsc-3.4.5

and add to your startup file ``$HOME/.bashrc``

.. code-block:: bash

    $ export PETSC_DIR=$HOME/petsc-3.4.5

LAPACK
------

BOUT++ comes with linear solvers for tridiagonal and band-diagonal
systems, but these are not particularly optimised and are in any case
descended from Numerical Recipes code (hence NOT covered by LGPL
license).

To replace these routines, BOUT++ can use the LAPACK library. This is
however written in FORTRAN 77, which can cause linking headaches. To
enable these routines use

.. code-block:: bash

    $ ./configure --with-lapack

and to specify a non-standard path

.. code-block:: bash

    $ ./configure --with-lapack=/path/to/lapack

.. _sec-mumps:

MUMPS
-----

This is still experimental, but does work on at least some systems at
York. The PETSc library can be used to call MUMPS for directly solving
matrices (e.g. for Laplacian inversions), or MUMPS can be used directly.
To enable MUMPS, configure with

.. code-block:: bash

    $ ./configure --with-mumps

MUMPS has many dependencies, including ScaLapack and ParMetis, which the
configuration script assumes are in the same place as MUMPS. The easiest
way to get MUMPS installed is to install PETSc with MUMPS, as the
configuration script will check the PETSc directory.

MPI compilers
-------------

These are usually called something like mpicc and mpiCC (or mpicxx), and
the configure script will look for several common names. If your
compilers aren’t recognised then set them using

.. code-block:: bash

    $ ./configure MPICC=<your C compiler> MPICXX=<your C++ compiler>

NOTES:

-  On LLNL’s Grendel, mpicxx is broken. Use mpiCC instead by passing
   “MPICXX=mpiCC” to configure. Also need to specify this to NetCDF
   library by passing “CXX=mpiCC” to NetCDF configure.

Issues
------

Wrong install script
~~~~~~~~~~~~~~~~~~~~

Before installing, make sure the correct version of ``install`` is being
used by running

.. code-block:: bash

     $ which install

This should point to a system directory like ``/usr/bin/install``.
Sometimes when IDL has been installed, this points to the IDL install
(e.g. something like ``/usr/common/usg/idl/idl70/bin/install`` on
Franklin). A quick way to fix this is to create a link from your local
bin to the system install:

.. code-block:: bash

     $ ln -s /usr/bin/install $HOME/local/bin/

“which install” should now print the install in your local bin
directory.

Compiling cvode.cxx fails
~~~~~~~~~~~~~~~~~~~~~~~~~

Occasionally compiling the CVODE solver interface will fail with an
error similar to:

.. code-block:: bash

    cvode.cxx: In member function ‘virtual int CvodeSolver::init(rhsfunc, bool, int, BoutR...
    cvode.cxx:234:56: error: invalid conversion from ‘int (*)(CVINT...
    ...

This is caused by different sizes of ints used in different versions of
the CVODE library. The configure script tries to determine the correct
type to use, but may fail in unusual circumstances. To fix, edit
``src/solver/impls/cvode/cvode.cxx``, and change line 48 from

.. code-block:: cpp

    typedef int CVODEINT;

to

.. code-block:: cpp

    typedef long CVODEINT;

.. _sec-running:

Running BOUT++
==============

The ``examples/`` directory contains some test cases for a variety of
fluid models. The ones starting ``test-`` are short tests, which often
just run a part of the code rather than a complete simulation. The
simplest example to start with is ``examples/conduction/``. This solves
a single equation for a 3D scalar field :math:`T`:

.. math::

   \frac{\partial T}{\partial t} = \nabla_{||}(\chi\partial_{||} T)

There are several files involved:

-  ``conduction.cxx`` contains the source code which specifies the
   equation to solve

-  ``conduct_grid.nc`` is the grid file, which in this case just
   specifies the number of grid points in :math:`X` and :math:`Y`
   (``nx`` & ``ny``) with everything else being left as the default
   (e.g. grid spacings dx and dy are :math:`1`, the metric tensor is the
   identity matrix). For details of the grid file format, see
   :ref:`sec-gridgen`.

-  ``generate.py`` is a Python script to create the grid file. In this
   case it just writes nx and ny

-  ``data/BOUT.inp`` is the settings file, specifying how many output
   timesteps to take, differencing schemes to use, and many other
   things. In this case it’s mostly empty so the defaults are used.

First you need to compile the example:

.. code-block:: bash

    $ gmake

which should print out something along the lines of

.. code-block:: bash

      Compiling  conduction.cxx
      Linking conduction

If you get an error, most likely during the linking stage, you may need
to go back and make sure the libraries are all set up correctly. A
common problem is mixing MPI implementations, for example compiling
NetCDF using Open MPI and then BOUT++ with MPICH2. Unfortunately the
solution is to recompile everything with the same compiler.

Then try running the example. If you’re running on a standalone server,
desktop or laptop then try:

.. code-block:: bash

    $ mpirun -np 2 ./conduction

If you’re running on a cluster or supercomputer, you should find out how
to submit jobs. This varies, but usually on these bigger machines there
will be a queueing system and you’ll need to use ``qsub``, ``msub``,
``llsubmit`` or similar to submit jobs.

When the example runs, it should print a lot of output. This is
recording all the settings being used by the code, and is also written
to log files for future reference. The test should take a few seconds to
run, and produce a bunch of files in the ``data/`` subdirectory.

-  ``BOUT.log.*`` contains a log from each process, so because we ran
   with “-np 2” there should be 2 logs. The one from processor :math:`0`
   will be the same as what was printed to the screen. This is mainly
   useful because if one process crashes it may only put an error
   message into its own log.

-  ``BOUT.restart.*.nc`` are the restart files for the last time point.
   Currently each processor saves its own state in a separate file, but
   there is experimental support for parallel I/O. For the settings, see
   :ref:`sec-iooptions`.

-  ``BOUT.dmp.*.nc`` contain the output data, including time history. As
   with the restart files, each processor currently outputs a separate
   file.

Restart files allow the run to be restarted from where they left off:

.. code-block:: bash

     $ mpirun -np 2 ./conduction restart

This will delete the output data ``BOUT.dmp.*.nc`` files, and start
again. If you want to keep the output from the first run, add “append”

.. code-block:: bash

     $ mpirun -np 2 ./conduction restart append

which will then append any new outputs to the end of the old data files.
For more information on restarting, see :ref:`sec-restarting`.

To analyse the output of the simulation, cd into the ``data``
subdirectory and start IDL. If you don’t have IDL, don’t panic as all
this is also possible in Python and discussed in
:ref:`sec-pythonroutines`. First, list the variables in one of the
data files:

.. code-block:: idl

    IDL> print, file_list("BOUT.dmp.0.nc")
    iteration MXSUB MYSUB MXG MYG MZ NXPE NYPE BOUT_VERSION t_array ZMAX ZMIN T

All of these except ’\ ``T``\ ’ are in all output files, and they
contain information about the layout of the mesh so that the data can be
put in the correct place. The most useful variable is ’\ ``t_array``\ ’
which is a 1D array of simulation output times. To read this, we can use
the ``collect`` function:

.. code-block:: idl

    IDL> time = collect(var="t_array")
    IDL> print, time
          1.10000      1.20000      1.30000      1.40000      1.50000 ...

The number of variables in an output file depends on the model being
solved, which in this case consists of a single scalar field
’\ ``T``\ ’. To read this into IDL, again use ``collect``:

.. code-block:: idl

    IDL> T = collect(var="T")
    IDL> help, T
    T               FLOAT     = Array[5, 64, 1, 20]

This is a 4D variable, arranged as ``[x, y, z, t]``. The :math:`x`
direction has 5 points, consisting of 2 points either side for the
boundaries and one point in the middle which is evolving. This case is
only solving a 1D problem in :math:`y` with 64 points so to display an
animation of this

.. code-block:: idl

    IDL> showdata, T[2,*,0,*]

which selects the only evolving :math:`x` point, all :math:`y`, the only
:math:`z` point, and all time points. If given 3D variables, showdata
will display an animated surface

.. code-block:: idl

    IDL> showdata, T[*,*,0,*]

and to make this a coloured contour plot

.. code-block:: idl

    IDL> showdata, T[*,*,0,*], /cont

The equivalent commands in Python are as follows. To print a list of
variables in a file:

.. code-block:: pycon

    >>> from boututils.datafile import DataFile
    >>> DataFile("BOUT.dmp.0.nc").list()

To collect a variable,

.. code-block:: pycon

    >>> from boutdata.collect import collect
    >>> T = collect("T")
    >>> T.shape

Note that the order of the indices is different in Python and IDL: In
Python, 4D variables are arranged as ``[t, x, y, z]``. To show an
animation

.. code-block:: pycon

    >>> from boututils.showdata import showdata
    >>> showdata(T[:,:,:,0])

The next example to look at is ``test-wave``, which is solving a wave
equation using

.. math::

   \frac{\partial f}{\partial t} = \partial_{||} g \qquad \frac{\partial g}{\partial t} = \partial_{||} f

using two different methods. Other examples contain two scripts: One
for running the example and then an IDL script to plot the results:

.. code-block:: bash

    ./runcase.sh
    idl runidl.pro

Assuming these examples work (which they should), looking through the
scripts and code may give you an idea of how BOUT++ works. More
information on setting up and running BOUT++ is given in
:ref:`sec-running`, and details of analysing the results using IDL
are given in :ref:`sec-output`.

Alternatively, one can run BOUT++ with the python wrapper
``bout_runners`` , as explained in section [sec:bout\_runners]. Examples
of using ``bout_runners`` can be found in
``examples/bout_runners_example``.

When things go wrong
--------------------

BOUT++ is still under development, and so occasionally you may be lucky
enough to discover a new bug. This is particularly likely if you’re
modifying the physics module source code (see :ref:`sec-equations`)
when you need a way to debug your code too.

-  Check the end of each processor’s log file (tail data/BOUT.log.\*).
   When BOUT++ exits before it should, what is printed to screen is just
   the output from processor 0. If an error occurred on another
   processor then the error message will be written to it’s log file
   instead.

-  By default when an error occurs a kind of stack trace is printed
   which shows which functions were being run (most recent first). This
   should give a good indication of where an error occurred. If this
   stack isn’t printed, make sure checking is set to level 2 or higher
   (``./configure –with-checks=2``)

-  If the error is a segmentation fault, you can try a debugger such as
   totalview

-  If the error is due to non-finite numbers, increase the checking
   level (``./configure –with-checks=3``) to perform more checking of
   values and (hopefully) find an error as soon as possible after it
   occurs.

Startup output
--------------

When BOUT++ is run, it produces a lot of output initially, mainly
listing the options which have been used so you can check that it’s
doing what you think it should be. It’s generally a good idea to scan
over this see if there are any important warnings or errors. Each
processor outputs its own log file ``BOUT.log.#`` and the log from
processor 0 is also sent to the screen. This output may look a little
different if it’s out of date, but the general layout will probably be
the same.

First comes the introductory blurb:

.. code-block:: bash

    BOUT++ version 1.0
    Revision: c8794400adc256480f72c651dcf186fb6ea1da49
    MD5 checksum: 8419adb752f9c23b90eb50ea2261963c
    Code compiled on May 11 2011 at 18:22:37

    B.Dudson (University of York), M.Umansky (LLNL) 2007
    Based on BOUT by Xueqiao Xu, 1999

The version number (1.0 here) gets increased occasionally after some
major feature has been added. To help match simulations to code
versions, the Git revision of the core BOUT++ code and the date and time
it was compiled is recorded. Because code could be modified from the
revision, an MD5 checksum of all the code is also calculated. This
information makes it possible to verify precisely which version of the
code was used for any given run.

Next comes the compile-time options, which depend on how BOUT++ was
configured (see :ref:`sec-installbout`)

.. code-block:: bash

    Compile-time options:
            Checking enabled, level 2
            Signal handling enabled
            netCDF support enabled
            Parallel NetCDF support disabled

This says that some run-time checking of values is enabled, that the
code will try to catch segmentation faults to print a useful error, that
NetCDF files are supported, but that the parallel flavour isn’t.

The processor number comes next:

.. code-block:: bash

    Processor number: 0 of 1

This will always be processor number ’0’ on screen as only the output
from processor ’0’ is sent to the terminal. After this the core BOUT++
code reads some options:

.. code-block:: bash

            Option /nout = 50 (data/BOUT.inp)
            Option /timestep = 100 (data/BOUT.inp)
            Option /grid = slab.6b5.r1.cdl (data/BOUT.inp)
            Option /dump_float = true   (default)
            Option /non_uniform = false (data/BOUT.inp)
            Option /restart = false  (default)
            Option /append = false  (default)
            Option /dump_format = nc (data/BOUT.inp)
            Option /StaggerGrids = false  (default)

This lists each option and the value it has been assigned. For every
option the source of the value being used is also given. If a value had
been given on the command line then ``(command line)`` would appear
after the option.

.. code-block:: bash

    Setting X differencing methods
            First       :  Second order central (C2)
            Second      :  Second order central (C2)
            Upwind      :  Third order WENO (W3)
            Flux        :  Split into upwind and central (SPLIT)
    Setting Y differencing methods
            First       :  Fourth order central (C4)
            Second      :  Fourth order central (C4)
            Upwind      :  Third order WENO (W3)
            Flux        :  Split into upwind and central (SPLIT)
    Setting Z differencing methods
            First       :  FFT (FFT)
            Second      :  FFT (FFT)
            Upwind      :  Third order WENO (W3)
            Flux        :  Split into upwind and central (SPLIT)

This is a list of the differential methods for each direction. These are
set in the BOUT.inp file (``[ddx]``, ``[ddy]`` and ``[ddz]`` sections),
but can be overridden for individual operators. For each direction,
numerical methods can be specified for first and second central
difference terms, upwinding terms of the form
:math:`{{\frac{\partial f}{\partial t}}} = {{\boldsymbol{v}}}\cdot\nabla f`,
and flux terms of the form
:math:`{{\frac{\partial f}{\partial t}}} = \nabla\cdot({{\boldsymbol{v}}}f)`.
By default the flux terms are just split into a central and an upwinding
term.

In brackets are the code used to specify the method in BOUT.inp. A list
of available methods is given in :ref:`sec-diffmethod`.

.. code-block:: bash

    Setting grid format
            Option /grid_format =  (default)
            Using NetCDF format for file 'slab.6b5.r1.cdl'
    Loading mesh
            Grid size: 10 by 64
            Option /mxg = 2 (data/BOUT.inp)
            Option /myg = 2 (data/BOUT.inp)
            Option /NXPE = 1 (default)
            Option /mz = 65 (data/BOUT.inp)
            Option /twistshift = false (data/BOUT.inp)
            Option /TwistOrder = 0 (default)
            Option /ShiftOrder = 0 (default)
            Option /shiftxderivs = false (data/BOUT.inp)
            Option /IncIntShear = false  (default)
            Option /BoundaryOnCell = false  (default)
            Option /StaggerGrids = false  (default)
            Option /periodicX = false  (default)
            Option /async_send = false  (default)
            Option /zmin = 0 (data/BOUT.inp)
            Option /zmax = 0.0028505 (data/BOUT.inp)

.. code-block:: bash

    WARNING: Number of inner y points 'ny_inner' not found. Setting to 32

Optional quantities (such as ``ny_inner`` in this case) which are not
specified are given a default (best-guess) value, and a warning is
printed.

.. code-block:: bash

            EQUILIBRIUM IS SINGLE NULL (SND)
            MYPE_IN_CORE = 0
            DXS = 0, DIN = -1. DOUT = -1
            UXS = 0, UIN = -1. UOUT = -1
            XIN = -1, XOUT = -1
            Twist-shift:

At this point, BOUT++ reads the grid file, and works out the topology of
the grid, and connections between processors. BOUT++ then tries to read
the metric coefficients from the grid file:

.. code-block:: bash

            WARNING: Could not read 'g11' from grid. Setting to 1.000000e+00
            WARNING: Could not read 'g22' from grid. Setting to 1.000000e+00
            WARNING: Could not read 'g33' from grid. Setting to 1.000000e+00
            WARNING: Could not read 'g12' from grid. Setting to 0.000000e+00
            WARNING: Could not read 'g13' from grid. Setting to 0.000000e+00
            WARNING: Could not read 'g23' from grid. Setting to 0.000000e+00

These warnings are printed because the coefficients have not been
specified in the grid file, and so the metric tensor is set to the
default identity matrix.

.. code-block:: bash

            WARNING: Could not read 'zShift' from grid. Setting to 0.000000e+00
            WARNING: Z shift for radial derivatives not found

To get radial derivatives, the quasi-ballooning coordinate method is
used . The upshot of this is that to get radial derivatives,
interpolation in Z is needed. This should also always be set to FFT.

.. code-block:: bash

            WARNING: Twist-shift angle 'ShiftAngle' not found. Setting from zShift
            Option /twistshift_pf = false  (default)

.. code-block:: bash

            Maximum error in diagonal inversion is 0.000000e+00
            Maximum error in off-diagonal inversion is 0.000000e+00

If only the contravariant components (``g11`` etc.) of the metric tensor
are specified, the covariant components (``g_11`` etc.) are calculated
by inverting the metric tensor matrix. Error estimates are then
calculated by calculating :math:`g_{ij}g^{jk}` as a check. Since no
metrics were specified in the input, the metric tensor was set to the
identity matrix, making inversion easy and the error tiny.

.. code-block:: bash

            WARNING: Could not read 'J' from grid. Setting to 0.000000e+00
            WARNING: Jacobian 'J' not found. Calculating from metric tensor

.. code-block:: bash

            Maximum difference in Bxy is 1.444077e-02
    Calculating differential geometry terms
            Communicating connection terms
    Boundary regions in this processor: core, sol, target, target,
            done

.. code-block:: bash

    Setting file formats
            Using NetCDF format for file 'data/BOUT.dmp.0.nc'

The laplacian inversion code is initialised, and prints out the options
used.

.. code-block:: bash

    Initialising Laplacian inversion routines
            Option comms/async = true   (default)
            Option laplace/filter = 0.2 (default)
            Option laplace/low_mem = false  (default)
            Option laplace/use_pdd = false  (default)
            Option laplace/all_terms = false  (default)
            Option laplace/laplace_nonuniform = false  (default)
            Using serial algorithm
            Option laplace/max_mode = 26 (default)

After this comes the physics module-specific output:

.. code-block:: bash

    Initialising physics module
            Option solver/type =  (default)
            .
            .
            .

This typically lists the options used, and useful/important
normalisation factors etc.

Finally, once the physics module has been initialised, and the current
values loaded, the solver can be started

.. code-block:: bash

    Initialising solver
            Option /archive = -1 (default)
            Option /dump_format = nc (data/BOUT.inp)
            Option /restart_format = nc (default)
            Using NetCDF format for file 'nc'

.. code-block:: bash

    Initialising PVODE solver
            Boundary region inner X
            Boundary region outer X
            3d fields = 2, 2d fields = 0 neq=84992, local_N=84992

This last line gives the number of equations being evolved (in this case
84992), and the number of these on this processor (here 84992).

.. code-block:: bash

            Option solver/mudq = 16 (default)
            Option solver/mldq = 16 (default)
            Option solver/mukeep = 0 (default)
            Option solver/mlkeep = 0 (default)

The absolute and relative tolerances come next:

.. code-block:: bash

            Option solver/atol = 1e-10 (data/BOUT.inp)
            Option solver/rtol = 1e-05 (data/BOUT.inp)

.. code-block:: bash

            Option solver/use_precon = false  (default)
            Option solver/precon_dimens = 50 (default)
            Option solver/precon_tol = 0.0001 (default)
            Option solver/mxstep = 500 (default)

.. code-block:: bash

            Option fft/fft_measure = false  (default)

This next option specifies the maximum number of internal timesteps
which CVODE will take between outputs.

.. code-block:: bash

            Option fft/fft_measure = false  (default)
    Running simulation

    Run started at  : Wed May 11 18:23:20 2011

            Option /wall_limit = -1 (default)

Per-timestep output
-------------------

At the beginning of a run, just after the last line in the previous
section, a header is printed out as a guide

.. code-block:: bash

    Sim Time  |  RHS evals  | Wall Time |  Calc    Inv   Comm    I/O   SOLVER

Each timestep (the one specified in BOUT.inp, not the internal
timestep), BOUT++ prints out something like

.. code-block:: bash

    1.001e+02         76       2.27e+02    87.1    5.3    1.0    0.0    6.6

This gives the simulation time; the number of times the time-derivatives
(RHS) were evaluated; the wall-time this took to run, and percentages
for the time spent in different parts of the code.

-  ``Calc`` is the time spent doing calculations such as
   multiplications, derivatives etc

-  ``Inv`` is the time spent in inversion code (i.e. inverting
   Laplacians), including any communication which may be needed to do
   the inversion.

-  ``Comm`` is the time spent communicating variables (outside the
   inversion routine)

-  ``I/O`` is the time spent writing dump and restart files to disk.
   Most of the time this should not be an issue

-  ``SOLVER`` is the time spent in the implicit solver code.

The output sent to the terminal (not the log files) also includes a run
time, and estimated remaining time.

.. _sec-restarting:

Restarting runs
---------------

Every output timestep, BOUT++ writes a set of files named
“BOUT.restart.#.nc” where ’#’ is the processor number (for parallel
output, a single file “BOUT.restart.nc” is used). To restart from where
the previous run finished, just add the keyword **restart** to the end
of the command, for example:

.. code-block:: bash

     $ mpirun -np 2 ./conduction restart

Equivalently, put “restart=true” near the top of the BOUT.inp input
file. Note that this will overwrite the existing data in the
“BOUT.dmp.\*.nc” files. If you want to append to them instead then add
the keyword append to the command, for example:

.. code-block:: bash

     $ mpirun -np 2 ./conduction restart append

or also put “append=true” near the top of the BOUT.inp input file.

When restarting simulations BOUT++ will by default output the initial
state, unless appending to existing data files when it will not output
until the first timestep is completed. To override this behaviour, you
can specify the option dump\_on\_restart manually. If dump\_on\_restart
is true then the initial state will always be written out, if false then
it never will be (regardless of the values of restart and append).

If you need to restart from a different point in your simulation, or the
BOUT.restart files become corrupted, you can either use archived restart
files, or create new restart files. Archived restart files have names
like “BOUT.restart\_0020.#.nc”, and are written every 20 outputs by
default. To change this, set “archive” in the BOUT.inp file. To use
these files, they must be renamed to “BOUT.restart.#.nc”. A useful tool
to do this is “rename”:

.. code-block:: bash

    $ rename 's/_0020//' *.nc

will strip out “\_0020” from any file names ending in “.nc”.

If you don’t have archived restarts, or want to start from a different
time-point, there are Python routines for creating new restart files. If
your PYTHONPATH environment variable is set up (see
:ref:`sec-configanalysis`) then you can use the
``boutdata.restart.create`` function in
``tools/pylib/boutdata/restart.py``:

.. code-block:: pycon

    >>> from boutdata.restart import create
    >>> create(final=10, path='data', output='.')

The above will take time point 10 from the BOUT.dmp.\* files in the
“data” directory. For each one, it will output a BOUT.restart file in
the output directory “.”.

Makefiles and compiling BOUT++
==============================

BOUT++ has its own makefile system. These can be used to

#. Write an example or executable (see :ref:`sec-executables`)

#. Add a feature to BOUT++ (see :ref:`sec-modules`)

In all makefiles, ``BOUT_TOP`` is required!

These makefiles are sufficient for most uses, but for more complicated,
an executable script ``bout-config`` can be used to get the compilation
flags (see section [sec:bout-config]).

.. _sec-executables:

Executables example
-------------------

If writing an example (or physics module that executes) then the
makefile is very simple:

.. code-block:: makefile

    BOUT_TOP        = ../..

    SOURCEC         = <filename>.cxx

    include $(BOUT_TOP)/make.config

where ``BOUT_TOP`` - refers to the relative (or absolute) location of
the BOUT directory (the one that includes ``/lib`` and ``/src``) and
``SOURCEC`` is the name of your file, e.g. ``gas_compress.cxx``.

Optionally, it is possible to specify ``TARGET`` which defines what the
executable should be called (e.g. if you have multiple source files).
That’s it!

Multiple subdirectories
~~~~~~~~~~~~~~~~~~~~~~~

Large physics modules can have many files, and it can be helpful to
organise these into subdirectories. An example of how to do this is in
``examples/make_subdir``.

In the top level, list the directories

.. code-block:: makefile

    DIRS = fuu bar

In the makefile in each subdirectory, specify

.. code-block:: makefile

    TARGET = sub

then specify the path to the top-level directory

.. code-block:: makefile

    MODULE_DIR = ..

and the name of the subdirectory that the makefile is in

.. code-block:: makefile

    SUB_NAME = fuu

.. _sec-modules:

Modules example
---------------

If you are writing a new module (or concrete implementation) to go into
the BOUT++ library, then it is again pretty simple

.. code-block:: makefile

    BOUT_TOP = ../..

    SOURCEC         = communicator.cxx difops.cxx geometry.cxx grid.cxx \
                      interpolation.cxx topology.cxx
    SOURCEH         = $(SOURCEC:%.cxx=%.h)
    TARGET          = lib

    include $(BOUT_TOP)/make.config

``TARGET`` - must be ``lib`` to signify you are adding to
``libbout++.a``.

The other variables should be pretty self explanatory.

Adding a new subdirectory to ’src’
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No worries, just make sure to edit ``src/makefile`` to add it to the
``DIRS`` variable.

bout-config script
------------------

The ``bout-config`` script is in the ``bin`` subdirectory of the BOUT++
distribution, and is generated by ``configure``. This script can be used
to get the compilers, flags and settings to compile BOUT++. To get a
list of available options:

.. code-block:: bash

    $ bout-config --help

so to get the library linking flags, for example

.. code-block:: bash

    $ bout-config --libs

This script can be used in makefiles to compile BOUT++ alongside other
libraries. An example is in ``examples/make-script``.

.. _sec-output:

Output and post-processing
==========================

The majority of the existing analysis and post-processing code is
written in IDL. The directory ``idllib`` contains many useful routines
for reading output files and analysing data. A summary of available IDL
routines is given in Appendix [apx:idl\_routines].

Post-processing using Python is also possible, and there are some
modules in the ``pylib`` directory, and a list of routines in
Appendix [apx:py\_routines]. This is a more recent addition, and so is
not yet as developed as the IDL support.

Reading BOUT++ output into IDL
------------------------------

There are several routines provided for reading data from BOUT++ output
into IDL. In the directory containing the BOUT++ output files (usually
``data/``), you can list the variables available using

.. code-block:: idl

    IDL> print, file_list("BOUT.dmp.0.nc")
    Ajpar Apar BOUT_VERSION MXG MXSUB MYG MYSUB MZ NXPE NYPE Ni Ni0 Ni_x Te0 Te_x
    Ti0 Ti_x ZMAX ZMIN iteration jpar phi rho rho_s t_array wci

The ``file_list`` procedure just returns an array, listing all the
variables in a given file.

One thing new users can find confusing is that different simulations may
have very different outputs. This is because **BOUT++ is not a single
physics model**: the variables evolved and written to file are
determined by the model, and will be very different between (for
example) full MHD and reduced Braginskii models. There are however some
variables which all BOUT++ output files contain:

-  ``BOUT_VERSION``, which gives the version number of BOUT++ which
   produced the file. This is mainly to help output processing codes
   handle changes to the output file format. For example, BOUT++ version
   0.30 introduced 2D domain decomposition which needs to be handled
   when collecting data.

-  ``MXG``,\ ``MYG``. These are the sizes of the X and Y guard cells

-  ``MXSUB``, the number of X grid points in each processor. This does
   not include the guard cells, so the total X size of each field will
   be ``MXSUB + 2*MXG``.

-  ``MYSUB``, the number of Y grid points per processor (like MXSUB)

-  ``MZ``, the number of Z points

-  ``NXPE, NYPE``, the number of processors in the X and Y directions.
   ``NXPE * MXSUB + 2*MXG= NX``, ``NYPE * MYSUB = NY``

-  ``ZMIN``, ``ZMAX``, the range of Z in fractions of :math:`2\pi`.

-  ``iteration``, the last timestep in the file

-  ``t_array``, an array of times

Most of these - particularly those concerned with grid size and
processor layout - are used by post-processing routines such as
``collect``, and are seldom needed directly. To read a single variable
from a file, there is the ``file_read`` function:

.. code-block:: idl

    IDL> wci = file_read("BOUT.dmp.0.nc", "wci")
    IDL> print, wci
      9.58000e+06

To read in all the variables in a file into a structure, use the
``file_import`` function:

.. code-block:: idl

    IDL> d = file_import("BOUT.dmp.0.nc")
    IDL> print, d.wci
      9.58000e+06

This is often used to read in the entire grid file at once. Doing this
for output data files can take a long time and use a lot of memory.

Reading from individual files is fine for scalar quantities and time
arrays, but reading arrays which are spread across processors (i.e.
evolving variables) is tedious to do manually. Instead, there is the
``collect`` function to automate this:

.. code-block:: idl

    IDL> ni = collect(var="ni")
    Variable 'ni' not found
    -> Variables are case-sensitive: Using 'Ni'
    Reading from .//BOUT.dmp.0.nc: [0-35][2-6] -> [0-35][0-4]

This function takes care of the case, so that reading “ni” is
automatically corrected to “Ni”. The result is a 4D variable:

.. code-block:: idl

    IDL> help, ni
    NI              FLOAT     = Array[36, 5, 64, 400]

with the indices ``[X, Y, Z, T]``. Note that in the output files, these
variables are stored in ``[T, X, Y, Z]`` format instead but this is
changed by ``collect``. Sometimes you don’t want to read in the entire
array (which may be very large). To read in only a subset, there are
several optional keywords with ``[min,max]`` ranges:

.. code-block:: idl

    IDL> ni = collect(var="Ni", xind=[10,20], yind=[2,2], zind=[0,31],
    tind=[300,399])
    Reading from .//BOUT.dmp.0.nc: [10-20][4-4] -> [10-20][2-2]
    IDL> help, ni
    NI              FLOAT     = Array[11, 1, 32, 100]

Summary of IDL file routines
----------------------------

Functions file\_ can currently only read/write NetCDF files. HDF5 is not
supported yet.

Open a NetCDF file:

.. code-block:: idl

    handle = file_open("filename", /write, /create)

Array of variable names:

.. code-block:: idl

    list = file_list(handle)
    list = file_list("filename")

Number of dimensions:

.. code-block:: idl

    nd = file_ndims(handle, "variable")
    nd = file_ndims("filename", "variable")

Read a variable from file. Inds = [xmin, xmax, ymin, ymax, ...]

.. code-block:: idl

    data = file_read(handle, "variable", inds=inds)
    data = file_read("filename", "variable", inds=inds)

Write a variable to file. For NetCDF it tries to match up dimensions,
and defines new dimensions when needed

.. code-block:: idl

    status = file_write(handle, "variable", data)

Close a file after use

.. code-block:: idl

    file_close, handle

To read in all the data in a file into a structure:

.. code-block:: idl

    data = file_import("filename")

and to write a structure to file:

.. code-block:: idl

    status = file_export("filename", data)

IDL analysis routines
---------------------

Now that the BOUT++ results have been read into IDL, all the usual
analysis and plotting routines can be used. In addition, there are many
useful routines included in the ``idllib`` subdirectory. There is a
``README`` file which describes what each of these routines, but some of
the most useful ones are listed here. All these examples assume there is
a variable ``P`` which has been read into IDL as a 4D [x,y,z,t]
variable:

-  ``fft_deriv`` and ``fft_integrate`` which differentiate and integrate
   periodic functions.

-  ``get_integer``, ``get_float``, and ``get_yesno`` request integers,
   floats and a yes/no answer from the user respectively.

-  ``showdata`` animates 1 or 2-dimensional variables. Useful for
   quickly displaying results in different ways. This is useful for
   taking a quick look at the data, but can also produce bitmap outputs
   for turning into a movie for presentation. To show an animated
   surface plot at a particular poloidal location (32 here):

   .. code-block:: idl

       IDL> showdata, p[*,32,*,*]

   To turn this into a contour plot,

   .. code-block:: idl

       IDL> showdata, p[*,32,*,*], /cont

   To show a slice through this at a particular toroidal location (0
   here):

   .. code-block:: idl

       IDL> showdata, p[*,32,0,*]

   There are a few other options, and ways to show data using this code;
   see the README file, or comments in ``showdata.pro``. Instead of
   plotting to screen, showdata can produce a series of numbered bitmap
   images by using the ``bmp`` option

   .. code-block:: idl

       IDL> showdata, p[*,32,*,*], /cont, bmp="result_"

   which will produce images called ``result_0000.bmp``,
   ``result_0001.bmp`` and so on. Note that the plotting should not be
   obscured or minimised, since this works by plotting to screen, then
   grabbing an image of the resulting plot.

-  ``moment_xyzt`` takes a 4D variable (such as those from ``collect``),
   and calculates RMS, DC and AC components in the Z direction.

-  ``safe_colors`` A general routine for IDL which arranges the color
   table so that colors are numbered 1 (black), 2 (red), 3 (green), 4
   (blue). Useful for plotting, and used by many other routines in this
   library.

There are many other useful routines in the ``idllib`` directory. See
the ``idllib/README`` file for a short description of each one.

.. _sec-pythonroutines:

Python routines
---------------

There are several modules available for reading NetCDF files, so to
provide a consistent interface, file access is wrapped into a class
DataFile. This provides a simple interface for reading and writing files
from any of the following modules: ``netCDF4``;
``Scientific.IO.NetCDF``; and ``scipy.io.netcdf``. The DataFile class
also provides allows access to HDF5 files through the same interface,
using the ``h5py`` module. To open a file using DataFile:

.. code-block:: python

    from boututils.datafile import DataFile

    f = DataFile("file.nc")  # Open the file
    var = f.read("variable") # Read a variable from the file
    f.close()                # Close the file

or similarly for an HDF5 file

.. code-block:: python

    from boututils.datafile import DataFile

    f = DataFile("file.hdf5")  # Open the file
    var = f.read("variable")   # Read a variable from the file
    f.close()                  # Close the file

A more robust way to read from DataFiles is to use the context manager
syntax:

.. code-block:: python

    from boututils.datafile import DataFile

    with DataFile("file.hdf5") as f: # Open the file
        var = f.read("variable")     # Read a variable from the file

This way the DataFile is automatically closed at the end of the ``with``
block, even if there is an error in ``f.read``. To list the variables in
a file e.g.

.. code-block:: pycon

    >>> f = DataFile("test_io.grd.nc")
    >>> print(f.list())
    ['f3d', 'f2d', 'nx', 'ny', 'rvar', 'ivar']

and to list the names of the dimensions

.. code-block:: pycon

    >>> print(f.dimensions("f3d"))
    ('x', 'y', 'z')

or to get the sizes of the dimensions

.. code-block:: pycon

    >>> print(f.size("f3d"))
    [12, 12, 5]

To read in all variables in a file into a dictionary there is the
``file_import`` function

.. code-block:: python

    from boututils.file_import import file_import

    grid = file_import("grid.nc")

As for IDL, there is a ``collect`` routine which reads gathers together
the data from multiple processors

.. code-block:: python

    from boutdata.collect import collect

    Ni = collect("Ni")  # Collect the variable "Ni"

Matlab routines
---------------

These are Matlab routines for collecting data, showing animation and
performing some basic analysis. To use these routines, either you may
copy these routines (from **tools/matlablib**) directly to your present
working directory or a path to **tools/matlablib** should be added
before analysis.

.. code-block:: matlab

    >> addpath <full_path_BOUT_directory>/tools/matlablib/

Now, the first routine to collect data and import it to Matlab for
further analysis is

.. code-block:: matlab

    >> var = import_dmp(path,var_name);

Here, *path* is the path where the output data in netcdf format has been
dumped. *var\_name* is the name of variable which user want to load for
further analysis. For example, to load “P” variable from present working
directory:

.. code-block:: matlab

    >> P = import_dmp('.','P');

Variable “P” can be any of [X,Y,Z,T]/[X,Y,Z]/[X,Y]/Constant formats. If
we are going to Import a large data set with [X,Y,Z,T] format. Normally
such data files are of very big size and Matlab goes out of memory/ or
may take too much time to load data for all time steps. To resolve this
limitation of above routine *import\_dmp*, another routine
*import\_data\_netcdf* is being provided. It serves all purposes the
routine *import\_dmp* does but also gives user freedom to import data at
only few/specific time steps.

.. code-block:: matlab

    >> var = import_data_netcdf(path,var_name,nt,ntsp);

Here, *path* and *var\_name* are same variables as described before.
*nt* is the number of time steps user wish to load data. *ntsp* is the
steps at which one wish to write data of of total simulation times the
data written.

.. code-block:: matlab

    >> P = import_data_netcdf('.','P',5,100);

Variable “P” has been imported from present working directory for 5 time
steps. As the original netcdf data contains time information of 500
steps (assume NT=500 in BOUT++ simulations), user will pick only 5 time
steps at steps of *ntsp* i.e. 100 here. Details of other Matlab routines
provided with BOUT++ package can be looked in to README.txt of
**tools/matlablib** directory. The Matlab users can develop their own
routines using ***ncread, ncinfo, ncwrite, ncdisp, netcdf etc.***
functions provided in Matlab package.

Mathematica routines
--------------------

A package to read BOUT++ output data into Mathematica is in
``tools/mathematicalib``. To read data into Mathematica, first add this
directory to Mathematica’s path by putting

.. code-block:: mathematica

       AppendTo[$Path,"/full/path/to/BOUT>/tools/mathematicalib"]

in your Mathematica startup file (usually
``\$HOME/.Mathematica/Kernel/init.m`` ). To use the package, call

.. code-block:: mathematica

       Import["BoutCollect.m"]

from inside Mathematica. Then you can use e.g.

.. code-block:: mathematica

       f=BoutCollect[variable,path->"data"]

or

.. code-block:: mathematica

       f=BoutCollect[variable,path->"data"]

’ ``bc``\ ’ is a shorthand for ’\ ``BoutCollect`` ’. All options
supported by the Python ``collect()`` function are included, though Info
does nothing yet.

Octave routines
---------------

There is minimal support for reading data into Octave, which has been
tested on Octave 3.2. It requires the ``octcdf`` library to access
NetCDF files.

.. code-block:: octave

    f = bcollect()  # optional path argument is "." by default

    f = bsetxrange(f, 1, 10) # Set ranges
    # Same for y, z, and t (NOTE: indexing from 1!)

    u = bread(f, "U")  # Finally read the variable

.. _sec-options:

BOUT++ options
==============

The inputs to BOUT++ are a text file containing options, and for complex
grids a binary grid file in NetCDF or HDF5 format. Generating input
grids for tokamaks is described in :ref:`sec-gridgen`. The grid file
describes the size and topology of the X-Y domain, metric tensor
components and usually some initial profiles. The option file specifies
the size of the domain in the symmetric direction (Z), and controls how
the equations are evolved e.g. differencing schemes to use, and boundary
conditions. In most situations, the grid file will be used in many
different simulations, but the options may be changed frequently.

The text input file ``BOUT.inp`` is always in a subdirectory called
``data`` for all examples. The files include comments (starting with
either ’;’ or ’#’) and should be fairly self-explanatory. The format is
the same as a windows INI file, consisting of ``name = value`` pairs.
Comments are started with a hash (#) or semi-colon, which comments out
the rest of the line. values can be:

-  Integers

-  Real values

-  Booleans

-  Strings

Options are also divided into sections, which start with the section
name in square brackets.

.. code-block:: cfg

    [section1]
    something = 132         # an integer
    another = 5.131         # a real value
    yetanother = true       # a boolean
    finally = "some text"   # a string

Subsections can also be used, separated by colons ’:’, e.g.

.. code-block:: cfg

    [section:subsection]

Have a look through the examples to see how the options are used.

Command line options
--------------------

All options can be set on the command line, and will override those set
in BOUT.inp. The most commonly used are “restart” and “append”,
described in :ref:`sec-running`. If values are not given for
command-line arguments, then the value is set to ``true`` , so putting
``restart`` is equivalent to ``restart=true`` .

Values can be specified on the command line for other settings, such as
the fraction of a torus to simulate (ZPERIOD):

.. code-block:: bash

     ./command zperiod=10

Remember **no** spaces around the ’=’ sign. Like the BOUT.inp file,
setting names are not case sensitive.

Sections are separated by colons ’:’, so to set the solver type
(:ref:`sec-timeoptions`) you can either put this in BOUT.inp:

.. code-block:: cfg

    [solver]
    type = rk4

or put ``solver:type=rk4`` on the command line. This capability is used
in many test suite cases to change the parameters for each run.

General options
---------------

At the top of the BOUT.inp file (before any section headers), options
which affect the core code are listed. These are common to all physics
models, and the most useful of them are:

.. code-block:: bash

    NOUT = 100       # number of time-points output
    TIMESTEP = 1.0   # time between outputs

which set the number of outputs, and the time step between them. Note
that this has nothing to do with the internal timestep used to advance
the equations, which is adjusted automatically. What time-step to use
depends on many factors, but for high-\ :math:`\beta` reduced MHD ELM
simulations reasonable choices are ``1.0`` for the first part of a run
(to handle initial transients), then around ``10.0`` for the linear
phase. Once non-linear effects become important, you will have to reduce
the timestep to around ``0.1``.

Most large clusters or supercomputers have a limit on how long a job can
run for called “wall time”, because it’s the time taken according to a
clock on the wall, as opposed to the CPU time actually used. If this is
the case, you can use the option

.. code-block:: bash

    wall_limit = 10 # wall clock limit (in hours)

BOUT++ will then try to quit cleanly before this time runs out. Setting
a negative value (default is -1) means no limit.

Often it’s useful to be able to restart a simulation from a chosen
point, either to reproduce a previous run, or to modify the settings and
re-run. A restart file is output every timestep, but this is overwritten
each time, and so the simulation can only be continued from the end of
the last simulation. Whilst it is possible to create a restart file from
the output data afterwards, it’s much easier if you have the restart
files. Using the option

.. code-block:: bash

    archive = 20

saves a copy of the restart files every 20 timesteps, which can then be
used as a starting point.

The X and Y size of the computational grid is set by the grid file, but
the number of points in the Z (axisymmetric) direction is specified in
the options file:

.. code-block:: bash

    MZ = 33

This must be :math:`\texttt{MZ} = 2^n + 1`, and can be
:math:`2,3,5,9,\ldots`. The power of 2 is so that FFTs can be used in
this direction; the :math:`+1` is for historical reasons (inherited from
BOUT) and is going to be removed at some point.

Since the Z dimension is periodic, the domain size is specified as
multiples or fractions of :math:`2\pi`. To specify a fraction of
:math:`2\pi`, use

.. code-block:: bash

    ZPERIOD = 10

This specifies a Z range from :math:`0` to
:math:`2\pi / {\texttt{ZPERIOD}}`, and is useful for simulation of
tokamaks to make sure that the domain is an integer fraction of a torus.
If instead you want to specify the Z range directly (for example if Z is
not an angle), there are the options

.. code-block:: bash

    ZMIN = 0.0
    ZMAX = 0.1

which specify the range in multiples of :math:`2\pi`.

In BOUT++, grids can be split between processors in both X and Y
directions. By default only Y decomposition is used, and to use X
decomposition you must specify the number of processors in the X
direction:

.. code-block:: bash

    NXPE = 1  # Set number of X processors

The grid file to use is specified relative to the root directory where
the simulation is run (i.e. running “``ls ./data/BOUT.inp``” gives the
options file)

.. code-block:: bash

    grid = "data/cbm18_8_y064_x260.nc"

Communications
--------------

The communication system has a section ``[comms]``, with a true/false
option ``async``. This determines whether asynchronous MPI sends are
used; which method is faster varies (though not by much) with machine
and problem.

.. _sec-diffmethodoptions:

Differencing methods
--------------------

Differencing methods are specified in three section (``[ddx]``,
``[ddy]`` and ``[ddz]``), one for each dimension.

-  ``first``, the method used for first derivatives

-  ``second``, method for second derivatives

-  ``upwind``, method for upwinding terms

-  ``flux``, for conservation law terms

The methods which can be specified are U1, U4, C2, C4, W2, W3, FFT Apart
from FFT, the first letter gives the type of method (U = upwind, C =
central, W = WENO), and the number gives the order.

Model-specific options
----------------------

The options which affect a specific physics model vary, since they are
defined in the physics module itself (see :ref:`sec-inputopts`). They
should have a separate section, for example the high-\ :math:`\beta`
reduced MHD code uses options in a section called ``[highbeta]``.

There are three places to look for these options: the BOUT.inp file; the
physics model C++ code, and the output logs. The physics module author
should ideally have an example input file, with commented options
explaining what they do; alternately they may have put comments in the
C++ code for the module. Another way is to look at the output logs: when
BOUT++ is run, (nearly) all options used are printed out with their
default values. This won’t provide much explanation of what they do, but
may be useful anyway. See :ref:`sec-output` for more details.

.. _sec-iooptions:

Input and Output
----------------

The format of the output (dump) files can be controlled, if support for
more than one output format has been configured, by setting the
top-level option **dump\_format** to one of the recognised file
extensions: ‘nc’ for NetCDF; ‘hdf5’, ‘hdf’ or ‘h5’ for HDF5. For example
to select HDF5 instead of the default NetCDF format put

.. code-block:: cfg

    dump_format = hdf5

before any section headers. The output (dump) files with time-history
are controlled by settings in a section called “output”. Restart files
contain a single time-slice, and are controlled by a section called
“restart”. The options available are listed in table [tab:outputopts].

+-------------+----------------------------------------------------+--------------+
| Option      | Description                                        | Default      |
+-------------+----------------------------------------------------+--------------+
|             |                                                    | value        |
+-------------+----------------------------------------------------+--------------+
| enabled     | Writing is enabled                                 | true         |
+-------------+----------------------------------------------------+--------------+
| floats      | Write floats rather than doubles                   | true (dmp)   |
+-------------+----------------------------------------------------+--------------+
| flush       | Flush the file to disk after each write            | true         |
+-------------+----------------------------------------------------+--------------+
| guards      | Output guard cells                                 | true         |
+-------------+----------------------------------------------------+--------------+
| openclose   | Re-open the file for each write, and close after   | true         |
+-------------+----------------------------------------------------+--------------+
| parallel    | Use parallel I/O                                   | false        |
+-------------+----------------------------------------------------+--------------+

Table: Output file options

**enabled** is useful mainly for doing performance or scaling tests,
where you want to exclude I/O from the timings. **floats** is used to
reduce the size of the output files: restart files are stored as double
by default (since these will be used to restart a simulation), but
output dump files are set to floats by default.

To enable parallel I/O for either output or restart files, set

::

    parallel = true

in the output or restart section. If you have compiled BOUT++ with a
parallel I/O library such as pnetcdf (see
:ref:`sec-advancedinstall`), then rather than outputting one file per
processor, all processors will output to the same file. For restart
files this is particularly useful, as it means that you can restart a
job with a different number of processors. Note that this feature is
still experimental, and incomplete: output dump files are not yet
supported by the collect routines.

Implementation
--------------

To control the behaviour of BOUT++ a set of options is used, with
options organised into sections which can be nested. To represent this
tree structure there is the ``Options`` class defined in
``bout++/include/options.hxx``

::

    class Options {
     public:
      // Setting options
      void set(const string &key,const int &val,const string &source="");
      ...
      // Testing if set
      bool isSet(const string &key);
      // Getting options
      void get(const string &key,int &val,const int &def,bool log=true);
      ...
      // Get a subsection. Creates if doesn't exist
      Options* getSection(const string &name);
    };

To access the options, there is a static function (singleton)

::

      Options *options = Options::getRoot();

which gives the top-level (root) options class. Setting options is done
using the ``set()`` methods which are currently defined for ``int``,
``BoutReal``, ``bool`` and ``string`` . For example:

::

      options->set("nout", 10);      // Set an integer
      options->set("restart", true); // A bool

Often it’s useful to see where an option setting has come from e.g. the
name of the options file or “command line”. To specify a source, pass it
as a third argument:

::

      options->set("nout", 10, "manual");

To create a section, just use ``getSection`` : if it doesn’t exist it
will be created.

::

      Options *section = options->getSection("mysection");
      section->set("myswitch", true);

To get options, use the ``get()`` method which take the name of the
option, the variable to set, and the default value.

::

      int nout;
      options->get("nout", nout, 1);

Internally, ``Options`` converts all types to strings and does type
conversion when needed, so the following code would work:

::

      Options *options = Options::getRoot();
      options->set("test", "123");
      int val;
      options->get("test", val, 1);

This is because often the type of the option is not known at the time
when it’s set, but only when it’s requested.

By default, the ``get`` methods output a message to the log files giving
the value used and the source of that value. To suppress this, set the
``log`` argument to ``false`` :

::

      options->get("test", val, 1, false);

Reading options
---------------

To allow different input file formats, each file parser implements the
``OptionParser`` interface defined in
``bout++/src/sys/options/optionparser.hxx``

::

    class OptionParser {
     public:
      virtual void read(Options *options, const string &filename) = 0;
     private:
    };

and so just needs to implement a single function which reads a given
file name and inserts the options into the given ``Options`` object.

To use these parsers and read in a file, there is the ``OptionsReader``
class defined in ``bout++/include/optionsreader.hxx``

::

    class OptionsReader {
     public:
     void read(Options *options, const char *file, ...);
     void parseCommandLine(Options *options, int argc, char **argv);
    };

This is a singleton object which is accessed using

::

      OptionsReader *reader = OptionsReader::getInstance();

so to read a file ``BOUT.inp`` in a directory given in a variable
``data_dir`` the following code is used in ``bout++.cxx``:

::

      Options *options = Options::getRoot();
      OptionsReader *reader = OptionsReader::getInstance();
      reader->read(options, "%s/BOUT.inp", data_dir);

To parse command line arguments as options, the ``OptionsReader`` class
has a method:

::

      reader->parseCommandLine(options, argc, argv);

This is currently quite rudimentary and needs improving.

Variable initialisation
=======================

Variables in BOUT++ are not initialised automatically, but must be
explicitly given a value. For example the following code declares a
``Field3D`` variable then attempts to access a particular element:

::

    Field3D f;    // Declare a variable
    f(0,0,0) = 1.0;  // Error!

This results in an error because the data array to store values in ``f``
has not been allocated. Allocating data can be done in several ways:

#. Initialise with a value

   ::

           Field3D f = 0.0; // Allocates memory, fills with zeros
           f(0,0,0) = 1.0; // ok
         

   That this cannot be done at a global scope, since it requires the
   mesh to already exist and have a defined size.

#. Set to a scalar value

   ::

           Field3D f;
           f = 0.0; // Allocates memory, fills with zeros
           f(0,0,0) = 1.0; // ok
         

   Note that setting a field equal to another field has the effect of
   making both fields share the same underlying data. This behaviour is
   similar to how NumPy arrays behave in Python.

   ::

           Field3D g = 0.0;  // Allocates memory, fills with zeros
           Field3D f = g; // f now shares memory with g
           
           f(0,0,0) = 1.0; // g also modified 
         

   To ensure that a field has a unique underlying memory array call the
   ``allocate`` method before writing to individual indices.

#. Use ``allocate()`` to allocate memory

   ::

           Field3D f;
           f.allocate(); // Allocates memory, values undefined
           f(0,0,0) = 1.0; // ok
         

In a BOUT++ simulation some variables are typically evolved in time. The
initialisation of these variables is handled by the time integration
solver.

Initialisation of time evolved variables
----------------------------------------

Each variable being evolved has its own section, with the same name as
the output data. For example, the high-\ :math:`\beta` model has
variables “P”, “jpar”, and “U”, and so has sections ``[P]``, ``[jpar]``,
``[U]`` (not case sensitive).

.. _sec-expressions:

Expressions
~~~~~~~~~~~

The recommended way to initialise a variable is to use the ``function``
option for each variable:

.. code-block:: cfg

    [p]
    function = 1 + gauss(x-0.5)*gauss(y)*sin(z)

This evaluates an analytic expression to initialise the :math:`P`
variable. Expressions can include the usual operators
(``+``,\ ``-``,\ ``*``,\ ``/``), including ``^`` for exponents. The
following values are also already defined:

+--------+------------------------------------------------------------------------------------+
| Name   | Description                                                                        |
+========+====================================================================================+
| x      | :math:`x` position between :math:`0` and :math:`1`                                 |
+--------+------------------------------------------------------------------------------------+
| y      | :math:`y` position between :math:`0` and :math:`2\pi` (excluding the last point)   |
+--------+------------------------------------------------------------------------------------+
| z      | :math:`z` position between :math:`0` and :math:`2\pi` (excluding the last point)   |
+--------+------------------------------------------------------------------------------------+
| pi     | :math:`3.1415\ldots`                                                               |
+--------+------------------------------------------------------------------------------------+

Table: Initialisation expression values

By default, :math:`x` is defined as ``i / (nx - 2*MXG)``, where ``MXG``
is the width of the boundary region, by default 2. Hence :math:`x`
actually goes from 0 on the leftmost point to ``(nx-1)/(nx-4)`` on the
rightmost point. This is not a particularly good definition, but for
most cases its sufficient to create some initial profiles. For some
problems like island reconnection simulations, it’s useful to define
:math:`x` in a particular way which is more symmetric than the default.
To do this, set in BOUT.inp

.. code-block:: cfg

      [mesh]
      symmetricGlobalX = true

This will change the definition of :math:`x` to ``i / (nx - 1)``, so
:math:`x` is then between :math:`0` and :math:`1` everywhere.

The functions in table [tab:initexprfunc] are also available in
expressions.

+----------------------------------------+------------------------------------------------------------------------------+
| Name                                   | Description                                                                  |
+========================================+==============================================================================+
| abs(x)                                 | Absolute value :math:`|x|`                                                   |
+----------------------------------------+------------------------------------------------------------------------------+
| asin(x), acos(x), atan(x), atan(y,x)   | Inverse trigonometric functions                                              |
+----------------------------------------+------------------------------------------------------------------------------+
| ballooning(x)                          | Ballooning transform (:eq:`ballooning_transform`, fig [fig:ballooning])      |
+----------------------------------------+------------------------------------------------------------------------------+
| ballooning(x,n)                        | Ballooning transform, using :math:`n` terms (default 3)                      |
+----------------------------------------+------------------------------------------------------------------------------+
| cos(x)                                 | Cosine                                                                       |
+----------------------------------------+------------------------------------------------------------------------------+
| cosh(x)                                | Hyperbolic cosine                                                            |
+----------------------------------------+------------------------------------------------------------------------------+
| exp(x)                                 | Exponential                                                                  |
+----------------------------------------+------------------------------------------------------------------------------+
| tanh(x)                                | Hyperbolic tangent                                                           |
+----------------------------------------+------------------------------------------------------------------------------+
| gauss(x)                               | Gaussian :math:`\exp(-x^2/2) / \sqrt{2\pi}`                                  |
+----------------------------------------+------------------------------------------------------------------------------+
| gauss(x, w)                            | Gaussian :math:`\exp[-x^2/(2w^2)] /                                          |
|                                        | (w\sqrt{2\pi})`                                                              |
+----------------------------------------+------------------------------------------------------------------------------+
| H(x)                                   | Heaviside function: :math:`1` if :math:`x > 0` otherwise :math:`0`           |
+----------------------------------------+------------------------------------------------------------------------------+
| log(x)                                 | Natural logarithm                                                            |
+----------------------------------------+------------------------------------------------------------------------------+
| max(x,y,...)                           | Maximum (variable arguments)                                                 |
+----------------------------------------+------------------------------------------------------------------------------+
| min(x,y,...)                           | Minimum (variable arguments)                                                 |
+----------------------------------------+------------------------------------------------------------------------------+
| mixmode(x)                             | A mixture of Fourier modes                                                   |
+----------------------------------------+------------------------------------------------------------------------------+
| mixmode(x, seed)                       | seed determines random phase (default 0.5)                                   |
+----------------------------------------+------------------------------------------------------------------------------+
| power(x,y)                             | Exponent :math:`x^y`                                                         |
+----------------------------------------+------------------------------------------------------------------------------+
| sin(x)                                 | Sine                                                                         |
+----------------------------------------+------------------------------------------------------------------------------+
| sinh(x)                                | Hyperbolic sine                                                              |
+----------------------------------------+------------------------------------------------------------------------------+
| sqrt(x)                                | :math:`\sqrt{x}`                                                             |
+----------------------------------------+------------------------------------------------------------------------------+
| tan(x)                                 | Tangent                                                                      |
+----------------------------------------+------------------------------------------------------------------------------+
| erf(x)                                 | The error function                                                           |
+----------------------------------------+------------------------------------------------------------------------------+
| TanhHat(x, width, centre, steepness)   | The hat function                                                             |
|                                        | :math:`\frac{1}{2}(\tanh[s (x-[c-\frac{w}{2}])]`                             |
|                                        | :math:`- \tanh[s (x-[c+\frac{w}{2}])] )`                                     |
+----------------------------------------+------------------------------------------------------------------------------+

Table: Initialisation expression functions

For field-aligned tokamak simulations, the Y direction is along the
field and in the core this will have a discontinuity at the twist-shift
location where field-lines are matched onto each other. To handle this,
the ``ballooning`` function applies a truncated Ballooning
transformation to construct a smooth initial perturbation:

.. math::
   :label: ballooning_transform

   U_0^{balloon} = \sum_{i=-N}^N F(x)G(y + 2\pi i)H(z + q2\pi i)

.. figure:: figs/init_balloon.*
   :alt: Initial profiles
   :width: 48.0%

   Initial profiles in twist-shifted grid. **Left**: Without ballooning
   transform, showing discontinuity at the matching location **Right**:
   with ballooning transform

There is an example code ``test-ballooning`` which compares methods of
setting initial conditions with the ballooning transform.

The ``mixmode(x)`` function is a mixture of Fourier modes of the form:

.. math::

   \mathrm{mixmode}(x) = \sum_{i=1}^{14} \frac{1}{(1 +
   |i-4|)^2}\cos[ix + \phi(i, \mathrm{seed})]

where :math:`\phi` is a random phase between :math:`-\pi` and
:math:`+\pi`, which depends on the seed. The factor in front of each
term is chosen so that the 4th harmonic (:math:`i=4`) has the highest
amplitude. This is useful mainly for initialising turbulence
simulations, where a mixture of mode numbers is desired.

Initalising variables with the ``FieldFactory`` class
-----------------------------------------------------

This class provides a way to generate a field with a specified form. For
example to create a variable ``var`` from options we could write

::

    FieldFactory f(mesh);
    Field2D var = f.create2D("var");

This will look for an option called “var”, and use that expression to
initialise the variable ``var``. This could then be set in the BOUT.inp
file or on the command line.

::

    var = gauss(x-0.5,0.2)*gauss(y)*sin(3*z)

To do this, ``FieldFactory`` implements a recursive descent parser to
turn a string containing something like
``"gauss(x-0.5,0.2)*gauss(y)*sin(3*z)"`` into values in a ``Field3D`` or
``Field2D`` object. Examples are given in the ``test-fieldfactory``
example:

::

    FieldFactory f(mesh);
    Field2D b = f.create2D("1 - x");
    Field3D d = f.create3D("gauss(x-0.5,0.2)*gauss(y)*sin(z)");

This is done by creating a tree of ``FieldGenerator`` objects which then
generate the field values:

::

    class FieldGenerator {
     public:
      virtual ~FieldGenerator() { }
      virtual FieldGenerator* clone(const list<FieldGenerator*> args) {return NULL;}
      virtual BoutReal generate(int x, int y, int z) = 0;
    };

All classes inheriting from ``FieldGenerator`` must implement a
``generate`` function, which returns the value at the given ``(x,y,z)``
position. Classes should also implement a ``clone`` function, which
takes a list of arguments and creates a new instance of its class. This
takes as input a list of other ``FieldGenerator`` objects, allowing a
variable number of arguments.

The simplest generator is a fixed numerical value, which is represented
by a ``FieldValue`` object:

::

    class FieldValue : public FieldGenerator {
     public:
      FieldValue(BoutReal val) : value(val) {}
      BoutReal generate(int x, int y, int z) { return value; }
     private:
      BoutReal value;
    };

Adding a new function
---------------------

To add a new function to the FieldFactory, a new ``FieldGenerator``
class must be defined. Here we will use the example of the ``sinh``
function, implemented using a class ``FieldSinh`` . This takes a single
argument as input, but ``FieldPI`` takes no arguments, and
``FieldGaussian`` takes either one or two. Study these after reading
this to see how these are handled.

First, edit ``src/field/fieldgenerators.hxx`` and add a class
definition:

::

    class FieldSinh : public FieldGenerator {
     public:
      FieldSinh(FieldGenerator* g) : gen(g) {}
      ~FieldSinh() {if(gen) delete gen;}

      FieldGenerator* clone(const list<FieldGenerator*> args);
      BoutReal generate(int x, int y, int z);
     private:
      FieldGenerator *gen;
    };

The ``gen`` member is used to store the input argument, and to make sure
it’s deleted properly we add some code to the destructor. The
constructor takes a single input, the ``FieldGenerator`` argument to the
``sinh`` function, which is stored in the member ``gen`` .

Next edit ``src/field/fieldgenerators.cxx`` and add the implementation
of the ``clone`` and ``generate`` functions:

::

    FieldGenerator* FieldSinh::clone(const list<FieldGenerator*> args) {
      if(args.size() != 1) {
        throw ParseException("Incorrect number of arguments to sinh function. Expecting 1, got %d", args.size());
      }

      return new FieldSinh(args.front());
    }

    BoutReal FieldSinh::generate(double x, double y, double z, double t) {
      return sinh(gen->generate(x,y,z,t));
    }

The ``clone`` function first checks the number of arguments using
``args.size()`` . This is used in ``FieldGaussian`` to handle different
numbers of input, but in this case we throw a ``ParseException`` if the
number of inputs isn’t one. ``clone`` then creates a new ``FieldSinh``
object, passing the first argument ( ``args.front()`` ) to the
constructor (which then gets stored in the ``gen`` member variable).

The ``generate`` function for ``sinh`` just gets the value of the input
by calling ``gen->generate(x,y,z)``, calculates ``sinh`` of it and
returns the result.

The ``clone`` function means that the parsing code can make copies of
any ``FieldGenerator`` class if it’s given a single instance to start
with. The final step is therefore to give the ``FieldFactory`` class an
instance of this new generator. Edit the ``FieldFactory`` constructor
``FieldFactory::FieldFactory()`` in ``src/field/field_factory.cxx`` and
add the line:

::

    addGenerator("sinh", new FieldSinh(NULL));

That’s it! This line associates the string ``"sinh"`` with a
``FieldGenerator`` . Even though ``FieldFactory`` doesn’t know what type
of ``FieldGenerator`` it is, it can make more copies by calling the
``clone`` member function. This is a useful technique for polymorphic
objects in C++ called the “Virtual Constructor” idiom.

Parser internals
----------------

When a ``FieldGenerator`` is added using the ``addGenerator`` function,
it is entered into a ``std::map`` which maps strings to
``FieldGenerator`` objects (``include/field_factory.hxx``):

::

    map<string, FieldGenerator*> gen;

Parsing a string into a tree of ``FieldGenerator`` objects is done by
first splitting the string up into separate tokens like operators like
’\*’, brackets ’(’, names like ’sinh’ and so on, then recognising
patterns in the stream of tokens. Recognising tokens is done in
``src/field/field_factory.cxx``:

::

    char FieldFactory::nextToken() {
     ...

This returns the next token, and setting the variable ``char curtok`` to
the same value. This can be one of:

-  -1 if the next token is a number. The variable ``BoutReal curval`` is
   set to the value of the token

-  -2 for a string (e.g. “sinh”, “x” or “pi”). This includes anything
   which starts with a letter, and contains only letters, numbers, and
   underscores. The string is stored in the variable ``string curident``
   .

-  0 to mean end of input

-  The character if none of the above. Since letters and numbers are
   taken care of (see above), this includes brackets and operators like
   ’+’ and ’-’.

The parsing stage turns these tokens into a tree of ``FieldGenerator``
objects, starting with the ``parse()`` function

::

    FieldGenerator* FieldFactory::parse(const string &input) {
       ...

which puts the input string into a stream so that ``nextToken()`` can
use it, then calls the ``parseExpression()`` function to do the actual
parsing:

::

    FieldGenerator* FieldFactory::parseExpression() {
       ...

This breaks down expressions in stages, starting with writing every
expression as

::

    expression := primary [ op primary ]

i.e. a primary expression, and optionally an operator and another
primary expression. Primary expressions are handled by the
``parsePrimary()`` function, so first ``parsePrimary()`` is called, and
then ``parseBinOpRHS`` which checks if there is an operator, and if so
calls ``parsePrimary()`` to parse it. This code also takes care of
operator precedence by keeping track of the precedence of the current
operator. Primary expressions are then further broken down and can
consist of either a number, a name (identifier), a minus sign and a
primary expression, or brackets around an expression:

::

    primary := number
            := identifier
            := '-' primary
            := '(' expression ')'
            := '[' expression ']'

The minus sign case is needed to handle the unary minus e.g. ``"-x"`` .
Identifiers are handled in ``parseIdentifierExpr()`` which handles
either variable names, or functions

::

    identifier := name
               := name '(' expression [ ',' expression [ ',' ... ] ] ')'

i.e. a name, optionally followed by brackets containing one or more
expressions separated by commas. names without brackets are treated the
same as those with empty brackets, so ``"x"`` is the same as ``"x()"``.
A list of inputs (``list<FieldGenerator*> args;`` ) is created, the
``gen`` map is searched to find the ``FieldGenerator`` object
corresponding to the name, and the list of inputs is passed to the
object’s ``clone`` function.

Time integration
================

.. _sec-timeoptions:

Options
-------

BOUT++ can be compiled with several different time-integration solvers ,
and at minimum should have Runge-Kutta (RK4) and PVODE (BDF/Adams)
solvers available.

The solver library used is set using the ``solver:type`` option, so
either in BOUT.inp:

.. code-block:: cfg

    [solver]
    type = rk4  # Set the solver to use

or on the command line by adding ``solver:type=pvode`` for example:

.. code-block:: bash

    mpirun -np 4 ./2fluid solver:type=rk4

**NB**: Make sure there are no spaces around the “=” sign:
``solver:type =pvode`` won’t work (probably). Table [tab:solvers] gives
a list of time integration solvers, along with any compile-time options
needed to make the solver available.

+---------------+-----------------------------------------+--------------------+
| Name          | Description                             | Compile options    |
+===============+=========================================+====================+
| euler         | Euler explicit method                   | Always available   |
+---------------+-----------------------------------------+--------------------+
| rk4           | Runge-Kutta 4th-order explicit method   | Always available   |
+---------------+-----------------------------------------+--------------------+
| karniadakis   | Karniadakis explicit method             | Always available   |
+---------------+-----------------------------------------+--------------------+
| pvode         | 1998 PVODE with BDF method              | Always available   |
+---------------+-----------------------------------------+--------------------+
| cvode         | SUNDIALS CVODE. BDF and Adams methods   | –with-cvode        |
+---------------+-----------------------------------------+--------------------+
| ida           | SUNDIALS IDA. DAE solver                | –with-ida          |
+---------------+-----------------------------------------+--------------------+
| petsc         | PETSc TS methods                        | –with-petsc        |
+---------------+-----------------------------------------+--------------------+
| imexbdf2      | IMEX-BDF2 scheme                        | –with-petsc        |
+---------------+-----------------------------------------+--------------------+

Table: Available time integration solvers

Each solver can have its own settings which work in slightly different
ways, but some common settings and which solvers they are used in are
given in table [tab:solveropts].

+------------------+--------------------------------------------+-------------------------------------+
| Option           | Description                                | Solvers used                        |
+==================+============================================+=====================================+
| atol             | Absolute tolerance                         | rk4, pvode, cvode, ida              |
+------------------+--------------------------------------------+-------------------------------------+
| rtol             | Relative tolerance                         | rk4, pvode, cvode, ida              |
+------------------+--------------------------------------------+-------------------------------------+
| mxstep           | Maximum internal steps                     | rk4                                 |
|                  | per output step                            |                                     |
+------------------+--------------------------------------------+-------------------------------------+
| max\_timestep    | Maximum timestep                           | rk4, cvode                          |
+------------------+--------------------------------------------+-------------------------------------+
| timestep         | Starting timestep                          | rk4, karniadakis, euler, imexbdf2   |
+------------------+--------------------------------------------+-------------------------------------+
| adaptive         | Adapt timestep? (Y/N)                      | rk4                                 |
+------------------+--------------------------------------------+-------------------------------------+
| use\_precon      | Use a preconditioner? (Y/N)                | pvode, cvode, ida                   |
+------------------+--------------------------------------------+-------------------------------------+
| mudq, mldq       | BBD preconditioner settings                | pvode, cvode, ida                   |
+------------------+--------------------------------------------+-------------------------------------+
| mukeep, mlkeep   |                                            |                                     |
+------------------+--------------------------------------------+-------------------------------------+
| maxl             |                                            |                                     |
+------------------+--------------------------------------------+-------------------------------------+
| use\_jacobian    | Use user-supplied Jacobian? (Y/N)          | cvode                               |
+------------------+--------------------------------------------+-------------------------------------+
| adams\_moulton   | Use Adams-Moulton method                   | cvode                               |
|                  | rather than BDF                            |                                     |
+------------------+--------------------------------------------+-------------------------------------+
| diagnose         | Collect and print additional diagnostics   | cvode                               |
+------------------+--------------------------------------------+-------------------------------------+

Table: Time integration solver options

The most commonly changed options are the absolute and relative solver
tolerances, ``ATOL`` and ``RTOL`` which should be varied to check
convergence.

ODE integration
---------------

The Solver class can be used to solve systems of ODEs inside a physics
model: Multiple Solver objects can exist besides the main one used for
time integration. Example code is in ``examples/test-integrate``.

To use this feature, systems of ODEs must be represented by a class
derived from ``PhysicsModel`` (see :ref:`sec-newapi`).

::

    class MyFunction : public PhysicsModel {
     public:
      int init(bool restarting) {
        // Initialise ODE
        // Add variables to solver as usual
        solver->add(result, "result");
        ...
      }

      int rhs(BoutReal time) {
        // Specify derivatives of fields as usual
        ddt(result) = ...
      }
     private:
      Field3D result;
    };

To solve this ODE, create a new Solver object:

::

    Solver* ode = Solver::create(Options::getRoot()->getSection("ode"));

This will look in the section ``[ode]`` in the options file.
**Important:** To prevent this solver overwriting the main restart files
with its own restart files, either disable restart files:

.. code-block:: cfg

    [ode]
    enablerestart = false

or specify a different directory to put the restart files:

.. code-block:: cfg

    [ode]
    restartdir = ode  # Restart files ode/BOUT.restart.0.nc, ...

Create a model object, and pass it to the solver:

::

    MyFunction* model = new MyFunction();
    ode->setModel(model);

Finally tell the solver to perform the integration:

::

    ode->solve(5, 0.1);

The first argument is the number of steps to take, and the second is the
size of each step. These can also be specified in the options, so
calling

::

    ode->solve();

will cause ode to look in the input for ``nout`` and ``timestep``
options:

.. code-block:: cfg

    [ode]
    nout = 5
    timestep = 0.1

Finally, delete the model and solver when finished:

::

    delete model;
    delete solver;

**Note:** If an ODE needs to be solved multiple times, at the moment it
is recommended to delete the solver, and create a new one each time.

Preconditioning
---------------

At every time step, an implicit scheme such as BDF has to solve a
non-linear problem to find the next solution. This is usually done using
Newton’s method, each step of which involves solving a linear (matrix)
problem. For :math:`N` evolving variables is an :math:`N\times N` matrix
and so can be very large. By default matrix-free methods are used, in
which the Jacobian :math:`\mathcal{J}` is approximated by finite
differences (see next subsection), and so this matrix never needs to be
explicitly calculated. Finding a solution to this matrix can still be
difficult, particularly as :math:`\delta t` gets large compared with
some time-scales in the system (i.e. a stiff problem).

A preconditioner is a function which quickly finds an approximate
solution to this matrix, speeding up convergence to a solution. A
preconditioner does not need to include all the terms in the problem
being solved, as the preconditioner only affects the convergence rate
and not the final solution. A good preconditioner can therefore
concentrate on solving the parts of the problem with the fastest
time-scales.

A simple example  [3]_ is a coupled wave equation, solved in the
``test-precon`` example code:

.. math::

   \frac{\partial u}{\partial t} = \partial_{||}v \qquad \frac{\partial
   v}{\partial t} = \partial_{||} u

First, calculate the Jacobian of this set of equations by taking
partial derivatives of the time-derivatives with respect to each of the
evolving variables

.. math::

   \mathcal{J} = (\begin{array}{cc}
   \frac{\partial}{\partial u}\frac{\partial u}{\partial t} &
   \frac{\partial}{\partial v}\frac{\partial u}{\partial t}\\
   \frac{\partial}{\partial u}\frac{\partial v}{\partial t} &
   \frac{\partial}{\partial v}\frac{\partial v}{\partial t}
   \end{array}
   ) = (\begin{array}{cc}
   0 & \partial_{||} \\
   \partial_{||} & 0
   \end{array}
   )

In this case :math:`\frac{\partial u}{\partial t}` doesn’t depend on
:math:`u` nor :math:`\frac{\partial v}{\partial t}` on :math:`v`, so the
diagonal is empty. Since the equations are linear, the Jacobian doesn’t
depend on :math:`u` or :math:`v` and so

.. math::

   \frac{\partial}{\partial t}(\begin{array}{c} u \\
   v \end{array}) = \mathcal{J} (\begin{array}{c} u \\
   v \end{array} )

In general for non-linear functions :math:`\mathcal{J}` gives the
change in time-derivatives in response to changes in the state variables
:math:`u` and :math:`v`.

In implicit time stepping, the preconditioner needs to solve an equation

.. math::

   \mathcal{I} - \gamma \mathcal{J}

where :math:`\mathcal{I}` is the identity matrix, and :math:`\gamma`
depends on the time step and method (e.g. :math:`\gamma = \delta t` for
backwards Euler method). For the simple wave equation problem, this is

.. math::

   \mathcal{I} - \gamma \mathcal{J} = (\begin{array}{cc}
   1 & -\gamma\partial_{||} \\
   -\gamma\partial_{||} & 1
   \end{array}
   )

This matrix can be block inverted using Schur factorisation  [4]_

.. math::

   (\begin{array}{cc}
     {\mathbf{E}} & {\mathbf{U}} \\
     {\mathbf{L}} & {\mathbf{D}}
   \end{array})^{-1}
    = (\begin{array}{cc}
     {\mathbf{I}} & -{\mathbf{E}}^{-1}{\mathbf{U}} \\
     0 & {\mathbf{I}}
   \end{array}
   )(\begin{array}{cc}
     {\mathbf{E}}^{-1} & 0 \\
     0 & {\mathbf{P}}_{Schur}^{-1}
   \end{array}
   )(\begin{array}{cc}
     {\mathbf{I}} & 0 \\
     -{\mathbf{L}}{\mathbf{E}}^{-1} & {\mathbf{I}}
   \end{array}
   )

where
:math:`{\mathbf{P}}_{Schur} = {\mathbf{D}} - {\mathbf{L}}{\mathbf{E}}^{-1}{\mathbf{U}}`
Using this, the wave problem becomes:

.. math::
   :label: precon

   (\begin{array}{cc} 1 & -\gamma\partial_{||} \\
   -\gamma\partial_{||} & 1 \end{array})^{-1} = (\begin{array}{cc} 1 & \gamma\partial_{||}\\
   0 & 1 \end{array} )(\begin{array}{cc} 1 & 0 \\
   0 & (1 -\gamma^2\partial^2_{||})^{-1} \end{array} )(\begin{array}{cc} 1 & 0\\
   \gamma\partial_{||} & 1 \end{array} )

The preconditioner is implemented by defining a function of the form

::

    int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
      ...
    }

which takes as input the current time, the :math:`\gamma` factor
appearing above, and :math:`\delta` which is only important for
constrained problems (not discussed here... yet). The current state of
the system is stored in the state variables (here ``u`` and ``v`` ),
whilst the vector to be preconditioned is stored in the time derivatives
(here ``ddt(u)`` and ``ddt(v)`` ). At the end of the preconditioner the
result should be in the time derivatives. A preconditioner which is just
the identity matrix and so does nothing is therefore:

::

    int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
    }

To implement the preconditioner in equation :eq:`precon`, first apply the
rightmost matrix to the given vector:

.. math::

   (\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   ) = (\begin{array}{cc}
   1 & 0 \\
   \gamma\partial_{||} & 1
   \end{array}
   )(\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   )

::

    int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
      mesh->communicate(ddt(u));
      //ddt(u) = ddt(u);
      ddt(v) = gamma*Grad_par(ddt(u)) + ddt(v);

note that since the preconditioner is linear, it doesn’t depend on
:math:`u` or :math:`v`. As in the RHS function, since we are taking a
differential of ``ddt(u)``, it first needs to be communicated to
exchange guard cell values.

The second matrix

.. math::

   (\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   ) arrow (\begin{array}{cc}
   1 & 0 \\
   0 & (1 - \gamma^2\partial^2_{||})^{-1}
   \end{array}
   )(\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   )

doesn’t alter :math:`u`, but solves a parabolic equation in the
parallel direction. There is a solver class to do this called
``InvertPar`` which solves the equation
:math:`(A + B\partial_{||}^2)x = b` where :math:`A` and :math:`B`
are ``Field2D`` or constants  [5]_. In ``physics_init`` we create one of
these solvers:

::

    InvertPar *inv; // Parallel inversion class
    int physics_init(bool restarting) {
       ...
       inv = InvertPar::Create();
       inv->setCoefA(1.0);
       ...
    }

In the preconditioner we then use this solver to update :math:`v`:

::

      inv->setCoefB(-SQ(gamma));
      ddt(v) = inv->solve(ddt(v));

which solves
:math:`ddt(v) arrow (1 - \gamma^2\partial_{||}^2)^{-1} ddt(v)`.
The final matrix just updates :math:`u` using this new solution for
:math:`v`

.. math::

   (\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   ) arrow (\begin{array}{cc}
   1 & \gamma\partial_{||} \\
   0 & 1
   \end{array}
   )(\begin{array}{c}
   \texttt{ddt(u)} \\
   \texttt{ddt(v)}
   \end{array}
   )

::

      mesh->communicate(ddt(v));
      ddt(u) = ddt(u) + gamma*Grad_par(ddt(v));

Finally, boundary conditions need to be imposed, which should be
consistent with the conditions used in the RHS

::

      ddt(u).applyBoundary("dirichlet");
      ddt(v).applyBoundary("dirichlet");

To use the preconditioner, pass the function to the solver in
``physics_init``

::

    int physics_init(bool restarting) {
      solver->setPrecon(precon);
      ...
    }

then in the ``BOUT.inp`` settings file switch on the preconditioner

.. code-block:: bash

    [solver]
    type = cvode          # Need CVODE or PETSc
    use_precon = true     # Use preconditioner
    rightprec = false     # Use Right preconditioner (default left)

Jacobian function
-----------------

DAE constraint equations
------------------------

Using the IDA or IMEX-BDF2 solvers, BOUT++ can solve Differential
Algebraic Equations (DAEs), in which algebraic constraints are used for
some variables. Examples of how this is used are in the
``examples/constraints`` subdirectory.

First the variable to be constrained is added to the solver, in a
similar way to time integrated variables. For example

::

    Field3D phi;
    ...
    solver->constraint(phi, ddt(phi), "phi");

The first argument is the variable to be solved for (constrained). The
second argument is the field to contain the residual (error). In this
example the time derivative field ``ddt(phi)`` is used, but it could be
another ``Field3D`` variable. The solver will attempt to find a solution
to the first argument (``phi`` here) such that the second argument
(``ddt(phi)``) is zero to within tolerances.

In the RHS function the residual should be calculated. In this example
(``examples/constraints/drift-wave-constraint``) we have:

::

    ddt(phi) = Delp2(phi) - Vort;

so the time integration solver includes the algebraic constraint
``Delp2(phi) = Vort`` i.e. (:math:`\nabla_\perp^2\phi = \omega`).

IMEX-BDF2
---------

This is an implicit-explicit multistep method, which uses the PETSc
library for the SNES nonlinear solver. To use this solver, BOUT++ must
have been configured with PETSc support, and the solver type set to
``imexbdf2``

::

    [solver]
    type = imexbdf2

For examples of using IMEX-BDF2, see the ``examples/IMEX/``
subdirectory, in particular the ``diffusion-nl``, ``drift-wave`` and
``drift-wave-constrain`` examples.

The time step is currently fixed (not adaptive), and defaults to the
output timestep. To set a smaller internal timestep, the
``solver:timestep`` option can be set. If the timestep is too large,
then the explicit part of the problem may become unstable, or the
implicit part may fail to converge.

The implicit part of the problem can be solved matrix-free, in which
case the Jacobian-vector product is approximated using finite
differences. This is currently the default, and can be set on the
command-line using the options:

::

     solver:matrix_free=true  -snes_mf

Note the ``-snes_mf`` flag which is passed to PETSc. When using a matrix
free solver, the Jacobian is not calculated and so the amount of memory
used is minimal. However, since the Jacobian is not known, many standard
preconditioning methods cannot be used, and so in many cases a custom
preconditioner is needed to obtain good convergence.

An experimental feature uses PETSc’s ability to calculate the Jacobian
using finite differences. This can then speed up the linear solve, and
allows more options for preconditioning. To enable this option:

::

     solver:matrix_free=false

There are two ways to calculate the Jacobian: A brute force method which
is set up by this call to PETSc which is generally very slow, and a
“coloring” scheme which can be quite fast and is the default. Coloring
uses knowledge of where the non-zero values are in the Jacobian, to work
out which rows can be calculated simultaneously. The coloring code in
IMEX-BDF2 currently assumes that every field is coupled to every other
field in a star pattern: one cell on each side, a 7 point stencil for 3D
fields. If this is not the case for your problem, then the solver may
not converge.

The brute force method can be useful for comparing the Jacobian
structure, so to turn off coloring:

::

     solver:use_coloring=false

Using MatView calls, or the ``-mat_view`` PETSc options, the non-zero
structure of the Jacobian can be plotted or printed.

Monitoring the simulation output
--------------------------------

Monitoring of the solution can be done at two levels: output monitoring,
and timestep monitoring. Output monitoring occurs only when data is
written to file, whereas timestep monitoring is every timestep and so
(usually) much more frequent. Examples of both are in
``examples/monitor`` and ``examples/monitor-newapi``.

**Output monitoring**: At every output timestep the solver calls a
monitor function, which writes the output dump file, calculates and
prints timing information and estimated time remaining. If you want to
run additional code or write data to a different file, you can add
monitor function(s).

You can call your output monitor function whatever you like, but it must
have 4 inputs and return an int:

::

    int my_output_monitor(Solver *solver, BoutReal simtime, int iter, int NOUT) {
      ...
    }

The first input is the solver object, the second is the current
simulation time, the third is the output number, and the last is the
total number of outputs requested. To get the solver to call this
function every output time, put in your ``physics_init`` code:

::

      solver->addMonitor(my_output_monitor);

If you want to later remove a monitor, you can do so with

::

      solver->removeMonitor(my_output_monitor);

A simple example using this monitor is:

::

    int my_output_monitor(Solver *solver, BoutReal simtime, int iter, int NOUT) {
      output.write("My monitor, time = %e, dt = %e\n",
          simtime, solver->getCurrentTimestep());
    }

    int physics_init(bool restarting) {
      solver->addMonitor(my_monitor);
    }

See the monitor example (``examples/monitor``) for full code.

**Timestep monitoring**: This works in the same way as output
monitoring. First define a monitor function:

::

    int my_timestep_monitor(Solver *solver, BoutReal simtime, BoutReal lastdt) {
      ...
    }

where ``simtime`` will again contain the current simulation time, and
``lastdt`` the last timestep taken. Add this function to the solver:

::

      solver->addTimestepMonitor(my_timestep_monitor);

Timestep monitoring is disabled by default, unlike output monitoring. To
enable timestep monitoring, set in the options file (BOUT.inp):

::

    [solver]
    monitor_timestep = true

or put on the command line ``solver:monitor_timestep=true`` . When this
is enabled, it will change how solvers like CVODE and PVODE (the default
solvers) are used. Rather than being run in NORMAL mode, they will
instead be run in SINGLE\_STEP mode (see the SUNDIALS notes
here:\ http://computation.llnl.gov/casc/sundials/support/notes.html).
This may in some cases be less efficient.

.. _sec-bndryopts:

Boundary conditions
===================

Like the variable initialisation, boundary conditions can be set for
each variable in individual sections, with default values in a section
``[All]``. Boundary conditions are specified for each variable, being
applied to variable itself during initialisation, and the
time-derivatives at each timestep. They are a combination of a basic
boundary condition, and optional modifiers.

When finding the boundary condition for a variable ``var`` on a boundary
region, the options are checked in order from most to least specific:

-  Section ``var``, ``bndry_`` + region name. Depending on the mesh
   file, regions of the grid are given labels. Currently these are
   ``core``, ``sol``, ``pf`` and ``target`` which are intended for
   tokamak edge simulations. Hence the variables checked are
   ``bndry_core``, ``bndry_pf`` etc.

-  Section ``var``, ``bndry_`` + boundary side. These names are ``xin``,
   ``xout``, ``yup`` and ``ydown``.

-  Section ``var``, variable ``bndry_all``

-  The same settings again except in section ``All``.

The default setting for everything is therefore ``bndry_all`` in the
``All`` section.

Boundary conditions are given names, with optional arguments in
brackets. Currently implemented boundary conditions are:

-  ``dirichlet`` - Set to zero

-  ``dirichlet(<number>)`` - Set to some number e.g. ``dirichlet(1)``
   sets the boundary to :math:`1.0`

-  ``neumann`` - Zero gradient

-  ``robin`` - A combination of zero-gradient and zero-value
   :math:`a f + b{{\frac{\partial f}{\partial x}}} = g` where the
   syntax is ``robin(a, b, g)``.

-  ``constgradient`` - Constant gradient across boundary

-  ``zerolaplace`` - Laplacian = 0, decaying solution (X boundaries
   only)

-  ``zerolaplace2`` - Laplacian = 0, using coefficients from the
   Laplacian inversion and Delp2 operator.

-  ``constlaplace`` - Laplacian = const, decaying solution (X boundaries
   only)

The zero- or constant-Laplacian boundary conditions works as follows:

.. math::

   \nabla_\perp^2 f =& 0 \\ &\simeq& g^{xx}\frac{\partial^2 f}{\partial x^2} +
       g^{zz}\frac{\partial^2 f}{\partial z^2}

 which when Fourier transformed in :math:`z` becomes:

.. math::

   g^{xx}\frac{\partial^2 \hat{f}}{\partial x^2} - g^{zz}k_z^2 \hat{f} = 0

 which has the solution

.. math::

   \hat{f} = Ae^{xk_z\sqrt{g^{zz}/g^{xx}}} + Be^{-xk_z\sqrt{g^{zz}/g^{xx}}}

Assuming that the solution should decay away from the domain, on the
inner :math:`x` boundary :math:`B = 0`, and on the outer boundary
:math:`A = 0`. Boundary modifiers change the behaviour of boundary
conditions, and more than one modifier can be used. Currently the
following are available:

-  ``relax`` - Relaxing boundaries. Evolve the variable towards the
   given boundary condition at a given rate

-  ``shifted`` - Apply boundary conditions in orthogonal X-Z
   coordinates, rather than field-aligned

-  ``width`` - Modifies the width of the region over which the boundary
   condition is applied

These are described in the following subsections.

Relaxing boundaries
-------------------

All boundaries can be modified to be “relaxing” which are a combination
of zero-gradient time-derivative, and whatever boundary condition they
are applied to. The idea is that this prevents sharp discontinuities at
boundaries during transients, whilst maintaining the desired boundary
condition on longer time-scales. In some cases this can improve the
numerical stability and timestep.

For example, ``relax(dirichlet)`` will make a field :math:`f` at point
:math:`i` in the boundary follow a point :math:`i-1` in the domain:

.. math::

   .{{\frac{\partial f}{\partial t}}}|_i = .{{\frac{\partial f}{\partial t}}}|_{i-1}  - f_i / \tau

where :math:`\tau` is a time-scale for the boundary (currently set to
0.1, but will be a global option). When the time-derivatives are slow
close to the boundary, the boundary relaxes to the desired condition
(Dirichlet in this case), but when the time-derivatives are large then
the boundary approaches Neumann to reduce discontinuities.

By default, the relaxation rate is set to :math:`10` (i.e. a time-scale
of :math:`\tau=0.1`). To change this, give the rate as the second
argument e.g. ``relax(dirichlet, 2)`` would relax to a Dirichlet
boundary condition at a rate of :math:`2`.

Shifted boundaries
------------------

By default boundary conditions are applied in field-aligned coordinates,
where :math:`y` is along field-lines but :math:`x` has a discontinuity
at the twist-shift location. If radial derivatives are being done in
shifted coordinates where :math:`x` and :math:`z` are orthogonal, then
boundary conditions should also be applied in shifted coordinates. To do
this, the ``shifted`` boundary modifier applies a :math:`z` shift,
applies the boundary condition, then shifts back. For example:

::

    bndry_core = shifted( neumann )

would ensure that radial derivatives were zero in shifted coordinates on
the core boundary.

Changing the width of boundaries
--------------------------------

To change the width of a boundary region, the ``width`` modifier changes
the width of a boundary region before applying the boundary condition,
then changes the width back afterwards. To use, specify the boundary
condition and the width, for example

::

    bndry_core = width( neumann , 4 )

would apply a Neumann boundary condition on the innermost 4 cells in the
core, rather than the usual 2. When combining with other boundary
modifiers, this should be applied first e.g.

::

    bndry_sol = width( relax( dirichlet ), 3)

would relax the last 3 cells towards zero, whereas

::

    bndry_sol = relax( width( dirichlet, 3) )

would only apply to the usual 2, since relax didn’t use the updated
width.

Limitations:

#. Because it modifies then restores a globally-used BoundaryRegion,
   this code is not thread safe.

#. Boundary conditions can’t be applied across processors, and no checks
   are done that the width asked for fits within a single processor.

Examples
--------

This example is taken from the UEDGE benchmark test (in
``examples/uedge-benchmark``):

.. code-block:: cfg

    [All]
    bndry_all = neumann # Default for all variables, boundaries

    [Ni]
    bndry_target = neumann
    bndry_core = relax(dirichlet(1.))   # 1e13 cm^-3 on core boundary
    bndry_all  = relax(dirichlet(0.1))  # 1e12 cm^-3 on other boundaries

    [Vi]
    bndry_ydown = relax(dirichlet(-1.41648))   # -3.095e4/Vi_x
    bndry_yup   = relax(dirichlet( 1.41648))

The variable ``Ni`` (density) is set to a Neumann boundary condition on
the targets (yup and ydown), relaxes towards :math:`1` on the core
boundary, and relaxes to :math:`0.1` on all other boundaries. Note that
the ``bndry_target = neumann`` needs to be in the ``Ni`` section: If we
just had

.. code-block:: cfg

    [All]
    bndry_all = neumann # Default for all variables, boundaries

    [Ni]
    bndry_core = relax(dirichlet(1.))   # 1e13 cm^-3 on core boundary
    bndry_all  = relax(dirichlet(0.1))  # 1e12 cm^-3 on other boundaries

then the “target” boundary condition for ``Ni`` would first search in
the ``[Ni]`` section for ``bndry_target``, then for ``bndry_all`` in the
``[Ni]`` section. This is set to ``relax(dirichlet(0.1))``, not the
Neumann condition desired.

Boundary regions
----------------

The boundary condition code needs ways to loop over the boundary
regions, without needing to know the details of the mesh.

At the moment two mechanisms are provided: A RangeIterator over upper
and lower Y boundaries, and a vector of BoundaryRegion objects.

::

    // Boundary region iteration
    virtual const RangeIterator iterateBndryLowerY() const = 0;
    virtual const RangeIterator iterateBndryUpperY() const = 0;

    bool hasBndryLowerY();
    bool hasBndryUpperY();

    bool BoundaryOnCell; // NB: DOESN'T REALLY BELONG HERE

The ``RangeIterator`` class is an iterator which allows looping over a
set of indices. For example, in ``src/solver/solver.cxx`` to loop over
the upper Y boundary of a 2D variable ``var``:

::

    for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
      ...
    }

The ``BoundaryRegion`` class is defined in
``include/boundary_region.hxx``

Boundary regions
----------------

Different regions of the boundary such as “core”, “sol” etc. are
labelled by the ``Mesh`` class (i.e. ``BoutMesh``), which implements a
member function defined in ``mesh.hxx``:

::

      // Boundary regions
      virtual vector<BoundaryRegion*> getBoundaries() = 0;

This returns a vector of pointers to ``BoundaryRegion`` objects, each of
which describes a boundary region with a label, a ``BndryLoc`` location
(i.e. inner x, outer x, lower y, upper y or all), and iterator functions
for looping over the points. This class is defined in
``boundary_region.hxx``:

::

    /// Describes a region of the boundary, and a means of iterating over it
    class BoundaryRegion {
      public:
      BoundaryRegion();
      BoundaryRegion(const string &name, int xd, int yd);
      virtual ~BoundaryRegion();

      string label; // Label for this boundary region

      BndryLoc location; // Which side of the domain is it on?

      int x,y; // Indices of the point in the boundary
      int bx, by; // Direction of the boundary [x+dx][y+dy] is going outwards

      virtual void first() = 0;
      virtual void next() = 0; // Loop over every element from inside out (in X or
    Y first)
      virtual void nextX() = 0; // Just loop over X
      virtual void nextY() = 0; // Just loop over Y
      virtual bool isDone() = 0; // Returns true if outside domain. Can use this
    with nested nextX, nextY
    };

**Example:** To loop over all points in ``BoundaryRegion *bndry`` , use

::

      for(bndry->first(); !bndry->isDone(); bndry->next()) {
        ...
      }

Inside the loop, ``bndry->x`` and ``bndry->y`` are the indices of the
point, whilst ``bndry->bx`` and ``bndry->by`` are unit vectors out of
the domain. The loop is over all the points from the domain outwards
i.e. the point ``[bndry->x - bndry->bx][bndry->y - bndry->by]`` will
always be defined.

Sometimes it’s useful to be able to loop over just one direction along
the boundary. To do this, it is possible to use ``nextX()`` or
``nextY()`` rather than ``next()``. It is also possible to loop over
both dimensions using:

::

      for(bndry->first(); !bndry->isDone(); bndry->nextX())
        for(; !bndry->isDone(); bndry->nextY()) {
          ...
        }

Boundary operations
-------------------

On each boundary, conditions must be specified for each variable. The
different conditions are imposed by ``BoundaryOp`` objects. These set
the values in the boundary region such that they obey e.g. Dirichlet or
Neumann conditions. The ``BoundaryOp`` class is defined in
``boundary_op.hxx``:

::

    /// An operation on a boundary
    class BoundaryOp {
     public:
      BoundaryOp() {bndry = NULL;}
      BoundaryOp(BoundaryRegion *region)

      // Note: All methods must implement clone, except for modifiers (see below)
      virtual BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);

      /// Apply a boundary condition on field f
      virtual void apply(Field2D &f) = 0;
      virtual void apply(Field3D &f) = 0;

      virtual void apply(Vector2D &f);

      virtual void apply(Vector3D &f);

      /// Apply a boundary condition on ddt(f)
      virtual void apply_ddt(Field2D &f);
      virtual void apply_ddt(Field3D &f);
      virtual void apply_ddt(Vector2D &f);
      virtual void apply_ddt(Vector3D &f);

      BoundaryRegion *bndry;
    };

(where the implementations have been removed for clarity). Which has a
pointer to a ``BoundaryRegion`` object specifying which region this
boundary is operating on.

Boundary conditions need to be imposed on the initial conditions (after
``physics_init()``), and on the time-derivatives (after
``physics_run()``). The ``apply()`` functions are therefore called
during initialisation and given the evolving variables, whilst the
``apply_ddt`` functions are passed the time-derivatives.

To implement a boundary operation, as a minimum the ``apply(Field2D)``,
``apply(Field2D)`` and ``clone()`` need to be implemented: By default
the ``apply(Vector)`` will call the ``apply(Field)`` functions on each
component individually, and the ``apply_ddt()`` functions just call the
``apply()`` functions.

**Example**: Neumann boundary conditions are defined in
``boundary_standard.hxx``:

::

    /// Neumann (zero-gradient) boundary condition
    class BoundaryNeumann : public BoundaryOp {
     public:
      BoundaryNeumann() {}
     BoundaryNeumann(BoundaryRegion *region):BoundaryOp(region) { }
      BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
      void apply(Field2D &f);
      void apply(Field3D &f);
    };

and implemented in ``boundary_standard.cxx``

::

    void BoundaryNeumann::apply(Field2D &f) {
      // Loop over all elements and set equal to the next point in
      for(bndry->first(); !bndry->isDone(); bndry->next())
        f[bndry->x][bndry->y] = f[bndry->x - bndry->bx][bndry->y - bndry->by];
    }

    void BoundaryNeumann::apply(Field3D &f) {
      for(bndry->first(); !bndry->isDone(); bndry->next())
        for(int z=0;z<mesh->ngz;z++)
          f[bndry->x][bndry->y][z] = f[bndry->x - bndry->bx][bndry->y -
    bndry->by][z];
    }

This is all that’s needed in this case since there’s no difference
between applying Neumann conditions to a variable and to its
time-derivative, and Neumann conditions for vectors are just Neumann
conditions on each vector component.

To create a boundary condition, we need to give it a boundary region to
operate over:

::

    BoundaryRegion *bndry = ...
    BoundaryOp op = new BoundaryOp(bndry);

The ``clone`` function is used to create boundary operations given a
single object as a template in ``BoundaryFactory``. This can take
additional arguments as a vector of strings - see explanation in
:ref:`sec-BoundaryFactory`.

Boundary modifiers
------------------

To create more complicated boundary conditions from simple ones (such as
Neumann conditions above), boundary operations can be modified by
wrapping them up in a ``BoundaryModifier`` object, defined in
``boundary_op.hxx``:

::

    class BoundaryModifier : public BoundaryOp {
     public:
      virtual BoundaryOp* clone(BoundaryOp *op, const list<string> &args) = 0;
     protected:
      BoundaryOp *op;
    };

Since ``BoundaryModifier`` inherits from ``BoundaryOp``, modified
boundary operations are just a different boundary operation and can be
treated the same (Decorator pattern). Boundary modifiers could also be
nested inside each other to create even more complicated boundary
operations. Note that the ``clone`` function is different to the
``BoundaryOp`` one: instead of a ``BoundaryRegion`` to operate on,
modifiers are passed a ``BoundaryOp`` to modify.

Currently the only modifier is ``BoundaryRelax``, defined in
``boundary_standard.hxx``:

::

    /// Convert a boundary condition to a relaxing one
    class BoundaryRelax : public BoundaryModifier {
     public:
      BoundaryRelax(BoutReal rate) {r = fabs(rate);}
      BoundaryOp* clone(BoundaryOp *op, const list<string> &args);

      void apply(Field2D &f);
      void apply(Field3D &f);

      void apply_ddt(Field2D &f);
      void apply_ddt(Field3D &f);
     private:
      BoundaryRelax() {} // Must be initialised with a rate
      BoutReal r;
    };

.. _sec-BoundaryFactory:

Boundary factory
----------------

The boundary factory creates new boundary operations from input strings,
for example turning “relax(dirichlet)” into a relaxing Dirichlet
boundary operation on a given region. It is defined in
``boundary_factory.hxx`` as a Singleton, so to get a pointer to the
boundary factory use

::

      BoundaryFactory *bfact = BoundaryFactory::getInstance();

and to delete this singleton, free memory and clean-up at the end use:

::

      BoundaryFactory::cleanup();

Because users should be able to add new boundary conditions during
``physics_init()``, boundary conditions are not hard-wired into
``BoundaryFactory``. Instead, boundary conditions must be registered
with the factory, passing an instance which can later be cloned. This is
done in ``bout++.cxx`` for the standard boundary conditions:

::

      BoundaryFactory* bndry = BoundaryFactory::getInstance();
      bndry->add(new BoundaryDirichlet(), "dirichlet");
      ...
      bndry->addMod(new BoundaryRelax(10.), "relax");

where the ``add`` function adds BoundaryOp objects, whereas ``addMod``
adds ``BoundaryModifier`` objects. **Note**: The objects passed to
``BoundaryFactory`` will be deleted when ``cleanup()`` is called.

When a boundary operation is added, it is given a name such as
“dirichlet”, and similarly for the modifiers (“relax” above). These
labels and object pointers are stored internally in ``BoundaryFactory``
in maps defined in ``boundary_factory.hxx``:

::

      // Database of available boundary conditions and modifiers
      map<string, BoundaryOp*> opmap;
      map<string, BoundaryModifier*> modmap;

These are then used by ``BoundaryFactory::create()``:

::

      /// Create a boundary operation object
      BoundaryOp* create(const string &name, BoundaryRegion *region);
      BoundaryOp* create(const char* name, BoundaryRegion *region);

to turn a string such as “relax(dirichlet)” and a ``BoundaryRegion``
pointer into a ``BoundaryOp`` object. These functions are implemented in
``boundary_factory.cxx``, starting around line 42. The parsing is done
recursively by matching the input string to one of:

-  ``modifier(<expression>, arg1, ...)``

-  ``modifier(<expression>)``

-  ``operation(arg1, ...)``

-  ``operation``

the ``<expression>`` variable is then resolved into a BoundaryOp object
by calling ``create(<expression, region)``.

When an operator or modifier is found, it is created from the pointer
stored in the ``opmap`` or ``modmap`` maps using the ``clone`` method,
passing a ``list<string>`` reference containing any arguments. It’s up
to the operation implementation to ensure that the correct number of
arguments are passed, and to parse them into floats or other types.

**Example**: The Dirichlet boundary condition can take an optional
argument to change the value the boundary’s set to. In
``boundary_standard.cxx``:

::

    BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion *region, const list<string>
    &args) {
      if(!args.empty()) {
        // First argument should be a value
        stringstream ss;
        ss << args.front();

        BoutReal val;
        ss >> val;
        return new BoundaryDirichlet(region, val);
      }
      return new BoundaryDirichlet(region);
    }

If no arguments are passed i.e. the string was “dirichlet” or
“dirichlet()” then the ``args`` list is empty, and the default value
(0.0) is used. If one or more arguments is used then the first argument
is parsed into a ``BoutReal`` type and used to create a new
``BoundaryDirichlet`` object. If more arguments are passed then these
are just ignored; probably a warning should be printed.

To set boundary conditions on a field, ``FieldData`` methods are defined
in ``field_data.hxx``:

::

    // Boundary conditions
      void setBoundary(const string &name); ///< Set the boundary conditions
      void setBoundary(const string &region, BoundaryOp *op); ///< Manually set
      virtual void applyBoundary() {}
      virtual void applyTDerivBoundary() {};
     protected:
      vector<BoundaryOp*> bndry_op; // Boundary conditions

The ``setBoundary(const string &name)`` method is implemented in
``field_data.cxx``. It first gets a vector of pointers to
``BoundaryRegion``\ s from the mesh, then loops over these calling
``BoundaryFactory::createFromOptions`` for each one and adding the
resulting boundary operator to the ``bndry_op`` vector.

Iterating over fields
=====================

In BOUT++ 4.0.0, we now have the ability to use C++ range-based
for-loops. This means that it is possible to iterate over a whole field
using a single loop:

::

    Field3D f(0.0);
    for (auto i : f) {
       f[i] = a[i] + b[i];
    }

The iterator provides access to the x, y, z indices:

::

    Field3D f(0.0);
    for (auto i : f) {
       f[i] = i.x + i.y + i.z;
    }

It is also possible to specify regions to iterate over using this
syntax:

::

    Field3D f(0.0);
    for (auto i : f.region(RGN_NOBNDRY) {
       f[i] = 1.0;
    }

Available regions are:

-  ``RGN_ALL``, which is the whole mesh;

-  ``RGN_NOBNDRY``, which skips all boundaries;

-  ``RGN_NOX``, which skips the x boundaries

-  ``RGN_NOY``, which skips the y boundaries

.. _sec-gridgen:

Generating input grids
======================

The simulation mesh describes the number and topology of grid points,
the spacing between them, and the coordinate system. For many problems,
a simple mesh can be created using options.

.. code-block:: bash

    [mesh]
    nx = 260  # X grid size
    ny = 256  # Y grid size

    dx = 0.1  # X mesh spacing
    dy = 0.1  # Y mesh spacing

The above options will create a :math:`260\times 256` mesh in X and Y
(MZ option sets Z resolution), with mesh spacing of :math:`0.1` in both
directions. By default the coordinate system is Cartesian (metric tensor
is the identity matrix), but this can be changed by specifying the
metric tensor components.

Integer quantities such as ``nx`` must be numbers (like “260”), not
expressions (like “256 + 2\*MXG”). Real (floating-point) values can be
expressions, allowing quite complicated analytic inputs. For example in
the example ``test-griddata``:

.. code-block:: bash

    # Screw pinch

    rwidth = 0.4

    Rxy = 0.1 + rwidth*x  # Radius from axis     [m]
    L   = 10              # Length of the device [m]

    dy = L/ny
    hthe = 1.0

    Zxy = L * y / (2*pi)

    Bpxy = 1.0      # Axial field [T]
    Btxy = 0.1*Rxy  # Azimuthal field [T]
    Bxy = sqrt(Btxy^2 + Bpxy^2)

    dr = rwidth / nx
    dx = dr * Bpxy * Rxy

These expressions use the same mechanism as used for variable
initialisation (:ref:`sec-expressions`): ``x`` is a variable from
:math:`0` to :math:`1` in the domain which is uniform in index space;
``y`` and ``z`` go from :math:`0` to :math:`2\pi`. As with variable
initialisation, common trigonometric and mathematical functions can be
used. In the above example, some variables depend on each other, for
example ``dy`` depends on ``L`` and ``ny``. The order in which these
variables are defined doesn’t matter, so ``L`` could be defined below
``dy``, but circular dependencies are not allowed. If the variables are
defined in the same section (as ``dy`` and ``L``) then no section prefix
is required. To refer to a variable in a different section, prefix the
variable with the section name e.g. “``section:variable``”.

More complex meshes can be created by supplying an input grid file to
describe the grid points, geometry, and starting profiles. Currently
BOUT++ supports either NetCDF, HDF5 format binary files. During startup,
BOUT++ looks in the grid file for the following variables. If any are
not found, a warning will be printed and the default values used.

-  X and Y grid sizes (integers) ``nx`` and ``ny`` **REQUIRED**

-  Differencing quantities in 2D arrays ``dx[nx][ny]`` and
   ``dy[nx][ny]``. If these are not found they will be set to 1.

-  Diagonal terms of the metric tensor :math:`g^{ij}` ``g11[nx][ny]``,
   ``g22[nx][ny]``, and ``g33[nx][ny]``. If not found, these will be set
   to 1.

-  Off-diagonal metric tensor :math:`g^{ij}` elements ``g12[nx][ny]``,
   ``g13[nx][ny]``, and ``g23[nx][ny]``. If not found, these will be set
   to 0.

-  Z shift for interpolation between field-aligned coordinates and
   shifted coordinates (see ``manual/coordinates.pdf``). Perpendicular
   differential operators are calculated in shifted coordinates when
   ``ShiftXderivs`` in ``mesh/mesh.hxx`` is enabled. ``ShiftXderivs``
   can be set in the root section of ``BOUT.inp`` as
   ``ShiftXderivs = true``. The shifts must be provided in the gridfile
   in a field ``zshift[nx][ny]``. If not found, ``zshift`` is set to
   zero.

The remaining quantities determine the topology of the grid. These are
based on tokamak single/double-null configurations, but can be adapted
to many other situations.

-  Separatrix locations ``ixseps1``, and ``ixseps2`` If neither is
   given, both are set to nx (i.e. all points in closed “core” region).
   If only ``ixseps1`` is found, ``ixseps2`` is set to nx, and if only
   ixseps2 is found, ixseps1 is set to -1.

-  Branch-cut locations ``jyseps1_1``, ``jyseps1_2``, ``jyseps2_1``, and
   ``jyseps2_2``

-  Twist-shift matching condition ``ShiftAngle[nx]`` for field aligned
   coordinates. This is applied in the “core” region between indices
   ``jyseps2_2``, and ``jyseps1_1 + 1``, if either ``TwistShift = True``
   enabled in the options file or in general the ``TwistShift`` flag in
   ``mesh/impls/bout/boutmesh.hxx`` is enabled by other means. BOUT++
   automatically reads the twist shifts in the gridfile if the shifts
   are stored in a field in a field ShiftAngle[nx]. If not given, this
   is set to zero.

The only quantities which are required are the sizes of the grid. If
these are the only quantities specified, then the coordinates revert to
Cartesian.

This section describes how to generate inputs for tokamak equilibria. If
you’re not interested in tokamaks then you can skip to the next section.

The directory ``tokamak_grids`` contains code to generate input grid
files for tokamaks. These can be used by the ``2fluid`` and
``highbeta_reduced`` modules, and are (mostly) compatible with inputs to
the BOUT-06 code.

Figure [fig:gridgen] shows the routines and file formats used in taking
output from different codes and converting into input to BOUT++.

BOUT++ Topology
---------------

Basic
~~~~~

In order to handle tokamak geometry BOUT++ contains an internal topology
which is determined by the branch-cut locations (``jyseps1_1``,
``jyseps1_2``, ``jyseps2_1``, and ``jyseps2_2``) and separatrix
locations (``ixseps1`` and ``ixseps2``).

The separatrix locations, ``ixseps1`` and ``ixseps2``, give the indices
in the ``x`` domain where the first and second separatrices are located.

If ``ixseps1 == ixseps2`` then there is a single separatrix representing
the boundary between the core region and the SOL region and the grid is
a connected double null configuration. If ``ixseps1 > ixseps2`` then
there are two separatrices and the inner separatrix is ``ixseps2`` so
the tokamak is an upper double null. If ``ixseps1 < ixseps2`` then there
are two separatrices and the inner separatrix is ``ixseps1`` so the
tokamak is a lower double null.

In other words: Let us for illustrative purposes say that
``ixseps1 > ixseps2`` (see figure [fig:topology\_cross\_section]). Let
us say that we have a field ``f(x,y,z)`` with a global ``x``-index which
includes ghost points. ``f(x<=xseps1,y,z)``) will then be periodic in
the ``y``-direction, ``f(xspes1<x<=xseps2,y,z)``) will have boundary
condition in the ``y``-direction set by the lowermost ``ydown`` and
``yup``. If ``f(xspes2<x,y,z)``) the boundary condition in the
``y``-direction will be set by the uppermost ``ydown`` and ``yup``. As
for now, there is no difference between the two sets of upper and lower
``ydown`` and ``yup`` boundary conditions (unless manually specified,
see section [sec:custom\_BC]).

These values are set either in the grid file or in ``BOUT.inp``. Figure
[fig:topology\_cross\_section] shows schematically how ``ixseps`` is
used.

The branch cut locations, ``jyseps1_1``, ``jyseps1_2``, ``jyseps2_1``,
and ``jyseps2_2``, split the ``y`` domain into logical regions defining
the SOL, the PFR (private flux region) and the core of the tokamak. This
is illustrated also in figure [fig:topology\_cross\_section]. If
``jyseps1_2 == jyseps2_1`` then the grid is a single null configuration,
otherwise the grid is a double null configuration.

Advanced
~~~~~~~~

The internal domain in BOUT++ is deconstructed into a series of
logically rectangular sub-domains with boundaries determined by the
``ixseps`` and ``jyseps`` parameters. The boundaries coincide with
processor boundaries so the number of grid points within each sub-domain
must be an integer multiple of ``ny/nypes`` where ``ny`` is the number
of grid points in ``y`` and ``nypes`` is the number of processors used
to split the y domain. Processor communication across the domain
boundaries is then handled internally. Figure [fig:topology\_schematic]
shows schematically how the different regions of a double null tokamak
with ``ixseps1 = ixseps2`` are connected together via communications.

Implementations
~~~~~~~~~~~~~~~

In BOUT++ each processor has a logically rectangular domain, so any
branch cuts needed for X-point geometry (see
figure [fig:topology\_schematic]) must be at processor boundaries.

In the standard “bout” mesh (``src/mesh/impls/bout/``), the
communication is controlled by the variables

::

    int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
    int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
    int IDATA_DEST, ODATA_DEST;

These control the behavior of the communications as shown in
figure [fig:boutmesh-comms].

In the Y direction, each boundary region (**U**\ p and **D**\ own in Y)
can be split into two, with ``0 <= x < UDATA_XSPLIT`` going to the
processor index ``UDATA_INDEST``, and ``UDATA_INDEST <= x < ngx`` going
to ``UDATA_OUTDEST``. Similarly for the Down boundary. Since there are
no branch-cuts in the X direction, there is just one destination for the
**I**\ nner and **O**\ uter boundaries. In all cases a negative
processor number means that there’s a domain boundary so no
communication is needed.

The communication control variables are set in the ``topology()``
function, in ``src/mesh/impls/bout/boutmesh.cxx`` starting around line
2056. First the function ``default_connections()`` sets the topology to
be a rectangle

To change the topology, the function ``set_connection`` checks that the
requested branch cut is on a processor boundary, and changes the
communications consistently so that communications are two-way and there
are no “dangling” communications.

3D variables
------------

BOUT++ was originally designed for tokamak simulations where the input
equilibrium varies only in X-Y, and Z is used as the axisymmetric
toroidal angle direction. In those cases, it is often convenient to have
input grids which are only 2D, and allow the Z dimension to be specified
independently, such as in the options file. The problem then is how to
store 3D variables in the grid file?

Two representations are now supported for 3D variables:

#. A Fourier representation. If the size of the toroidal domain is not
   specified in the grid file (``nz`` is not defined), then 3D fields
   are stored as Fourier components. In the Z dimension the coefficients
   must be stored as

   .. math::

      [n = 0, n = 1 (\textrm{real}), n = 1 (\textrm{imag}), n = 2
      (\textrm{real}), n = 2 (\textrm{imag}), \ldots ]

   where :math:`n` is the toroidal mode number. The size of the array
   must therefore be odd in the Z dimension, to contain a constant
   (:math:`n=0`) component followed by real/imaginary pairs for the
   non-axisymmetric components.

   If you are using IDL to create a grid file, there is a routine in
   ``tools/idllib/bout3dvar.pro`` for converting between BOUT++’s real
   and Fourier representation.

#. Real space, as values on grid points. If ``nz`` is set in the grid
   file, then 3D variables in the grid file must have size
   ``nx``\ :math:`\times`\ ``ny``\ :math:`\times`\ ``nz``. These are
   then read in directly into ``Field3D`` variables as required.

From EFIT files
---------------

An IDL code called “Hypnotoad” has been developed to create BOUT++ input
files from R-Z equilibria. This can read EFIT ’g’ files, find flux
surfaces, and calculate metric coefficients. The code is in
``tools/tokamak_grids/gridgen``, and has its own manual under the
``doc`` subdirectory.

From ELITE and GATO files
-------------------------

Currently conversions exist for ELITE ``.eqin`` and GATO ``dskgato``
equilibrium files. Conversion of these into BOUT++ input grids is in two
stages: In the first, both these input files are converted into a common
NetCDF format which describes the Grad-Shafranov equilibrium. These
intermediate files are then converted to BOUT++ grids using an
interactive IDL script.

Generating equilibria
---------------------

The directory ``tokamak_grids/shifted_circle`` contains IDL code to
generate shifted circle (large aspect ratio) Grad-Shafranov equilibria.

.. _sec-laplacian:

Laplacian inversion
===================

A common problem in plasma models is to solve an equation of the form

.. math::
   :label: full_laplace_inv

   d\nabla^2_\perp x + \frac{1}{c_1}(\nabla_\perp c_2)\cdot\nabla_\perp x +
   a x = b

For example,

.. math::

   \nabla_\perp^2 x + a x = b

appears in reduced MHD for the vorticity inversion and :math:`j_{||}`.

Alternative formulations and ways to invert equation
(:eq:`full_laplace_inv`) can be found in section :ref:`sec-LaplaceXY` and
:ref:`sec-LaplaceXZ`

Usage of the laplacian inversion
--------------------------------

| In BOUT++, equation (:eq:`full_laplace_inv`) can be solved in two
  ways. The first method Fourier transforms in the :math:`z`-direction,
  whilst the other is solving the full two dimensional problem by matrix
  inversion. The derivation of :math:`\nabla_\perp^2f` for a general
  coordinate system can be found in the ``coordinates`` manual. What is
  important, is to note that if :math:`g_{xy}` and :math:`g_{yz}` are
  non-zero, BOUT++ is neglecting the :math:`y`-parallel derivatives when
  using the solvers ``Laplacian`` and ``LaplaceXZ``.
|  
| By neglecting the :math:`y`-derivatives (or if
  :math:`g_{xy}=g_{yz}=0`), one can solve equation
  (:eq:`full_laplace_inv`) :math:`y` plane by :math:`y` plane.

The first approach utilizes that it is possible Fourier transform the
equation in :math:`z` (using some assumptions described in section
[sec:num\_laplace]), and solve a tridiagonal system for each
mode. These inversion problems are band-diagonal (tri-diagonal in the
case of 2nd-order differencing) and so inversions can be very
efficient: :math:`O(n_z \log n_z)` for the FFTs,
:math:`O(n_x)` for tridiagonal inversion using the Thomas
algorithm [1]_, where :math:`n_x` and :math:`n_z` are the number of
grid-points in the :math:`x` and :math:`z` directions respectively.

.. [1] Numerical recipes in C. The art of scientific computing, Press, W H and Teukolsky, S A and Vetterling, W T and Flannery, B P

In the second approach, the full :math:`2`\ -D system is being solved.
This requires PETSc to be built with BOUT++.

The ``Laplacian`` class is defined in ``invert_laplace.hxx`` and solves
problems formulated like equation (:eq:`full_laplace_inv`) To use
this class, first create an instance of it:

::

    Laplacian *lap = Laplacian::create();

By default, this will use the options in a section called “laplace”, but
can be given a different section as an argument. By default
:math:`d = 1`, :math:`a = 0`, and the :math:`c=1`. To set the values of
these coefficients, there are the ``setCoefA()``, ``setCoefC()``, and
``setCoefD()`` methods:

::

    Field2D a = ...;
    lap->setCoefA(a);
    lap->setCoefC(0.5);

arguments can be ``Field2D``, ``Field3D`` , or real values.

Settings for the inversion can be set in the input file under the
section ``laplace`` (default) or whichever settings section name was
specified when the ``Laplacian`` class was created. Commonly used
settings are listed in tables [tab:laplacesettings] to
[tab:laplaceflags].

In particular boundary conditions on the :math:`x` boundaries can be set
using the and ``outer_boundary_flags`` variables, as detailed in table
[tab:laplaceBCflags]. Note that DC (‘direct-current’) refers to
:math:`k = 0` Fourier component, AC (‘alternating-current’) refers to
:math:`k
\neq 0` Fourier components. Non-Fourier solvers use AC options (and
ignore DC ones). Multiple boundary conditions can be selected by adding
together the required boundary condition flag values together. For
example, ``inner_boundary_flags = 3`` will set a Neumann boundary
condition on both AC and DC components.

| It is pertinent to note here that the boundary in BOUT++ is defined by
  default to be located half way between the first guard point and first
  point inside the domain. For example, when a Dirichlet boundary
  condition is set, using ``inner_boundary_flags = 0`` , ``16``, or
  ``32``, then the first guard point, :math:`f_{-}` will be set to
  :math:`f_{-} = 2v - f_+`, where :math:`f_+` is the first grid point
  inside the domain, and :math:`v` is the value to which the boundary is
  being set to.
|  
| The ``global_flags``, ``inner_boundary_flags``,
  ``outer_boundary_flags`` and ``flags`` values can also be set from
  within the physics module using ``setGlobalFlags``,
  ``setInnerBoundaryFlags`` , ``setOuterBoundaryFlags`` and ``setFlags``
  .

::

    lap->setGlobalFlags(Global_Flags_Value);
    lap->setInnerBoundaryFlags(Inner_Flags_Value);
    lap->setOuterBoundaryFlags(Outer_Flags_Value);
    lap->setFlags(Flags_Value);

+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| Name                     | Meaning                                                              | Default value                          |
+==========================+======================================================================+========================================+
| ``type``                 | Which implementation to use                                          | ``tri`` (serial), ``spt`` (parallel)   |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``filter``               | Filter out modes above :math:`(1-`\ ``filter``\                      | 0                                      |
|                          | :math:`)\times k_{max}`, if using Fourier solver                     |                                        |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``maxmode``              | Filter modes with :math:`n >`\ ``maxmode``                           | ``MZ``/2                               |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``all_terms``            | Include first derivative terms                                       | ``true``                               |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``global_flags``         | Sets global inversion options See table                              | ``0``                                  |
|                          | :ref:`Laplace global flags<tab-laplaceglobalflags>`                  |                                        |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``inner_boundary_flags`` | Sets boundary conditions on inner boundary. See table                | ``0``                                  |
|                          | :ref:`Laplace boundary flags<tab-laplaceBCflags>`                    |                                        |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``outer_boundary_flags`` | Sets boundary conditions on outer boundary. See table                | ``0``                                  |
|                          | :ref:`Laplace boundary flags<tab-laplaceBCflags>`                    |                                        |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``flags``                | DEPRECATED. Sets global solver options and boundary                  | ``0``                                  |
|                          | conditions. See :ref:`Laplace flags<tab-laplaceflags>` or            |                                        |
|                          | :doc:`invert_laplace.cxx<_breathe_autogen/file/invert__laplace_8cxx>`|                                        |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+
| ``include_yguards``      | Perform inversion in :math:`y`\ -boundary guard cells                | ``true``                               |
+--------------------------+----------------------------------------------------------------------+----------------------------------------+

Table: Laplacian inversion options

.. _tab-laplaceglobalflags:

+--------+--------------------------------------------------------------------------------+-----------------------------+
| Flag   | Meaning                                                                        | Code variable               |
+========+================================================================================+=============================+
| 0      | No global option set                                                           | :math:`-`                   |
+--------+--------------------------------------------------------------------------------+-----------------------------+
| 1      | zero DC component (Fourier solvers)                                            | ``INVERT_ZERO_DC``          |
+--------+--------------------------------------------------------------------------------+-----------------------------+
| 2      | set initial guess to 0 (iterative solvers)                                     | ``INVERT_START_NEW``        |
+--------+--------------------------------------------------------------------------------+-----------------------------+
| 4      | equivalent to                                                                  | ``INVERT_BOTH_BNDRY_ONE``   |
|        | ``outer_boundary_flags = 128``,                                                |                             |
|        | ``inner_boundary_flags = 128``                                                 |                             |
+--------+--------------------------------------------------------------------------------+-----------------------------+
| 8      | Use 4th order differencing (Apparently not actually implemented anywhere!!!)   | ``INVERT_4TH_ORDER``        |
+--------+--------------------------------------------------------------------------------+-----------------------------+
| 16     | Set constant component (:math:`k_x = k_z = 0`) to zero                         | ``INVERT_KX_ZERO``          |
+--------+--------------------------------------------------------------------------------+-----------------------------+

Table: Laplacian inversion ``global_flags`` values: add the required
quantities together.

.. _tab-laplaceBCflags:

+--------+----------------------------------------------------------------------+----------------------------+
| Flag   | Meaning                                                              | Code variable              |
+========+======================================================================+============================+
| 0      | Dirichlet (Set boundary to 0)                                        | :math:`-`                  |
+--------+----------------------------------------------------------------------+----------------------------+
| 1      | Neumann on DC component (set gradient to 0)                          | ``INVERT_DC_GRAD``         |
+--------+----------------------------------------------------------------------+----------------------------+
| 2      | Neumann on AC component (set gradient to 0)                          | ``INVERT_AC_GRAD``         |
+--------+----------------------------------------------------------------------+----------------------------+
| 4      | Zero or decaying Laplacian on AC components (                        | ``INVERT_AC_LAP``          |
|        | :math:`\frac{\partial^2}{\partial x^2}+k_z^2` vanishes/decays)       |                            |
+--------+----------------------------------------------------------------------+----------------------------+
| 8      | Use symmetry to enforce zero value or gradient (redundant for 2nd    | ``INVERT_SYM``             |
|        | order now)                                                           |                            |
+--------+----------------------------------------------------------------------+----------------------------+
| 16     | Set boundary condition to values in boundary guard cells of second   | ``INVERT_SET``             |
|        | argument, ``x0``, of ``Laplacian::solve(const Field3D &b, const      |                            |
|        | Field3D &x0)`` . May be combined with any combination of 0, 1 and 2, |                            |
|        | i.e. a Dirichlet or Neumann boundary condition set to values which   |                            |
|        | are :math:`\neq 0` or :math:`f(y)`                                   |                            |
+--------+----------------------------------------------------------------------+----------------------------+
| 32     | Set boundary condition to values in boundary guard cells of RHS,     | ``INVERT_RHS``             |
|        | ``b`` in ``Laplacian::solve(const Field3D &b, const Field3D &x0)``   |                            |
|        | . May be combined with any combination of 0, 1 and 2, i.e. a         |                            |
|        | Dirichlet or Neumann boundary condition set to values which are      |                            |
|        | :math:`\neq 0` or :math:`f(y)`                                       |                            |
+--------+----------------------------------------------------------------------+----------------------------+
| 64     | Zero or decaying Laplacian on DC components                          | ``INVERT_DC_LAP``          |
|        | (:math:`\frac{\partial^2}{\partial x^2}` vanishes/decays)            |                            |
+--------+----------------------------------------------------------------------+----------------------------+
| 128    | Assert that there is only one guard cell in the :math:`x`-boundary   | ``INVERT_BNDRY_ONE``       |
+--------+----------------------------------------------------------------------+----------------------------+
| 256    | DC value is set to parallel gradient, :math:`\nabla_\parallel f`     | ``INVERT_DC_GRADPAR``      |
+--------+----------------------------------------------------------------------+----------------------------+
| 512    | DC value is set to inverse of parallel gradient                      | ``INVERT_DC_GRADPARINV``   |
|        | :math:`1/\nabla_\parallel f`                                         |                            |
+--------+----------------------------------------------------------------------+----------------------------+
| 1024   | Boundary condition for inner ‘boundary’ of cylinder                  | ``INVERT_IN_CYLINDER``     |
+--------+----------------------------------------------------------------------+----------------------------+

Table: Laplacian inversion ``outer_boundary_flags`` or
``inner_boundary_flags`` values: add the required quantities together.

.. _tab-laplaceflags:

+--------+------------------------------------------------------------------------------------------+
| Flag   | Meaning                                                                                  |
+========+==========================================================================================+
| 1      | Zero-gradient DC on inner (X) boundary. Default is zero-value                            |
+--------+------------------------------------------------------------------------------------------+
| 2      | Zero-gradient AC on inner boundary                                                       |
+--------+------------------------------------------------------------------------------------------+
| 4      | Zero-gradient DC on outer boundary                                                       |
+--------+------------------------------------------------------------------------------------------+
| 8      | Zero-gradient AC on outer boundary                                                       |
+--------+------------------------------------------------------------------------------------------+
| 16     | Zero DC component everywhere                                                             |
+--------+------------------------------------------------------------------------------------------+
| 32     | Not used currently                                                                       |
+--------+------------------------------------------------------------------------------------------+
| 64     | Set width of boundary to 1 (default is ``MXG``)                                          |
+--------+------------------------------------------------------------------------------------------+
| 128    | Use 4\ :math:`^{th}`-order band solver (default is 2\ :math:`^{nd}` order tridiagonal)   |
+--------+------------------------------------------------------------------------------------------+
| 256    | Attempt to set zero laplacian AC component on inner boundary by combining                |
|        | 2nd and 4th-order differencing at the boundary.                                          |
|        | Ignored if tridiagonal solver used.                                                      |
+--------+------------------------------------------------------------------------------------------+
| 512    | Zero laplacian AC on outer boundary                                                      |
+--------+------------------------------------------------------------------------------------------+
| 1024   | Symmetric boundary condition on inner boundary                                           |
+--------+------------------------------------------------------------------------------------------+
| 2048   | Symmetric outer boundary condition                                                       |
+--------+------------------------------------------------------------------------------------------+

Table: Laplacian inversion ``flags`` values (DEPRECATED!): add the
required quantities together.

To perform the inversion, there’s the ``solve`` method

::

    x = lap->solve(b);

If you prefer, there are functions compatible with older versions of the
BOUT++ code:

::

    Field2D a, c, d;
    invert_laplace(b, x, flags, &a, &c, &d);

and

::

    x = invert_laplace(b, flags, &a, &c, &d);

The input ``b`` and output ``x`` are 3D fields, and the coefficients
``a``, ``c``, and ``d`` are pointers to 2D fields. To omit any of the
three coefficients, set them to NULL.

Numerical implementation
------------------------

We will here go through the implementation of the laplacian inversion
algorithm, as it is performed in BOUT++. We would like to solve the
following equation for :math:`f`

.. math::
   :label: to_invert

   d\nabla_\perp^2f + \frac{1}{c_1}(\nabla_\perp c_2)\cdot\nabla_\perp f + af = b

BOUT++ is neglecting the :math:`y`-parallel derivatives if
:math:`g_{xy}` and :math:`g_{yz}` are no-zero when using the solvers
``Laplacian`` and ``LaplaceXZ``. For these two solvers, equation
(:eq:`to_invert`) becomes (see ``coordinates`` manual for derivation)

.. math::
   :label: invert_expanded

   \, &d (g^{xx} \partial_x^2 + G^x \partial_x + g^{zz} \partial_z^2 +
   G^z \partial_z + 2g^{xz} \partial_x \partial_z ) f \\
   +& \frac{1}{c_1}( {{\boldsymbol{e}}}^x \partial_x +
   {\boldsymbol{e}}^z \partial_z ) c_2 \cdot ({\boldsymbol{e}}^x
   \partial_x + {\boldsymbol{e}}^z \partial_z ) f \\ +& af = b


Using tridiagonal solvers
~~~~~~~~~~~~~~~~~~~~~~~~~

When using the tridiagonal solvers, :math:`c_1 = c_2` in equation
(:eq:`to_invert`), hence, it is rather solving

.. math::
   :label: to_invert_tri

   d\nabla_\perp^2f + \frac{1}{c}(\nabla_\perp c)\cdot\nabla_\perp f + af = b

Since there are no parallel :math:`y`-derivatives if
:math:`g_{xy}=g_{yz}=0` (or if they are neglected), equation
(:eq:`to_invert_tri`) will only contain derivatives of :math:`x` and
:math:`z` for the dependent variable. The hope is that the modes in the
periodic :math:`z` direction will decouple, so that we in the end only
have to invert for the :math:`x` coordinate.

If the modes decouples when Fourier transforming equation
(:eq:`invert_expanded`), we can use a tridiagonal solver to solve the
equation for each Fourier mode.

Using the discrete Fourier transform

.. math::

   F(x,y)_{k} = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp(\frac{-2\pi i k
   Z}{N})

we see that the modes will not decouple if a term consist of a product
of two terms which depends on :math:`z`, as this would give terms like

.. math::

   \frac{1}{N}\sum_{Z=0}^{N-1} a(x,y)_Z f(x,y)_Z \exp(\frac{-2\pi i k
   Z}{N})

Thus, in order to use a tridiagonal solver, :math:`a`, :math:`c` and
:math:`d` cannot be functions of :math:`z`. Because of this, the
:math:`{{\boldsymbol{e}}}^z \partial_z c` term in equation
(:eq:`invert_expanded`) is zero. In principle the modes would still
decouple if the :math:`{{\boldsymbol{e}}}^z \partial_z f`
part of equation (:eq:`invert_expanded`) was kept, but currently this
part is also neglected in solvers using a tridiagonal matrix. Thus the
tridiagonal solvers are solving equations on the form

.. math::

   \, &d(x,y) ( g^{xx}(x,y) \partial_x^2 + G^x(x,y) \partial_x +
       g^{zz}(x,y) \partial_z^2 + G^z(x,y) \partial_z + 2g^{xz}(x,y)
       \partial_x \partial_z ) f(x,y,z) \\
     +& \frac{1}{c(x,y)}({{\boldsymbol{e}}}^x \partial_x ) c(x,y) \cdot (
       {{\boldsymbol{e}}}^x \partial_x ) f(x,y,z) \\
     +& a(x,y)f(x,y,z) = b(x,y,z)

after using the discrete Fourier transform (see section
[sec:deriv\_of\_FT]), we get

.. math::

   \, &d (    g^{xx} \partial_x^2F_z + G^x \partial_xF_z + g^{zz} [i k]^2F_z
        + G^z [i k]F_z + 2g^{xz} \partial_x[i k]F_z ) \\
     +& \frac{1}{c}( {{\boldsymbol{e}}}^x \partial_x ) c \cdot ( {{\boldsymbol{e}}}^x
        \partial_xF_z ) \\
     +& aF_z = B_z

which gives

.. math::
   :label: FT_laplace_inversion

   \, &d ( g^{xx} \partial_x^2 + G^x \partial_x - k^2 g^{zz} + i
   kG^z + i k2g^{xz} \partial_x )F_z \\
   +& \frac{g^{xx}}{c} (\partial_x c ) \partial_xF_z \\
   +& aF_z = B_z

As nothing in equation (:eq:`FT_laplace_inversion`) couples points in
:math:`y` together (since we neglected the :math:`y`-derivatives if
:math:`g_{xy}` and :math:`g_{yz}` were non-zero). Also, as the modes are
decoupled, we may solve equation (:eq:`FT_laplace_inversion`) :math:`k`
mode by :math:`k` mode in addition to :math:`y`\ -plane by
:math:`y`\ -plane.

The second order centred approximation of the first and second
derivatives in :math:`x` reads

.. math::

       &&\partial_x f \simeq \frac{-f_{n-1} + f_{n+1}}{2\text{d}x}&&
       &&\partial_x^2 f \simeq \frac{f_{n-1} - f_{n} + f_{n+1}}{\text{d}x^2}&&

This gives

.. math::

       \, &d (    g^{xx} \frac{F_{z,n-1} - 2F_{z,n} + F_{z, n+1}}{\text{d}x^2} +
       G^x \frac{-F_{z,n-1} + F_{z,n+1}}{2\text{d}x} - k^2 g^{zz}F_{z,n} .\\
       &\quad.  + i kG^zF_{z,n} + i k2g^{xz} \frac{-F_{z,n-1} +
   F_{z,n+1}}{2\text{d}x} ) \\
       +& \frac{g^{xx}}{c} ( \frac{-c_{n-1} + c_{n+1}}{2\text{d}x} )
   \frac{-F_{z,n-1} + F_{z,n+1}}{2\text{d}x} \\
       +& aF_{z,n} = B_{z,n}

collecting point by point

.. math::
   :label: discretized_laplace

       &( \frac{dg^{xx}}{\text{d}x^2} - \frac{dG^x}{2\text{d}x} -
       \frac{g^{xx}}{c_{n}} \frac{-c_{n-1} + c_{n+1}}{4\text{d}x^2} - i\frac{d
       k2g^{xz}}{2\text{d}x} ) F_{z,n-1} \\
           +&( - \frac{ dg^{xx} }{\text{d}x^2} - dk^2 g^{zz} + a + idkG^z )
       F_{z,n} \\
           +&( \frac{dg^{xx}}{\text{d}x^2} + \frac{dG^x}{2\text{d}x} +
       \frac{g^{xx}}{c_{n}} \frac{-c_{n-1} + c_{n+1}}{4\text{d}x^2} +
       i\frac{dk2g^{xz}}{2\text{d}x} ) F_{z, n+1} \\
        =& B_{z,n}

We now introduce

.. math::

       &c_1 = \frac{dg^{xx}}{\text{d}x^2}& &c_2 = dg^{zz}& &c_3 =
       \frac{2dg^{xz}}{2\text{d}x}& && \\ &c_4 = \frac{dG^x + g^{xx}\frac{-c_{n-1}
       + c_{n+1}}{2c_n\text{d}x}}{2\text{d}x}& &c_5 = dG^z& &&

which inserted in equation (:eq:`discretized_laplace`) gives

.. math::

       &( c_1 - c_4 -ikc_3 ) F_{z,n-1} \\
           +&( -2c_1 - k^2c_2 +ikc_5 + a ) F_{z,n} \\
           +&( c_1 + c_4 + ikc_3 ) F_{z, n+1} \\
        =& B_{z,n}

This can be formulated as the matrix equation

.. math::

   AF_z=B_z

where the matrix :math:`A` is tridiagonal. The boundary conditions are
set by setting the first and last rows in :math:`A` and :math:`B_z`.

Using PETSc solvers
~~~~~~~~~~~~~~~~~~~

When using PETSc, all terms of equation (:eq:`invert_expanded`) is being
used when inverting to find :math:`f`. Note that when using PETSc, we
are not Fourier decomposing in the :math:`z`-direction, so it may take
substantially longer time to find the solution. As with the tridiagonal
solver, the fields are being sliced in the :math:`y`-direction, and a
solution is being found for one :math:`y` plane at the time.

Before solving, equation (:eq:`invert_expanded`) is rewritten to the
form
:math:`A{{\boldsymbol{x}}} ={{\boldsymbol{b}}}`
(however, the full :math:`A` is not expanded in memory). To do this, a
row :math:`i` in the matrix :math:`A` is indexed from bottom left of the
two dimensional field :math:`= (0,0) = 0` to top right
:math:`= (\texttt{meshx}-1,
\texttt{meshz}-1) = \texttt{meshx}\cdot\texttt{meshz}-1` of the two
dimensional field. This is done in such a way so that a row :math:`i` in
:math:`A` increments by :math:`1` for an increase of :math:`1` in the
:math:`z-`\ direction, and by :math:`\texttt{meshz}` for an increase of
:math:`1` in the :math:`x-`\ direction, where the variables
:math:`\texttt{meshx}` and :math:`\texttt{meshz}` represents the highest
value of the field in the given direction.

| Similarly to equation (:eq:`discretized_laplace`), the discretised
  version of equation (:eq:`invert_expanded`) can be written. Doing the
  same for the full two dimensional case yields

0.45 Second order approximation

.. math::

       \; & c_{i,j} f_{i,j} \\
           &+ c_{i-1,j-1} f_{i-1,j-1} + c_{i-1,j} f_{i-1,j} \\
           &+ c_{i-1,j+1} f_{i-1,j+1} + c_{i,j-1} f_{i,j-1} \\
           &+ c_{i,j+1} f_{i,j+1} + c_{i+1,j-1} f_{i+1,j-1} \\
           &+ c_{i+1,j} f_{i+1,j} + c_{i+1,j+1} f_{i+1,j+1} \\
       =& b_{i,j}

0.45 Fourth order approximation

.. math::

       \; & c_{i,j} f_{i,j} \\
           &+ c_{i-2,j-2} f_{i-2,j-2} + c_{i-2,j-1} f_{i-2,j-1} \\
           &+ c_{i-2,j} f_{i-2,j} + c_{i-2,j+1} f_{i-2,j+1} \\
           &+ c_{i-2,j+2} f_{i-2,j+2} + c_{i-1,j-2} f_{i-1,j-2} \\
           &+ c_{i-1,j-1} f_{i-1,j-1} + c_{i-1,j} f_{i-1,j} \\
           &+ c_{i-1,j-1} f_{i-1,j-1} + c_{i-1,j} f_{i-1,j} \\
           &+ c_{i-1,j+1} f_{i-1,j+1} + c_{i-1,j+2} f_{i-1,j+2} \\
           &+ c_{i,j-2} f_{i,j-2} + c_{i,j-1} f_{i,j-1} \\
           &+ c_{i,j+1} f_{i,j+1} + c_{i,j+2} f_{i,j+2} \\
           &+ c_{i+1,j-2} f_{i+1,j-2} + c_{i+1,j-1} f_{i+1,j-1} \\
           &+ c_{i+1,j} f_{i+1,j} + c_{i+1,j+1} f_{i+1,j+1} \\
           &+ c_{i+1,j+2} f_{i+1,j+2} + c_{i+2,j-2} f_{i+2,j-2} \\
           &+ c_{i+2,j-1} f_{i+2,j-1} + c_{i+2,j} f_{i+2,j} \\
           &+ c_{i+2,j+1} f_{i+2,j+1} + c_{i+2,j+2} f_{i+2,j+2} \\
       =& b_{i,j}

| 
| To determine the coefficient for each node point, it is convenient to
  introduce some quantities

  .. math::

         &A_0 = a(x,y_{\text{current}},z)& &A_1 = dg^{xx}&\\ &A_2 = dg^{zz}& &A_3 =
         2dg^{xz}&

   In addition, we have

0.45 Second order approximation (5-point stencil)

.. math::

       \texttt{ddx\_c} = \frac{\texttt{c2}_{x+1} - \texttt{c2}_{x-1} }{
       2\texttt{c1}\text{d}x}\\
           \texttt{ddz\_c} = \frac{\texttt{c2}_{z+1} - \texttt{c2}_{z-1} }{
       2\texttt{c1}\text{d}z}

0.45 Fourth order approximation (9-point stencil)

.. math::

       \texttt{ddx\_c} = \frac{-\texttt{c2}_{x+2} + 8\texttt{c2}_{x+1} -
       8\texttt{c2}_{x-1} + \texttt{c2}_{x-1} }{ 12\texttt{c1}\text{d}x}\\
           \texttt{ddz\_c} = \frac{-\texttt{c2}_{z+2} + 8\texttt{c2}_{z+1} -
       8\texttt{c2}_{z-1} + \texttt{c2}_{z-1} }{ 12\texttt{c1}\text{d}z}

| 
| This gives

  .. math::

         &A_4 = dG^x + g^{xx}\texttt{ddx\_c} + g^{xz}\texttt{ddz\_c}& &A_5 = dG^z +
         g^{xz}\texttt{ddx\_c} + g^{xx}\texttt{ddz\_c}&

  The coefficients :math:`c_{i+m,j+n}` are finally being set according
  to the appropriate order of discretisation. The coefficients can be
  found in the file ``petsc_laplace.cxx``.

Example: The 5-point stencil
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now consider the 5-point stencil for a mesh with :math:`3` inner
points in the :math:`x`-direction, and :math:`3` inner points in the
:math:`z`-direction. The :math:`z` direction will be periodic, and the
:math:`x` direction will have the boundaries half between the grid-point
and the first ghost point (see figure [fig:lapl\_inv\_mesh]).

Applying the :math:`5`-point stencil to point :math:`f_{22}` this mesh
will result in figure [fig:lapl\_inv\_mesh\_w\_stencil].

We want to solve a problem on the form
:math:`A{{\boldsymbol{x}}}={{\boldsymbol{b}}}`. We
will order :math:`{{\boldsymbol{x}}}` in a row-major order
(so that :math:`z` is varying faster than :math:`x`). Further, we put
the inner :math:`x` boundary points first in
:math:`{{\boldsymbol{x}}}`, and the outer :math:`x` boundary
points last in :math:`{{\boldsymbol{x}}}`. The matrix problem
for our mesh can then be written like in figure [fig:lapl\_inv\_matrix].

As we are using a row-major implementation, the global indices of the
matrix will be as in figure [fig:lapl\_inv\_global]

.. _sec-equations:

Fluid equations
===============

Once you have tried some example codes, and generally got the hang of
running BOUT++ and analysing the results, there will probably come a
time when you want to change the equations being solved. This section
uses the ideal MHD equations as an example, demonstrating how a BOUT++
physics module is put together. It assumes you have a working knowledge
of C or C++, but you don’t need to be an expert - most of the messy code
is hidden away from the physics module. There are several good books on
C and C++, but I’d recommend online tutorials over books because there
are a lot more of them, they’re quicker to scan through, and they’re
cheaper.

When going through this section, it may help to refer to the finished
code, which is given in the file ``mhd.cxx`` in the BOUT++ examples
directory. The equations to be solved are:

.. math::

   {{\frac{\partial \rho}{\partial t}}} =& -\mathbf{v}\cdot\nabla\rho - \rho\nabla\cdot\mathbf{v} \\
       {{\frac{\partial p}{\partial t}}} =& -\mathbf{v}\cdot\nabla p - \gamma p\nabla\cdot\mathbf{v} \\
       {{\frac{\partial \mathbf{v}}{\partial t}}} =& -\mathbf{v}\cdot\nabla\mathbf{v} +
       \frac{1}{\rho}(-\nabla p +
       (\nabla\times\mathbf{B})\times\mathbf{B}) \\ {{\frac{\partial \mathbf{B}}{\partial t}}} =&
       \nabla\times(\mathbf{v}\times\mathbf{B})

There are two ways to specify a set of equations to solve in BOUT++.
For advanced users, an object-oriented interface is available and
described in :ref:`sec-newapi`. The simplest way to start is to use a
C-like interface and define two functions:

::

    int physics_init(bool restarting) {
      return 0;
    }

    int physics_run(BoutReal t) {
      return 0;
    }

The first of these is called once at the start of the simulation, and
should set up the problem, specifying which variables are to be evolved.
The argument ``restarting`` is false the first time a problem is run,
and true if loading the state from a restart file.

The second function ``physics_run`` is called every time-step, and
should calculate the time-derivatives for a given state. In both cases
returning non-zero tells BOUT++ that an error occurred.

Variables
---------

We need to define the variables to evolve as global variables (so they
can be used in ``physics_init`` and ``physics_run``.

For ideal MHD, we need two 3D scalar fields density :math:`\rho` and
pressure :math:`p`, and two 3D vector fields velocity :math:`v`, and
magnetic field :math:`B`:

::

    Field3D rho, p; // 3D scalar fields
    Vector3D v, B;  // 3D vector fields

    int physics_init(bool restarting) {
    }

Scalar and vector fields behave much as you would expect: ``Field3D``
objects can be added, subtracted, multiplied, divided and exponentiated,
so the following examples are all valid operations:

::

    Field3D a, b, c;
    BoutReal r;

    a = b + c; a = b - c;
    a = b * c; a = r * b;
    a = b / c; a = b / r; a = r / b;
    a = b ^ c; a = b ^ r; a = r ^ b;

Similarly, vector objects can be added/subtracted from each other,
multiplied/divided by scalar fields and real numbers, for example:

::

    Vector3D a, b, c;
    Field3D f;
    BoutReal r;

    a = b + c; a = b - c;
    a = b * f; a = b * r;
    a = b / f; a = b / r;

In addition the dot and cross products are represented by ``*`` and
:math:`\wedge` \ symbols:

::

    Vector3D a, b, c;
    Field3D f;

    f = a * b // Dot-product
    a = b ^ c // Cross-product

For both scalar and vector field operations, so long as the result of an
operation is of the correct type, the usual C/C++ shorthand notation can
be used:

::

    Field3D a, b;
    Vector3D v, w;

    a += b; v *= a; v -= w; v ^= w; // valid
    v *= w; // NOT valid: result of dot-product is a scalar

Evolution equations
-------------------

At this point we can tell BOUT++ which variables to evolve, and where
the state and time-derivatives will be stored. This is done using the
``bout_solve(variable, name)`` function in ``physics_init``:

::

    int physics_init(bool restarting) {
      bout_solve(rho, "density");
      bout_solve(p,   "pressure");
      bout_solve(v,   "v");
      bout_solve(B,   "B");

      return 0;
    }

The name given to this function will be used in the output and restart
data files. These will be automatically read and written depending on
input options (see :ref:`sec-options`). Input options based on these
names are also used to initialise the variables.

If the name of the variable in the output file is the same as the
variable name, you can use a shorthand macro. In this case, we could use
this shorthand for ``v`` and ``B``:

::

    SOLVE_FOR(v);
    SOLVE_FOR(B);

To make this even shorter, we can use macros ``SOLVE_FOR2``,
``SOLVE_FOR3``, ..., ``SOLVE_FOR6`` to shorten our initialisation code
to

::

    int physics_init(bool restarting) {
      bout_solve(rho, "density");
      bout_solve(p,   "pressure");
      SOLVE_FOR2(v, B);

      return 0;
    }

The equations to be solved can now be written in the ``physics_run``
function. The value passed to the function (``BoutReal t``) is the
simulation time - only needed if your equations contain time-dependent
sources or similar terms. To refer to the time-derivative of a variable
``var``, use ``ddt(var)``. The ideal MHD equations can be written as:

::

    int physics_run(BoutReal t) {
      ddt(rho) = -V_dot_Grad(v, rho) - rho*Div(v);
      ddt(p) = -V_dot_Grad(v, p) - gamma*p*Div(v);
      ddt(v) = -V_dot_Grad(v, v) + ( (Curl(B)^B) - Grad(p) ) / rho;
      ddt(B) = Curl(v^B);
    }

Where the differential operators ``vector = Grad(scalar)``,
``scalar = Div(vector)``, and ``vector = Curl(vector)`` are used. For
the density and pressure equations, the
:math:`\mathbf{v}\cdot\nabla\rho` term could be written as
``v*Grad(rho)``, but this would then use central differencing in the
Grad operator. Instead, the function ``V_dot_Grad`` uses upwinding
methods for these advection terms. In addition, the ``Grad`` function
will not operate on vector objects (since result is neither scalar nor
vector), so the :math:`\mathbf{v}\cdot\nabla\mathbf{v}` term CANNOT be
written as ``v*Grad(v)``.

.. _sec-inputopts:

Input options
-------------

Note that in the above equations the extra parameter ``gamma`` has been
used. To enable this to be set in the input options file (see
:ref:`sec-options`), we use the ``options`` object in the
initialisation function:

::

    BoutReal gamma;

    int physics_init(bool restarting) {
      Options *globalOptions = Options::getRoot();
      Options *options = globalOptions->getSection("mhd");

      options->get("gamma", gamma, 5.0/3.0);

This specifies that an option called “gamma” in a section called “mhd”
should be put into the variable ``gamma``. If the option could not be
found, or was of the wrong type, the variable should be set to a default
value of :math:`5/3`. The value used will be printed to the output file,
so if gamma is not set in the input file the following line will appear:

::

          Option mhd / gamma = 1.66667 (default)

This function can be used to get integers and booleans. To get strings,
there is the function (``char* options.getString(section, name)``. To
separate options specific to the physics model, these options should be
put in a separate section, for example here the “mhd” section has been
specified. To save having to write the section name for every option,
there is the ``setSection`` function:

::

    BoutReal gamma;
    int someint;

    int physics_init(bool restarting) {
      Options *globalOptions = Options::getRoot();
      Options *options = globalOptions->getSection("mhd");

      options->get("gamma", gamma, 5.0/3.0);
      options->get("someint", someint, 0);

Most of the time, the name of the variable (e.g. ``gamma``) will be the
same as the identifier in the options file (“gamma”). In this case,
there is the macro

::

    OPTION(options, gamma, 5.0/3.0);

which is equivalent to

::

    options->get("gamma", gamma, 5.0/3.0);

See :ref:`sec-options` for more details of how to use the input
options.

Communication
-------------

If you plan to run BOUT++ on more than one processor, any operations
involving y derivatives will require knowledge of data stored on other
processors. To handle the necessary parallel communication, there is the
``mesh->communicate`` function. This takes care of where the data needs
to go to/from, and only needs to be told which variables to transfer.

If you only need to communicate a small number (up to 5 currently) of
variables then just call the ``mesh->communicate`` function directly.
For the MHD code, we need to communicate the variables ``rho,p,v,B`` at
the beginning of the ``physics_run`` function before any derivatives are
calculated:

::

    int physics_run(BoutReal t) {
      mesh->communicate(rho, p, v, B);

If you need to communicate lots of variables, or want to change at
run-time which variables are evolved (e.g. depending on input options),
then you can create a group of variables and communicate them later. To
do this, first create a ``FieldGroup`` object , in this case called
``comms`` , then use the add method. This method does no communication,
but records which variables to transfer when the communication is done
later.

::

    FieldGroup comms;

    int physics_init() {
      .
      .
      .
      comms.add(rho);
      comms.add(p);
      comms.add(v);
      comms.add(B);

      return 0;
    }

The ``comms.add()`` routine can be given up to 6 variables at once
(there’s no practical limit on the total number of variables which are
added to a ``FieldGroup`` ), so this can be shortened to

::

    FieldGroup comms;

    int physics_init() {
      .
      .
      .
      comms.add(rho, p, v, B);

      return 0;
    }

To perform the actual communication, call the ``mesh->communicate``
function with the group. In this case we need to communicate all these
variables before performing any calculations, so call this function at
the start of the ``physics_run`` routine:

::

    int physics_run(BoutReal t) {
      mesh->communicate(comms);
      .
      .
      .

In many situations there may be several groups of variables which can be
communicated at different times. The function ``mesh->communicate``
consists of a call to ``mesh->send`` followed by ``mesh->wait`` which
can be done separately to interleave calculations and communications.
This will speed up the code if parallel communication bandwidth is a
problem for your simulation.

In our MHD example, the calculation of ``ddt(rho)`` and ``ddt(p)`` does
not require ``B``, so we could first communicate ``rho``, ``p``, and
``v``, send ``B`` and do some calculations whilst communications are
performed:

::

    int physics_run(BoutReal t) {
      mesh->communicate(rho, p, v); // sends and receives rho, p and v
      comm_handle ch = mesh->send(B);// only send B

      ddt(rho) = ...
      ddt(p) = ...

      mesh->wait(ch); // now wait for B to arrive

      ddt(v) = ...
      ddt(B) = ...

      return 0;
    }

This scheme is not used in ``mhd.cxx``, partly for clarity, and partly
because currently communications are not a significant bottleneck (too
much inefficiency elsewhere!).

When a differential is calculated, points on neighbouring cells are
assumed to be in the guard cells. There is no way to calculate the
result of the differential in the guard cells, and so after every
differential operator the values in the guard cells are invalid.
Therefore, if you take the output of one differential operator and use
it as input to another differential operator, you must perform
communications (and set boundary conditions) first. See
:ref:`sec-diffops`.

Boundary conditions
-------------------

All evolving variables have boundary conditions applied automatically
after the ``physics_run`` has finished. Which condition is applied
depends on the options file settings (see :ref:`sec-bndryopts`). If
you want to disable this and apply your own boundary conditions then set
boundary condition to ``none`` in the ``BOUT.inp`` options file.

In addition to evolving variables, it’s sometimes necessary to impose
boundary conditions on other quantities which are not explicitly
evolved.

The simplest way to set a boundary condition is to specify it as text,
so to apply a Dirichlet boundary condition:

::

      Field3D var;
      ...
      var.applyBoundary("dirichlet");

The format is exactly the same as in the options file. Each time this is
called it must parse the text, create and destroy boundary objects. To
avoid this overhead and have different boundary conditions for each
region, it’s better to set the boundary conditions you want to use first
in ``physics_init``, then just apply them every time:

::

    Field3D var;

    int physics_init() {
      ...
      var.setBoundary("myVar");
      ...
    }

    int physics_run(BoutReal t) {
      ...
      var.applyBoundary();
      ...
    }

This will look in the options file for a section called ``[myvar]``
(upper or lower case doesn’t matter) in the same way that evolving
variables are handled. In fact this is precisely what is done: inside
``bout_solve`` (or ``SOLVE_FOR``) the ``setBoundary`` method is called,
and then after ``physics_run`` the applyBoundary() method is called on
each evolving variable. This method also gives you the flexibility to
apply different boundary conditions on different boundary regions (e.g.
radial boundaries and target plates); the first method just applies the
same boundary condition to all boundaries.

Another way to set the boundaries is to copy them from another variable:

::

    Field3D a, b;
      ...
      a.setBoundaryTo(b); // Copy b's boundaries into a
      ...

Custom boundary conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The boundary conditions supplied with the BOUT++ library cover the most
common situations, but cannot cover all of them. If the boundary
condition you need isn’t available, then it’s quite straightforward to
write your own. First you need to make sure that your boundary condition
isn’t going to be overwritten. To do this, set the boundary condition to
“none” in the BOUT.inp options file, and BOUT++ will leave that boundary
alone. For example:

::

    [P]
    bndry_all = dirichlet
    bndry_xin = none
    bndry_xout = none

would set all boundaries for the variable “P” to zero value, except for
the X inner and outer boundaries which will be left alone for you to
modify.

To set an X boundary condition, it’s necessary to test if the processor
is at the left boundary (first in X), or right boundary (last in X).
Note that it might be both if ``NXPE = 1``, or neither if ``NXPE > 2``.

::

      Field3D f;
      ...
      if(mesh->firstX()) {
        // At the left of the X domain
        // set f[0:1][*][*] i.e. first two points in X, all Y and all Z
        for(int x=0; x < 2; x++)
          for(int y=0; y < mesh->ngy; y++)
            for(int z=0; z < mesh->ngz; z++) {
              f[x][y][z] = ...
            }
      }
      if(mesh->lastX()) {
        // At the right of the X domain
        // Set last two points in X
        for(int x=mesh->ngx-2; x < mesh->ngx; x++)
          for(int y=0; y < mesh->ngy; y++)
            for(int z=0; z < mesh->ngz; z++) {
              f[x][y][z] = ...
            }
      }

note the size of the local mesh including guard cells is given by
``mesh->ngx``, ``mesh->ngy``, and ``mesh->ngz``. The functions
``mesh->firstX()`` and ``mesh->lastX()`` return true only if the current
processor is on the left or right of the X domain respectively.

Setting custom Y boundaries is slightly more complicated than X
boundaries, because target or limiter plates could cover only part of
the domain. Rather than use a ``for`` loop to iterate over the points in
the boundary, we need to use a more general iterator:

::

      Field3D f;
      ...
      RangeIterator it = mesh->iterateBndryLowerY();
      for(it.first(); !it.isDone(); it++) {
        // it.ind contains the x index
        for(int y=2;y>=0;y--)  // Boundary width 3 points
          for(int z=0;z<mesh->ngz;z++) {
            ddt(f)[it.ind][y][z] = 0.;  // Set time-derivative to zero in boundary
          }
      }

This would set the time-derivative of f to zero in a boundary of width 3
in Y (from 0 to 2 inclusive). In the same way
``mesh->iterateBndryUpperY()`` can be used to iterate over the upper
boundary:

::

      RangeIterator it = mesh->iterateBndryUpperY();
      for(it.first(); !it.isDone(); it++) {
        // it.ind contains the x index
        for(int y=mesh->ngy-3;y<mesh->ngy;y--)  // Boundary width 3 points
          for(int z=0;z<mesh->ngz;z++) {
            ddt(f)[it.ind][y][z] = 0.;  // Set time-derivative to zero in boundary
          }
      }

Initial profiles
----------------

Up to this point the code is evolving total density, pressure etc. This
has advantages for clarity, but has problems numerically: For small
perturbations, rounding error and tolerances in the time-integration
mean that linear dispersion relations are not calculated correctly. The
solution to this is to write all equations in terms of an initial
“background” quantity and a time-evolving perturbation, for example
:math:`\rho(t) arrow \rho_0 +
\tilde{\rho}(t)`. For this reason, **the initialisation of all
variables passed to the ``bout_solve`` function is a combination of
small-amplitude gaussians and waves; the user is expected to have
performed this separation into background and perturbed quantities.**

To read in a quantity from a grid file, there is the ``grid.get``
function:

::

    Field2D Ni0; // Background density

    int physics_init(bool restarting) {
      ...
      mesh->get(Ni0, "Ni0");
      ...
    }

As with the input options, most of the time the name of the variable in
the physics code will be the same as the name in the grid file to avoid
confusion. In this case, you can just use

::

    GRID_LOAD(Ni0);

which is equivalent to

::

    mesh->get(Ni0, "Ni0");

Output variables
----------------

BOUT++ always writes the evolving variables to file, but often it’s
useful to add other variables to the output. For convenience you might
want to write the normalised starting profiles or other non-evolving
values to file. For example:

::

      Field2D Ni0;
      ...
      GRID_LOAD(Ni0);
      dump.add(Ni0, "Ni0", 0);

where the ’0’ at the end means the variable should only be written to
file once at the start of the simulation. For convenience there are some
macros e.g.

::

      SAVE_ONCE(Ni0);

is equivalent to

::

      dump.add(Ni0, "Ni0", 0);

In some situations you might also want to write some data to a different
file. To do this, create a Datafile object:

::

    Datafile mydata;

in physics\_init, you then:

#. (optional) Initialise the file, passing it the options to use. If you
   skip this step, default (sane) options will be used. This just allows
   you to enable/disable, use parallel I/O, set whether files are opened
   and closed every time etc.

   ::

       mydata = Datafile(Options::getRoot()->getSection("mydata"));

   which would use options in a section [mydata] in BOUT.inp

#. Open the file for writing

   ::

       mydata.openw("mydata.nc")

   By default this only specifies the file name; actual opening of the
   file happens later when the data is written. If you are not using
   parallel I/O, the processor number is also inserted into the file
   name before the last “.”, so mydata.nc” becomes “mydata.0.nc”,
   “mydata.1.nc” etc. The file format used depends on the extension, so
   “.nc” will open NetCDF, and “.hdf5” or “.h5” an HDF5 file.

   (see e.g. src/fileio/datafile.cxx line 139, which calls
   src/fileio/dataformat.cxx line 23, which then calls the file format
   interface e.g. src/fileio/impls/netcdf/nc\_format.cxx line 172).

#. Add variables to the file

   ::

       mydata.add(variable, "name") ;  // Not evolving. Every time the file is
       written, this will be overwritten
       mydata.add(variable2, "name2", 1); // Evolving. Will output a sequence of values

Whenever you want to write values to the file, for example in
physics\_run or a monitor, just call

::

    mydata.write();

To collect the data afterwards, you can specify the prefix to collect.
In Python:

::

    >>> var = collect("name", prefix="mydata")

or in IDL:

::

    IDL> var = collect(var="name", prefix="mydata")

By default the prefix is “BOUT.dmp”.

Fluid equations 2: reduced MHD
==============================

The MHD example presented previously covered some of the functions
available in BOUT++, which can be used for a wide variety of models.
There are however several other significant functions and classes which
are commonly used, which will be illustrated using the
``reconnect-2field`` example. This is solving equations for
:math:`A_{||}` and vorticity :math:`U`

.. math::

   {{\frac{\partial U}{\partial t}}} =& -\frac{1}{B}\mathbf{b}_0\times\nabla\phi\cdot\nabla U + B^2
       \nabla_{||}(j_{||} / B) \\ {{\frac{\partial A_{||}}{\partial t}}} =&
       -\frac{1}{\hat{\beta}}\nabla_{||}\phi - \eta\frac{1}{\hat{\beta}} j_{||}

with :math:`\phi` and :math:`j_{||}` given by

.. math::

   U =& \frac{1}{B}\nabla_\perp^2\phi \\ j_{||} =& -\nabla_\perp^2 A_{||}

First create the variables which are going to be evolved, ensure
they’re communicated

::

    Field3D U, Apar; // Evolving variables

    int physics_init(bool restarting) {

      SOLVE_FOR2(U, Apar);
    }

    int physics_run(BoutReal t) {
      mesh->communicate(U, Apar);

    }

In order to calculate the time derivatives, we need the auxiliary
variables :math:`\phi` and :math:`j_{||}`. Calculating :math:`j_{||}`
from :math:`A_{||}` is a straightforward differential operation, but
getting :math:`\phi` from :math:`U` means inverting a Laplacian.

::

    Field3D U, Apar;
    Field3D phi, jpar; // Auxilliary variables

    int physics_init(bool restarting) {
      SOLVE_FOR2(U, Apar);
      SAVE_REPEAT2(phi, jpar); // Save variables in output file
      return 0;
    }

    int physics_run(BoutReal t) {
      phi = invert_laplace(mesh->Bxy*U, phi_flags); // Solve for phi
      mesh->communicate(U, Apar, phi);  // Communicate phi
      jpar = -Delp2(Apar);     // Calculate jpar
      mesh->communicate(jpar); // Communicate jpar
      return 0;
    }

Note that the Laplacian inversion code takes care of boundary regions,
so ``U`` doesn’t need to be communicated first. The differential
operator ``Delp2`` , like all differential operators, needs the values
in the guard cells and so ``Apar`` needs to be communicated before
calculating ``jpar`` . Since we will need to take derivatives of
``jpar`` later, this needs to be communicated as well.

::

    int physics_run(BoutReal t) {
      ...
      mesh->communicate(jpar);

      ddt(U) = -b0xGrad_dot_Grad(phi, U) + SQ(mesh->Bxy)*Grad_par(Jpar / mesh->Bxy)
      ddt(Apar) = -Grad_par(phi) / beta_hat - eta*jpar / beta_hat; }

.. _sec-printing:

Printing messages/warnings
--------------------------

In order to print to screen and/or a log file, the object ``output`` is
provided. This provides two different ways to write output: the C
(``printf``) way, and the C++ stream way. This is because each method
can be clearer in different circumstances, and people have different
tastes in these matters.

The C-like way (which is the dominant way in BOUT++) is to use the
``write`` function, which works just like ``printf``, and takes all the
same codes (it uses ``sprintf`` internally).

::

    output.write(const char *format, ...)

For example:

::

    output.write("This is an integer: %d, and this a real: %e\n", 5, 2.0)

For those who prefer the C++ way of doing things, a completely
equivalent way is to treat ``output`` as you would ``cout``:

::

    output << "This is an integer: " << 5 << ", and this a real: " << 2.0 << endl;

which will produce exactly the same result as the ``output.write`` call
above.

On all processors, anything sent to ``output`` will be written to a log
file called ``BOUT.log.#`` with # replaced by the processor number. On
processor 0, anything written to the output will be written to screen
(stdout), in addition to the log file. Unless there is a really good
reason not to, please use this ``output`` object when writing text
output.

Error handling
--------------

Finding where bugs have occurred in a (fairly large) parallel code is a
difficult problem. This is more of a concern for developers of BOUT++
(see the developers manual), but it is still useful for the user to be
able to hunt down bug in their own code, or help narrow down where a bug
could be occurring.

If you have a bug which is easily reproduceable i.e. it occurs almost
immediately every time you run the code, then the easiest way to hunt
down the bug is to insert lots of ``output.write`` statements (see
:ref:`sec-printing`). Things get harder when a bug only occurs after
a long time of running, and/or only occasionally. For this type of
problem, a useful tool can be the message stack. At the start of a
section of code, put a message onto the stack:

::

       msg_stack.push("Some message here");

which can also take arguments in ``printf`` format, as with
``output.write``. At the end of the section of code, take the message
off the stack again:

::

       msg_stack.pop();

If an error occurs, the message stack is printed out, and this can then
help track down where the error originated.

.. _sec-newapi:

Object-orientated interface
===========================

| If you prefer to create classes rather than global variables and C
  functions for your physics model, this can be done using a (somewhat
  experimental) interface. To see the difference, compare
  ``examples/advect1d/gas_compress.cxx`` with
| ``examples/advect1d-newapi/gas_compress.cxx``. The disadvantage of
  this interface is that it’s marginally more complicated to set up, but
  it has several advantages: It makes splitting the model into multiple
  files easier (sharing global variables is a pain), models can be
  combined together to enable coupling of models, and BOUT++ can be more
  easily used alongside other libraries. For large models, it’s
  recommended to use this method. Converting C-style interface to a
  class is also quite straightforward, and discussed below.

In a header file (e.g. ``examples/advect1d-newapi/gas_compress.hxx``),
first put

::

    #include <bout/physicsmodel.hxx>

(do NOT include ``boutmain.hxx``, as that defines the C-like interface
and a ``main()`` function).

Next define a class which inherits from ``PhysicsModel``

::

    class GasCompress : public PhysicsModel {
    protected:
      int init(bool restarting);
      int rhs(BoutReal t);
    private:
      // Evolving variables, parameters etc. here
    };

As a minimum, you need to define the initialisation function ``init``
(it’s a pure virtual member of PhysicsModel, so if you don’t you’ll get
a compile-time error). Any variables being evolved should now be members
of this class. If you are converting a C-style model, just move all the
global variables into the ``private`` section.

Next create a source file (e.g.
``examples/advect1d-newapi/gas_compress.cxx``, which includes your
header file

::

    #include "gas_compress.hxx"

Then implement the init and rhs functions:

::

    int GasCompress::init(bool restarting) {
      ...
    }

    int GasCompress::rhs(BoutReal t) {
      ...
    }

To convert simple physics models, just rename ``physics_init`` to
``YourModel::init`` , and ``physics_run`` to ``YourModel::run`` .

Finally, you need to create a ``main()`` function for your code. The
easiest way to do this is to use the macro ``BOUTMAIN`` :

::

    BOUTMAIN(GasCompress);

This is defined in ``include/bout/physicsmodel.hxx``, and expands to

::

      int main(int argc, char **argv) {
        BoutInitialise(argc, argv); // Initialise BOUT++

        GasCompress *model = new GasCompress(); // Create a model

        Solver *solver = Solver::create(); // Create a solver
        solver->setModel(model); // Specify the model to solve
        solver->addMonitor(bout_monitor); // Monitor the solver

        solver->solve(); // Run the solver

        delete model;
        delete solver;
        BoutFinalise(); // Finished with BOUT++
        return 0;
      }

If you like, you can define your own ``main()`` function, making it
easier to combine BOUT++ with other libraries.

.. _sec-diffops:

Differential operators
======================

There are a huge number of possible ways to perform differencing in
computational fluid dynamics, and BOUT++ is intended to be able to
implement a large number of them. This means that the way differentials
are handled internally is quite involved; see the developer’s manual for
full gory details. Much of the time this detail is not all that
important, and certainly not while learning to use BOUT++. Default
options are therefore set which work most of the time, so you can start
using the code without getting bogged down in these details.

In order to handle many different differencing methods and operations,
many layers are used, each of which handles just part of the problem.
The main division is between differencing methods (such as 4th-order
central differencing), and differential operators (such as
:math:`\nabla_{||}`).

.. _sec-diffmethod:

Differencing methods
--------------------

Methods are implemented on 5-point stencils, and are divided into three
categories:

-  Central-differencing methods, for diffusion operators
   :math:`\frac{df}{dx}`, :math:`\frac{d^2f}{dx^2}`. Each method has a
   short code, and currently include

   -  ``C2``: 2\ :math:`^{nd}` order :math:`f_{-1} - 2f_0 + f_1`

   -  ``C4``: 4\ :math:`^{th}` order
      :math:`(-f_{-2} + 16f_{-1} - 30f_0 + 16f_1 - f_2)/12`

   -  ``W2``: 2\ :math:`^{nd}` order CWENO

   -  ``W3``: 3\ :math:`^{rd}` order CWENO

   -  ``FFT``: Fourier Transform method in Z (axisymmetric) direction
      only

-  Upwinding methods for advection operators :math:`v_x\frac{df}{dx}`

   -  ``U1``: 1\ :math:`^{st}` order upwinding

   -  ``U4``: 4\ :math:`^{th}` order upwinding

   -  ``W3``: 3\ :math:`^{rd}` order Weighted Essentially
      Non-Oscillatory (WENO):raw-latex:`\cite{jiang-1997}`

-  Flux conserving and limiting methods for terms of the form
   :math:`\frac{d}{dx}(v_x f)`

   -  ``SPLIT``: split into upwind and central terms
      :math:`\frac{d}{dx}(v_x f) = v_x\frac{df}{dx} + f\frac{dv_x}{dx}`

   -  ``NND``: Non-oscillatory, containing No free parameters and
      Dissipative (NND) scheme:raw-latex:`\cite{nnd-2010}`

Both of these methods avoid overshoots (Gibbs phenomena) at sharp
gradients such as shocks, but the simple 1st-order method has very large
artificial diffusion. WENO schemes are a development of the ENO
reconstruction schemes which combine good handling of sharp-gradient
regions with high accuracy in smooth regions.

To use these differencing operators directly, add the following to the
top of your physics module

::

    #include <derivs.hxx>

+--------------+-----------------------------------------------+
| Function     | Formula                                       |
+==============+===============================================+
| DDX(f)       | :math:`\partial f / \partial x`               |
+--------------+-----------------------------------------------+
| DDY(f)       | :math:`\partial f / \partial y`               |
+--------------+-----------------------------------------------+
| DDZ(f)       | :math:`\partial f / \partial z`               |
+--------------+-----------------------------------------------+
| D2DX2(f)     | :math:`\partial^2 f / \partial x^2`           |
+--------------+-----------------------------------------------+
| D2DY2(f)     | :math:`\partial^2 f / \partial y^2`           |
+--------------+-----------------------------------------------+
| D2DZ2(f)     | :math:`\partial^2 f / \partial z^2`           |
+--------------+-----------------------------------------------+
| D2DX4(f)     | :math:`\partial^4 f / \partial x^4`           |
+--------------+-----------------------------------------------+
| D2DY4(f)     | :math:`\partial^4 f / \partial y^4`           |
+--------------+-----------------------------------------------+
| D2DZ4(f)     | :math:`\partial^4 f / \partial z^4`           |
+--------------+-----------------------------------------------+
| D2DXDZ(f)    | :math:`\partial^2 f / \partial x\partial z`   |
+--------------+-----------------------------------------------+
| D2DYDZ(f)    | :math:`\partial^2 f / \partial y\partial z`   |
+--------------+-----------------------------------------------+
| VDDX(f, g)   | :math:`f \partial g / \partial x`             |
+--------------+-----------------------------------------------+
| VDDY(f, g)   | :math:`f \partial g / \partial y`             |
+--------------+-----------------------------------------------+
| VDDZ(f, g)   | :math:`f \partial g / \partial z`             |
+--------------+-----------------------------------------------+
| FDDX(f, g)   | :math:`\partial/\partial x( f * g )`          |
+--------------+-----------------------------------------------+
| FDDY(f, g)   | :math:`\partial/\partial x( f * g )`          |
+--------------+-----------------------------------------------+
| FDDZ(f, g)   | :math:`\partial/\partial x( f * g )`          |
+--------------+-----------------------------------------------+

Table: Coordinate derivatives

By default the method used will be the one specified in the options
input file (see :ref:`sec-diffmethodoptions`), but most of these
methods can take an optional ``DIFF\_METHOD`` argument, specifying
exactly which method to use.

Non-uniform meshes
------------------

**examples/test-nonuniform seems to not work?** Setting
``non_uniform = true`` in the BOUT.inp options file enables corrections
to second derivatives in :math:`X` and :math:`Y`. This correction is
given by writing derivatives as:

.. math::

   {{\frac{\partial f}{\partial x}}} \simeq \frac{1}{\Delta x} {{\frac{\partial f}{\partial i}}}

where :math:`i` is the cell index number. The second derivative is
therefore given by

.. math::

   \frac{\partial^2 f}{\partial x^2} \simeq \frac{1}{\Delta x^2}\frac{\partial^2
   f}{\partial i^2} + \frac{1}{\Delta x}{{\frac{\partial f}{\partial x}}} \cdot
   {{\frac{\partial }{\partial i}}}(\frac{1}{\Delta x})

The correction factor :math:`\partial/\partial i(1/\Delta x)` can
be calculated automatically, but you can also specify ``d2x`` in the
grid file which is

.. math::

   \texttt{d2x} = {{\frac{\partial \Delta x}{\partial i}}} = \frac{\partial^2 x}{\partial i^2}

The correction factor is then calculated from ``d2x`` using

.. math::

   {{\frac{\partial }{\partial i}}}(\frac{1}{\Delta x}) = -\frac{1}{\Delta x^2} {{\frac{\partial \Delta x}{\partial i}}}

General operators
-----------------

These are differential operators which are for a general coordinate
system.

.. math::

   \begin{array}{rclrcl}
   \mathbf{v} =& \nabla f &\qquad {\texttt{Vector}} =& {\texttt{Grad(Field)}} \\
   f =& \nabla\cdot\mathbf{a} &\qquad {\texttt{Field}} =& {\texttt{Div(Vector)}} \\
   \mathbf{v} =& \nabla\times\mathbf{a} &\qquad {\texttt{Vector}} =&
   {\texttt{Curl(Vector)}} \\
   f =& \mathbf{v}\cdot\nabla g &\qquad {\texttt{Field}} =& {\texttt{V\_dot\_Grad(Vector,
   Field)}} \\
   \mathbf{v} =& \mathbf{a}\cdot\nabla\mathbf{c} &\qquad {\texttt{Vector}} =&
   {\texttt{V\_dot\_Grad(Vector, Vector)}} \\
   f =& \nabla^2 f &\qquad {\texttt{Field}} =& {\texttt{Laplace(Field)}}
   \end{array}

.. math::

   \nabla\phi =& {{\frac{\partial \phi}{\partial u^i}}}\nabla u^i arrow (\nabla\phi)_i =
       {{\frac{\partial \phi}{\partial u^i}}} \\ \nabla\cdot A =& =
       \frac{1}{J}{{\frac{\partial }{\partial u^i}}}(Jg^{ij}A_j) \\ \nabla^2\phi =&
       G^j{{\frac{\partial \phi}{\partial u^i}}} + g^{ij}\frac{\partial^2\phi}{\partial u^i\partial
       u^j}

where we have defined

.. math::

   G^j =& \frac{1}{J}{{\frac{\partial }{\partial u^i}}}(Jg^{ij})

**not** to be confused with the Christoffel symbol of the second kind
(see the coordinates manual for more details).

Clebsch operators
-----------------

Another set of operators assume that the equilibrium magnetic field is
written in Clebsch form as

.. math::

   \mathbf{B}_0 = \nabla z\times\nabla x \qquad B_0 = \frac{\sqrt{g_{yy}}}{J}

 where

.. math::

   \mathbf{B}_0 = |\mathbf{B}_0|\mathbf{b}_0 = B_0 \mathbf{b}_0

 is the background *equilibrium* magnetic field.

| l c Function & Formula
| ``Grad_par`` &
  :math:`\displaystyle\partial^0_{||} = \mathbf{b}_0\cdot\nabla =
  \frac{1}{\sqrt{g_{yy}}}{{\frac{\partial }{\partial y}}}`
| ``Div_par`` & :math:`\displaystyle \nabla^0_{||}f =
  B_0\partial^0_{||}(\frac{f}{B_0})`
| ``Grad2_par2`` & :math:`\displaystyle \partial^2_{||}\phi =
  \partial^0_{||}(\partial^0_{||}\phi) =
  \frac{1}{\sqrt{g_{yy}}}{{\frac{\partial }{\partial y}}}(\frac{1}{\sqrt{g_{yy}}}){{\frac{\partial 
  \phi}{\partial y}}} + \frac{1}{g_{yy}}\frac{\partial^2\phi}{\partial y^2}`
| ``Laplace_par`` & :math:`\displaystyle \nabla_{||}^2\phi =
  \nabla\cdot\mathbf{b}_0\mathbf{b}_0\cdot\nabla\phi =
  \frac{1}{J}{{\frac{\partial }{\partial y}}}(\frac{J}{g_{yy}}{{\frac{\partial \phi}{\partial y}}})`
| ``Laplace_perp`` &
  :math:`\displaystyle \nabla_\perp^2 = \nabla^2 - \nabla_{||}^2`
| ``Delp2`` & Perpendicular Laplacian, neglecting all :math:`y`
  derivatives
| & The ``Laplacian`` solver performs the inverse operation
| ``brackets`` & Poisson brackets
| & The Arakawa option, neglects the parallel :math:`y` derivatives if
  :math:`g_{xy}` and :math:`g_{yz}` are non-zero

We have that

.. math::

   \mathbf{b}_0\cdot\nabla\phi\times\nabla A =&
       \frac{1}{J\sqrt{g_{yy}}}[(g_{yy}{{\frac{\partial \phi}{\partial z}}} -
       g_{yz}{{\frac{\partial \phi}{\partial y}}}){{\frac{\partial A}{\partial x}}} + (g_{yz}{{\frac{\partial \phi}{\partial x}}} -
   g_{xy}{{\frac{\partial \phi}{\partial z}}}){{\frac{\partial A}{\partial y}}} + (g_{xy}{{\frac{\partial \phi}{\partial y}}} -
   g_{yy}{{\frac{\partial \phi}{\partial x}}}){{\frac{\partial A}{\partial z}}}]

.. math::

   \nabla_\perp \equiv \nabla - {{\boldsymbol{b}}}({{\boldsymbol{b}}}\cdot\nabla) \qquad
   {{\boldsymbol{b}}}\cdot\nabla = \frac{1}{JB}\frac{\partial}{\partial y}

.. math::

   {{\boldsymbol{b}}} = \frac{1}{JB}{{\boldsymbol{e}}}_y = \frac{1}{JB}[g_{xy}\nabla x + g_{yy}\nabla y
   + g_{yz}\nabla z]

In a Clebsch coordinate system
:math:`{{\boldsymbol{B}}} = \nabla z \times \nabla x = \frac{1}{J}{{\boldsymbol{e}}}_y`,
:math:`g_{yy} = {{\boldsymbol{e}}}_y\cdot{{\boldsymbol{e}}}_y = J^2B^2`,
and so the :math:`\nabla y` term cancels out:

.. math::

   \nabla_\perp =& \nabla x({{\frac{\partial }{\partial x}}} -
       \frac{g_{xy}}{(JB)^2}{{\frac{\partial }{\partial y}}}) + \nabla z({{\frac{\partial }{\partial z}}} -
       \frac{g_{yz}}{(JB)^2}{{\frac{\partial }{\partial y}}})

The bracket operators
---------------------

| The bracket operator ``brackets(phi, f, method)`` aims to
  differentiate equations on the form

  .. math::

         -\frac{\nabla\phi\times{{\boldsymbol{b}}}}{B}\cdot\nabla f

| Notice that when we use the Arakawa scheme, :math:`y`-derivatives are
  neglected if :math:`g_{xy}` and :math:`g_{yz}` are non-zero. An
  example of usage of the brackets can be found in for example
  ``examples/MMS/advection`` or ``examples/blob2d``.

Setting differencing method
---------------------------

Staggered grids
===============

Until now all quantities have been cell-centred i.e. both velocities and
conserved quantities were defined at the same locations. This is because
these methods are simple and this was the scheme used in the original
BOUT. This class of methods can however be susceptible to grid-grid
oscillations, and so most shock-capturing schemes involve densities and
velocities (for example) which are not defined at the same location:
their grids are staggered.

By default BOUT++ runs with all quantities at cell centre. To enable
staggered grids, set

::

    StaggerGrids = true

in the top section of the ``BOUT.inp`` file. The **test-staggered**
example illustrates how to use staggered grids in BOUT++.

There are four possible locations in a grid cell where a quantity can be
defined in BOUT++: centre, lower X, lower Y, and lower Z. These are
illustrated in figure [fig:stagLocations].

To specify the location of a variable, use the method ``setLocation()``
with one of the locations ``CELL\_CENTRE``, ``CELL\_XLOW``,
``CELL\_YLOW`` , or ``CELL\_ZLOW`` .

The key lines in the **test-staggered** example which specify the
locations of the evolving variables are

::

    Field3D n, v;

    int physics_init(bool restart) {
      v.setLocation(CELL_YLOW); // Staggered relative to n
      SOLVE_FOR2(n, v);
      ...

which makes the velocity ``v`` staggered to the lower side of the cell
in Y, whilst the density :math:`n` remains cell centred.

Arithmetic operations between staggered quantities are handled by
interpolating them to the same location according to the algorithm in
figure [fig:stagArith].

If performing an operation between variables defined at two different
locations, the order of the variables matter: the result will be defined
at the locations of the **left** variable. For example, ``n*v`` would be
``CELL_CENTRE`` because this is the location of ``n`` , whilst ``v*n``
would be ``CELL_YLOW`` . Relying on this behaviour could lead to
trouble, to make your code clearer it’s probably best to use the
interpolation routines. Include the header file

::

    #include <interpolation.hxx>

then use the ``interp_to(field, location)`` function. Using this,
``interp_to(n, CELL_YLOW)*v`` would be ``CELL_YLOW`` as ``n`` would be
interpolated.

Differential operators by default return fields which are defined at the
same location as their inputs, so here ``Grad_par(v)`` would be
``CELL_YLOW`` . If this is not what is wanted, give the location of the
result as an additional argument: ``Grad_par(v, CELL_CENTRE)`` uses
staggered differencing to produce a result which is defined at the cell
centres. As with the arithmetic operators, if you ask for the result to
be staggered in a different direction from the input then the
differencing will be to cell centre and then be interpolated. For
example ``Grad_par(v, CELL_XLOW)`` would first perform staggered
differencing from ``CELL_YLOW`` to get a result at ``CELL_CENTRE`` , and
then interpolate the result to ``CELL_XLOW`` .

Advection operators which take two arguments return a result which is
defined at the location of the field being advected. For example
``Vpar_Grad_par(v, f)`` calculates :math:`v \nabla_{||} f` and returns a
result at the same location as ``f``. If ``v`` and ``f`` are defined at
the same locations then centred differencing is used, if one is centred
and the other staggered then staggered differencing is used, and if both
are staggered to different locations then the behaviour is less well
defined (don’t do it). As with other differential operators, the
required location of the result can be given as an optional argument.

Advanced methods
================

This section describes the more advanced methods which can be used to
speed up simulations using implicit time stepping schemes. At the time
of writing (Dec ’12), they can be used with either the SUNDIALS CVODE or
PETSc solvers.

Global field gather / scatter
-----------------------------

In BOUT++ each processor performs calculations on a sub-set of the mesh,
and communicates with other processors primarily through exchange of
guard cells (the ``mesh->commmunicate`` function). If you need to gather
data from the entire mesh onto a single processor, then this can be done
using either 2D or 3D ``GlobalFields`` .

First include the header file

::

    #include <bout/globalfield.hxx>

which defines both ``GlobalField2D`` and ``GlobalField3D`` . To create a
3D global field, pass it the mesh pointer:

::

      GlobalField3D g3d(mesh);

By default all data will be gathered onto processor 0. To change this,
specify which processor the data should go to as the second input

::

      GlobalField3D g3d(mesh, processor);

Gather and scatter methods are defined:

::

      Field3D localData;
      // Set local data to some value

      g3d.gather(localData);  // Gathers all data onto one processor

      localData = g3d.scatter(); // Scatter data back

**Note:** Boundary guard cells are **not** handled by the scatter step,
as this would mean handling branch-cuts etc. To obtain valid data in the
guard and Y boundary cells, you will need to communicate and set Y
boundaries.

**Note:** Gather and Scatter are global operations, so all processors
must call these functions.

Once data has been gathered, it can be used on one processor. To check
if the data is available, call the method ``dataIsLocal()``, which will
return ``true`` only on one processor

::

      if(g3d.dataIsLocal()) {
        // Data is available on this processor

      }

The sizes of the global array are available through ``xSize()``,
``ySize()`` and ``zSize()`` methods. The data itself can be accessed
indirectly using ``(x,y,z)`` operators:

::

      for(int x=0; x<g3d.xSize(); x++)
        for(int y=0; y<g3d.ySize(); y++)
          for(int z=0; z<g3d.zSize(); z++)
            output.write("Value at (%d,%d,%d) is %e\n",
            x,y,z,
            g3d(x,y,z) );

or by getting a pointer to the underlying data, which is stored as a 1D
array:

::

      BoutReal *data = g3d.getData();
      nx = g3d.xSize();
      ny = g3d.ySize();
      nz = g3d.zSize();

      data[x*ny*nz + y*nz + z]; // Value at g3d(x,y,z)

See the example ``examples/test-globalfield`` for more examples.

.. _sec-LaplaceXY:

LaplaceXY
---------

Perpendicular Laplacian solver in X-Y.

.. math::
   :label: nabl_perp_f

   \nabla_\perp f =& \nabla f - \mathbf{b}\left(\mathbf{b}\cdot\nabla\right)
       \nonumber \\ =& \left(\frac{\partial f}{\partial x} -
   \frac{g_{xy}}{g_{yy}}\frac{\partial f}{\partial y}\right)\nabla x +
   \left(\frac{\partial f}{\partial z} - \frac{g_{yz}}{g_{yy}}\frac{\partial
   f}{\partial y}\right)\nabla z

In 2D (X-Y), the :math:`g_{xy}` component can be dropped since this
depends on integrated shear :math:`I` which will cancel with the
:math:`g_{xz}` component. The :math:`z` derivative is zero and so this
simplifies to

.. math::

   \nabla_\perp f = \frac{\partial f}{\partial x}\nabla x -
   \frac{g_{yz}}{g_{yy}}\frac{\partial f}{\partial y}\nabla z

The divergence operator in conservative form is

.. math::

   \nabla\cdot\mathbf{A} = \frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}A_j\right)

and so the perpendicular Laplacian in X-Y is

.. math::

   \nabla_\perp^2f = \frac{1}{J}\frac{\partial}{\partial
   x}\left(Jg^{xx}\frac{\partial f}{\partial x}\right) -
   \frac{1}{J}\frac{\partial}{\partial
   y}\left(Jg^{yz}\frac{g_{yz}}{g_{yy}}\frac{\partial f}{\partial y}\right)

In field-aligned coordinates, the metrics in the :math:`y` derivative
term become:

.. math::

   g^{yz}\frac{g_{yz}}{g_{yy}} = \frac{B_{tor}^2}{B^2}\frac{1}{h_\theta^2}

In the LaplaceXY operator this is implemented in terms of fluxes at
cell faces.

.. math::

   \frac{1}{J}\frac{\partial}{\partial x}\left(Jg^{xx}\frac{\partial f}{\partial
   x}\right) &\rightarrow&
           \frac{1}{J_i\mathrm{dx_i}}\left[J_{i+1/2}g^{xx}_{i+1/2}\left(\frac{f_{i+1}
               - f_{i}}{\mathrm{dx}_{i+1/2}}\right) -
               J_{i-1/2}g^{xx}_{i-1/2}\left(\frac{f_{i} -
           f_{i-1}}{\mathrm{dx}_{i-1/2}}\right)\right]

Notes:

-  The ShiftXderivs option must be true for this to work, since it
   assumes that :math:`g^{xz} = 0`

.. _sec-LaplaceXZ:

LaplaceXZ
---------

This is a Laplacian inversion code in X-Z, similar to the ``Laplacian``
solver described in :ref:`sec-laplacian`. The difference is in the
form of the Laplacian equation solved, and the approach used to derive
the finite difference formulae. The equation solved is:

.. math::

     \nabla\cdot\left( A \nabla_\perp f \right) + Bf = b

where :math:`A` and :math:`B` are coefficients, :math:`b` is the known
RHS vector (e.g. vorticity), and :math:`f` is the unknown quantity to be
calculated (e.g. potential), and :math:`\nabla_\perp f` is the same as
equation (:eq:`nabl_perp_f`), but with negligible :math:`y`-parallel
derivatives if :math:`g_{xy}`, :math:`g_{yz}` and :math:`g_{xz}` is
non-vanishing. The Laplacian is written in conservative form like the
``LaplaceXY`` solver, and discretised in terms of fluxes through cell
faces.

.. math::

     \frac{1}{J}\frac{\partial}{\partial x}\left(J A g^{xx}\frac{\partial
     f}{\partial x}\right) + \frac{1}{J}\frac{\partial}{\partial z}\left(J A
     g^{zz}\frac{\partial f}{\partial z}\right) + B f = b

The header file is ``include/bout/invert/laplacexz.hxx``. The solver is
constructed by using the ``LaplaceXZ::create`` function:

::

      LaplaceXZ *lap = LaplaceXZ::create(mesh);

Note that a pointer to a ``Mesh`` object must be given, which for now is
the global variable ``mesh`` . By default the options section
``laplacexz`` is used, so to set the type of solver created, set in the
options

.. code-block:: cfg

      [laplacexz]
      type = petsc  # Set LaplaceXZ type

or on the command-line ``laplacexz:type=petsc`` .

The coefficients must be set using ``setCoefs`` . All coefficients must
be set at the same time:

::

      lap->setCoefs(1.0, 0.0);

Constants, ``Field2D`` or ``Field3D`` values can be passed. If the
implementation doesn’t support ``Field3D`` values then the average over
:math:`z` will be used as a ``Field2D`` value.

To perform the inversion, call the ``solve`` function:

::

      Field3D vort = ...;

      Field3D phi = lap->solve(vort, 0.0);

The second input to ``solve`` is an initial guess for the solution,
which can be used by iterative schemes e.g. using PETSc.

Implementations
~~~~~~~~~~~~~~~

The currently available implementations are:

-  ``cyclic``: This implementation assumes coefficients are constant in
   :math:`Z`, and uses FFTs in :math:`z` and a complex tridiagonal
   solver in :math:`x` for each :math:`z` mode (the ``CyclicReduction``
   solver). Code in ``src/invert/laplacexz/impls/cyclic/``.

-  ``petsc``: This uses the PETSc KSP interface to solve a matrix with
   coefficients varying in both :math:`x` and :math:`z`. To improve
   efficiency of direct solves, a different matrix is used for
   preconditioning. When the coefficients are updated the preconditioner
   matrix is not usually updated. This means that LU factorisations of
   the preconditioner can be re-used. Since this factorisation is a
   large part of the cost of direct solves, this should greatly reduce
   the run-time.

Test case
~~~~~~~~~

The code in ``examples/test-laplacexz`` is a simple test case for
``LaplaceXZ`` . First it creates a ``LaplaceXZ`` object:

::

      LaplaceXZ *inv = LaplaceXZ::create(mesh);

For this test the ``petsc`` implementation is the default:

.. code-block:: cfg

      [laplacexz]
      type = petsc
      ksptype = gmres # Iterative method
      pctype  = lu  # Preconditioner

By default the LU preconditioner is used. PETSc’s built-in factorisation
only works in serial, so for parallel solves a different package is
needed. This is set using:

::

      factor_package = superlu_dist

This setting can be “petsc” for the built-in (serial) code, or one of
“superlu”, “superlu\_dist”, “mumps”, or “cusparse”.

Then we set the coefficients:

::

      inv->setCoefs(Field3D(1.0),Field3D(0.0));

Note that the scalars need to be cast to fields (Field2D or Field3D)
otherwise the call is ambiguous. Using the PETSc command-line flag
``-mat_view ::ascii_info`` information on the assembled matrix is
printed:

.. code-block:: bash

      $ mpirun -np 2 ./test-laplacexz -mat_view ::ascii_info
      ...
      Matrix Object: 2 MPI processes
      type: mpiaij
      rows=1088, cols=1088
      total: nonzeros=5248, allocated nonzeros=5248
      total number of mallocs used during MatSetValues calls =0
        not using I-node (on process 0) routines
      ...

which confirms that the matrix element pre-allocation is setting the
correct number of non-zero elements, since no additional memory
allocation was needed.

A field to invert is created using FieldFactory:

::

      Field3D rhs = FieldFactory::get()->create3D("rhs",
                                                  Options::getRoot(),
                                                  mesh);

which is currently set to a simple function in the options:

::

      rhs = sin(x - z)

and then the system is solved:

::

      Field3D x = inv->solve(rhs, 0.0);

Using the PETSc command-line flags ``-ksp_monitor`` to monitor the
iterative solve, and ``-mat_superlu_dist_statprint`` to monitor
SuperLU\_dist we get:

.. code-block:: bash

            Nonzeros in L       19984
            Nonzeros in U       19984
            nonzeros in L+U     38880
            nonzeros in LSUB    11900
            NUMfact space (MB) sum(procs):  L\U     0.45    all     0.61
            Total highmark (MB):  All       0.62    Avg     0.31    Max     0.36
            Mat conversion(PETSc->SuperLU_DIST) time (max/min/avg):
                                  4.69685e-05 / 4.69685e-05 / 4.69685e-05
            EQUIL time             0.00
            ROWPERM time           0.00
            COLPERM time           0.00
            SYMBFACT time          0.00
            DISTRIBUTE time        0.00
            FACTOR time            0.00
            Factor flops    1.073774e+06    Mflops    222.08
            SOLVE time             0.00
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     28.67
      0 KSP Residual norm 5.169560044060e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     60.50
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     49.86
      1 KSP Residual norm 1.359142853145e-12

So after the initial setup and factorisation, the system is solved in
one iteration using the LU direct solve.

As a test of re-using the preconditioner, the coefficients are then
modified:

::

      inv->setCoefs(Field3D(2.0),Field3D(0.1));

and solved again:

::

            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     84.15
      0 KSP Residual norm 5.169560044060e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     90.42
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     98.51
      1 KSP Residual norm 2.813291076609e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     94.88
      2 KSP Residual norm 1.688683980433e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     87.27
      3 KSP Residual norm 7.436784980024e+01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     88.77
      4 KSP Residual norm 1.835640800835e+01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     89.55
      5 KSP Residual norm 2.431147365563e+00
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     88.00
      6 KSP Residual norm 5.386963293959e-01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     93.50
      7 KSP Residual norm 2.093714782067e-01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     91.91
      8 KSP Residual norm 1.306701698197e-02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     89.44
      9 KSP Residual norm 5.838501185134e-04
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     81.47

Note that this time there is no factorisation step, but the direct solve
is still very effective.

Blob2d comparison
~~~~~~~~~~~~~~~~~

The example ``examples/blob2d-laplacexz`` is the same as
``examples/blob2d`` but with ``LaplaceXZ`` rather than ``Laplacian``.

Tests on one processor: Using Boussinesq approximation, so that the
matrix elements are not changed, the cyclic solver produces output

::

    1.000e+02        125       8.28e-01    71.8    8.2    0.4    0.6   18.9
    2.000e+02         44       3.00e-01    69.4    8.1    0.4    2.1   20.0

whilst the PETSc solver with LU preconditioner outputs

::

    1.000e+02        146       1.15e+00    61.9   20.5    0.5    0.9   16.2
    2.000e+02         42       3.30e-01    58.2   20.2    0.4    3.7   17.5

so the PETSc direct solver seems to take only slightly longer than the
cyclic solver. For comparison, GMRES with Jacobi preconditioning gives:

::

    1.000e+02        130       2.66e+00    24.1   68.3    0.2    0.8    6.6
    2.000e+02         78       1.16e+00    33.8   54.9    0.3    1.1    9.9

and with SOR preconditioner

::

    1.000e+02        124       1.54e+00    38.6   50.2    0.3    0.4   10.5
    2.000e+02         45       4.51e-01    46.8   37.8    0.3    1.7   13.4

When the Boussinesq approximation is not used, the PETSc solver with LU
preconditioning, re-setting the preconditioner every 100 solves gives:

::

    1.000e+02        142       3.06e+00    23.0   70.7    0.2    0.2    6.0
    2.000e+02         41       9.47e-01    21.0   72.1    0.3    0.6    6.1

i.e. around three times slower than the Boussinesq case. When using
jacobi preconditioner:

::

    1.000e+02        128       2.59e+00    22.9   70.8    0.2    0.2    5.9
    2.000e+02         68       1.18e+00    26.5   64.6    0.2    0.6    8.1

For comparison, the ``Laplacian`` solver using the tridiagonal solver as
preconditioner gives:

::

    1.000e+02        222       5.70e+00    17.4   77.9    0.1    0.1    4.5
    2.000e+02        172       3.84e+00    20.2   74.2    0.2    0.2    5.2

or with Jacobi preconditioner:

::

    1.000e+02        107       3.13e+00    15.8   79.5    0.1    0.2    4.3
    2.000e+02        110       2.14e+00    23.5   69.2    0.2    0.3    6.7

The ``LaplaceXZ`` solver does not appear to be dramatically faster **in
serial** than the ``Laplacian`` solver when the matrix coefficients are
modified every solve. When matrix elements are not modified then the
solve time is competitive with the tridiagonal solver.

As a test, timing only the ``setCoefs`` call for the non-Boussinesq case
gives

::

    1.000e+02        142       1.86e+00    83.3    9.5    0.2    0.3    6.7
    2.000e+02         41       5.04e-01    83.1    8.0    0.3    1.2    7.3

so around 9% of the run-time is in setting the coefficients, and the
remaining :math:`\sim 60`\ % in the solve itself.

Eigenvalue solver
=================

By using the SLEPc library, BOUT++ can be used as an eigenvalue solver
to find the eigenvectors and eigenvalues of sets of equations.

Configuring with SLEPc
----------------------

The BOUT++ interface has been tested with SLEPc version 3.4.3, itself
compiled with PETSc 3.4.2. SLEPc version 3.4 should work, but other
versions will not yet.

SLEPc options
-------------

Time derivatives can be taken directly from the RHS function, or by
advancing the simulation in time by a relatively large increment. This
second method acts to damp high frequency components

Examples
--------

Wave in a box
~~~~~~~~~~~~~

``examples/eigen-box``

Testing
=======

Two types of tests are currently used in BOUT++ to catch bugs as early
as possible: Unit tests, which check a small piece of the code
separately, and a test suite which runs the entire code on a short
problem. Unit tests can be run using the ``src/unit_tests`` Python
script. This searches through the directories looking for an executable
script called ``unit_test``, runs them, and collates the results. Not
many tests are currently available as much of the code is too tightly
coupled. If done correctly, the unit tests should describe and check the
behaviour of each part of the code, and hopefully the number of these
will increase over time. The test suite is in the ``examples``
directory, and is run using the ``test_suite`` python script. At the top
of this file is a list of the subdirectories to run (e.g. ``test-io``,
``test-laplace``, and ``interchange-instability``). In each of those
subdirectories the script ``runtest`` is executed, and the return value
used to determine if the test passed or failed.

All tests should be short, otherwise it discourages people from running
the tests before committing changes. A few minutes or less on a typical
desktop, and ideally only a few seconds. If you have a large simulation
which you want to stop anyone breaking, find starting parameters which
are as sensitive as possible so that the simulation can be run quickly.

Method of Manufactured Solutions
--------------------------------

The Method of Manufactured solutions (MMS) is a rigorous way to check
that a numerical algorithm is implemented correctly. A known solution is
specified (manufactured), and it is possible to check that the code
output converges to this solution at the expected rate.

To enable testing by MMS, switch an input option “mms” to true:

.. code-block:: cfg

    [solver]
    mms = true

This will have the following effect:

#. For each evolving variable, the solution will be used to initialise
   and to calculate the error

Choosing manufactured solutions
-------------------------------

Manufactured solutions must be continuous and have continuous
derivatives. Common mistakes:

-  Don’t use terms multiplying coordinates together e.g. ``x * z`` or
   ``y * z``. These are not periodic in :math:`y` and/or :math:`z`, so
   will give strange answers and usually no convergence. Instead use
   ``x * sin(z)`` or similar, which are periodic.

Timing
------

To time parts of the code, and calculate the percentage of time spent in
communications, file I/O, etc. there is the ``Timer`` class defined in
``include/bout/sys/timer.hxx``. To use it, just create a ``Timer``
object at the beginning of the function you want to time:

::

    #include <bout/sys/timer.hxx>

    void someFunction() {
      Timer timer("test")
      ...
    }

Creating the object starts the timer, and since the object is destroyed
when the function returns (since it goes out of scope) the destructor
stops the timer.

::

    class Timer {
    public:
      Timer();
      Timer(const std::string &label);
      ~Timer();

      double getTime();
      double resetTime();
    };

The empty constructor is equivalent to setting ``label = ""`` .
Constructors call a private function ``getInfo()`` , which looks up the
``timer_info`` structure corresponding to the label in a
``map<string, timer_info*>`` . If no such structure exists, then one is
created. This structure is defined as:

::

    struct timer_info {
      double time;    ///< Total time
      bool running;   ///< Is the timer currently running?
      double started; ///< Start time
    };

Since each timer can only have one entry in the map, creating two timers
with the same label at the same time will lead to trouble. Hence this
code is **not** thread-safe.

The member functions ``getTime()`` and ``resetTime()`` both return the
current time. Whereas ``getTime()`` only returns the time without
modifying the timer, ``resetTime()`` also resets the timer to zero.

If you don’t have the object, you can still get and reset the time using
static methods:

::

    double Timer::getTime(const std::string &label);
    double Timer::resetTime(const std::string &label);

These look up the ``timer_info`` structure, and perform the same task as
their non-static namesakes. These functions are used by the monitor
function in ``bout++.cxx`` to print the percentage timing information.

Examples
========

The code and input files in the ``examples/`` subdirectory are for
research, demonstrating BOUT++, and to check for broken functionality.
Some proper unit tests have been implemented, but this is something
which needs improving. The examples which were published in
:raw-latex:`\cite{Dudson2009,dudson-2008-arxiv}` were
``drift-instability``, ``interchange-instability`` and ``orszag-tang``.

advect1d
--------

The model in ``gas_compress.cxx`` solves the compressible gas dynamics
equations for the density :math:`n`, velocity :math:`\mathbf{V}`, and
pressure :math:`P`:

drift-instability
-----------------

The physics code ``2fluid.cxx`` implements a set of reduced Braginskii
2-fluid equations, similar to those solved by the original BOUT code.
This evolves 6 variables: Density, electron and ion temperatures,
parallel ion velocity, parallel current density and vorticity.

Input grid files are the same as the original BOUT code, but the output
format is different.

em-drift
--------

gyro-gem
--------

interchange-instability
-----------------------

.. figure:: figs/interchange_inst_test.*
   :alt: Interchange instability test

   Interchange instability test. Solid lines are from analytic theory,
   symbols from BOUT++ simulations, and the RMS density is averaged over
   :math:`z`. Vertical dashed line marks the reference point, where
   analytic and simulation results are set equal

jorek-compare
-------------

lapd-drift
----------

orszag-tang
-----------

The file ``mhd.cxx`` solves the full MHD equations for the full values
(perturbation + initial), whilst the file ``mhd_perturb.cxx`` solves for
a perturbation about the equilibrium.

shear-alfven-wave
-----------------

sod-shock
---------

.. figure:: figs/sod_result.*
   :alt: Sod shock-tube problem for testing shock-handling methods
   :width: 48.0%

   Sod shock-tube problem for testing shock-handling methods

uedge-benchmark
---------------

Notes
=====

Compile options
---------------

Compiling with ``-DCHECK`` enables a lot of checks of operations
performed by the field objects. This is very useful for debugging a
code, and can be omitted once bugs have been removed.

For (sometimes) more useful error messages, there is the ``-DTRACK``
option. This keeps track of the names of variables and includes these in
error messages.

Adaptive grids
--------------

Two types of adaptive grids can be used in BOUT++: Moving meshes, and
changing resolution.

Moving meshes
~~~~~~~~~~~~~

During either the initialisation, or the simulation itself, the metric
tensors can be modified. This could be used to make the coordinate
system time-dependent. Since currently the metric tensors are 2D fields,
this would only allow axisymmetric motion. Changing the tensors to be 3D
objects is however possible with fairly small modification to the code.

Whenever one of the metrics :math:`g^{ij}` are changed, a call to
``geometry()`` must be made.

Changing resolution
~~~~~~~~~~~~~~~~~~~

Since all 2D and 3D fields/vectors are located internally in global
lists, the resolution of the grid can be changed when required by
interpolation. **This requires a new, more efficient implementation of
the Fields classes**.

Machine-specific installation
=============================

Archer
------

As of 30th April 2014, the following configuration should work

.. code-block:: bash

    $ module swap PrgEnv-cray PrgEnv-gnu/5.1.29
    $ module load fftw
    $ module load netcdf/4.1.3

Compiling and running under AIX
===============================

Most development and running of BOUT++ is done under Linux, with the
occasional FreeBSD and OSX. The configuration scripts are therefore
heavily tested on these architectures. IBM’s POWER architecture however
runs AIX, which has some crucial differences which make compiling a
pain.

-  Under Linux/BSD, it’s usual for a Fortran routine ``foo`` to appear
   under C as ``foo_``, whilst under AIX the name is unchanged

-  MPI compiler scripts are usually given the names ``mpicc`` and either
   ``mpiCC`` or ``mpicxx``. AIX uses ``mpcc`` and ``mpCC``.

-  Like BSD, the ``make`` command isn’t compatible with GNU make, so you
   have to run ``gmake`` to compile everything.

-  The POWER architecture is big-endian, different to the little endian
   Intel and AMD chips. This can cause problems with binary file
   formats.

.. _sec-sundials:

SUNDIALS
--------

To compile SUNDIALS, use

.. code-block:: bash

    $ export CC=cc
    $ export CXX=xlC
    $ export F77=xlf
    $ export OBJECT_MODE=64
    $ ./configure --prefix=$HOME/local/ --with-mpicc=mpcc --with-mpif77=mpxlf CFLAGS=-maix64

You may get an error message like:

.. code-block:: bash

    make: Not a recognized flag: w

This is because the AIX ``make`` is being used, rather than ``gmake``.
The easiest way to fix this is to make a link to ``gmake`` in your local
bin directory:

.. code-block:: bash

    $ ln -s /usr/bin/gmake $HOME/local/bin/make

Running ``which make`` should now point to this ``local/bin/make``, and
if not then you need to make sure that your bin directory appears first
in the ``PATH``:

.. code-block:: bash

    $ export PATH=$HOME/local/bin:$PATH

If you see an error like this:

.. code-block:: bash

    ar: 0707-126 ../../src/sundials/sundials_math.o is not valid with the current object file mode.
            Use the -X option to specify the desired object mode.

then you need to set the environment variable ``OBJECT_MODE``

.. code-block:: bash

    $ export OBJECT_MODE=64

Configuring BOUT++, you may get the error:

.. code-block:: bash

    configure: error: C compiler cannot create executables

In that case, you can try using:

.. code-block:: bash

    $ ./configure CFLAGS="-maix64"

When compiling, you may see warnings

.. code-block:: bash

    xlC_r: 1501-216 (W) command option -64 is not recognized - passed to ld

At this point, the main BOUT++ library should compile, and you can try
compiling one of the examples.

.. code-block:: bash

    ld: 0711-317 ERROR: Undefined symbol: .NcError::NcError(NcError::Behavior)
    ld: 0711-317 ERROR: Undefined symbol: .NcFile::is_valid() const
    ld: 0711-317 ERROR: Undefined symbol: .NcError::~NcError()
    ld: 0711-317 ERROR: Undefined symbol: .NcFile::get_dim(const char*) const

This is probably because the NetCDF libraries are 32-bit, whilst BOUT++
has been compiled as 64-bit. You can try compiling BOUT++ as 32-bit:

.. code-block:: bash

    $ export OBJECT_MODE=32
    $ ./configure CFLAGS="-maix32"
    $ gmake

If you still get undefined symbols, then go back to 64-bit, and edit
make.config, replacing ``-lnetcdf_c++`` with -lnetcdf64\_c++, and
``-lnetcdf`` with -lnetcdf64. This can be done by running:

.. code-block:: bash

    $ sed 's/netcdf/netcdf64/g' make.config > make.config.new
    $ mv make.config.new make.config

BOUT++ functions (alphabetical)
===============================

This is a list of functions which can be called by users writing a
physics module. For a full list of functions, see the Reference manual,
DOxygen documentation, and source code.

-  ``Field = abs(Field | Vector)``

-  | ``(Communicator).add(Field | Vector)``
   | Add a variable to a communicator object.

-  ``apply_boundary(Field. name)``

-  ``Field = b0xGrad_dot_Grad(Field, Field, CELL_LOC)``

-  ``bout_solve(Field, Field, name)``

-  ``bout_solve(Vector, Vector, name)``

-  | ``(Communicator).clear()``
   | Remove all variables from a Communicator object

-  ``Field = cos(Field)``

-  ``Field = cosh(Field)``

-  ``Vector = Curl(Vector)``

-  | ``Field = Delp2(Field)``
   | :math:`\nabla_\perp^2` operator

-  | ``Field = Div(Vector)``
   | Divergence of a vector

-  | ``Field = Div_par(Field f)``
   | Parallel divergence :math:`B_0\mathbf{b}\cdot\nabla(f / B_0)`

-  ``dump.add(Field, name, 1/0)``

-  ``Field = filter(Field, modenr)``

-  | ``geometry_derivs()``
   | Calculates useful quantities from the metric tensor. Call this
     every time the metric tensor is changed.

-  ``Vector = Grad(Field)``

-  ``Field = Grad_par(Field)``

-  ``Field = Grad2_par2(Field)``

-  | ``grid_load(BoutReal, name)``
   | Load a scalar real from the grid file

-  | ``grid_load2d(Field2D, name)``
   | Load a 2D scalar field from the grid file

-  | ``grid_load3d(Field3D, name)``
   | Load a 3D scalar field from the grid file

-  ``invert_laplace(Field input, Field output, flags, Field2D *A)``

-  | ``Field = invert_parderiv(Field2D|BoutReal A, Field2D|BoutReal B, Field3D r)``
   | Inverts an equation ``A*x + B*Grad2_par2(x) = r``

-  ``Field = Laplacian(Field)``

-  ``Field3D = low_pass(Field3D, max_modenr)``

-  ``BoutReal = max(Field)``

-  ``BoutReal = min(Field)``

-  | ``msg_stack.pop( |int)``
   | Remove a message from the top of the stack. If a message ID is
     passed, removes all messages back to that point.

-  | ``int = msg_stack.push(format, ...)``
   | Put a message onto the stack. Works like ``printf`` (and
     ``output.write``).

-  | ``options.get(name, variable, default)``
   | Get an integer, real or boolean value from the options file. If not
     in the file, the default value is used. The value used is printed
     to log file.

-  ``options.setSection(name)`` Set the section name in the input file

-  | ``output < < values``
   | Behaves like cout for stream output

-  | ``output.write(format, ...)``
   | Behaves like printf for formatted output

-  | ``(Communicator).receive()``
   | Receive data from other processors. Must be preceded by a ``send``
     call.

-  | ``(Communicator).run()``
   | Sends and receives data.

-  | ``(Communicator).send()``
   | Sends data to other processors (and posts receives). This must be
     followed by a call to ``receive()`` before calling send again, or
     adding new variables.

-  ``(Field3D).setLocation(CELL_LOC)``

-  ``(Field3D).ShiftZ(bool)``

-  ``Field = sin(Field)``

-  ``Field = sinh(Field)``

-  | ``solver.setPrecon(PhysicsPrecon)``
   | Set a preconditioner function

-  ``Field = sqrt(Field)``

-  ``Field = tan(Field)``

-  ``Field = tanh(Field)``

-  | ``Field = V_dot_Grad(Vector v, Field f)``
   | Calculates an advection term :math:`\mathbf{v}\cdot\nabla f`

-  | ``Vector = V_dot_Grad(Vector v, Vector u)``
   | Advection term :math:`\mathbf{v}\cdot\nabla\mathbf{u}`

-  ``Field = Vpar_Grad_par(Field v, Field f)``

-  | ``Field3D = where(Field2D test, Field|BoutReal gt0, Field|BoutReal lt0)``
   | Chooses between two values, depending on sign of ``test``.

IDL routines
============

List of IDL routines available in idllib. There are broadly three
categories of routine:

-  Completely general routines which could be useful outside BOUT++ work

   -  Data plotting and animation: **contour2** and **showdata**

   -  File reading and writing: **file\_open**, **file\_read** etc.

   -  User input and output: **get\_float**, **get\_integer**,
      **get\_yesno** and **str**

   -  FFT routines for integrating, differentiating and filtering:
      **fft\_integrate**, **fft\_deriv**, **fft\_filter**

-  Routines for BOUT++, but not specific to any application

   -  Modifying restart files: **expand\_restarts**, **scale\_restarts**
      and **split\_restarts**

   -  Processing 3D variables for input grid: **bout3dvar**

-  Routines specifically for tokamak simulations

   -  Reading A- and G-EQDSK format files into IDL: **read\_aeqdsk** and
      **read\_neqdsk**

   -  Plotting results: **polslice**, **plotpolslice**

Here the format is

**name**, arguments, [optional arguments]

-  | var = **bout3dvar** ( var )
   | Converts 3D variables to and from BOUT++’s Fourier representation
     which is used for input grids. By default converts from [x,y,z] to
     [x,y,f]

   -  **/reverse** Convert from [x,y,f] to [x,y,z]

   -  **nf**\ =nf Set number of frequencies in the result

   -  **nz**\ =nz When using /reverse, set number of Z points in the
      result

-  | var = **collect**\ ()
   | Read in data from a set of BOUT++ dump files

   -  **var** = “name of variable”

   -  **path** = “path/to/variable/”

   -  **xind**, **yind**, **zind**, **tind** = [min, max] index pairs

   -  **t\_array** = Output 1D array of times

-  | **contour2**, data [, x, y]
   | This is a replacement for the IDL contour which includes a scale
     color bar.

   -  **data** can be either 2D (x,y) or 3D (x,y,t). If data is 3D then
      the color is scaled to the entire range.

   -  **x** is an optional 2D (x,y) array of X coordinates

   -  **y** is an optional 2D (x,y) array of Y coordinates

   -  **t**\ =t is a time index for 3D data

   -  **nlev**\ =nlev

   -  **centre**\ =centre Make zero the middle of the color range (white
      if redblue)

   -  **redblue**\ =redblue Use a blue-white-red color scheme

   -  **revcolor**\ =revcolor Reverse color scheme

-  | **expand\_restarts**, newz
   | Increases the number of Z points in restart files. Together with
     scale\_restarts and split\_restarts, this makes it easier to modify
     a linear simulation as a start for non-linear runs.

   -  **newz** is the new value of NZ

   -  **path**\ =path Input path

   -  **output**\ =output Output path

   -  **format**\ =format File extension of output

-  | result = **fft\_deriv** ( var1d )
   | Calculates the derivative of a variable on a periodic domain.

-  result = **fft\_filter** (var, nf) Fourier filter a variable on a
   periodic domain. Arguments are a 1D variable and the number of
   Fourier components to keep

-  result = **fft\_integrate** ( var1d ) Integrates a variable on a
   periodic domain.

   -  **loop**\ =loop The loop integral is returned in this variable

-  | **file\_close**, handle
   | Close a file opened using file\_open()

-  | list = **file\_list** ( handle )
   | Return a list of variable names in the file

-  | integer = **file\_ndims** ( handle , “variable” )
   | Get the number of dimensions of a variable

-  | handle = **file\_open** ( “file” )
   | Open a NetCDF file.

   -  **/write** Open file for writing (default is read only)

   -  **/create** Create a new file, over-writing if already exists

-  var = **file\_read** ( handle, “variable” )

   -  **inds** = [xmin, xmax, ymin, ymax, ... ]

-  | float = **get\_float** ( “prompt” )
   | Ask the user for a float, using the given prompt

-  | integer = **get\_integer** ( “prompt” )
   | Ask the user for an integer

-  | integer = **get\_yesno** ( “prompt” )
   | Ask for a yes (1) or no (0) answer

-  | result = **gmres** ( x0, operator, b )
   | General Minimal Residual (GMRES)

   -  **x0** is the starting guess at the solution

   -  **operator**

   -  **b**

   Optional arguments

   -  **restart**\ =restart

   -  **max\_iter**\ =max\_iter

   -  **tol**\ =tol

   -  **stats**\ =stats

   -  **show**\ =show

   -  **output**\ =output

-  | result = **int\_func** ( [x,] f )
   | Integrate a function, always using the maximum number of
     grid-points possible for highest accuracy

-  | bool = **is\_pow2** ( value )
   | Returns 1 (true) if the given number is a power of 2, 0 (false)
     otherwise

-  | **plotpolslice**, var3d, grid
   | Takes a slice through a field-aligned tokamak domain, showing a
     poloidal cross-section.

   -  **var3d** is a 3D (x,y,z) variable to plot. Needs all of the
      points to work properly.

   -  **grid** is a structure from importing a grid file

   Optional arguments:

   -  **period**\ =period

   -  **zangle**\ =zangle

   -  **nlev**\ =nlev

   -  **yr**\ =yr

   -  **profile**\ =profile

   -  **output**\ =output

   -  **lines**\ =lines

   -  **linecol**\ =linecol

   -  **filter**\ =filter

-  | **polslice**, data, gridfile
   | Plots a 2D poloidal contour for single or double-null
     configurations, including color bar.

   -  **xstart**\ =xstart X index where the data begins. Useful if only
      part of the domain has been collected

   -  **ystart**\ =ystart Y index where data begins

-  | struct = **read\_aeqdsk** ( “filename” )
   | Reads an A-EQDSK file. Format is specified here:
     https://fusion.gat.com/THEORY/efit/a_eqdsk.html

-  | struct = **read\_neqdsk** ( “filename” )
   | Reads in an ’neqdsk’ or G-EQDSK formatted tokamak equilibrium file.
     Format of G-EQDSK file is specified here:
     https://fusion.gat.com/THEORY/efit/g_eqdsk.html

-  | stringarray = **regex\_extract** ( line, pattern )
   | Extract all matches to Regular Expression pattern contained in
     line. Useful for extracting numbers from FORTRAN-formatted text
     files.

   -  **line** Input string

   -  **pattern** Regular expression pattern to match

   -  **nmatch**\ =nmatch

-  | var = **reverse\_inds** ( var )
   | Reverse array indices e.g. ``arr[t,z,y,x] -> arr[x,y,z,t]``. Works
     on up to 5 dimensional variables

-  | **safe\_colors**
   | Sets the color table to useful values for plotting.

   -  **/first** Sets the first 10 colors to specific values, otherwise
      sets last 7

-  **scale\_restarts**, factor

   -  **path**\ =path Path to the restart files (default is current
      directory ’.’)

   -  **format**\ =format Specify what the file format is, otherwise
      goes on the file name

-  | **showdata**, data
   | Display animations of 1D,2D and 3D data. Defaults:

   -  2D data Animate a line plot

   -  3D data Animate a surface plot

   -  4D data Animate a poloidal cross-section (tokamaks only)

   Optional arguments:

   -  **/addsym** For 2D data (1D plots), add symbols to mark data
      points

   -  **az**\ =angle Rotate surface plots

   -  **/bw** Make contour plots grey scale

   -  **chars**\ =size character size

   -  **/contour** For 3D input, show color contour plot

   -  **delay**\ =time Time delay between plots (default 0.2 seconds)

   -  **/noscale** By default, all plots are on the same scale. This
      changes the scale for each plot’s range

   -  **profile**\ =array Background profile. Data is 3D: profile is 1D
      (X). Data is 4D -> profile is 2D (X,Y)

   -  **yr**\ =[min,max] Y range

-  | result = **sign** ( var )
   | This returns +1 if the variable is :math:`> 0`, -1 otherwise

-  **spectrum**

-  | **split\_restarts**, [nxpe], nype
   | split restart files between a different number of processors

   -  **nxpe** is an optional argument giving the number of processors
      in the X direction

   -  **nype** is the number of processors in the Y direction

   -  **path**\ =path Input path

   -  **output**\ =output Output path

   -  **format**\ =format File extension of output

-  | string = **str** ( value )
   | Convert a value to a string with whitespace trimmed. Arrays are
     converted to a comma-separated list in brackets.

-  | result = **zfamp** ( var4d )
   | Given a 4D variable [x,y,z,t], returns the Fourier amplitudes in
     [x,y,f,t]

-  | var = **zshift** ( var, shift )
   | Shifts a variable in the Z direction, useful for mapping between
     field-aligned and orthogonal coordinates.

   -  **period**\ =period How many domains fit in :math:`2\pi`. Default
      is 1 (full torus)

Python routines (alphabetical)
==============================

boututils
---------

-  ``class Datafile`` provides a convenient way to read and write NetCDF
   or HDF5 files. There are many different NetCDF libraries available
   for Python, so this class tries to provide a consistent interface to
   many of them, as well as to h5py.

-  ``deriv()``

-  ``determineNumberOfCPUs()``

-  ``file_import()`` reads the contents of a NetCDF file into a
   dictionary

-  ``integrate()``

-  ``launch()``

-  ``linear_regression()``

boutdata
--------

-  ``collect()`` provides an interface to read BOUT++ data outputs,
   returning NumPy arrays of data. It deals with the processor layout,
   working out which file contains each part of the domain.

   .. code-block:: python

           from boutdata.collect import collect

           t = collect("t_array")  # Collect the time values
         

-  ``pol_slice()`` takes a 3 or 4-D data set for a toroidal equilibrium,
   and calculates a slice through it at fixed toroidal angle.

-  ``gen_surface()`` is a generator for iterating over flux surfaces

bout\_runners
-------------

``bout_runners`` contains classes which gives an alternative way of
running BOUT++ simulations either normally using the class
``basic_runner`` , or on a cluster through a generated Portable Batch
System (PBS) script using the child class ``PBS_runner`` . Examples can
be found in ``examples/bout_runners_example/``.

``bout_runners`` is especially useful if one needs to make several runs
with only small changes in the options (which is normally written in
``BOUT.inp`` or in the command-line), as is the case when performing a
parameter scan, or when performing a MMS test.

Instead of making several runs with several different input files with
only small changes in the option, one can with ``bout_runners`` specify
the changes as member data of an instance of the appropriate
``bout_runners`` class. One way to do this is to write a *driver* in the
same directory as the executable. The *driver* is just a python script
which imports ``bout_runners`` , creates an instance, specifies the
running option as member data of that instance and finally calls the
member function ``self.execute_runs()`` .

In addition, the ``bout_runners`` provides a way to run any python
post-processing script after finished simulations (as long as it accept
at least one parameter containing the folder name(s) of the run(s)). If
the simulations have been performed using the ``PBS_runner`` , the
post-processing function will be submitted to the cluster (although it
is possible to submit it to a different queue, using a different amount
of nodes etc.).

When the function ``self.execute_runs()`` is executed, a folder
structure like the one presented in figure [fig:folder\_tree] is
created. ``BOUT.inp`` is copied to the folder of execution, where the
``BOUT.*.dmp`` files are stored. Secondly a list of combination of the
options specified in the driver is made. Eventually unset options are
obtained from ``BOUT.inp`` or given a default value if the option is
nowhere to be found.

.. figure:: figs/folder_tree.*
   :alt: Longest possible folder tree

   Longest possible folder tree made by the ``self.execute_runs()``
   function.

Derivatives of the Fourier transform
====================================

By using the definition of the Fourier transformed, we have

.. math::

   F(x,y,\xi) = {\int_{-\infty}^{\infty} {f(x,y,z)\exp(-2\pi iz\xi)} \; \text{d} {z}}

this gives

.. math::
   :label: f_derivative

   &{\int_{-\infty}^{\infty} {(\partial_zf[x,y,z])\exp(-2\pi iz\xi)} \; \text{d} {z}}\\
   =& {\int_{-\infty}^{\infty} {\partial_z(f[x,y,z]\exp[-2\pi iz\xi])} \; \text{d} {z}}
   - {\int_{-\infty}^{\infty} {f(x,y,z)\partial_z\exp(-2\pi iz\xi)} \; \text{d} {z}}\\
   =& (f[x,y,z]\exp[-2\pi iz\xi])\bigg|_{-\infty}^{\infty} - (-2\pi
   i\xi){\int_{-\infty}^{\infty} {f(x,y,z)\exp(-2\pi iz\xi)} \; \text{d} {z}}\\
   =& 2\pi i\xi F(x,y,\xi)

where we have used that :math:`f(x,y,\pm\infty)=0` in order to have a
well defined Fourier transform. This means that

.. math::

   \partial_z^n F(x,y,\xi) = (2\pi i \xi)^n F(x,y,\xi)

In our case, we are dealing with periodic boundary conditions. Strictly
speaking, the Fourier transform does not exist in such cases, but it is
possible to define a Fourier transform in the limit which in the end
lead to the Fourier series  [6]_ By discretising the spatial domain, it
is no longer possible to represent the infinite amount of Fourier modes,
but only :math:`N+1` number of modes, where :math:`N` is the number of
points (this includes the modes with negative frequencies, and the
zeroth offset mode). For the discrete Fourier transform, we have

.. math::
   :label: DFT

   F(x,y)_{k} = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp(\frac{-2\pi i k Z}{N})

where :math:`k` is the mode number, :math:`N` is the number of points
in :math:`z`. If we call the sampling points of :math:`z` for
:math:`z_Z`, where :math:`Z = 0, 1 \ldots N-1`, we have that
:math:`z_Z = Z \text{d}z`. As our domain goes from :math:`[0, 2\pi[`,
we have that (since we have one less line segment than point)
:math:`\text{d}z (N-1) = L_z = 2\pi - \text{d}z`, which gives
:math:`\text{d}z = \frac{2\pi}{N}`.  Inserting this is equation
(:eq:`DFT`) yields

.. math::

   F(x,y)_{k} = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp( - i k
   Z\text{d}z) = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp( - i k z_Z)

The discrete version of equation (:eq:`f_derivative`) thus gives

.. math::

   \partial_z^n F(x,y)_k = (i k)^n F(x,y)_k

.. [3]
   Taken from a talk by L.Chacon available here
   https://bout2011.llnl.gov/pdf/talks/Chacon_bout2011.pdf

.. [4]
   See paper http://arxiv.org/abs/1209.2054 for an application to
   2-fluid equations

.. [5]
   This ``InvertPar`` class can handle cases with closed field-lines and
   twist-shift boundary conditions for tokamak simulations

.. [6]
   For more detail see Bracewell, R. N. - The Fourier Transform and Its
   Applications 3rd Edition chapter :math:`10`
