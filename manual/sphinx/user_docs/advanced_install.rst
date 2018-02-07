.. Use bash as the default language for syntax highlighting in this file
.. highlight:: console

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

.. _NetCDF-4: https://www.unidata.ucar.edu/software/netcdf/
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

OpenMP
------

BOUT++ can make use of Single-Instruction Multiple-Data (SIMD)
parallelism through OpenMP. To enable OpenMP, use the
``--enable-openmp`` flag to configure::

    ./configure --enable-openmp

.. note::
    If you want to use OpenMP with Clang, you will need Clang 3.7+,
    and either ``libomp`` or ``libiomp``.

    You will be able to compile BOUT++ with OpenMP with lower versions
    of Clang, or using the GNU OpenMP library ``libgomp``, but it will
    only run with a single thread.

    
.. note::
    By default PVODE is built without OpenMP support. To enable this
    add ``--enable-pvode-openmp`` to the configure command.

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

.. note:: SUNDIALS is only downloadable from the home page, as submitting your
   name and e-mail is required for the download. As for the date of this
   typing, SUNDIALS version :math:`3.0.0` is the newest. In order for a
   smooth install it is recommended to install SUNDIALS from an install
   directory. The full installation guide is found in the downloaded
   ``.tar.gz``, but we will provide a step-by-step guide to install it
   and make it compatible with BOUT++ here

::

     $ cd ~
     $ mkdir -p local/examples
     $ mkdir -p install/sundials-install
     $ cd install/sundials-install
     $ # Move the downloaded sundials-3.0.0.tar.gz to sundials-install
     $ tar -xzvf sundials-3.0.0.tar.gz
     $ mkdir build
     $ cd build

     $ cmake \
       -DCMAKE_INSTALL_PREFIX=$HOME/local \
       -DEXAMPLES_INSTALL_PATH=$HOME/local/examples \
       -DCMAKE_LINKER=$HOME/local/lib \
       -DLAPACK_ENABLE=ON \
       -DOPENMP_ENABLE=ON \
       -DMPI_ENABLE=ON \
       ../sundials-3.0.0

     $ make
     $ make install

The SUNDIALS IDA solver is a Differential-Algebraic Equation (DAE)
solver, which evolves a system of the form
:math:`\mathbf{f}(\mathbf{u},\dot{\mathbf{u}},t) = 0`. This allows
algebraic constraints on variables to be specified.

To configure BOUT++ with SUNDIALS only (see section :ref:`sec-PETSc-install` on how
to build PETSc with SUNDIALS), go to the root directory of BOUT++ and
type::

    $ ./configure --with-sundials

SUNDIALS will allow you to select at run-time which solver to use. See
:ref:`sec-timeoptions` for more details on how to do this.

.. _sec-PETSc-install:

PETSc
-----

BOUT++ can use PETSc https://www.mcs.anl.gov/petsc/ for time-integration
and for solving elliptic problems, such as inverting Poisson and
Helmholtz equations.

Currently, BOUT++ supports PETSc versions 3.1, 3.2, 3.3 and 3.4
(support for newer versions are planned for the future). To install
PETSc version 3.4.5, use the following steps::

    $ cd ~
    $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.4.5.tar.gz
    $ tar -xzvf petsc-3.4.5.tar.gz
    $ # Optional
    $ # rm petsc-3.4.5.tar.gz
    $ cd petsc-3.4.5

To build PETSc without SUNDIALS, configure with::

    $ ./configure \
      --with-clanguage=cxx \
      --with-mpi=yes \
      --with-precision=double \
      --with-scalar-type=real \
      --with-shared-libraries=0

Add ``--with-debugging=yes`` to ``./configure`` in order to allow
debugging.

.. note:: To build PETSc with SUNDIALS, install SUNDIALS as explained
          in section :ref:`sec-sundials`, and append ``./configure``
          with ``--with-sundials-dir=$HOME/local``

.. note:: It is also possible to get PETSc to download and install
          MUMPS (see :ref:`sec-MUMPS`), by adding::

              --download-mumps \
              --download-scalapack \
              --download-blacs \
              --download-f-blas-lapack=1 \
              --download-parmetis \
              --download-ptscotch \
              --download-metis

          to ``./configure``.

To make PETSc, type::

    $ make PETSC_DIR=$HOME/petsc-3.4.5 PETSC_ARCH=arch-linux2-cxx-debug all

Should BLAS, LAPACK, or any other packages be missing, you will get an
error, and a suggestion that you can append
``--download-name-of-package`` to the ``./configure`` line. You may want
to test that everything is configured properly. To do this, type::

    $ make PETSC_DIR=$HOME/petsc-3.4.5 PETSC_ARCH=arch-linux2-cxx-debug test

To use PETSc, you have to define the ``PETSC_DIR`` and ``PETSC_ARCH``
environment variables to match how PETSc was built::

    $ export PETSC_DIR=$HOME/petsc-3.4.5
    $ export PETSC_ARCH=arch-linux2-cxx-debug

and add to your startup file ``$HOME/.bashrc``::

    export PETSC_DIR=$HOME/petsc-3.4.5
    export PETSC_ARCH=arch-linux2-cxx-debug

To configure BOUT++ with PETSc, go to the BOUT++ root directory, and
type::

    $ ./configure --with-petsc

You can configure BOUT++ against different PETSc installations either
through the ``PETSC_DIR/ARCH`` variables as above, or by giving
``--with-petsc`` a path::

  $ ./configure --with-petsc=/path/to/other/petsc

To configure BOUT++ with PETSc and SUNDIALS, type instead::

    $ ./configure --with-petsc --with-sundials


LAPACK
------

BOUT++ comes with linear solvers for tridiagonal and band-diagonal
systems, but these are not particularly optimised and are in any case
descended from Numerical Recipes code (hence NOT covered by LGPL
license).

To replace these routines, BOUT++ can use the LAPACK library. This is
however written in FORTRAN 77, which can cause linking headaches. To
enable these routines use::

    $ ./configure --with-lapack

and to specify a non-standard path::

    $ ./configure --with-lapack=/path/to/lapack

.. _sec-mumps:

MUMPS
-----

This is still experimental, but does work on at least some systems at
York. The PETSc library can be used to call MUMPS for directly solving
matrices (e.g. for Laplacian inversions), or MUMPS can be used directly.
To enable MUMPS, configure with::

    $ ./configure --with-mumps

MUMPS has many dependencies, including ScaLapack and ParMetis, which the
configuration script assumes are in the same place as MUMPS. The easiest
way to get MUMPS installed is to install PETSc with MUMPS, as the
configuration script will check the PETSc directory.

MPI compilers
-------------

These are usually called something like mpicc and mpiCC (or mpicxx), and
the configure script will look for several common names. If your
compilers aren’t recognised then set them using::

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
used by running::

     $ which install

This should point to a system directory like ``/usr/bin/install``.
Sometimes when IDL has been installed, this points to the IDL install
(e.g. something like ``/usr/common/usg/idl/idl70/bin/install`` on
Franklin). A quick way to fix this is to create a link from your local
bin to the system install::

     $ ln -s /usr/bin/install $HOME/local/bin/

“which install” should now print the install in your local bin
directory.

Compiling cvode.cxx fails
~~~~~~~~~~~~~~~~~~~~~~~~~

Occasionally compiling the CVODE solver interface will fail with an
error similar to::

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
