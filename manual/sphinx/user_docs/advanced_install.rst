.. Use bash as the default language for syntax highlighting in this file
.. highlight:: console

.. _sec-advancedinstall:

Advanced installation options
=============================

This section describes some common issues encountered when configuring
and compiling BOUT++, how to manually install dependencies if they are
not available, and how to configure optional libraries like
SUNDIALS and PETSc.

Optimisation and run-time checking
----------------------------------

Configure with ``-DCHECK=3`` enables a lot of checks of
operations performed by the field objects. This is very useful for
debugging a code, and can be omitted once bugs have been removed.
``-DCHECK=2`` enables less checking, especially the
computationally rather expensive ones, while ``-DCHECK=0``
disables most checks.

For (sometimes) more useful error messages, there is the
``-DBOUT_ENABLE_TRACK=ON`` option. This keeps track of the names of
variables and includes these in error messages.

To get a backtrace, you can set the environment variable
``BOUT_SHOW_BACKTRACE`` in order for the exception to include the
backtrace.

To enable optimization, configure with appropriate flags for your
compiler, e.g. with ``-DCMAKE_CXX_FLAGS=" -O3 "`` for a gnu compiler.


Install dependencies:
---------------------

BOUT++ provides a way to install some (optional) dependencies that are
not always found on HPC systems. To do this, run from your BOUT++
source directory:

.. code-block:: bash

    bin/bout-build-deps.sh
    # or without any checks:
    CHECK=no bin/bout-build-deps.sh
    # or with openmp - not tested, maybe not good to add it to FFTW
    PETSCFLAGS=--with-openmp=1 FFTWFLAGS="--enable-avx512 --enable-avx-128-fma --with-openmp --enable-threads" bin/bout-build-deps.sh
    # and add "-DBOUT_ENABLE_OPENMP=ON" to cmake configure line

Infos about options and further info can be obtained by running:

.. code-block:: bash

    bin/bout-build-deps.sh --help

If the script fails, it might be fixed by removing the folders that
are used for compiling and installing, and start again.

.. _sec-machine-specific:

Machine-specific installation
-----------------------------

These are some configurations which have been found to work on
particular machines. There is also the repo
https://github.com/boutproject/BOUT-configs which provides scripts for one or
two line compilation, with dependencies, of known-good versions on several
machines (different machines are in different branches).

Archer
~~~~~~

As of 20th April 2018, the following configuration should work

.. code-block:: bash

    $ module swap PrgEnv-cray PrgEnv-gnu/5.1.29
    $ module load fftw
    $ module load archer-netcdf/4.1.3

When using CMake on Cray systems like Archer, you need to pass
``-DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment`` so that the Cray compiler
wrappers are detected properly.

Cori
~~~~

First set up the environment by loading the correct modules. For Bash shell use:

.. code-block:: bash

   source config/cori/setup-env-cgpu.sh

and for C shell:

.. code-block:: csh

   source config/cori/setup-env-cgpu.sh

Then configure BOUT++ by running a script which calls CMake. Under bash:

.. code-block:: bash

   ./config/cori/config-bout-cgpu.sh

and C shell:

.. code-block:: csh

   ./config/cori/config-bout-cgpu.csh

At the time of writing, Hypre linking is not working with CUDA. If you come across
errors with the above configuration, try turning off Hypre support:

.. code-block:: bash

   ./config/cori/config-bout-cgpu-nohypre.sh

or

.. code-block:: csh

   ./config/cori/config-bout-cgpu-nohypre.csh

See section :ref:`sec-gpusupport` for details of compiling and running
on GPU machines, including Cori. Note that in order to access GPU
nodes a request must be made through `NERSC services
<https://nersc.servicenowservices.com/>`_.

MacOS / Apple Darwin
~~~~~~~~~~~~~~~~~~~~

Compiling with Apple Clang 12, the following configuration has been known to work

.. code-block:: tcsh

   cmake . -B <build-directory> -DBOUT_ENABLE_BACKTRACE=Off -DBUILD_SHARED_LIBS=Off -DBOUT_USE_NLS=Off -DBOUT_USE_UUID_SYSTEM_GENERATOR=Off
   cd <build-directory>
   cmake --build <build-directory>

where ``<build-directory>`` is the path to the build directory


MPCDF HPC Systems
~~~~~~~~~~~~~~~~~
After cloning BOUT-dev and checking out the branch you want (e.g. db-outer), run:
.. code-block:: bash

    module purge # or at least onload intel
    module load gcc/13 anaconda/3/2021.11 impi/2021.9 hdf5-serial/1.12.2 mkl/2022.0 netcdf-serial/4.8.1 fftw-mpi/3.3.10
    BUILD=/ptmp/$USER/bout-deps NO_HDF5=1 NO_NETCDF=1 NO_FFTW=1 bin/bout-build-deps.sh

and follow the instructions for configuring BOUT++. To enable openMP
for a production run use:

.. code-block:: bash

    module load bout-dep
    cmake .. -DBOUT_USE_NETCDF=ON -DnetCDFCxx_ROOT=$BOUT_DEP \
      -DBOUT_USE_PETSC=ON -DPETSC_DIR=$BOUT_DEP \
      -DBOUT_USE_FFTW=ON \
      -DBOUT_USE_SUNDIALS=ON -DSUNDIALS_ROOT=$BOUT_DEP \
      -DBOUT_ENABLE_OPENMP=OFF \
      -DCMAKE_BUILD_TYPE=Release


File formats
------------

BOUT++ can currently use the NetCDF-4_ file format and the ADIOS2 library
for high-performance parallel output.

NetCDF is a widely used format and
has many tools for viewing and manipulating files.

.. _NetCDF-4: https://www.unidata.ucar.edu/software/netcdf/

BOUT++ will look for ``ncxx4-config`` or ``nc-config`` in your
``$PATH``. If it cannot find the libraries, or finds a different
version than the one you want, you can point it at the correct version
using::

   cmake -S .. -B . -DBOUT_USE_NETCDF=ON -DnetCDFCxx_ROOT=/path/to/ncxx4-config

where ``/path/to/ncxx4-config`` is the location of the
``ncxx4-config`` tool (``nc-config`` will also work, but
``ncxx4-config`` is preferred).


.. _sec-netcdf-from-source:

Installing NetCDF from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The latest versions of NetCDF have separated out the C++ API from the
main C library. As a result, you will need to download and install both.
Download the latest versions of the NetCDF-C and NetCDF-4 C++ libraries
from https://www.unidata.ucar.edu/downloads/netcdf. As of
September 2020, these are versions 4.7.4 and 4.3.1 respectively.

Untar the file and ’cd’ into the resulting directory::

    $ tar -xzvf netcdf-4.7.4.tar.gz
    $ cd netcdf-4.7.4

Then run ``configure``, ``make`` and ``make install``::

    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

Sometimes configure can fail, in which case try disabling Fortran::

    $ ./configure --prefix=$HOME/local --disable-fortran
    $ make
    $ make install

Similarly for the C++ API::

    $ tar -xzvf netcdf-cxx4-4.3.1.tar.gz
    $ cd netcdf-cxx4-4.3.1
    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

You may need to set a couple of environment variables as well::

    $ export PATH=$HOME/local/bin:$PATH
    $ export LD_LIBRARY_PATH=$HOME/local/lib:$LD_LIBRARY_PATH

You should check where NetCDF actually installed its libraries. On some
systems this will be ``$HOME/local/lib``, but on others it may be, e.g.
``$HOME/local/lib64``. Check which it is, and set ``$LD_LIBRARY_PATH``
appropriately.

OpenMP
------

BOUT++ can make use of OpenMP parallelism. To enable OpenMP, use the
``-DBOUT_ENABLE_OPENMP=ON`` flag to configure::

    cmake -S .. -B . -DBOUT_ENABLE_OPENMP=ON

OpenMP can be used to parallelise in more directions than can be
achieved with MPI alone. For example, it is currently difficult to
parallelise in X using pure MPI if FCI is used, and impossible to
parallelise at all in Z with pure MPI.

OpenMP is in a large number of places now, such that a decent speed-up
can be achieved with OpenMP alone. Hybrid parallelisation with both
MPI and OpenMP can lead to more significant speed-ups, but it
sometimes requires some fine tuning of numerical parameters in order
to achieve this. This greatly depends on the details not just of your
system, but also your particular problem. We have tried to choose
"sensible" defaults that will work well for the most common cases, but
this is not always possible. You may need to perform some testing
yourself to find e.g. the optimum split of OpenMP threads and MPI
ranks.

One such parameter that can potentially have a significant effect (for
some problem sizes on some machines) is setting the OpenMP schedule
used in some of the OpenMP loops (specifically those using
`BOUT_FOR`). This can be set using::

    cmake . -DBOUT_ENABLE_OPENMP=ON -DBOUT_OPENMP_SCHEDULE=<schedule>

with ``<schedule>`` being one of: ``static`` (the default),
``dynamic``, ``guided``, ``auto`` or ``runtime``.


.. note::
    If you want to use OpenMP with Clang, you will need Clang 3.7+,
    and either ``libomp`` or ``libiomp``.

    You will be able to compile BOUT++ with OpenMP with lower versions
    of Clang, or using the GNU OpenMP library ``libgomp``, but it will
    only run with a single thread.


.. note::
    By default PVODE is built without OpenMP support. To enable this
    add ``--enable-pvode-openmp`` to the configure command.


.. note::
    OpenMP will attempt to use all available threads by default. This
    can cause oversubscription problems on certain systems. You can
    limit the number of threads OpenMP uses with the
    ``OMP_NUM_THREADS`` environment variable. See your system
    documentation for more details.

.. _sec-sundials:

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

.. note:: BOUT++ currently supports SUNDIALS > 2.6, up to 6.7.0 as of
          January 2024. It is advisable to use the highest possible
          version. Support for SUNDIALS versions < 4 will be removed
          in the next release.

The full installation guide is found in the downloaded ``.tar.gz``,
but we will provide a step-by-step guide to install it and make it
compatible with BOUT++ here::

     $ tar -xzvf sundials-5.4.0.tar.gz
     $ cd sundials-5.4.0
     $ mkdir build && cd build

     $ cmake .. \
       -DCMAKE_INSTALL_PREFIX=$HOME/local \
       -DLAPACK_ENABLE=ON \
       -DOPENMP_ENABLE=ON \
       -DMPI_ENABLE=ON \
       -DCMAKE_C_COMPILER=$(which mpicc) \
       -DCMAKE_CXX_COMPILER=$(which mpicxx) \

     $ make
     $ make test
     $ make install

The SUNDIALS IDA solver is a Differential-Algebraic Equation (DAE)
solver, which evolves a system of the form
:math:`\mathbf{f}(\mathbf{u},\dot{\mathbf{u}},t) = 0`. This allows
algebraic constraints on variables to be specified.

Use the ``-DBOUT_USE_SUNDIALS=ON -DSUNDIALS_ROOT=`` option to configure BOUT++ with SUNDIALS::

    $ cmake . -DBOUT_USE_SUNDIALS=ON -DSUNDIALS_ROOT=/path/to/sundials/install

SUNDIALS will allow you to select at run-time which solver to use. See
:ref:`sec-timeoptions` for more details on how to do this.

Notes:

* If compiling SUNDIALS, make sure that it is configured with MPI (``MPI_ENABLE=ON``)

.. _sec-PETSc-install:

PETSc
-----

BOUT++ can use PETSc https://www.mcs.anl.gov/petsc/ for time-integration
and for solving elliptic problems, such as inverting Poisson and
Helmholtz equations.

Currently, BOUT++ supports PETSc versions 3.7 - 3.23. More recent versions may
well work, but the PETSc API does sometimes change in backward-incompatible
ways, so this is not guaranteed. To install PETSc version 3.19, use the
following steps::

    $ cd ~
    $ wget https://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.19.1.tar.gz
    $ tar -xzvf petsc-3.19.1.tar.gz
    $ cd petsc-3.19.1

Use the following configure options to ensure PETSc is compatible with BOUT++::

    $ ./configure \
      --with-mpi=yes \
      --with-precision=double \
      --with-scalar-type=real \
      --with-shared-libraries=1 \
      --with-debugging=0 \
      {C,CXX,F}OPTFLAGS="-O3 -march=native" \
      --prefix=$HOME/local/petsc-version-options

You may also wish to change to ``--with-debugging=yes`` in the
arguments to ``./configure``, in order to allow debugging of PETSc.
The optimisation flags need changing for cross compiling or non gcc 
compilers. Set a different prefix to change the place PETSc will be
installed to.

.. note:: If you build BOUT++ using a standalone version of SUNDIALS,
          it is advisable to not also build PETSc with SUNDIALS.

.. note:: It is also possible to get PETSc to download and install
          MUMPS, by adding::

              --download-mumps \
              --download-scalapack \
              --download-blacs \
              --download-fblaslapack=1 \
              --download-parmetis \
              --download-ptscotch \
              --download-metis

          to ``./configure``.

To make PETSc, type what is shown in the terminal output after the configure
step, something like::

    $ make PETSC_DIR=$HOME/petsc-3.19.1 PETSC_ARCH=arch-linux2-cxx-debug all

Should BLAS, LAPACK, or any other packages be missing, you will get an
error, and a suggestion that you can append
``--download-name-of-package`` to the ``./configure`` line.

You may want to test that everything is configured properly. To do this replace
``all`` with ``test`` in the make command. It should be something like::

    $ make PETSC_DIR=$HOME/petsc-3.19.1 PETSC_ARCH=arch-linux2-cxx-debug test

To install PETSc, replace ``test``/``all`` with ``install`` and run
something like::

    $ make PETSC_DIR=$HOME/petsc-3.19.1 PETSC_ARCH=arch-linux2-cxx-debug install

To configure BOUT++ with PETSc, add to the cmake configure command::

    -DBOUT_USE_PETSC=ON -DPETSC_DIR=$HOME/local/petsc-version-options

For example like this::

    $ cmake -S . -B <build-directory> -DBOUT_USE_PETSC=ON -DPETSC_DIR=$HOME/local/petsc-version-options

BOUT++ can also work with PETSc if it has not been installed. In this
case ensure that ``PETSC_DIR`` and ``PETSC_ARCH`` are set, for example
like this::

    $ PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-linux2-cxx-debug cmake -DBOUT_USE_PETSC=ON

.. _sec-lapack:

LAPACK
------

BOUT++ comes with linear solvers for tridiagonal and band-diagonal
systems. Some implementations of these solvers (for example Laplacian
inversion, section :ref:`sec-laplacian`) use LAPACK for efficient
serial performance. This does not add new features, but may be faster
in some cases. LAPACK is however written in FORTRAN 77, which can
cause linking headaches. To enable these routines use::

    $ cmake -S . -B <build-directory> -DBOUT_USE_LAPACK=ON

and to specify a non-standard path::

    $ cmake -S . -B <build-directory> -DBOUT_USE_LAPACK=ON -DLAPACK_ROOT=/path/to/lapack


MPI compilers
-------------

These are usually called something like mpicc and mpiCC (or mpicxx), and
the configure script will look for several common names. If your
compilers aren’t recognised then check the `cmake documentation for MPI <https://cmake.org/cmake/help/latest/module/FindMPI.html#variables-for-locating-mpi>`_

NOTES:

-  On LLNL’s Grendel, mpicxx is broken. Use mpiCC instead by passing
   “MPICXX=mpiCC” to configure. Also need to specify this to NetCDF
   library by passing “CXX=mpiCC” to NetCDF configure.

.. _sec-fftw-from-source:

Installing FFTW from source
---------------------------

If you haven’t already, create directories “install” and “local”
in your home directory::

    $ cd
    $ mkdir install
    $ mkdir local

Download the latest stable version from
http://www.fftw.org/download.html into the “install” directory. At the
time of writing, this was called ``fftw-3.3.2.tar.gz``. Untar this file,
and ’cd’ into the resulting directory. As with the MPI compiler,
configure and install the FFTW library into ``$HOME/local`` by running::

    $ ./configure --prefix=$HOME/local
    $ make
    $ make install


Compiling and running under AIX
-------------------------------

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

SUNDIALS under AIX
~~~~~~~~~~~~~~~~~~

To compile SUNDIALS, use:

.. code-block:: bash

    export CC=cc
    export CXX=xlC
    export F77=xlf
    export OBJECT_MODE=64
    ./configure --prefix=$HOME/local/ --with-mpicc=mpcc --with-mpif77=mpxlf CFLAGS=-maix64

You may get an error message like

.. code-block:: bash

    make: Not a recognized flag: w

This is because the AIX ``make`` is being used, rather than ``gmake``.
The easiest way to fix this is to make a link to ``gmake`` in your local
bin directory

.. code-block:: bash

    ln -s /usr/bin/gmake $HOME/local/bin/make

Running ``which make`` should now point to this ``local/bin/make``, and
if not then you need to make sure that your bin directory appears first
in the ``PATH``

.. code-block:: bash

    export PATH=$HOME/local/bin:$PATH

If you see an error like this

.. code-block:: bash

    ar: 0707-126 ../../src/sundials/sundials_math.o is not valid with the current object file mode.
            Use the -X option to specify the desired object mode.


then you need to set the environment variable ``OBJECT_MODE``

.. code-block:: bash

    export OBJECT_MODE=64

Configuring BOUT++, you may get the error

.. code-block:: bash

    configure: error: C compiler cannot create executables

In that case, you can try using:

.. code-block:: bash

    ./configure CFLAGS="-maix64"

When compiling, you may see warnings:

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
has been compiled as 64-bit. You can try compiling BOUT++ as 32-bit

.. code-block:: bash

    export OBJECT_MODE=32
    ./configure CFLAGS="-maix32"
    gmake

If you still get undefined symbols, then go back to 64-bit, and edit
make.config, replacing ``-lnetcdf_c++`` with -lnetcdf64\_c++, and
``-lnetcdf`` with -lnetcdf64. This can be done by running

.. code-block:: bash

     sed 's/netcdf/netcdf64/g' make.config > make.config.new
     mv make.config.new make.config

Compiling on Windows
~~~~~~~~~~~~~~~~~~~~

It is possible to compile BOUT++ on Windows using the CMake
interface. Support is currently very experimental, and some features do
not work. Testing has been done with MSVC 19.24 and Visual Studio 16.4,
although previous versions may still work.

The main difficulty of using BOUT++ on Windows is getting the
dependencies sorted. The easiest way to install dependencies on Windows
is using `vcpkg <https://github.com/microsoft/vcpkg/>`_. You may need to
set the CMake toolchain file if calling ``cmake`` from PowerShell, or on
older versions of Visual Studio. This will be a file somewhere like
``C:/vcpkg/scripts/buildsystems/vcpkg.cmake``

The minimal required CMake options are as follows:

.. code-block:: bash

    -DBOUT_ENABLE_BACKTRACE=OFF \
    -DCMAKE_CXX_FLAGS="/permissive- /EHsc /bigobj" \
    -DBUILD_SHARED_LIBS=OFF

``ENABLE_BACKTRACE`` must be turned off due to the currently required
``addr2line`` executable not being available on Windows.

The following flags for the MSVC compiler are required:

- ``/permissive-`` for standards compliance, such as treating the binary
  operator alternative tokens (``and``, ``or``, etc) as tokens
- ``/EHsc`` for standard C++ exception handling, and to assume that
  ``extern "C"`` functions never throw
- ``/bigobj`` to increase the number of sections in the .obj file,
  required for the template-heavy derivatives machinery

No modification to the source has been done to export the correct
symbols for shared libraries on Windows, so you must either specifiy
``-DBUILD_SHARED_LIBS=OFF`` to only build static libraries, or, if you
really want shared libraries, ``-DCMAKE_WINDOWS_EXPORT_ALL_SYMBOLS=ON``.
The latter is untested, use at your own risk!

The unit tests should all pass, but most of the integrated tests will
not run work out of the box yet as Windows doesn't understand
shabangs. That is, without a file extension, it doesn't know what
program to use to run ``runtest``. The majority of the tests can be
run manually with ``python.exe runtest``. You will stil need to set
``PYTHONPATH`` and have a suitable Python environment.

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

Compiling with IBM xlC compiler fails
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When using the ``xlC`` compiler, an error may occur::

  variant.hpp(1568) parameter pack "Ts" was referenced but not expanded


The workaround is to change line 428 of  ``externalpackages/mpark.variant/include/mpark/lib.hpp`` from::

  #ifdef MPARK_TYPE_PACK_ELEMENT

to::

  #ifdef CAUSES_ERROR // MPARK_TYPE_PACK_ELEMENT

This will force an alternate implementation of type_pack_element to be defined.
See also https://software.intel.com/en-us/forums/intel-c-compiler/topic/501502


Compiling fails after changing branch
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If compiling fails after changing branch, for example from ``master``
to ``next``, with an error like the following::

   $ make
   Downloading mpark.variant
   You need to run this command from the toplevel of the working tree.
   make[2]: *** [BOUT-dev/externalpackages/mpark.variant/include/mpark/variant.hpp] Error 1
   make[1]: *** [field] Error 2
   make: *** [src] Error 2

it's possible something has gone wrong with the submodules. To fix,
just run::

  $ git submodule update --init --recursive  ./externalpackages/*

If you regularly work on two different branches and need to run the above command a lot, you may consider telling git to automatically
update the submodules::

  git config submodule.recurse=true

This requires ``git >= 2.14``.
