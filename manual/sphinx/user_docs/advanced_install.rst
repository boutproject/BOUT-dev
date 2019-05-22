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

Configure with ``--enable-checks=3`` enables a lot of checks of
operations performed by the field objects. This is very useful for
debugging a code, and can be omitted once bugs have been removed.
``--enable=checks=2`` enables less checking, especially the
computationally rather expensive ones, while ``--enable-checks=0``
disables most checks.

To get most checking, both from BOUT++ and from the compiler
``--enable-debug`` can be used. That enables checks of level 3, as
well as debug flags, e.g. ``-g`` for gcc.

For (sometimes) more useful error messages, there is the
``--enable-track`` option. This keeps track of the names of variables
and includes these in error messages.

To enable optimization, configure with ``--enable-optimize=3``.
This will try to set appropriate flags, but may not set the best ones.
This should work well for gcc. Similar to checks, different levels can
be specified, where 3 is high, and 0 means disabling all
optimization. ``--enable-optimize=fast`` will set the ``-Ofast`` flag
for gcc which enables optimizations that are not standard conforming, so
proceed at own risk.

Manually set compilation flags
------------------------------

You can set the following environment variables if you need more
control over how BOUT++ is built:

- ``LDFLAGS``: extra flags for linking, e.g. ``-L<library dir>``

- ``LIBS``: extra libraries for linking, e.g. ``-l<library>``

- ``CPPFLAGS``: preprocessor flags, e.g. ``-I<include dir>``

- ``CXXFLAGS``: compiler flags, e.g. ``-Wall``

- ``SUNDIALS_EXTRA_LIBS`` specifies additional libraries for linking
  to SUNDIALS, which are put at the end of the link command.

It is possible to change flags for BOUT++ after running configure, by
editing the ``make.config`` file. Note that this is not recommended,
as e.g. PVODE will not be built with these flags.

.. _sec-machine-specific:

Machine-specific installation
-----------------------------

These are some configurations which have been found to work on
particular machines.

Archer
~~~~~~

As of 20th April 2018, the following configuration should work

.. code-block:: bash

    $ module swap PrgEnv-cray PrgEnv-gnu/5.1.29
    $ module load fftw
    $ module load archer-netcdf/4.1.3

KNL @ Archer
~~~~~~~~~~~~

To use the KNL system, configure BOUT++ as follows:

.. code-block:: bash

    ./configure MPICXX=CC --host=knl --with-netcdf --with-pnetcdf=no --with-hypre=no CXXFLAGS="-xMIC-AVX512 -D_GLIBCXX_USE_CXX11_ABI=0"

Atlas
~~~~~

.. code-block:: bash

   ./configure --with-netcdf=/usr/local/tools/hdf5-gnu-serial-1.8.1/lib --with-fftw=/usr/local --with-pdb=/usr/gapps/pact/new/lnx-2.5-ib/gnu

Cab
~~~

.. code-block:: bash

   ./configure --with-netcdf=/usr/local/tools/hdf5-gnu-serial-1.8.1/lib --with-fftw=/usr/local/tools/fftw3-3.2 --with-pdb=/usr/gapps/pact/new/lnx-2.5-ib/gnu

Edison
~~~~~~

.. code-block:: bash

   module swap PrgEnv-intel PrgEnv-gnu
   module load fftw
   ./configure MPICC=cc MPICXX=CC --with-netcdf=/global/u2/c/chma/PUBLIC/netcdf_edison/netcdf --with-fftw=/opt/fftw/3.3.0.1/x86_64

Hoffman2
~~~~~~~~

.. code-block:: bash

   ./configure --with-netcdf=/u/local/apps/netcdf/current --with-fftw=/u/local/apps/fftw3/current --with-cvode=/u/local/apps/sundials/2.4.0 --with-lapack=/u/local/apps/lapack/current

Hopper
~~~~~~

.. code-block:: bash

    module swap PrgEnv-pgi PrgEnv-gnu
    module load netcdf
    module swap netcdf netcdf/4.1.3
    module swap gcc gcc/4.6.3
    ./configure MPICC=cc MPICXX=CC --with-fftw=/opt/fftw/3.2.2.1 --with-pdb=/global/homes/u/umansky/PUBLIC/PACT_HOPP2/pact

Hyperion
~~~~~~~~

With the bash shell use

.. code-block:: bash

   export PETSC_DIR=~farley9/projects/petsc/petsc-3.2-p1
   export PETSC_ARCH=arch-c
   ./configure --with-netcdf=/usr/local/tools/netcdf-gnu-4.1 --with-fftw=/usr/local MPICXX=mpiCC EXTRA_LIBS=-lcurl --with-petsc --with-cvode=~farley9/local --with-ida=~farley9/local

With the tcsh shell use

.. code-block:: tcsh

   setenv PETSC_DIR ~farley9/projects/petsc/petsc-3.2-p1
   setenv PETSC_ARCH arch-c
   ./configure --with-netcdf=/usr/local/tools/netcdf-gnu-4.1 --with-fftw=/usr/local MPICXX=mpiCC EXTRA_LIBS=-lcurl --with-petsc --with-cvode=~farley9/local --with-ida=~farley9/local

Marconi
~~~~~~~

.. code-block:: bash

   module load intel intelmpi fftw lapack
   module load szip zlib/1.2.8--gnu--6.1.0
   module load hdf5/1.8.17--intel--pe-xe-2017--binary
   module load netcdf-cxx4
   module load python

To compile for the SKL partition, configure with

.. code-block:: bash

   ./configure --enable-checks=0 CPPFLAGS="-Ofast -funroll-loops -xCORE-AVX512 -mtune=skylake" --host skl

to enable AVX512 vectorization.

.. note:: As of 20/04/2018, an issue with the netcdf and netcdf-cxx4
          modules means that you will need to remove ``-lnetcdf`` from
          ``EXTRA_LIBS`` in ``make.config`` after running
          ``./configure`` and before running ``make``. ``-lnetcdf``
          needs also to be removed from ``bin/bout-config`` to allow a
          successful build of the python interface. Recreation of
          ``boutcore.pyx`` needs to be manually triggered, if
          ``boutcore.pyx`` has already been created.

Ubgl
~~~~

.. code-block:: bash

   ./configure --with-netcdf CXXFLAGS=-DMPICH_IGNORE_CXX_SEEK CFLAGS=-DMPICH_IGNORE_CXX_SEEK --with-pdb=/usr/gapps/pact/new_s/lnx-2.5-ib --with-netcdf=/usr/local/tools/netcdf/netcdf-4.1_c++


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

.. _sec-netcdf-from-source:

Installing NetCDF from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The latest versions of NetCDF have separated out the C++ API from the
main C library. As a result, you will need to download and install both.
Download the latest versions of the NetCDF-C and NetCDF-4 C++ libraries
from https://www.unidata.ucar.edu/downloads/netcdf. As of
January 2017, these are versions 4.4.1.1 and 4.3.0 respectively.

Untar the file and ’cd’ into the resulting directory::

    $ tar -xzvf netcdf-4.4.1.1.tar.gz
    $ cd netcdf-4.4.1.1

Then run ``configure``, ``make`` and ``make install``::

    $ ./configure --prefix=$HOME/local
    $ make
    $ make install

Sometimes configure can fail, in which case try disabling Fortran::

    $ ./configure --prefix=$HOME/local --disable-fortran
    $ make
    $ make install

Similarly for the C++ API::

    $ tar -xzvf netcdf-cxx4-4.3.0.tar.gz
    $ cd netcdf-cxx4-4.3.0
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
``--enable-openmp`` flag to configure::

    ./configure --enable-openmp

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

    ./configure --enable-openmp --with-openmp-schedule=<schedule>

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

.. note:: BOUT++ currently supports SUNDIALS > 2.6, up to 4.1.0 as of
          March 2019. It is advisable to use the highest possible
          version

In order for a smooth install it is recommended to install SUNDIALS
from an install directory. The full installation guide is found in the
downloaded ``.tar.gz``, but we will provide a step-by-step guide to
install it and make it compatible with BOUT++ here::

     $ cd ~
     $ mkdir -p install/sundials-install
     $ cd install/sundials-install
     $ # Move the downloaded sundials-4.1.0.tar.gz to sundials-install
     $ tar -xzvf sundials-4.1.0.tar.gz
     $ mkdir build && cd build

     $ cmake \
       -DCMAKE_INSTALL_PREFIX=$HOME/local \
       -DLAPACK_ENABLE=ON \
       -DOPENMP_ENABLE=ON \
       -DMPI_ENABLE=ON \
       -DCMAKE_C_COMPILER=$(which mpicc) \
       -DCMAKE_CXX_COMPILER=$(which mpicxx) \
       ../sundials-4.1.0

     $ make
     $ make test
     $ make install

The SUNDIALS IDA solver is a Differential-Algebraic Equation (DAE)
solver, which evolves a system of the form
:math:`\mathbf{f}(\mathbf{u},\dot{\mathbf{u}},t) = 0`. This allows
algebraic constraints on variables to be specified.

To configure BOUT++ with SUNDIALS only (see section
:ref:`sec-PETSc-install` on how to build PETSc with SUNDIALS), go to
the root directory of BOUT++ and type::

    $ ./configure --with-sundials=/path/to/sundials/install

SUNDIALS will allow you to select at run-time which solver to use. See
:ref:`sec-timeoptions` for more details on how to do this.

.. _sec-PETSc-install:

PETSc
-----

BOUT++ can use PETSc https://www.mcs.anl.gov/petsc/ for time-integration
and for solving elliptic problems, such as inverting Poisson and
Helmholtz equations.

Currently, BOUT++ supports PETSc versions 3.4 - 3.9. To install PETSc
version 3.4.5, use the following steps::

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
              --download-fblas-lapack=1 \
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
through the ``PETSC_DIR/ARCH`` variables as above, or by specifying
them on the command line::

  $ ./configure --with-petsc PETSC_DIR=/path/to/other/petsc PETSC_ARCH=other-arch

.. note:: Unfortunately, there are a variety of ways PETSc can be
          installed on a system, and it is hard to automatically work
          out how to compile against a particular installation. In
          particular, there are two PETSc-supported ways of installing
          PETSc that are subtly different.

          The first way is as above, using ``PETSC_DIR`` and
          ``PETSC_ARCH``. A second way is to use the ``--prefix``
          argument to ``configure`` (much like the traditional GNU
          ``configure`` scripts) when building PETSc. In this case,
          ``PETSC_DIR`` will be the path passed to ``--prefix`` and
          ``PETSC_ARCH`` will be empty. When configuring BOUT++, one
          can use ``--with-petsc=$PETSC_DIR`` as a shortcut in this
          case. This will NOT work if PETSc was installed with a
          ``PETSC_ARCH``.

          However, there are at least some Linux distributions that
          install PETSc in yet another way and you may need to set
          ``PETSC_DIR/ARCH`` differently. For example, for Fedora, as
          of May 2018, you will need to configure and build BOUT++
          like so::

            $ ./configure --with-petsc=/usr/lib64/openmpi
            $ PETSC_DIR=/usr make

          Replace `openmpi` with the correct MPI implementation that
          you have installed.

.. _sec-lapack:

LAPACK
------

BOUT++ comes with linear solvers for tridiagonal and band-diagonal
systems. Some implementations of these solvers (for example Laplacian
inversion, section :ref:`sec-laplacian`) use LAPACK for efficient
serial performance. This does not add new features, but may be faster
in some cases. LAPACK is however written in FORTRAN 77, which can
cause linking headaches. To enable these routines use::

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

MUMPS has many dependencies, including ScaLapack and
ParMetis. Unfortunately, the exact dependencies and configuration of
MUMPS varies a lot from system to system. The easiest way to get MUMPS
installed is to install PETSc with MUMPS, or supply the ``CPPFLAGS``,
``LDFLAGS`` and ``LIBS`` environment variables to ``configure``::

   $ ./configure --with-mumps CPPFLAGS=-I/path/to/mumps/includes \
       LDFLAGS=-L/path/to/mumps/libs \
       LIBS="-ldmumps -lmumps_common -lother_libs_needed_for_mumps"

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

.. _sec-mpi-from-source:

Installing MPICH from source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In your home directory, create
two subdirectories: One called “install” where we’ll put the source
code, and one called “local” where we’ll install the MPI compiler::

    $ cd
    $ mkdir install
    $ mkdir local

Download the latest stable version of MPICH from https://www.mpich.org/ and put the
file in the “install” subdirectory created above. At the time of writing
(January 2018), the file was called ``mpich-3.2.1.tar.gz``. Untar the file::

    $ tar -xzvf mpich-3.2.1.tar.gz

which will create a directory containing the source code. ’cd’ into this
directory and run::

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
environment variable. To do this, run::

    $ export PATH=$PATH:$HOME/local/bin

and add this to the end of your startup file ``$HOME/.bashrc``. If
you’re using CSH rather than BASH, the command is::

    % setenv PATH ${PATH}:${HOME}/local/bin

and the startup file is ``$HOME/.cshrc``. You should now be able to run
``mpicc`` and so have a working MPI compiler.

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
