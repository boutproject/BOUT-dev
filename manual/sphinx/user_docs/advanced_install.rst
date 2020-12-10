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

To get a backtrace, you can set the environment variable
`BOUT_SHOW_BACKTRACE` in order for the exception to include the
backtrace.

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

When using CMake on Cray systems like Archer, you need to pass
``-DCMAKE_SYSTEM_NAME=CrayLinuxEnvironment`` so that the Cray compiler
wrappers are detected properly.

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

Marconi with gnu compilers
**************************

It is also possible to configure on Marconi using gnu compilers, which may give better performance. A set of modules which work as of 30/9/2020 is

.. code-block:: bash

    module load env-skl
    module load profile/advanced
    module load intel/pe-xe-2018--binary  # note need to keep the 'intel' module loaded in order for shared libraries needed by numpy/scipy to be available
    module load gnu/7.3.0
    module load openmpi/4.0.1--gnu--7.3.0
    module load mkl/2017--binary
    module load python/3.6.4
    module load szip/2.1--gnu--6.1.0 zlib/1.2.8--gnu--6.1.0

Then download source code for hdf5-1.12.0 (hdf5 is available in a module on
Marconi, but has issues linking OpenMPI), netCDF-c-4.7.4, netCDF-cxx4-4.3.1,
and FFTW-3.3.8. Optionally also SUNDIALS-5.1.0 or PETSc-3.13.0. Configure and
compile all of the downloaded packages. Make sure to install netCDF and
netCDF-cxx4 into the same directory (this is assumed by netCDF's linking
strategy, and makes netCDF configuration simpler).

The following configuration commands have been used successfully:

* hdf5-1.12.0::

    ./configure --prefix /directory/to/install/hdf5 --enable-build-mode=production
    make
    make install

* netCDF-4.7.4::

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/directory/to/install/netcdf -DCMAKE_BUILD_TYPE=Release ..
    make
    make install

* netCDF-cxx4-4.3.1::

    mkdir build
    cd build
    cmake -DCMAKE_INSTALL_PREFIX=/directory/to/install/netcdf -DCMAKE_BUILD_TYPE=Release ..
    make
    make install

* FFTW-3.3.8::

    ./configure --prefix /directory/to/install/fftw --enable-shared --enable-sse2 --enable-avx --enable-avx2 --enable-avx512 --enable-avx-128-fma
    make
    make install

* SUNDIALS-5.1.0::

    mkdir build
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/directory/to/install/sundials -DMPI_ENABLE=ON ..
    make
    make install

* PETSc-3.13.0::

    unset PETSC_DIR
    ./configure COPTFLAGS="-O3" CXXOPTFLAGS="-O3" FOPTFLAGS="-O3" --with-batch --known-mpi-shared-libraries=1 --with-mpi-dir=$OPENMPI_HOME --download-fblaslapack --known-64-bit-blas-indices=0 --download-hypre --with-debugging=0 --prefix=/directory/to/install/petsc

  then follow the instructions printed by PETSc at the end of each step to make, install and check the build.

Finally example configurations for BOUT++, where you should replace <...> by appropriate directories that you used to install the libraries:

* for an optimized build (some experimentation with optimisation flags would be welcome, please share the results if you do!)::

    ./configure --enable-optimize=3 --enable-checks=no --without-hdf5 --enable-static --with-netcdf=<...> --with-sundials=<...> --with-fftw=<...> --with-petsc=<...>

* for a debugging build::

    ./configure --enable-debug --without-hdf5 --enable-static --with-netcdf=<...> --with-sundials=<...> --with-fftw=<...> --with-petsc=<...>

Ubgl
~~~~

.. code-block:: bash

   ./configure --with-netcdf CXXFLAGS=-DMPICH_IGNORE_CXX_SEEK CFLAGS=-DMPICH_IGNORE_CXX_SEEK --with-pdb=/usr/gapps/pact/new_s/lnx-2.5-ib --with-netcdf=/usr/local/tools/netcdf/netcdf-4.1_c++


File formats
------------

BOUT++ can currently use two different file formats: NetCDF-4_, and
HDF5_ and experimental support for parallel flavours of both. NetCDF
is a widely used format and so has many more tools for viewing and
manipulating files. HDF5 is another widely used format. If you have
multiple libraries installed then BOUT++ can use them simultaneously,
for example reading in grid files in NetCDF format, but writing output
data in HDF5 format.

.. _NetCDF-4: https://www.unidata.ucar.edu/software/netcdf/
.. _HDF5: https://www.hdfgroup.org/HDF5/

BOUT++ will try to use NetCDF by default. It will look for
``ncxx4-config`` or ``nc-config`` in your ``$PATH``. If it cannot find
the libraries, or finds a different version than the one you want, you
can point it at the correct version using::

   ./configure --with-netcdf=/path/to/ncxx4-config

where ``/path/to/ncxx4-config`` is the location of the
``ncxx4-config`` tool (``nc-config`` will also work, but
``ncxx4-config`` is preferred).

To use HDF5, you will need to explicitly enable it::

   ./configure --with-hdf5

BOUT++ will look for ``h5cc`` in your ``$PATH``. Similar to NetCDF,
you can pass the location of the particular version you wish to use
with::

   ./configure --with-hdf5=/path/to/h5cc


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

.. note:: BOUT++ currently supports SUNDIALS > 2.6, up to 5.4.0 as of
          September 2020. It is advisable to use the highest possible
          version

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

Use the ``--with-sundials`` option to configure BOUT++ with SUNDIALS::

    $ ./configure --with-sundials=/path/to/sundials/install

SUNDIALS will allow you to select at run-time which solver to use. See
:ref:`sec-timeoptions` for more details on how to do this.

.. _sec-PETSc-install:

PETSc
-----

BOUT++ can use PETSc https://www.mcs.anl.gov/petsc/ for time-integration
and for solving elliptic problems, such as inverting Poisson and
Helmholtz equations.

Currently, BOUT++ supports PETSc versions 3.7 - 3.14. More recent versions may
well work, but the PETSc API does sometimes change in backward-incompatible
ways, so this is not guaranteed. To install PETSc version 3.13, use the
following steps::

    $ cd ~
    $ wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-3.13.4.tar.gz
    $ tar -xzvf petsc-3.13.4.tar.gz
    $ cd petsc-3.13.4

Use the following configure options to ensure PETSc is compatible with BOUT++::

    $ ./configure \
      --with-clanguage=cxx \
      --with-mpi=yes \
      --with-precision=double \
      --with-scalar-type=real \
      --with-shared-libraries=0

You may also wish to add ``--with-debugging=yes`` to ``./configure``
in order to allow debugging.

.. note:: If you build BOUT++ using a standalone version of SUNDIALS,
          it is advisable to not also build PETSc with SUNDIALS.

.. note:: It is also possible to get PETSc to download and install
          MUMPS, by adding::

              --download-mumps \
              --download-scalapack \
              --download-blacs \
              --download-fblas-lapack=1 \
              --download-parmetis \
              --download-ptscotch \
              --download-metis

          to ``./configure``.

To make PETSc, type::

    $ make PETSC_DIR=$HOME/petsc-3.13.4 PETSC_ARCH=arch-linux2-cxx-debug all

Should BLAS, LAPACK, or any other packages be missing, you will get an
error, and a suggestion that you can append
``--download-name-of-package`` to the ``./configure`` line. You may want
to test that everything is configured properly. To do this, type::

    $ make PETSC_DIR=$HOME/petsc-3.13.4 PETSC_ARCH=arch-linux2-cxx-debug test

To use PETSc, you have to define the ``PETSC_DIR`` and ``PETSC_ARCH``
environment variables to match how PETSc was built::

    $ export PETSC_DIR=$HOME/petsc-3.13.4
    $ export PETSC_ARCH=arch-linux2-cxx-debug

and add to your startup file ``$HOME/.bashrc``::

    export PETSC_DIR=$HOME/petsc-3.13.4
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
just run ``make submodules``::

  $ make submodules
  Downloading mpark.variant
  git submodule update --init --recursive /home/peter/Codes/BOUT-dev/externalpackages/mpark.variant
  Submodule path 'externalpackages/mpark.variant': checked out '0b488da9bebac980e7ba0e158a959c956a449676'

If you regularly work on two different branches and need to run ``make
submodules`` a lot, you may consider telling git to automatically
update the submodules::

  git config submodule.recurse=true

This requires ``git >= 2.14``.
