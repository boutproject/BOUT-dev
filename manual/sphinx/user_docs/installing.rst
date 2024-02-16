.. Use bash as the default language for syntax highlighting in this file
.. highlight:: console

.. _sec-install:

Getting started
===============

.. _sec-getting-started:

This section goes through the process of getting, installing, and
starting to run BOUT++.

The quickest way to get started is to use a pre-built binary. These
take care of all dependencies, configuration and compilation. See
section :ref:`sec-prebuiltinstall`.

The remainder of this section will go through the following steps to
manually install BOUT++. Only the basic functionality needed to use
BOUT++ is described here; the next section (:ref:`sec-advancedinstall`) goes
through more advanced options, configurations for particular machines,
and how to fix some common problems.

#. :ref:`Obtaining a copy of BOUT++ <sec-obtainbout>`

#. :ref:`Installing dependencies <sec-dependencies>`

#. :ref:`Configuring BOUT++ <sec-config-bout>`

#. :ref:`Configuring BOUT++ analysis codes <sec-configanalysis>`

   #. :ref:`Python <sec-config-python>`

   #. :ref:`IDL <sec-config-idl>`

#. :ref:`Compiling BOUT++ <sec-compile-bout>`

#. :ref:`Running the test suite <sec-runtestsuite>`

#. :ref:`Installing BOUT++ (experimental) <sec-install-bout>`

**Note**: In this manual commands to run in a BASH shell will begin with
’$’, and commands specific to CSH with a ’%’.

Pre-built binaries
------------------

.. _sec-prebuiltinstall:

Docker image
~~~~~~~~~~~~

`Docker <https://www.docker.com>`_ is a widely used container system,
which packages together the operating system environment, libraries
and other dependencies into an image. This image can be downloaded and
run reproducibly on a wide range of hosts, including Windows, Linux and OS X.
Here is the starting page for `instructions on installing Docker
<https://docs.docker.com/install/>`_.

The BOUT++ docker images are `hosted on dockerhub
<https://hub.docker.com/u/boutproject/>`_ for some releases and
snapshots. Check the `list of BOUT-next tags <https://hub.docker.com/r/boutproject/bout-next/tags/>`_
if you want a recent version of BOUT++ “next” (development) branch.
First download the image::

    $ sudo docker pull boutproject/boutproject/bout-next:9f4c663-petsc

then run::

    $ sudo docker run --rm -it boutproject/bout-next:9f4c663-petsc

This should give a terminal in a "boutuser" home directory, in which
there is "BOUT-next", containing BOUT++ configured and compiled with
NetCDF, SUNDIALS, PETSc and SLEPc. Python 3 is also installed, with
ipython, NumPy, Scipy and Matplotlib libaries. To plot to screen an
X11 display is needed. Alternatively a shared directory can be created
to pass files between the docker image and host. The following
commands both enable X11 and create a shared directory::

    $ mkdir shared
    $ sudo docker run --rm -it \
       -e DISPLAY -v $HOME/.Xauthority:/home/boutuser/.Xauthority --net=host \
       -v $PWD/shared:/home/boutuser/bout-img-shared \
       boutproject/bout-next:9f4c663-petsc

This should enable plotting from python, and files in the docker image
put in "/home/boutuser/bout-img-shared" should be visible on the host in
the "shared" directory.

If this is successful, then you can skip to section :ref:`sec-running`.

.. _sec-obtainbout:

Obtaining BOUT++
----------------

BOUT++ is hosted publicly on github at
https://github.com/boutproject/BOUT-dev. You can the latest stable
version from https://github.com/boutproject/BOUT-dev/releases. If you
want to develop BOUT++, you should use git to clone the repository. To
obtain a copy of the latest version, run::

    $ git clone https://github.com/boutproject/BOUT-dev.git


which will create a directory ``BOUT-dev`` containing the code::

    $ cd BOUT-dev

To get the latest changes later, go into the ``BOUT-dev`` directory and run::

    $ git pull

Development is done on the “next” branch, which you can checkout with::

    $ git checkout next

.. _sec-installmpi:

Installing dependencies
-----------------------

.. _sec-dependencies:

The bare-minimum requirements for compiling and running BOUT++ are:

#. A C++ compiler that supports C++17

#. An MPI compiler such as OpenMPI (`www.open-mpi.org/ <https://www.open-mpi.org/>`__),
   MPICH ( `https://www.mpich.org/ <https://www.mpich.org/>`__)

#. The NetCDF library (`https://www.unidata.ucar.edu/downloads/netcdf
   <https://www.unidata.ucar.edu/downloads/netcdf>`__)

The FFTW-3 library (`http://www.fftw.org/ <http://www.fftw.org/>`__)
is also strongly recommended. Fourier transforms are used for some
derivative methods, as well as the `ShiftedMetric` parallel transform
which is used in the majority of BOUT++ tokamak simulations. Without
FFTW-3, these options will not be available.

.. note::
   If you use an Intel compiler, you must also make sure that you have
   a version of GCC that supports C++17 (GCC 8+).

   On supercomputers, or in other environments that use a module
   system, you may need to load modules for both Intel and GCC.

On a cluster or supercomputer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are installing on a cluster or supercomputer then the MPI C++ compilers will
already be installed, and on Cray or IBM machines will probably be
called ``CC`` and ``xlC`` respectively.

On large facilities (e.g NERSC or Archer), the compilers and libraries
needed should already be installed, but you may need to load them to use them.
It is common to organise libraries using the ``modules`` system, so try typing::

   modules avail

to get a list of available modules. Some instructions for specific machines can
be found in :ref:`sec-machine-specific`. See your system’s
documentation on modules and which ones to load. If you don’t know, or
modules don’t work, you can still install libraries in your home
directory by following the instructions below for :ref:`FFTW <sec-fftw-from-source>`
and :ref:`NetCDF <sec-netcdf-from-source>`.


Ubuntu / Debian
~~~~~~~~~~~~~~~

On Ubuntu or Debian distributions if you have administrator rights then you can install
the basic dependencies with::

    $ sudo apt-get install libmpich-dev libfftw3-dev libnetcdf-c++4-dev git make

To additionally build the Python interface, you need some Python packages::

    $ sudo apt-get install python3 python3-distutils python3-pip python3-numpy python3-netcdf4 python3-scipy
    $ pip3 install --user Cython

Further, the encoding for python needs to be utf8 - it may be required
to set ``export LC_CTYPE=C.utf8``.

If you do not have administrator rights, so can't install packages, then
you need to install these libraries from source into your home directory.
See :ref:`sec-advancedinstall` for details on installing some of these.


Arch Linux
~~~~~~~~~~

::

   $ pacman -S openmpi fftw netcdf-cxx make gcc


Fedora
~~~~~~

On Fedora the required libraries can be installed by running::

   $ sudo dnf build-dep bout++

This will install all the dependencies that are used to install
BOUT++ for fedora. Feel free to install only a subset of the
suggested packages. For example, only mpich or openmpi is required.
To load an mpi implementation type::

   $ module load mpi

After that the mpi library is loaded.
Precompiled binaries are available for fedora as well.
To get precompiled BOUT++ run::

   $ # install the mpich version - openmpi is available as well
   $ sudo dnf install bout++-mpich-devel
   $ # get the python3 modules - python2 is available as well
   $ sudo dnf install python3-bout++

.. _sec-config-bout:

Configuring BOUT++
------------------

BOUT++ uses the `CMake <https://cmake.org/>`_ build system
generator. You will need CMake >= 3.17.

.. note::
   It is possible to get the latest version of CMake using ``pip``::

      $ pip install --user --upgrade cmake

   or ``conda``::

      $ conda install cmake

   You may need to put ``~/.local/bin`` in your ``$PATH``

CMake supports out-of-source builds by default, which are A Good Idea.
Basic configuration with CMake looks like::

  $ cmake -S . -B build

which creates a new directory ``build``. You can call this directory
anything you like, and you also put it anywhere you like, you just
need to specify the path to the BOUT++ source directory with the
``-S`` argument. This makes it very easy to keep two build directories
alongside one another, one with a debug build and one optimised, for
example.

After configuring the build directory, you can then compile BOUT++ with::

  # Build the library
  $ cmake --build build
  # Build the library with 8 threads
  $ cmake --build build -j 8
  # Build the "blob2d" example
  $ cmake --build build --target blob2d

By default, CMake will use ``makefiles``, and so it is possible to
also build BOUT++ with ``make`` from the build directory -- note that
you must still run ``cmake`` once first to configure BOUT++::

  $ cmake . -B build
  $ cd build
  $ make

.. note::
   You might see some instructions in the documentation using ``make``
   -- they should be run from the ``build`` directory.

You can see what build options are available with::

  $ cmake . -B build -LH
  ...
  // Enable backtrace
  BOUT_ENABLE_BACKTRACE:BOOL=ON

  // Output coloring
  BOUT_ENABLE_COLOR:BOOL=ON

  // Enable OpenMP support
  BOUT_ENABLE_OPENMP:BOOL=OFF

  // Enable support for PETSc time solvers and inversions
  BOUT_USE_PETSC:BOOL=OFF
  ...

CMake uses the ``-D<variable>=<choice>`` syntax to control these
variables. You can set ``<package>_ROOT`` to guide CMake in finding
the various optional third-party packages (except for PETSc/SLEPc,
which use ``_DIR``). Note that some packages have funny
captialisation, for example ``NetCDF_ROOT``! Use ``-LH`` to see the
form that each package expects.

CMake understands the usual environment variables for setting the
compiler, compiler/linking flags, as well as having built-in options
to control them and things like static vs shared libraries, etc. See
the `CMake documentation <https://cmake.org/documentation/>`_ for more
infomation.

A more complicated CMake configuration command
might look like::

  $ CC=mpicc CXX=mpic++ cmake . -B build \
      -DBOUT_USE_PETSC=ON -DPETSC_DIR=/path/to/petsc/ \
      -DBOUT_USE_SLEPC=ON -DSLEPC_DIR=/path/to/slepc/ \
      -DBOUT_USE_SUNDIALS=ON -DSUNDIALS_ROOT=/path/to/sundials \
      -DBOUT_USE_NETCDF=ON -DNetCDF_ROOT=/path/to/netcdf \
      -DBOUT_ENABLE_OPENMP=ON \
      -DBOUT_ENABLE_SIGFPE=OFF \
      -DCMAKE_BUILD_TYPE=Debug \
      -DBUILD_SHARED_LIBS=ON
      -DCMAKE_INSTALL_PREFIX=/path/to/install/BOUT++

If you wish to change the configuration after having built ``BOUT++``,
it's wise to delete the ``CMakeCache.txt`` file in the build
directory. The equivalent of ``make distclean`` with CMake is to just
delete the entire build directory and reconfigure.

If you need to debug a CMake build, you can see the compile and link
commands which are being issued by adding ``--verbose`` to the build
command::

  $ cmake --build build --verbose

Common CMake Options
~~~~~~~~~~~~~~~~~~~~

The default build configuration options try to be sensible for new
users and developers, but there are a few you probably want to set
manually for production runs or for debugging:

* ``CMAKE_BUILD_TYPE``: The default is ``RelWithDebInfo``, which
  builds an optimised executable with debug symbols included. Change
  this to ``Release`` to remove the debug symbols, or ``Debug`` for an
  unoptimised build, but better debug experience

* ``CHECK``: This sets the level of internal runtime checking done in
  the BOUT++ library, and ranges from 0 to 4 (inclusive). By default,
  this is 2, which aims to be a balance between useful checks and
  speed. Set this to 0 for faster production runs, and to 4 for more
  in-depth (and slower) checking.

* ``BOUT_UPDATE_GIT_SUBMODULE``: This is on by default, and ensures
  that the bundled git submodules are up-to-date. You should turn this
  off if you are using system versions, or if you run into problems
  updating the submodules.

* ``NetCDF_ROOT``: NetCDF is one of the few required, non-bundled
  dependencies. If CMake is having trouble finding netCDF, or the
  correct version, you should set this variable to the installed
  location of the netCDF C library.

* ``BOUT_BUILD_EXAMPLES``, ``BOUT_TESTS``: These two options are
  particularly useful for developers of the BOUT++ library, and for
  new users. You can turn them off to save some time configuring the
  library. By default, these are on, but the examples and tests are
  not built unless you specifically ask for them, using the targets
  ``build-all-examples`` and ``build-check`` respectively.


Downloading Dependencies
~~~~~~~~~~~~~~~~~~~~~~~~

If you don't have some dependencies installed, CMake can be used to download,
configure and compile them alongside BOUT++.

For NetCDF, use ``-DBOUT_DOWNLOAD_NETCDF_CXX4=ON``

For SUNDIALS, use ``-DBOUT_DOWNLOAD_SUNDIALS=ON``. If using ``ccmake`` this option
may not appear initially. This automatically sets ``BOUT_USE_SUNDIALS=ON``, and
configures SUNDIALS to use MPI.

For ADIOS2, use ``-DBOUT_DOWNLOAD_ADIOS=ON``. This will download and
configure `ADIOS2 <https://adios2.readthedocs.io/>`_, enabling BOUT++
to read and write this high-performance parallel file format.

Bundled Dependencies
~~~~~~~~~~~~~~~~~~~~

BOUT++ bundles some dependencies, currently `mpark.variant
<https://github.com/mpark/variant>`_, `fmt <https://fmt.dev>`_ and
`googletest <https://github.com/google/googletest>`_. If you wish to
use an existing installation of ``mpark.variant``, you can set
``-DBOUT_USE_SYSTEM_MPARK_VARIANT=ON``, and supply the installation
path using ``mpark_variant_ROOT`` via the command line or environment
variable if it is installed in a non standard loction. Similarly for
``fmt``, using ``-DBOUT_USE_SYSTEM_FMT=ON`` and ``fmt_ROOT``
respectively. To turn off both, you can set
``-DBOUT_USE_GIT_SUBMODULE=OFF``.

The recommended way to use ``googletest`` is to compile it at the same
time as your project, therefore there is no option to use an external
installation for that.


Working with an active ``conda`` environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When ``conda`` is used, it installs separate versions of several libraries. These
can cause warnings or even failures when linking BOUT++ executables. There are
several alternatives to deal with this problem:
* The simplest but least convenient option is to use ``conda deactivate`` before
  configuring, compiling, or running any BOUT++ program.
* You might sometimes want to link to the conda-installed libraries. This is
  probably not ideal for production runs on an HPC system (as conda downloads
  binary packages that will not be optimized for specific hardware), but can be
  a simple way to get packages for testing or on a personal computer. In this
  case just keep your ``conda`` environment active, and with luck the libraries
  should be picked up by the standard search mechanisms.
* In case you do want a fully optimized and as-stable-as-possible build for
  production runs, it is probably best not to depend on any conda packages for
  compiling or running BOUT++ executables (restrict ``conda`` to providing Python
  packages for post-processing, and their dependencies). Passing
  ``-DBOUT_IGNORE_CONDA_ENV=ON`` (default ``OFF``) excludes anything in the conda
  environment from CMake search paths. This should totally separate BOUT++ from
  the ``conda`` environment.

.. _sec-config-nls:

Natural Language Support
------------------------

BOUT++ has support for languages other than English, using GNU
gettext. If you are planning on installing BOUT++ (see
:ref:`sec-install-bout`) then this should work automatically, but if
you will be running BOUT++ from the directory you downloaded it into,
then configure with the option::

  cmake . -DCMAKE_INSTALL_LOCALEDIR=$PWD/locale

This will enable BOUT++ to find the translations.

See :ref:`sec-run-nls` for details of how to switch language when running
BOUT++ simulations.

.. _sec-configanalysis:

Configuring analysis routines
-----------------------------

The BOUT++ installation comes with a set of useful routines which can be
used to prepare inputs and analyse outputs. Most of this code is now in Python,
though IDL was used for many years. Python is useful In particular because the test suite
scripts and examples use Python, so to run these you’ll need python configured.

When the configure script finishes, it prints out the paths you need to
get IDL, Python, and Octave analysis routines working. If you
just want to compile BOUT++ then you can skip to the next section, but
make a note of what configure printed out.


.. _sec-config-python:

Python configuration
~~~~~~~~~~~~~~~~~~~~

To use Python, you will need the dependencies of the `boututils
<https://github.com/boutproject/boututils>`__ and `boutdata
<https://github.com/boutproject/boutdata>`__ libraries. The simplest way to get these is
to install the packages with pip::

    $ pip install --user boutdata

or conda::

    $ conda install boutdata

You can also install all the packages directly (see the documentation in the `boututils
<https://github.com/boutproject/boututils>`__ and `boutdata
<https://github.com/boutproject/boutdata>`__ repos for the most up to date list)
using pip::

    $ pip install --user numpy scipy matplotlib sympy netCDF4 future importlib-metadata

or conda::

    $ conda install numpy scipy matplotlib sympy netcdf4 future importlib-metadata

They may also be available from your Linux system's package manager. 

For example on Fedora::

    $ sudo dnf install python3-boututils python3-boutdata

To use the versions of ``boututils`` and ``boutdata`` provided by BOUT++,  the path to
``tools/pylib`` should be added to the ``PYTHONPATH`` environment variable. This is not
necessary if you have installed the ``boututils`` and ``boutdata`` packages.  Instructions
for doing this are printed at the end of the configure script, for example::

    Make sure that the tools/pylib directory is in your PYTHONPATH
    e.g. by adding to your ~/.bashrc file

       export PYTHONPATH=/home/ben/BOUT/tools/pylib/:$PYTHONPATH

To test if this command has worked, try running::

    $ python -c "import boutdata"

If this doesn’t produce any error messages then Python is configured correctly.

Note that ``boututils`` and ``boutdata`` are provided by BOUT++ as submodules, so versions
compatible with the checked out version of BOUT++ are downloaded into the
``externalpackages`` directory. These are the versions used by the tests run by ``make
check`` even if you have installed ``boututils`` and ``boutdata`` on your system.


.. _sec-config-idl:

IDL configuration
~~~~~~~~~~~~~~~~~

If you want to use `IDL <https://en.wikipedia.org/wiki/IDL_(programming_language)>`__ to analyse
BOUT++ outputs, then the ``IDL_PATH`` environment variable should include the
``tools/idllib/`` subdirectory included with BOUT++.
The required command (for Bash) is printed at the end of the BOUT++ configuration::

    $ export IDL_PATH=...

After running that command, check that ``idl`` can find the analysis routines by running::

    $ idl
    IDL> .r collect
    IDL> help, /source

You should see the function ``COLLECT`` in the ``BOUT/tools/idllib``
directory. If not, something is wrong with your ``IDL_PATH`` variable.
On some machines, modifying ``IDL_PATH`` causes problems, in which case
you can try modifying the path inside IDL by running::

    IDL> !path = !path + ":/path/to/BOUT-dev/tools/idllib"

where you should use the full path. You can get this by going to the
``tools/idllib`` directory and typing ``pwd``. Once this is done
you should be able to use ``collect`` and other routines.

.. _sec-compile-bout:

Compiling BOUT++
----------------

Once BOUT++ has been configured, you can compile the bulk of the code by
going to the ``BOUT-dev`` directory and running::

    $ cmake --build <build-directory>

where ``<build-directory>`` is the path to the build directory

At the end of this, you should see a file ``libbout++.so`` in the
``lib/`` subdirectory of the BOUT++ build directory. If you get an error,
please `create an issue on Github <https://github.com/boutproject/BOUT-dev/issues>`__
including:

-  Which machine you’re compiling on

-  The output from make, including full error message

-  The ``CMakeCache.txt`` file in the BOUT++ build directory

.. _sec-runtestsuite:

Running the test suite
----------------------

BOUT++ comes with three sets of test suites: unit tests, integrated
tests and method of manufactured solutions (MMS) tests. The easiest
way to run all of them is to simply do::

    $ cmake --build <build-directory> --target check

Alternatively, if you just want to run
one set of them individually, you can do::

    $ cmake --build <build-directory> --target check-unit-tests
    $ cmake --build <build-directory> --target check-integrated-tests
    $ cmake --build <build-directory> --target check-mms-tests

**Note:** The integrated and MMS test suites currently uses the ``mpirun``
command to launch the runs, so won’t work on machines which use a job
submission system like slurm or PBS.

These tests should all pass, but if not please `create an issue on Github <https://github.com/boutproject/BOUT-dev/issues>`__
containing:

-  Which machine you’re running on

-  The ``CMakeCache.txt`` file in the BOUT++ build directory

-  The ``run.log.*`` files in the directory of the test which failed

If the tests pass, congratulations! You have now got a working
installation of BOUT++. Unless you want to use some experimental
features of BOUT++, skip to section [sec-running] to start running the
code.

.. _sec-install-bout:

Installing BOUT++ (experimental)
--------------------------------

Most BOUT++ users install and develop their own copies in their home directory,
so do not need to install BOUT++ to a system directory.
As of version 4.1 (August 2017), it is possible to install BOUT++ but this is
not widely used and so should be considered experimental.

After configuring and compiling BOUT++ as above, BOUT++ can be installed
to system directories by running as superuser or ``sudo``::

   $ sudo cmake --build <build-directory> --target install

.. DANGER:: Do not do this unless you know what you're doing!

This will install the following files under ``/usr/local/``:

* ``/usr/local/bin/bout-config``  A script which can be used to query BOUT++ configuration and compile codes with BOUT++.

* ``/usr/local/include/bout++/...`` header files for BOUT++

* ``/usr/local/lib/libbout++.so``  The main BOUT++ library

* ``/usr/local/lib/libpvode.so`` and ``/usr/local/lib/libpvpre.so``, the PVODE library

* ``/usr/local/share/bout++/pylib/...`` Python analysis routines

* ``/usr/local/share/bout++/idllib/...`` IDL analysis routines


To install BOUT++ under a different directory, use the ``prefix=``
flag e.g. to install in your home directory::

   $ cmake --build <build-directory> --target install -DCMAKE_INSTALL_PREFIX=$HOME/local/

You can also specify this prefix when configuring, in the usual way
(see :ref:`sec-config-bout`)::

     $ cmake -S . -B <build-directory> -DCMAKE_INSTALL_PREFIX=$HOME/local/
     $ cmake --build <build-directory> -j 4
     $ cmake --build <build-directory> --target install

More control over where files are installed is possible by passing options to
``cmake``, following the GNU conventions:

* ``-DCMAKE_INSTALL_BINDIR=``  sets where ``bout-config`` will be installed ( default ``/usr/local/bin``)

* ``-DCMAKE_INSTALL_INCLUDEDIR=`` sets where the ``bout++/*.hxx`` header files wil be installed (default ``/usr/local/include``)

* ``-DCMAKE_INSTALL_LIBDIR=`` sets where the ``libbout++.so``, ``libpvode.so`` and ``libpvpre.so`` libraries are installed (default ``/usr/local/lib``)


After installing, that you can run ``bout-config`` e.g::

    $ bout-config --all

which should print out the list of configuration settings which ``bout-config`` can provide.
If this doesn't work, check that the directory containing ``bout-config`` is in your ``PATH``.

The python and IDL analysis scripts can be configured using
``bout-config`` rather than manually setting paths as in
:ref:`sec-configanalysis`. Add this line to your startup file
(e.g. ``$HOME/.bashrc``)::

   export PYTHONPATH=`bout-config --python`:$PYTHONPATH

note the back ticks around ``bout-config --python`` not
quotes. Similarly for IDL::

   export IDL_PATH=`bout-config --idl`:'<IDL_DEFAULT>':$IDL_PATH

More details on using bout-config are in the :ref:`section on makefiles <sec-bout-config>`.
