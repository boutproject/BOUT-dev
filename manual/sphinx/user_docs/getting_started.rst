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
