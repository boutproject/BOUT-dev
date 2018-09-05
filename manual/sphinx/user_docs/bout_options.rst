.. _sec-options:

BOUT++ options
==============

The inputs to BOUT++ are a text file containing options, command-line options,
and for complex grids a binary grid file in NetCDF or HDF5 format. Generating input
grids for tokamaks is described in :ref:`sec-gridgen`. The grid file
describes the size and topology of the X-Y domain, metric tensor
components and usually some initial profiles. The option file specifies
the size of the domain in the symmetric direction (Z), and controls how
the equations are evolved e.g. differencing schemes to use, and boundary
conditions. In most situations, the grid file will be used in many
different simulations, but the options may be changed frequently.

All options used in a simulation are saved to a ``BOUT.settings`` file.
This includes values which are not explicitly set in ``BOUT.inp``.

BOUT.inp input file
-------------------

The text input file ``BOUT.inp`` is always in a subdirectory called
``data`` for all examples. The files include comments (starting with
either ``;`` or ``#``) and should be fairly self-explanatory. The format is
the same as a windows INI file, consisting of ``name = value`` pairs.
Supported value types are:

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

Numerical quantities can be plain numbers or expressions:

.. code-block:: cfg

   short_pi = 3.145
   foo = 6 * 9

Variables can even reference other variables:

.. code-block:: cfg

   pressure = temperature * density
   temperature = 12
   density = 3

Note that variables can be used before their definition; all variables
are first read, and then processed afterwards.

All expressions are calculated in floating point and then converted to
an integer when read inside BOUT++. The conversion is done by rounding
to the nearest integer, but throws an error if the floating point
value is not within :math:`1e-3` of an integer. This is to minimise
unexpected behaviour. If you want to round any result to an integer,
use the ``round`` function:

.. code-block:: cfg

    bad_integer = 256.4
    ok_integer = round(256.4)

Note that it is still possible to read ``bad_integer`` as a real
number though.

Have a look through the examples to see how the options are used.

Command line options
--------------------

Command-line switches are:

==============  ============================================================
   Switch               Description
==============  ============================================================
-h, --help      Prints a help message and quits
-v, --verbose   Outputs more messages to BOUT.log files
-q, --quiet     Outputs fewer messages to log files
-d <directory>  Look in <directory> for input/output files (default "data")
-f <file>       Use OPTIONS given in <file>
-o <file>       Save used OPTIONS given to <file> (default BOUT.settings)
==============  ============================================================

In addition all options in the BOUT.inp file can be set on the command line,
and will override those set in BOUT.inp. The most commonly used are “restart” and “append”,
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

.. code-block:: cfg

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

.. code-block:: cfg

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

.. code-block:: cfg

    archive = 20

saves a copy of the restart files every 20 timesteps, which can then be
used as a starting point.

.. _sec-grid-options:

Grids
~~~~~~~~~

You can set the size of the computational grid in the ``mesh`` section
of the input file (see :ref:`sec-gridgen` for more information):

.. code-block:: cfg

    [mesh]
    nx = 16  # Number of points in X
    ny = 16  # Number of points in Y
    nz = 32  # Number of points in Z

It is recommended, but not necessary, that this be :math:`\texttt{nz}
= 2^n`, i.e.  :math:`1,2,4,8,\ldots`. This is because FFTs are usually
slightly faster with power-of-two length arrays, and FFTs are used
quite frequently in many models.

.. note:: In previous versions of BOUT++, ``nz`` was constrained to be
          a power-of-two, and had to be specified as a power-of-two
          plus one (i.e. a number of the form :math:`2^n + 1` like
          :math:`2, 3, 5, 9,\ldots`) in order to account for an
          additional, unused, point in Z. Both of these conditions
          were relaxed in BOUT++ 4.0. If you use an input file from a
          previous version, check that this superfluous point is not
          included in ``nz``.

Since the Z dimension is periodic, the domain size is specified as
multiples or fractions of :math:`2\pi`. To specify a fraction of
:math:`2\pi`, use

.. code-block:: cfg

    ZPERIOD = 10

This specifies a Z range from :math:`0` to
:math:`2\pi / {\texttt{ZPERIOD}}`, and is useful for simulation of
tokamaks to make sure that the domain is an integer fraction of a torus.
If instead you want to specify the Z range directly (for example if Z is
not an angle), there are the options

.. code-block:: cfg

    ZMIN = 0.0
    ZMAX = 0.1

which specify the range in multiples of :math:`2\pi`.

In BOUT++, grids can be split between processors in both X and Y
directions. By default BOUT++ automatically divides the grid in both X and Y,
finding the decomposition with domains closest to square, whilst satisfying
constraints. These constraints are:

- Every processor must have the same size and shape domain

- Branch cuts, mostly at X-points, must be on processor boundaries.
  This is because the connection between grid points is modified in BOUT++
  by changing which processors communicate.

To specify a splitting manually, the number of processors in the X
direction can be specified:

.. code-block:: cfg

    NXPE = 1  # Set number of X processors

If you need to specify complex input values, e.g. numerical values
from experiment, you may want to use a grid file. The grid file to use
is specified relative to the root directory where the simulation is
run (i.e. running “``ls ./data/BOUT.inp``” gives the options
file). You can use the global option ``grid``, or ``mesh:file``:

.. code-block:: cfg

    grid = "data/cbm18_8_y064_x260.nc"

    # Alternatively:
    [mesh]
    file = "data/cbm18_8_y064_x260.nc"


Communications
--------------

The communication system has a section ``[comms]``, with a true/false
option ``async``. This determines whether asynchronous MPI sends are
used; which method is faster varies (though not by much) with machine
and problem.

.. _sec-diffmethodoptions:

Differencing methods
--------------------

Differencing methods are specified in the section (``[mesh:ddx]``,
``[mesh:ddy]``, ``[mesh:ddz]`` and ``[mesh:diff]``), one for each
dimension. The ``[mesh:diff]`` section is only used if the section for
the dimension does not contain an option for the differencing method.
Note that ``[mesh]`` is the name of the section passed to the mesh
constructor, which is most often ``mesh`` - but could have another
name, e.g. if multiple meshes are used.

-  ``first``, the method used for first derivatives

-  ``second``, method for second derivatives

-  ``upwind``, method for upwinding terms

-  ``flux``, for conservation law terms

The methods which can be specified are U1, U4, C2, C4, W2, W3, FFT Apart
from FFT, the first letter gives the type of method (U = upwind, C =
central, W = WENO), and the number gives the order.

The staggered derivatives can be specified as ``FirstStag`` or if the
value is not set, then ``First`` is checked.
Note that for the staggered quantities, if the staggered quantity in a
dimension is not set, first the staggered quantity in the ``[mesh:diff]``
section is checked. This is useful, as the staggered quantities are
more restricted in the available choices than the non-staggered
differenciating operators.

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
“restart”. The options available are listed in table :numref:`tab-outputopts`.

.. _tab-outputopts:
.. table:: Output file options
	   
   +-------------+----------------------------------------------------+--------------+
   | Option      | Description                                        | Default      |
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

|

**enabled** is useful mainly for doing performance or scaling tests,
where you want to exclude I/O from the timings. **floats** is used to
reduce the size of the output files: restart files are stored as double
by default (since these will be used to restart a simulation), but
output dump files are set to floats by default.

To enable parallel I/O for either output or restart files, set

.. code-block:: cfg

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
tree structure there is the `Options` class defined in
``bout++/include/options.hxx``::

    class Options {
     public:
      // Setting options
      void set(const string &key, int value, const string &source="");
      ...
      // Testing if set
      bool isSet(const string &key);
      // Getting options
      void get(const string &key,int &val,const int default);
      ...
      // Get a subsection. Creates if doesn't exist
      Options* getSection(const string &name);
    };

To access the options, there is a static function (singleton)::

      Options *options = Options::getRoot();

which gives the top-level (root) options class. Setting options is done
using the ``set()`` methods which are currently defined for ``int``,
``BoutReal``, ``bool`` and ``string`` . For example::

      options->set("nout", 10);      // Set an integer
      options->set("restart", true); // A bool

Often it’s useful to see where an option setting has come from e.g. the
name of the options file or “command line”. To specify a source, pass it
as a third argument::

      options->set("nout", 10, "manual");

To create a section, just use ``getSection`` : if it doesn’t exist it
will be created::

      Options *section = options->getSection("mysection");
      section->set("myswitch", true);

To get options, use the ``get()`` method which take the name of the
option, the variable to set, and the default value::

      int nout;
      options->get("nout", nout, 1);

Internally, `Options` converts all types to strings and does type
conversion when needed, so the following code would work::

      Options *options = Options::getRoot();
      options->set("test", "123");
      int val;
      options->get("test", val, 1);

This is because often the type of the option is not known at the time
when it’s set, but only when it’s requested.


If the verbose flag is set (``-v`` on command line) , then ``get`` methods output a
message to the log files giving the value used and the source of that value.
This is controlled by the ``log`` member of Options and can be turned on and off
by calling ``setLogging`` e.g.::

   options->getRoot()->setLogging(false); // Turn off logging of options

Changes to logging propagate to all sub-sections, so setting the logging
of the root options object sets it for all sections. To see logs from
a particular subset of the options tree set the logging for that section::

   options->getRoot()->getSection("mesh")->setLogging(true); // Turn on some logging

so the above code turns on logging of options for the "mesh" section and all subsections of mesh.

Note: This logging only affects messages printed to screen and to ``BOUT.log`` files.
It does not affect the ``BOUT.settings`` file, which should always contain a full list
of options used and their values.


Reading options
---------------

To allow different input file formats, each file parser implements the
`OptionParser` interface defined in
``bout++/src/sys/options/optionparser.hxx``::

    class OptionParser {
     public:
      virtual void read(Options *options, const string &filename) = 0;
     private:
    };

and so just needs to implement a single function which reads a given
file name and inserts the options into the given `Options` object.

To use these parsers and read in a file, there is the `OptionsReader`
class defined in ``bout++/include/optionsreader.hxx``::

    class OptionsReader {
     public:
     void read(Options *options, const char *file, ...);
     void parseCommandLine(Options *options, int argc, char **argv);
    };

This is a singleton object which is accessed using::

      OptionsReader *reader = OptionsReader::getInstance();

so to read a file ``BOUT.inp`` in a directory given in a variable
``data_dir`` the following code is used in ``bout++.cxx``::

      Options *options = Options::getRoot();
      OptionsReader *reader = OptionsReader::getInstance();
      reader->read(options, "%s/BOUT.inp", data_dir);

To parse command line arguments as options, the `OptionsReader` class
has a method::

      reader->parseCommandLine(options, argc, argv);

This is currently quite rudimentary and needs improving.


FFT
---

There is one global option for Fourier transforms, ``fft_measure``
(default: ``false``). Setting this to true enables the
``FFTW_MEASURE`` mode when performing FFTs, otherwise
``FFTW_ESTIMATE`` is used:

.. code-block:: cfg

    [fft]
    fft_measure = true

In ``FFTW_MEASURE`` mode, FFTW runs and measures how long several
FFTs take, and tries to find the optimal method.

.. note:: Technically, ``FFTW_MEASURE`` is non-deterministic and
          enabling ``fft_measure`` may result in slightly different
          answers from run to run, or be dependent on the number of
          MPI processes. This may be important if you are trying to
          benchmark or measure performance of your code.

          See the `FFTW FAQ`_ for more information.


.. _FFTW FAQ: http://www.fftw.org/faq/section3.html#nondeterministic
