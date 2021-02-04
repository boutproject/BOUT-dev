.. _sec-options:

BOUT++ options
==============

The inputs to BOUT++ are a text file containing options, command-line options,
and for complex grids a binary grid file in NetCDF format. Generating input
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
Any type which can be read from a stream using the ``>>`` operator can
be stored in an option (see later for the implementation details).
Supported value types include:

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
    工作的 = true            # a boolean
    इनपुट = "some text"      # a string

Option names can contain almost any character except ’=’ and ’:’,
including unicode.  If they start with a number or ``.``, contain
arithmetic symbols (``+-*/^``), brackets (``(){}[]``), equality
(``=``), whitespace or comma ``,``, then these will need to be escaped
in expressions. See below for how this is done.

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
are first read, and then processed afterwards on demand.
The value ``pi`` is already defined, as is ``π``, and can be used in expressions.

Uses for expressions include initialising variables
:ref:`sec-expressions` and input sources, defining grids
:ref:`sec-gridgen` and MMS convergence tests :ref:`sec-mms`.

Expressions can include addition (``+``), subtraction (``-``),
multiplication (``*``), division (``/``) and exponentiation (``^``)
operators, with the usual precedence rules. In addition to ``π``,
expressions can use predefined variables ``x``, ``y``, ``z`` and ``t``
to refer to the spatial and time coordinates.
A number of functions are defined, listed in table
:numref:`tab-initexprfunc`. One slightly unusual feature (borrowed from `Julia <https://julialang.org/>`_)
is that if a number comes before a symbol or an opening bracket (``(``)
then a multiplication is assumed: ``2x+3y^2`` is the same as
``2*x + 3*y^2``, which with the usual precedence rules is the same as
``(2*x) + (3*(y^2))``. 

Expressions can span more than one line, which can make long expressions
easier to read:

.. code-block:: cfg

   pressure = temperature * ( density0 +
                              density1 )
   temperature = 12
   density0 = 3
   density1 = 1

The convention is the same as in `Python <https://www.python.org/>`_:
If brackets are not balanced (closed) then the expression continues on the next line.

All expressions are calculated in floating point and then converted to
an integer if needed when read inside BOUT++. The conversion is done by rounding
to the nearest integer, but throws an error if the floating point
value is not within :math:`1e-3` of an integer. This is to minimise
unexpected behaviour. If you want to round any result to an integer,
use the ``round`` function:

.. code-block:: cfg

    bad_integer = 256.4
    ok_integer = round(256.4)

Note that it is still possible to read ``bad_integer`` as a real
number, since the type is determined by how it is used.

Have a look through the examples to see how the options are used.

Special symbols in Option names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If option names start with numbers or ``.`` or contain symbols such as
``+`` and ``-`` then these symbols need to be escaped in expressions
or they will be treated as arithmetic operators like addition or
subtraction. To escape a single character 
``\`` (backslash) can be used, for example ``plasma\-density * 10``
would read the option ``plasma-density`` and multiply it
by 10 e.g

.. code-block:: cfg

    plasma-density = 1e19
    2ndvalue = 10
    value = plasma\-density * \2ndvalue

To escape multiple characters, ` (backquote) can be used:

.. code-block:: cfg

    plasma-density = 1e19
    2ndvalue = 10
    value = `plasma-density` * `2ndvalue`

The character ``:`` cannot be part of an option or section name, and cannot be escaped,
as it is always used to separate sections.

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
    nx = 20  # Number of points in X
    ny = 16  # Number of points in Y
    nz = 32  # Number of points in Z

Due to historical reasons, ``nx`` is defined differently to ``ny`` and ``nz``:

- ``nx`` is the number of points in X **including** the boundaries
- ``ny`` and ``nz`` are the number of points in Y and Z **not including** the
  boundaries

The default number of boundary points in X is 2, so taking into account the
boundary at each end of the domain, ``nx`` usually means "the number of interior
grid points in X plus four". In the example above, both X and Y have 16 interior
grid points.

It is recommended, but not necessary, that this be :math:`\texttt{nz}
= 2^n`, that is :math:`1,2,4,8,\ldots`. This is because FFTs are usually
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

Alternatively, the number in the Y direction can be specified (if both are
given, ``NXPE`` takes precedence and ``NYPE`` is ignored):

.. code-block:: cfg

    NYPE = 1  # Set number of Y processors

When choosing ``NXPE`` or ``NYPE``, they must also obey some constraints:

- ``NXPE`` must be a factor of the number of grid points in the x-direction

  - That is, ``(nx - 4) / NXPE`` must be an integer, assuming the usual two
    boundary points

- ``NYPE`` must be a factor of the number of grid points in the y-direction

  - That is, ``ny / NYPE`` must be an integer

- For more general topologies, the number of points per processor ``ny / NYPE``
  must also be a factor of the number of points in each region. For example, in
  the usual tokamak topologies:

  - in single-null there are two divertor leg and one core regions
  - in double-null there are four divertor leg, one inner core and one outer
    core regions

Please note that here "core" means "core and adjacent SOL". See
:ref:`sec-bout-topology` for a more detailed explanation of these regions.

When BOUT++ automatically chooses ``NXPE`` and ``NYPE`` it finds all valid pairs
which give ``total number of processors == NPES = NXPE * NYPE`` and also satisfy
the constraints above. It then chooses the pair that makes the grid on each
processor as close to square as possible (technically it chooses the pair that
minimises ``abs(sqrt(NPES * (nx - 4) / ny) - NXPE)``).

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

-  ``fourth``, method for fourth derivatives

-  ``upwind``, method for upwinding terms

-  ``flux``, for conservation law terms

The methods which can be specified include U1, U4, C2, C4, W2, W3, FFT Apart
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

The output (dump) files with time-history are controlled by settings
in a section called “output”. Restart files contain a single
time-slice, and are controlled by a section called “restart”. The
options available are listed in table :numref:`tab-outputopts`.

.. _tab-outputopts:
.. table:: Output file options
	   
   +-------------+----------------------------------------------------+--------------+
   | Option      | Description                                        | Default      |
   |             |                                                    | value        |
   +-------------+----------------------------------------------------+--------------+
   | enabled     | Writing is enabled                                 | true         |
   +-------------+----------------------------------------------------+--------------+
   | floats      | Write floats rather than doubles                   | false        |
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

**enabled** is useful mainly for doing performance or scaling tests, where you
want to exclude I/O from the timings. **floats** can be used to reduce the size
of the output files: files are stored as double by default, but setting
**floats = true** changes the output to single-precision floats.

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
``bout++/include/options.hxx``. 

To access the options, there is a static function (singleton)::

    auto& options = Options::root();

which returns a reference (type ``Options&``). Note that without
the ``&`` the options tree will be copied, so any changes made will not
be retained in the global tree. Options can be set by
assigning, treating options as a map or dictionary::

    options["nout"] = 10;    // Integer
    options["restart"] = true;  // bool
    
Internally these values are stored in a variant type, which supports commonly
used types including strings, integers, real numbers and fields (2D and
3D). Since strings can be stored, any type can be assigned, so long as it can be
streamed to a string (using ``<<`` operator and a ``std::stringstream``).

Often it’s useful to see where an option setting has come from e.g. the
name of the options file or “command line”. To specify a source, use
the ``assign`` function to assign values::

    options["nout"].assign(10, "manual");

A value cannot be assigned more than once with different values and
the same source ("manual" in this example). This is to catch a common
error in which a setting is inconsistently specified in an input
file. To force a value to change, overwriting the existing value (if
any)::

    options["nout"].force(20, "manual");

Sub-sections are created as they are accessed, so a value in a
sub-section could be set using::

    auto& section = options["mysection"];
    section["myswitch"] = true;

or just::

    options["mysection"]["myswitch"] = true;

Names including sections, subsections, etc. can be specified using ``":"`` as a
separator, e.g.::
    options["mysection:mysubsection:myswitch"] = true;

To get options, they can be assigned to a variable::

    int nout = options["nout"];

If the option is not found then a ``BoutException`` will be thrown. A
default value can be given, which will be used if the option has not
been set::

    int nout = options["nout"].withDefault(1);

If ``options`` is not ``const``, then the given default value will be
cached. If a default value has already been cached for this option,
then the default values must be consistent: A ``BoutException`` is
thrown if inconsistent default values are detected.

The default can also be set from another option. This may be useful if two or
more options should usually be changed together::

    BoutReal value2 = options["value2"].withDefault(options["value1"]);

Note that if the result should be a real number (e.g. ``BoutReal``) then ``withDefault``
should be given a real. Otherwise it will convert the number to an integer::

  BoutReal value = options["value"].withDefault(42);  // Convert to integer

  BoutReal value = options["value"].withDefault(42.0); // ok

  auto value = options["value"].withDefault<BoutReal>(42); // ok

It is common for BOUT++ models to read in many settings which have the
same variable name as option setting (e.g. "nout" here). A convenient
macro reads options into an already-defined variable::

    int nout;
    OPTION(options, nout, 1);

where the first argument is a section, second argument is the variable
whose name will also be used as the option string, and third argument
is the default value.

Every time an option is accessed, a message is written
to ``output_info``. This message includes the value used and the
source of that value. By default this message is printed to the
terminal and saved in the log files, but this can be disabled by
changing the logging level: Add ``-q`` to the command line to reduce
logging level. See section :ref:`sec-logging` for more details about
logging.

The type to be returned can also be specified as a template argument::

    BoutReal nout = options["nout"].as<BoutReal>();

Any type can be used which can be streamed (operator ``>>``) from a
``stringstream``. There are special implementations for ``bool``,
``int`` and ``BoutReal`` which enable use of expressions in the input
file. The type can also be specified to ``withDefault``, or will be
inferred from the argument::

    BoutReal nout = options["nout"].withDefault<BoutReal>(1);

Documentation
~~~~~~~~~~~~~

Options can be given a ``doc`` attribute describing what they do. This documentation
will then be written to the ``BOUT.settings`` file at the end of a run::

  Te0 = options["Te0"].doc("Temperature in eV").withDefault(30.0);

The ``.doc()`` function returns a reference ``Options&`` so can be chained
with ``withDefault`` or ``as`` functions, or as part of an assignment::

  options["value"].doc("Useful setting info") = 42;

This string is stored in the attributes of the option::

  std::string docstring = options["value"].attributes["doc"];

Creating Options
~~~~~~~~~~~~~~~~

Options and subsections can be created by setting values, creating subsections as needed::

  Options options;
  options["value1"] = 42;
  options["subsection1"]["value2"] = "some string";
  options["subsection1"]["value3"] = 3.1415;

or using an initializer list::

  Options options {{"value1", 42},
                   {"subsection1", {{"value2", "some string"},
                                    {"value3", 3.1415}}}};

These are equivalent, but the initializer list method makes the tree structure clearer.
Note that the list can contain many of the types which ``Options`` can hold, including
``Field2D`` and ``Field3D`` objects.

Overriding library defaults
~~~~~~~~~~~~~~~~~~~~~~~~~~~

BOUT++ sets defaults for options controlling the mesh, etc. A physics model (or
other user code) can override these defaults by using the convenience macro
BOUT_OVERRIDE_DEFAULT_OPTION, for example if you want to change the default
value of ``mesh::staggergrids`` from false to true, put (outside any
class/function body)::

    BOUT_OVERRIDE_DEFAULT_OPTION("mesh:staggergrids", true);

Older interface
~~~~~~~~~~~~~~~

Some code in BOUT++ currently uses an older interface to ``Options``
which uses pointers rather than references. Both interfaces are
currently supported, but use of the newer interface above is
encouraged.

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

Reading and writing to NetCDF
-----------------------------

If NetCDF4 support is enabled, then the ``OptionsNetCDF`` class
provides an experimental way to read and write options. To use this class::

  #include "options_netcdf.hxx"
  using bout::experimental::OptionsNetCDF;

Examples are in integrated test ``tests/integrated/test-options-netcdf/``

To write the current ``Options`` tree (e.g. from ``BOUT.inp``) to a
NetCDF file::

  OptionsNetCDF("settings.nc").write(Options::root());

and to read it in again::

  Options data = OptionsNetCDF("settings.nc").read();

Fields can also be stored and written::

  Options fields;
  fields["f2d"] = Field2D(1.0);
  fields["f3d"] = Field3D(2.0);
  OptionsNetCDF("fields.nc").write(fields);

This should allow the input settings and evolving variables to be
combined into a single tree (see above on joining trees) and written
to the output dump or restart files.

Reading fields is a bit more difficult. Currently 1D data is read as
an ``Array<BoutReal>``, 2D as ``Matrix<BoutReal>`` and 3D as
``Tensor<BoutReal>``. These can be extracted directly from the
``Options`` tree, or converted to a Field::

  Options fields_in = OptionsNetCDF("fields.nc").read();
  Field2D f2d = fields_in["f2d"].as<Field2D>();
  Field3D f3d = fields_in["f3d"].as<Field3D>();

Note that by default reading as ``Field2D`` or ``Field3D`` will use the global
``bout::globals::mesh``. To use a different mesh, or different cell location,
pass a field which the result should be similar to::

  Field3D example = ... // Some existing field
  
  Field3D f3d = fields_in["f3d"].as<Field3D>(example);

Meta data like ``Mesh`` pointer, will be taken from ``example``. 

Currently converting from ``Matrix`` or ``Tensor`` types only works if
the data in the ``Matrix`` or ``Tensor`` is the same size as the
``Field``. In the case of grid files, the fields only needs a part of
the global values. Some kind of mapping from the global index to local
index is needed, probably defined by ``Mesh``. For now it should be
possible to be compatible with the current system, so that all
quantities from the grid file are accessed through Mesh::get.

Time dependence
~~~~~~~~~~~~~~~

When writing NetCDF files, some variables should have a time
dimension added, and then be added to each time they are written. This
has been implemented using an attribute: If variables in the ``Options``
tree have an attribute "time_dimension" then that is used as the name
of the time dimension in the output file. This allows multiple time
dimensions e.g. high frequency diagnostics and low frequency outputs,
to exist in the same file::

  Options data;
  data["scalar"] = 1.0;
  data["scalar"].attributes["time_dimension"] = "t";
  
  data["field"] = Field3D(2.0);
  data["field"].attributes["time_dimension"] = "t";
  
  OptionsNetCDF("time.nc").write(data);
  
  // Update time-dependent values. This can be done without `force` if the time_dimension
  // attribute is set
  data["scalar"] = 2.0;
  data["field"] = Field3D(3.0);
  
  // Append data to file
  OptionsNetCDF("time.nc", OptionsNetCDF::FileMode::append).write(data);

Some issues:

* Currently all variables in the Options tree are written when passed
  to ``OptionsNetCDF::write``. This means that the variables with
  different time dimensions should be stored in different Options
  trees, so they can be written at different times. One possibility is
  to have an optional argument to write, so that only variables with
  one specified time dimension are updated.


FFT
---

There is one option for Fourier transforms, ``fft_measurement_flag`` (default:
``estimate``). This can be used to control FFTW's measurement mode:
``estimate`` for ``FFTW_ESTIMATE``, ``measure`` for ``FFTW_MEASURE`` or
``exhaustive`` for ``FFTW_EXHAUSTIVE``:

.. code-block:: cfg

    [fft]
    fft_measurement_flag = measure

In ``FFTW_MEASURE`` mode, FFTW runs and measures how long several FFTs take,
and tries to find the optimal method; ``FFTW_EXHAUSTIVE`` tests even more
algorithms.

.. note:: Technically, ``FFTW_MEASURE`` and ``FFTW_EXHAUSTIVE`` are
          non-deterministic and enabling ``fft_measure`` may result in slightly
          different answers from run to run, or be dependent on the number of
          MPI processes. This may be important if you are trying to benchmark
          or measure performance of your code.

          See the `FFTW FAQ`_ for more information.


.. _FFTW FAQ: http://www.fftw.org/faq/section3.html#nondeterministic


Types for multi-valued options
------------------------------

An ``enum class`` can be a useful construct for options in a physics model. It
can have an arbitrary number of user-defined, named values (although the code
in ``include/bout/bout_enum_class.hxx`` needs extending for more than 10
values). The advantage over using a ``std::string`` for an option is that a
typo cannot produce an unexpected value: in C++ code it is a compile-time error
and reading from ``BOUT.inp`` it is a run-time exception. We provide a utility
macro ``BOUT_ENUM_CLASS`` to define an ``enum class`` with some extra
convenience methods. For example, after defining ``myoption`` like::

    BOUT_ENUM_TYPE(myoption, foo, bar, baz);

it is possible not only to test for a value, e.g.::

    myoption x = <something>;
    ...
    if (x == myoption::foo) {
      do a foo thing
    }

but also to convert the option to a string::

    std::string s = toString(x);

pass it to a stream::

    output << x;

or get an option like ``myinput=baz`` from an input file or the command line as
a ``myoption``::

    myoption y = Options::root()["myinput"].as<myoption>();

or with a default value::

    myoption y = Options::root()["myinput"].withDefault(myoption::bar);

Only strings exactly (but case-insensitively) matching the name of one of the
defined ``myoption`` values are allowed, anything else results in an exception
being thrown.
