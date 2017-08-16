.. _sec-running:

Running BOUT++
==============

Quick start
-----------

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

-  ``BOUT.settings`` contains all the options used in the code, including
   options which were not set and used the default values. It's in the same
   format as BOUT.inp, so can be renamed and used to re-run simulations
   if needed.
   
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

To see some of the other command-line options try "-h":

.. code-block:: bash

   $ ./conduction -h

and see the section on options (:ref:`sec-options`).

Analysing the output
--------------------

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

Stopping simulations
--------------------

If you need to stop a simulation early this can be done by Ctrl-C in a terminal,
but this will stop the simulation immediately without shutting down cleanly. Most
of the time this will be fine, but interrupting a simulation while it is writing
data to file could result in inconsistent or corrupted data.

Stop file
~~~~~~~~~

**Note** This method needs to be enabled before the simulation starts by setting
``stopCheck=true`` on the command line or input options:

.. code-block:: bash

    $ mpirun -np 4 ./conduction stopCheck=true

or in the top section of ``BOUT.inp`` set ``stopCheck=true``.

At every output time, the monitor checks for the existence of a file, by default called
``BOUT.stop``, in the same directory as the output data. If the file exists then
the monitor signals the time integration solver to quit. This should result in a clean
shutdown.

To stop a simulation using this method, just create an empty file in the output directory

.. code-block:: bash

    $ mpirun -np 4 ./conduction stopCheck=true
    ...
    $ touch data/BOUT.stop

just remember to delete the file afterwards.

Manipulating restart files
--------------------------

It is sometimes useful to change the number of processors used in a simulation,
or to modify restart files in various ways. For example, a 3D turbulence
simulation might start with a quick 2D simulation with diffusive transport to reach
a steady-state. The restart files can then be extended into 3D, noise added to seed
instabilities, and the files split over a more processors.

Routines to modify restart files are in ``tools/pylib/boutdata/restart.py``:

.. code-block:: pycon

    >>> from boutdata import restart
    >>> help(restart)

Changing number of processors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To change the number of processors use the ``redistribute`` function:

.. code-block:: pycon

    >>> from boutdata import restart
    >>> restart.redistribute(32, path="../oldrun", output=".")

where in this example ``32`` is the number of processors desired; ``path`` sets
the path to the existing restart files, and ``output`` is the path where
the new restart files should go.
**Note** Make sure that ``path`` and ``output`` are different.

If your simulation is divided in X and Y directions then you should also specify
the number of processors in the X direction, ``NXPE``:

.. code-block:: pycon

    >>> restart.redistribute(32, path="../oldrun", output=".", nxpe=8)

**Note** Currently this routine doesn't check that this split is consistent with
branch cuts, e.g. for X-point tokamak simulations. If an inconsistent choice is made
then the BOUT++ restart will fail.

**Note** It is a good idea to set ``nxpe`` in the ``BOUT.inp`` file to be consistent with
what you set here. If it is inconsistent then the restart will fail, but the error message may
not be particularly enlightening.
