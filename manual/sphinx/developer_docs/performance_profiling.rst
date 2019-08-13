.. Use bash as the default language for syntax highlighting in this file
.. highlight:: console

.. _sec-performanceprofiling:

Performance profiling
=====================

Analyzing code behaviour is vital for getting the best performance from BOUT++.
This is done by profiling the code, that is, building and running the code 
using tools that report the amount of time each processor spends in functions,
on communications, etc.

This section describes how to compile and run BOUT++ using the 
`Scorep <http://www.vi-hps.org/projects/score-p/>`_/`Scalasca <http://www.scalasca.org/>`_
and 
`Extrae <https://tools.bsc.es/extrae/>`_/`Paraver <https://tools.bsc.es/paraver/>`_
tool chains.
Both are suitable for analyzing code parallelized with MPI and/or OpenMP.
Scorep+Scalasca gives timings and call trees for each processor/thread,
while Extrae/Paraver produces visualizations showing what each processor/thread
is doing at a point in time.

Scorep/Scalasca profiling
-------------------------

Instrumentation
~~~~~~~~~~~~~~~

Scorep automatically reports the time spent in MPI communications and OpenMP
loops. However, to obtain information on the time spent in specific functions,
it is necessary to instrument the source code. The macros to do this are 
provided in ``scorepwrapper.hxx``.

To include a function in Scorep's timing, include the scorep wrapper in the 
source code

.. code-block:: c++

    #include <bout/scorepwrapper.hxx>

and then write the macro ``SCOREP0()`` at the top of the function, e.g.

.. code-block:: c++

    int Field::getNx() const{
      SCOREP0();
      return getMesh()->LocalNx;
    };

Regions of a function can also be timed by enclosing the region in braces and using the
``BOUT_SCOREP_REGION`` macro. For example,

.. code-block:: c++

    void Field2D::applyBoundary(BoutReal time) {
      SCOREP0();

      checkData(*this);

      {
      BOUT_SCOREP_REGION("display name");
        for (const auto& bndry : bndry_op) {
          bndry->apply(*this, time);
        }
      }
    };

Here, the ``SCOREP0`` macro ensures the whole ``applyBoundary`` function is timed. In
addition, the for loop is also timed and appears in the Scalasca profile as a region
inside ``applyBoundary`` with the name "display name". Any number of Scorep user regions
can be used in a function; user regions can also be nested.

**Caution** Instrumenting a function makes it execute more slowly. This can
result in misleading profiling information, particularly if 
fast-but-frequently-called functions are instrumented. Try to instrument 
significant functions only.

The profiling overhead in sensibly-instrumented code should be only a few
percent of runtime.

Configure and build
~~~~~~~~~~~~~~~~~~~

Configure with ``--with-scorep`` to enable Scorep instrumentation, then build
as normal.  This option can be combined with other options, but it is usually
desirable to profile the optimized code, configuring with the flags
``--enable-optimize=3 --enable-checks=0``. Build the code with ``make`` as
normal.

With CMake:

.. code-block:: bash

    $ SCOREP_WRAPPER=off cmake \
      -DCMAKE_C_COMPILER=scorep-mpicc \
      -DCMAKE_CXX_COMPILER=scorep-mpicxx \
      <other CMake options>

This will turn off the instrumentation during the configure
step. Please be aware that if you change ``CMakeLists.txt``, CMake
will try to automatically reconfigure the build, which the Score-P
wrappers interfere with. In this case you will need to restart the
configure step from scratch (i.e. remove the build directory and start
again).

Run and analysis
~~~~~~~~~~~~~~~~

When running the code, prepend the run command with ``scalasca -analyze``, e.g.

.. code-block:: bash

    $ scalasca -analyze mpirun -np 2 elm_pb

The run then produces an "archive" containing profiling data in a directory
called ``scorep_<exec_name>_<proc_info>_sum``.  To view the profiling 
information with the cube viewer, do

.. code-block:: bash

    $ cube scorep_<exec_name>_<proc_info>_sum/profile.cubex

Note that Scorep does not run if doing so would produce an archive with the 
same name as an existing archive. Therefore to rerun an executable on the same
number of processors, it is necessary to move or delete the first archive.

Machine-specific installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are some configurations which have been found to work on
particular machines.

Archer
^^^^^^

As of 23rd January 2019, the following configuration should work

.. code-block:: bash

    $ module swap PrgEnv-cray PrgEnv-gnu
    $ module load fftw
    $ module load archer-netcdf/4.1.3
    $ module load scalasca

Note that due to a bug in the ``CC`` compiler, it is necessary to modify 
``make.config`` after configuration if profiling OpenMP-parallelized code:

* add the flag ``-fopenmp`` to ``BOUT_FLAGS``
* add the flag ``--thread=omp:ancestry`` as an argument to ``scorep`` in ``CXX`` 


Extrae/Paraver profiling
------------------------

`Extrae <https://tools.bsc.es/extrae/>`_ is a powerful tool allowing visualization
of commumication and computation in parallel codes. It requires minimal 
instrumentation; however the trace files produced can be extremely large. 

Instrumentation, configure and build
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

No changes to the code are necessary. On some systems, environment variables
must be set before building.  Otherwise, compile and build as normal.

Run
~~~

To run, add a trace script into the normal run command, so that for example

.. code-block:: bash

    $ aprun -n 16 blob2d -d delta_1

becomes

.. code-block:: bash

    $ aprun -n 16 ./trace.sh blob2d -d delta_1

where ``trace.sh`` is the script file

.. code-block:: bash

    #!/bin/bash

    export EXTRAE_CONFIG_FILE=./extrae.xml
    export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitrace.so

    $*

The run directory must also contain the file ``extrae.xml``, which configures
which data Extrae collects. Example ``extrae.xml`` files may be found in
``${EXTRAE_HOME}/share/example/*/extrae.xml``

Running produces a file called ``TRACE.mpits``. To generate the ``.prv`` trace
file that can be read by Paraver, do

.. code-block:: bash

    TRACE_NAME=bout.prv
    ${EXTRAE_HOME}/bin/mpi2prv -f ${EXTRAE_WORK_DIR}/TRACE.mpits -o ${TRACE_NAME}

Analysis
~~~~~~~~

Open the trace file in `Paraver <https://tools.bsc.es/paraver/>`_ with

.. code-block:: bash

    $ wxparaver ${TRACE_NAME}

To view time traces, go to ``File -> Load Configuration``.  There are many
configurations to choose from!  Two useful configurations are:

* ``mpi/views/MPI_call.cfg`` to show when MPI calls are made
* ``General/views/useful_duration.cfg`` to show continuous bursts of computation

Reducing trace file size
^^^^^^^^^^^^^^^^^^^^^^^^

When trace files are very large, Paraver will prompt the user to filter or cut
the file to reduce its size.
Filtering removes some information from the trace, making it small enough to 
open and allow the user to select a region of interest.
Cutting crops the trace to a region of interest.
Both operations create new trace files, and never overwrite the original trace.

The following prescription should work for manipulating large trace files:

1. Open the large trace file in Paraver and click 'Yes' to filter it
2. Click on the tick box 'Filter'
3. Filter the trace file:
        a) select box for Events
        b) select box for Communications
        c) in 'Keep States' select box for 'Running'
        d) in 'Keep States' select box for 'IO'
        e) select a min duration of 1000
        f) click 'Apply' 
4. View 'useful duration' configuration and locate the region of interest
5. Zoom into the region of interest, and start and end the zoom on equivalent
   large sections of computation (blue/green) 
6. Right click -> Run -> Cutter
7. Change the 'Input' trace file to cut from the filtered to the original one.
8. Click cut.

This produces a trace file which has all the original profiling information, 
but is much smaller as it is limited in time to a region of interest.

Machine-specific installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are some configurations which have been found to work on
particular machines.

Archer
^^^^^^

As of 1st February 2019, the following configuration should work

.. code-block:: bash

    $ module swap PrgEnv-cray PrgEnv-gnu
    $ module load fftw
    $ module load archer-netcdf/4.1.3
    $ module load papi
    $ module load bsctools/extrae
    $
    $ export CRAYPE_LINK_TYPE=dynamic

Note that due to a bug in the ``CC`` compiler, it is necessary to modify 
``make.config`` after configuration to add the flag  ``-fopenmp`` to 
``BOUT_FLAGS``, when profiling OpenMP-parallelized code.
