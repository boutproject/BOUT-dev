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
`Scorep <http://www.vi-hps.org/projects/score-p/>`_+`Scalasca <http://www.scalasca.org/>`_
and 
`Extrae <https://tools.bsc.es/extrae/>`_+`Paraver <https://tools.bsc.es/paraver/>`_
tool chains.
Both are suitable for analyzing code parallelized with MPI and/or OpenMP.
Scorep+Scalasca gives timings and call trees for each processor/thread,
while Extrae/Paraver produces visualizations showing what each processor/thread
is doing at a point in time.

Scorep/Scalasca profiling
-------------------------

Instrumentation
~~~~~~~~~~~~~~~

Scorep automatically reports the time spend in MPI communications and OpenMP
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

Running
~~~~~~~

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

Running
~~~~~~~

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
