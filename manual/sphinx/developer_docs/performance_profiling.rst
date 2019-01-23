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
tool chain.

Scorep/Scalasca profiling
-------------------------

Instrumentation
~~~~~~~~~~~~~~~

Scorep automatically reports the time spend in MPI communications and OpenMP
loops. However, to obtain information on the time spent in specific functions,
it is necessary to instrument the source. The tools to do this are provided in
``scorepwrapper.hxx``.

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

This warning notwithstanding, it is reasonable to expect sensibly-instrumented
code to run ~50% to 100% slower than the uninstrumented code.

Configure and build
~~~~~~~~~~~~~~~~~~~

Configure with ``--with-scorep`` to enable Scorep instrumentation, then build
as normal.  This option can be combined with other options, but it is usually
desirable to profile the optimized code, configuring with the flags
``--enable-optimize=3 --enable-checks=0``. Build the code with ``make`` as
normal.

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

.. _sec-machine-specific:

Machine-specific installation
-----------------------------

These are some configurations which have been found to work on
particular machines.

Archer
~~~~~~

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



