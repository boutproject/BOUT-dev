.. _sec-debugging:

==================
 Debugging Models
==================

When developing a new physics model, or using an existing one in new
regimes, it's expected that things will occasionally go wrong and
you'll need to debug the program. While debuggers like ``gdb`` are
very powerful, using them in parallel can be difficult unless you have
access to a dedicated parallel debugger. BOUT++ has some utilities
that make it easier to debug issues that only arise in parallel and/or
long-running simulations.

Loggers
=======

The first of these is the standard "write to screen" with the
``output.write`` :ref:`family of logging functions <sec-logging>`. If
you have a bug which is easily reproducible and occurs almost
immediately every time you run the code, then this is probably the
easiest way to hunt it down.

The main downside of (most of) these loggers is that if you have a lot of
output they will slow down simulations. Even if you use the
``--quiet`` command line option to turn them off, they will still add
some overhead. The `output_debug` logger can be disabled entirely at
compile-time (so there will be no overhead at all), which means it's
well suited to adding in-depth diagnostic or debugging information
that can be kept permanently in the code and only enabled if needed.

To enable the ``output_debug`` messages, configure BOUT++ with a
``CHECK`` level ``>= 3``. To enable it at lower check levels,
configure BOUT++ with ``-DENABLE_OUTPUT_DEBUG``. When running BOUT++
add a ``-v -v`` flag to see ``output_debug`` messages.

Backtrace
=========

BOUT++ can also automatically print a backtrace in the event of a crash. This is
very useful to include if you ever need to report a bug to the developers! The
output looks something like this:

.. code:: text

    ...
    Error encountered: Stack trace (most recent call first):
    #0 (filtered)
    #1 (filtered)
    #2 (filtered)
    #3 in BoutMesh::createCommunicators()
       at BOUT-dev/src/mesh/impls/bout/boutmesh.cxx:709:64
         707:       }
         708:       // Unconditional exception for demo purposes
       > 709:       throw BoutException("Single null outer SOL not correct\n");
                                                                             ^
         710:       MPI_Group_free(&group);
         711:     }
    #4 in BoutMesh::load()
       at BOUT-dev/src/mesh/impls/bout/boutmesh.cxx:575:22
         573:   /// Communicator
         574:
       > 575:   createCommunicators();
                                   ^
         576:   output_debug << "Got communicators" << endl;
    #5 in BoutInitialise(int&, char**&)
       at BOUT-dev/src/bout++.cxx:201:28
         199:   bout::globals::mesh = Mesh::create();
         200:   // Load from sources. Required for Field initialisation
       > 201:   bout::globals::mesh->load();
                                         ^
         202:
         203:   // time_report options are used in BoutFinalise, i.e. after we
    #6 in main
       at BOUT-dev/examples/elm-pb/elm_pb.cxx:2161:1
        2159: };
        2160:
      > 2161: BOUTMAIN(ELMpb);
              ^
    #7 (filtered)
    #8 (filtered)
    #9 (filtered)

    ====== Exception thrown ======
    Single null outer SOL not correct


Debug symbols are required to get the filename/line number and code snippets. If
they are missing in either BOUT++ or the physics model, only the function name
and signature will be included for that part.

Including debug symbols is a configure time ``CMake`` option, set either:
``-DCMAKE_BUILD_TYPE=Debug`` or ``-DCMAKE_BUILD_TYPE=RelWithDebInfo`` (the
default).

The formatting of this backtrace is controlled in
`BoutException::getBacktrace()` using the `cpptrace
<https://github.com/jeremy-rifkin/cpptrace>`_ library.


Message Stack
=============

Another utility BOUT++ has to help debugging is the message stack using the
`TRACE` macro. This is very useful for when a bug only occurs after a long time
of running, and/or only occasionally. The ``TRACE`` macro can simply be dropped
in anywhere in the code::

    {
      TRACE("Some message here"); // message pushed

    } // Scope ends, message popped

This will push the message, then pop the message when the current
scope ends. If an error occurs or BOUT++ crashes, any un-popped
messages will be printed, along with the file name and line number, to
help find where an error occurred. For example, given this snippet::

  {
    TRACE("1. Outer-most scope");
    {
      TRACE("2. Middle scope");
      {
        TRACE("3. Scope not appearing in the output");
      }
      {
        TRACE("4. Inner-most scope");
        throw BoutException("Something went wrong");
      }
    }
  }

we would see something like the following output:

.. code:: text

    ====== Exception thrown ======
    Something went wrong

    === Additional information ===
    -> 4. Inner-most scope on line 58 of '/path/to/model.cxx'
    -> 2. Middle scope on line 53 of '/path/to/model.cxx'
    -> 1. Outer-most scope on line 51 of '/path/to/model.cxx'


The third ``TRACE`` message doesn't appear in the output because we've
left its scope and it's no longer relevant.

The run-time overhead of this should be small, but can be removed
entirely if the compile-time flag ``-DCHECK`` is not defined or set to
``0``. This turns off checking, and ``TRACE`` becomes an empty
macro. This means that ``TRACE`` macros can be left in your code
permanently, providing some simple diagnostics without compromising
performance, as well as demarcating separate sections with
user-friendly names.

If you need to capture runtime information in the message, you can use
the ``fmt`` syntax also used by the loggers::

    TRACE("Value of i={}, some arbitrary {}", i, "string");


Time evolution
==============

It can be useful to know what happened when the simulation failed.  The pvode
solver can dump the state of the simulation, at the time the solver
failed. This information includes the individual terms in the derivative. This
allows to identify which term is causing the issue.  Additionally, the
residuum is dumped. This identifies not only which term is causing the issue,
but also where in the domain the solver is struggling.  This can be enabled
with::

     solver:type=pvode solver:debug_on_failure=true

It is also possible to dump at a specific time using the euler solver.
This can be useful for tracking down what is causing differences between two
different versions. It can be used with::

      solver:type=euler solver:dump_at_time=0 input:error_on_unused_options=false
