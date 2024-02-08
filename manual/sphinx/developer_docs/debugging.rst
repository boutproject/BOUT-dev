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

Message Stack
=============

The second utility BOUT++ has to help debugging is the message stack
using the `TRACE` (and related `AUTO_TRACE`) macro. These are very
useful for when a bug only occurs after a long time of running, and/or
only occasionally. The ``TRACE`` macro can simply be dropped in
anywhere in the code::

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
       
    ====== Back trace ======
    -> 4. Inner-most scope on line 58 of '/path/to/model.cxx'
    -> 2. Middle scope on line 53 of '/path/to/model.cxx'
    -> 1. Outer-most scope on line 51 of '/path/to/model.cxx'

    ====== Exception thrown ======
    Something went wrong

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

There is also an ``AUTO_TRACE`` macro that automatically captures the
name of the function it's used in. This is used throughout the main
library, especially in functions where numerical issues are likely to
arise.

Backtrace
=========

Lastly, BOUT++ can also automatically print a backtrace in the event
of a crash. This is a compile-time option in the BOUT++ library
(``-DBOUT_ENABLE_BACKTRACE=ON``, the default, requires the
``addr2line`` program to be installed), and debug symbols to be turned
on (``-DCMAKE_BUILD_TYPE=Debug`` or ``=RelWithDebInfo``) in BOUT++
_and_ the physics model. If debug symbols are only present in part, the
backtrace will be missing names for the other part.

The output looks something like this:

.. code:: text

    ...
    Error encountered
    ====== Exception path ======
    [bt] #10 ./backtrace() [0x40a27e]
    _start at /home/abuild/rpmbuild/BUILD/glibc-2.33/csu/../sysdeps/x86_64/start.S:122
    [bt] #9 /lib64/libc.so.6(__libc_start_main+0xd5) [0x7fecbfa28b25]
    __libc_start_main at /usr/src/debug/glibc-2.33-4.1.x86_64/csu/../csu/libc-start.c:332
    [bt] #8 ./backtrace() [0x40a467]
    main at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:32 (discriminator 9)
    [bt] #7 /path/to/BOUT-dev/build/libbout++.so(_ZN6Solver8setModelEP12PhysicsModel+0xb5) [0x7fecc0ca2e93]
    Solver::setModel(PhysicsModel*) at /path/to/BOUT-dev/build/../src/solver/solver.cxx:94
    [bt] #6 /path/to/BOUT-dev/build/libbout++.so(_ZN12PhysicsModel10initialiseEP6Solver+0xc0) [0x7fecc0cad594]
    PhysicsModel::initialise(Solver*) at /path/to/BOUT-dev/build/../include/bout/physicsmodel.hxx:93 (discriminator 5)
    [bt] #5 ./backtrace() [0x40a986]
    Backtrace::init(bool) at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:27
    [bt] #4 ./backtrace() [0x40a3cf]
    f3() at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:19
    [bt] #3 ./backtrace() [0x40a3be]
    f2(int) at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:15
    [bt] #2 ./backtrace() [0x40a386]
    f1() at /path/to/BOUT-dev/build/../examples/backtrace/backtrace.cxx:13 (discriminator 2)
    [bt] #1 ./backtrace(_ZN13BoutExceptionC1IA19_cJEEERKT_DpRKT0_+0xba) [0x40ae16]
    BoutException::BoutException<char [19]>(char const (&) [19]) at /path/to/BOUT-dev/build/../include/bout/../boutexception.hxx:28 (discriminator 2)
              

This output tends to be much harder to read than the message stack
from ``TRACE`` macros, but the advantage is that it doesn't require
any modifications to the code to use, and can give you more precise
location information.
