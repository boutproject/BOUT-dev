Testing
=======

Two types of tests are currently used in BOUT++ to catch bugs as early
as possible: Unit tests, which check a small piece of the code
separately, and a test suite which runs the entire code on a short
problem. Unit tests can be run using the ``src/unit_tests`` Python
script. This searches through the directories looking for an executable
script called ``unit_test``, runs them, and collates the results. Not
many tests are currently available as much of the code is too tightly
coupled. If done correctly, the unit tests should describe and check the
behaviour of each part of the code, and hopefully the number of these
will increase over time. The test suite is in the ``examples``
directory, and is run using the ``test_suite`` python script. At the top
of this file is a list of the subdirectories to run (e.g. ``test-io``,
``test-laplace``, and ``interchange-instability``). In each of those
subdirectories the script ``runtest`` is executed, and the return value
used to determine if the test passed or failed.

All tests should be short, otherwise it discourages people from running
the tests before committing changes. A few minutes or less on a typical
desktop, and ideally only a few seconds. If you have a large simulation
which you want to stop anyone breaking, find starting parameters which
are as sensitive as possible so that the simulation can be run quickly.

.. _sec-mms:

Method of Manufactured Solutions
--------------------------------

The Method of Manufactured solutions (MMS) is a rigorous way to check
that a numerical algorithm is implemented correctly. A known solution is
specified (manufactured), and it is possible to check that the code
output converges to this solution at the expected rate.

To enable testing by MMS, switch an input option “mms” to true:

.. code-block:: cfg

    [solver]
    mms = true

This will have the following effect:

#. For each evolving variable, the solution will be used to initialise
   and to calculate the error

#. For each evolving variable, a source function will be read from the input file
   and added to the time derivative.

.. note:: The convergence behaviour of derivatives using FFTs is quite
          different to the finite difference methods: once the highest
          frequency in the manufactured solution is resolved, the
          accuracy will jump enormously, and after that, finer grids
          will not increase the accuracy. Whereas with finite
          difference methods, accuracy varies smoothly as the grid is
          refined.
   
Choosing manufactured solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Manufactured solutions must be continuous and have continuous
derivatives. Common mistakes:

-  Don’t use terms multiplying coordinates together e.g. ``x * z`` or
   ``y * z``. These are not periodic in :math:`y` and/or :math:`z`, so
   will give strange answers and usually no convergence. Instead use
   ``x * sin(z)`` or similar, which are periodic.

.. _sec-timerclass:

Timing
------

To time parts of the code, and calculate the percentage of time spent in
communications, file I/O, etc. there is the ``Timer`` class defined in
``include/bout/sys/timer.hxx``. To use it, just create a ``Timer``
object at the beginning of the function you want to time:

::

    #include <bout/sys/timer.hxx>

    void someFunction() {
      Timer timer("test")
      ...
    }

Creating the object starts the timer, and since the object is destroyed
when the function returns (since it goes out of scope) the destructor
stops the timer.

::

    class Timer {
    public:
      Timer();
      Timer(const std::string &label);
      ~Timer();

      double getTime();
      double resetTime();
    };

The empty constructor is equivalent to setting ``label = ""`` .
Constructors call a private function ``getInfo()`` , which looks up the
``timer_info`` structure corresponding to the label in a
``map<string, timer_info*>`` . If no such structure exists, then one is
created. This structure is defined as:

::

    struct timer_info {
      double time;    ///< Total time
      bool running;   ///< Is the timer currently running?
      double started; ///< Start time
    };

Since each timer can only have one entry in the map, creating two timers
with the same label at the same time will lead to trouble. Hence this
code is **not** thread-safe.

The member functions ``getTime()`` and ``resetTime()`` both return the
current time. Whereas ``getTime()`` only returns the time without
modifying the timer, ``resetTime()`` also resets the timer to zero.

If you don’t have the object, you can still get and reset the time using
static methods:

::

    double Timer::getTime(const std::string &label);
    double Timer::resetTime(const std::string &label);

These look up the ``timer_info`` structure, and perform the same task as
their non-static namesakes. These functions are used by the monitor
function in ``bout++.cxx`` to print the percentage timing information.


