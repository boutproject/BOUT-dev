Testing
=======

There are three types of test used in BOUT++, in order of complexity:
unit tests, integrated tests, and "method of manufactured solutions"
(MMS) tests. Unit tests are very short, quick tests that test a single
"unit" -- usually a single function or method. Integrated tests are
longer tests that range from tests that need a lot of set up and check
multiple conditions, to full physics model tests. MMS tests check the
numerical properties of operators, such as the error scaling of
derivatives.

There is a test suite that runs through all of the unit tests, and
selected integrated and MMS tests. The easiest way to run this is
with:

.. code-block:: console

   $ make check

We expect that any new feature or function implemented in BOUT++ also
has some corresponding tests, and *strongly* prefer unit tests.

.. _sec-automated-testing:

Automated tests and code coverage
---------------------------------

BOUT++ uses `Travis CI`_ to automatically run the test suite on every
push to the GitHub repository, as well as on every submitted Pull
Request. The Travis settings are in ``.travis.yml``. Pull requests
that fail the tests will not be merged.

We also gather information from how well the unit tests cover the
library using `CodeCov`_, the settings for which are stored in
``.codecov.yml``.

.. _Travis CI: https://travis-ci.org/boutproject/BOUT-dev/
.. _CodeCov: https://codecov.io/gh/boutproject/BOUT-dev


.. _sec-unit-tests:

Unit tests
----------

The unit test suits aims to be a comprehensive set of tests that run
*very* fast and ensure the basic functionality of BOUT++ is
correct. At the time of writing, we have around 500 tests that run in
less than a second. Because these tests run very quickly, they should
be run on every commit (or even more often!). For more information on
the unit tests, see ``tests/unit/README.md``.

You can run the unit tests with:

.. code-block:: console

   $ make check-unit-tests


.. _sec-integrated-tests:

Integrated tests
----------------

This set of tests are designed to test that different components of
the BOUT++ library work together. These tests are more expensive than
the unit tests, but are expected to be run on at least every pull
request, and the majority on every commit.

You can run the integrated tests with:

.. code-block:: console

   $ make check-integrated-tests

The test suite is in the ``tests/integrated`` directory, and is run
using the ``test_suite`` python
script. ``tests/integrated/test_suite_list`` contains a list of the
subdirectories to run (e.g. ``test-io``, ``test-laplace``,
``interchange-instability``). In each of those subdirectories the
script ``runtest`` is executed, and the return value used to determine
if the test passed or failed.

All tests should be short, otherwise it discourages people from running
the tests before committing changes. A few minutes or less on a typical
desktop, and ideally only a few seconds. If you have a large simulation
which you want to stop anyone breaking, find starting parameters which
are as sensitive as possible so that the simulation can be run quickly.

Custom test requirements
~~~~~~~~~~~~~~~~~~~~~~~~

Some tests require particular libraries or environments, so should be
skipped if these are not available. To do this, each ``runtest``
script can contain a line starting with ``#requires``, followed by a
python expression which evaluates to ``True`` or ``False``. For
example, a test which doesn't work if both ARKODE and PETSc are used:

.. code-block:: console

   #requires not (arkode and petsc)

or if there were a test which required PETSc to be available, it could
specify

.. code-block:: console

   #requires petsc
   
Currently the requirements which can be combined are ``travis``,
``netcdf``, ``pnetcdf``, ``hdf5``, ``pvode``, ``cvode``,
``ida``, ``lapack``, ``petsc``, ``slepc``, ``arkode``,
``openmp`` and ``make``. The ``make`` requirement is set to True when
the tests are being compiled (but not run), and False when the scripts
are run. It's used for tests which do not have a compilation stage.


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

To time parts of the code, and calculate the percentage of time spent
in communications, file I/O, etc. there is the `Timer` class defined
in ``include/bout/sys/timer.hxx``. To use it, just create a `Timer`
object at the beginning of the function you want to time::

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
created. This structure is defined as::

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
static methods::

    double Timer::getTime(const std::string &label);
    double Timer::resetTime(const std::string &label);

These look up the ``timer_info`` structure, and perform the same task as
their non-static namesakes. These functions are used by the monitor
function in ``bout++.cxx`` to print the percentage timing information.
