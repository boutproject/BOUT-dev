# BOUT++ Unit Tests

These unit tests are built around the [Googletest][gtest_readme] unit
testing framework. The overarching goal of this set of tests is that
we should be able to run it as often as possible, at least every
commit, and ideally more often! Therefore, each individual test should
be **very** fast -- preferably less than 10 ms per test. Each test
should aim to test a single function, but there may be multiple tests
per function.

## Running tests

There is a top-level make target `make check-unit-tests` that will run
the unit test suite. Currently, this just builds and runs a single
executable, `serial_tests`, that contains of the tests that need
to/can be run serially. To list the names of these tests, run:

    ./serial_tests --gtest_list_tests
    
You can then run a subset of tests using

    ./serial_tests --gtest_filter=POSTIVE_PATTERNS[-NEGATIVE_PATTERNS]
    
which will run only the tests whose name matches one of the positive
patterns but none of the negative patterns. '?' matches any single
character; '*' matches any substring; ':' separates two patterns.

Run

    ./serial_tests --help

to get a more complete list of options.

### Test coverage

In a perfect world, the unit tests would test every single line of
code in the BOUT++ library. You can see exactly which lines they do
cover using gcov. BOUT++ provides a configure flag and make target to
set this up for you. To get a coverage summary, run the following:

    # Set up the correct flags, etc.
    ./configure --enable-code-coverage <other configure flags>
    # Compile the main BOUT++ library
    make
    # Remove any old coverage information
    make code-coverage-clean
    # Run the test suite
    make check-unit-tests
    # Gather the code coverage information and produce a summary
    make code-coverage-capture

This will produce a set of HTML files with information about which
source lines are hit while running the tests. Once complete, you can
open [bout-coverage/index.html](../../bout-coverage/index.html) to see
a summary.

- Note: `make code-coverage-capture` will most likely produce lots of
  warnings like:
  
      geninfo: WARNING: no data found for /usr/include/c++/4.8/bits/shared_ptr_base.h

  These can be safely ignored, so you may wish to generate the
  coverage report with `make code-coverage-capture 2> /dev/null`

## Writing tests

We expect that any new feature or function implemented in BOUT++ also
has some corresponding unit tests. See
the [Googletest primer][gtest_primer] for an introduction to writing
tests using Googletest.

The tests are laid to (roughly) mirror the structure of the main
source code, with one test file per implementation file, in the
corresponding directory, e.g. for `src/field/field3d.cxx`, there is
`tests/unit/field/test_field3d.cxx` which contains the tests for the
`Field3D` class.

There are some helper functions and classes
in [test_extras.hxx][test_extras], notably `FakeMesh`. This is an
implementation of the `Mesh` class with the bare minimum needed to
work. It should as little work as possible -- no file I/O, no
communication, etc. To see it in action, look at
the [Field3D tests][test_field3d].


[gtest_readme]: ../../googletest/README.md
[gtest_primer]: ../../googletest/googletest/docs/Primer.md
[test_extras]: ./test_extras.hxx
[test_field3d]: ./field/test_field3d.cxx
