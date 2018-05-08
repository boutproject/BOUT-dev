# BOUT++ Tests

The BOUT++ test suite is split into three parts:

- [unit tests](./unit/README.md)
- [integrated tests](./integrated/README.md)
- [MMS tests](./MMS/README.md) (Method of Manufactured Solutions)

<!-- TODO: Pylib, etc. -->

## Unit tests

These aim to be a comprehensive set of tests that run *very* fast and
ensure the basic functionality of BOUT++ is correct. Because these
tests run very quickly, they should be run on every commit (or even
more often!).

## Integrated tests

These tests are more expensive than the unit tests and are designed to
make sure that the different parts of BOUT++ work together correctly.
The integrated tests range from checking that simulations can be
restarted to full-blown physics test cases.

Most of these tests should be run on every pull request and on most
commits.

## MMS tests

The Method of Manufactured Solutions (MMS) is a technique that aims to
verify the correct numerical implementation of mathematical
operators. These tests check that the numerical operators have the
expected error scalings from the discretisation used.

Some of the MMS tests may be very expensive and so we don't expect
them to be run for every commit or pull request.

## Test requirements

Some tests require a particular library to be available, or are known
to not work or take too long in a particular environment (looking
at you Travis). The integrated and MMS test scripts (called "runtest")
can contain lines starting with "#requires", which specify what that
test requires. If a requirement is not met, then the test is skipped.

The "requirements" subdirectory contains python code to handle these
requirements (in "__init__.py"), and some scripts which are used to
test for some requirements.




