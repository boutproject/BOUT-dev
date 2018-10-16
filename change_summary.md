# Change summary

This is a slightly more readable, and therefore incomplete, summary of
the changes from the full [changelog](CHANGELOG.md)

4.2.0 is a big feature release:
- Large number of optimisations (as much as 140% faster than v4.1.2!)
- OpenMP in many more places, enables parallelisation in Z (as well as
  X for FCI)
- Better support for OpenMP, including in Python tools
- Much more versatile region system, allowing arbitrary regions (can
  e.g. mask certain parts of the domain for most common operations)
- Specialised macro for looping over Fields, handles arbitrary
  regions, OpenMP parallelisation, while also supporting native
  vectorisation
- Add support for new region system to many functions
- Better support for staggered grids: many bugfixes and many more
  functions support setting the location
- `Coordinates` objects can be created at different locations, through
  the `Mesh::getCoordinates` and `Field::getCoordinates` methods
- Support for compiling as a shared library
- Experimental Python API via Cython module
- Arithmetic operators on fields are now generated using Jinja2
- Improved PETSc compatibility (better support out of the box,
  supports up to 3.9, drops support for versions before 3.4)
- New support classes for 2D/3D arrays (Matrix/Tensor)
- New interface for Options
- Divergence operators for FCI
- Support for attributes in NetCDF files
- Default Laplacian changed to cyclic
- Many C++ modernisation fixes
- New monotonic Hermite spline interpolator
- Better configure experience
- Zoidberg can produce curvilinear grids in all three directions
  (enables e.g. stellarator geometry. Current version of BOUT++ can't
  actually handle this yet -- upcoming version!)
- Many more tests, and a better testing framework for the integrated
  tests
- Some potential memory leaks and null pointer dereferences fixed

Deprecations:
- `DataIterator` is deprecated in favour of the new `Region` and
  `Ind2D/3D/Perp` family. This should not affect user code -- if it
  does, replacing `DataIterator` with `auto` should do the right thing
  in most cases
- `DataFile::writeVar`: use `DataFile::addOnce`
- `Field::setName` and `Field::getName`: just use `Field::name`
  directly instead
- `Field::error` and `bout_error`: use `BoutException` instead
- `rvector`/`rmatrix`/`rtensor` families of functions: use
  `Matrix`/`Tensor` instead
- `operator^(Vector2D/Vector3D)`: use `cross()` instead
- The derivative function overloads with this order of arguments:
  `DD?(..., DIFF_METHOD, CELL_LOC, REGION)`. Instead, use
  `DD?(..., CELL_LOC, DIFF_METHOD, REGION)`
- Vector derivative function overloads with three separate
  `outloc_[xyz]` arguments: use the versions with a single `outloc`
  argument instead
- `CyclicReduce::setCoefs` and `solve` overloads that take `T[]` or
  `T**`: use the version that takes `Array<T>` instead
- The `FCI` class constructors that take a `bool yperiodic` argument:
  this is no longer supported
- `Mesh::coordinates` is deprecated in favour of the more
  consistently-named `Mesh::getCoordinates`. There is also now
  `Field::getCoordinates` which may be more convenient

Removed functions:
- `PhysicsModel::addToRestart` and `Solver::addToRestart`: use
  `restart.add` directly instead
- `Solver::addMonitor(MonitorFunc)`: use the `Monitor*` overloads instead
- The `get/set` `array/data` methods in the `Field` classes: these
  methods are no longer supported

4.1.0 has some interesting new features:
- A way to cleanly stop simulations either through a stop-file, or via a KILL
  signal
- User-defined multiple monitors with different frequencies
- Four new boundary iterators for the inner and outer boundaries in
  double null configurations
- Better handling of 1D/2D Fields
- Dumping the actual settings used during a simulation to a file
- Colour output to highlight warnings/errors, along with options to
  increase/decrease verbosity
- Configure-time options to enable profiling with Score-P and coverage
  checking with gcov
- Removal of various functions deprecated in v4.0.0
- An overhauled testing framework with unit tests using GoogleTest

4.0.1 is purely a bug-fix release, with various fixes to some of the
python tools, patching several memory leaks, and some important fixes
for Intel compilers.
