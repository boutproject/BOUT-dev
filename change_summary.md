# Change summary

This is a slightly more readable, and therefore incomplete, summary of
the changes from the full [changelog](CHANGELOG.md)

4.3.2 is a bugfix release:
- Make downloading the submodules a bit nicer, including an option for
  using non-bundled versions when using `configure`
- Make `dz` in the Python API a property
- Make `Div_par_K_Grad_par` check the staggered location of its inputs
- Enable split-flux derivatives on staggered fields
- Fix `Grad2_par2` implementation in `InvertParCR`
- Fix an issue writing `FieldPerp`s
- Make it easier to compile the examples with different versions of
  BOUT++, plus fixes for the `tokamak-2fluid` example
- Fix `Div_par` when using more than one y-guard cell
- Fix an issue with text attributes in HDF5 files

4.3.1 is a bugfix release, with a few minor fixes to library code,
notably:
- Fix the creation of the `RGN_OUTER_X` region
- Several small bugs in the Python API
- Preserve restart files if there's a crash during initialisation
- Fix some segfaults in the PvodeSolver
- Fix some issues with Hypnotoad (see
  [\#1783](https://github.com/boutproject/BOUT-dev/pull/1783)
  ([friva000](https://github.com/friva000)))

Other changes are mostly housekeeping changes for the BOUT++ project.

4.3.0 is a big feature release:
- `Field`s are now "tagged" with their "y-direction": that is, whether
  they are in field-aligned space or not. This allows us to perform
  more internal checking, for example, that only field-aligned
  `Field3D`s are passed to `fromFieldAligned` and that calculations
  are done in the correct y-direction space. Users may need to call
  `f.setDirectionY(YDirectionType::Aligned)` for a `Field3D f` to set
  the direction tag
- Add `toFieldAligned` and `fromFieldAligned` free functions
- `bout::utils::is_Field` and variants provide simpler methods of
  checking that input types are `Field`s in templated code
- Many, many more functions and operators support staggering. We're
  now also much more consistent about checking function arguments have
  compatible staggered locations
- New `emptyFrom(f)` and `zeroFrom(f)` helper functions for creating
  `Field`s either allocated but not initialised, or allocated and
  initialised to `0.0` respectively, while ensuring the result is
  compatible with the `Field` `f` (same mesh, same staggering, etc.)
- Expressions used in input files now have support for unicode,
  escaping characters and implicit multiplication. See
  [\#1333](https://github.com/boutproject/BOUT-dev/pull/1333) for more
  details.
- Internationalisation support, including translations for French,
  German, Spanish, and Simplified and Traditional Chinese
- Keyword arguments for boundaries in input files, e.g. `dirichlet(1,
  width=3)`
- File-level attributes in output files
- Write more things to output files, including processor indices and
  parallel transform information
- Complete overhaul of the derivative operators:
    - derivative operators can now use native vectorisation where
      possible, as well parallelisation via OpenMP.
    - users can register their own derivative operators and choose
      them at runtime
    - more consistent handling of staggering in all directions and for
      all `Field` types
    - better handling of field-aligned `Field`s
- Various bug fixes for parallel derivative inversions
- Introduced `zstart` and `zend`, in preparation for introducing
  guard cells in the z-direction
- `Options` has several new features:
    - it can now store `Field`s. This uses an [implementation][mpark]
      of C++17's `std::variant` backported to work with
      C++11. Unfortunately, there are some compilers which have
      problems (see [installation issues][xlc] for help)
    - `Options::isSection` gained the ability to check whether an
      input argument is a subsection or not
    - now records the type of the option when used
    - gained the `.doc` method which allows documentation to be added
      as an attribute, which can be recorded in the `BOUT.settings`
      file post-run
- FFTW is now an optional dependency
- A non-Fourier implementation of `Delp2`
- A generic linear operator inversion class (see
  [\#1439](https://github.com/boutproject/BOUT-dev/pull/1439))
- `Array`, `Matrix` and `Tensor` all gained a `reallocate`
  method. This allows dynamic resizing of those objects, but
  invalidates the existing data
- `FieldFactory` now has separate parsing and generating stages, so
  functions can be parsed once and evaluated multiple times (e.g. for
  time-dependent functions)
- Enable communications for simulations with no core, only divertor
  legs
- `Coordinates` on staggered grids can now be read from grid files
- New `FDDX_U2` implementation
- Support for SUNDIALS versions 2.6 to 4.1.0
- `BoutInitialise` has been pulled apart into separate utility
  functions under the `bout::experimental` namespace so that they can
  be used individually
- `LaplaceCyclic` now accepts `C1` and `C2` coefficients which may be
  different
- `LaplaceNaulin` may use the DC parts of `C` for faster convergence
- `LaplaceNaulin` also gained an under-relaxation factor, which may
  improve convergence and robustness
- The `Laplace` solvers gained a `uses3DCoefs` method. This returns
  `true` if the solver can make use of `Field3D` coefficients rather
  than using the DC component of them
- A new time `Solver`: Runge-Kutta-Legendre stabilised explicit
  method, `splitrk`. See
  [\#1673](https://github.com/boutproject/BOUT-dev/pull/1673) for more
  details
- Experimental support for CMake
- Added `HeatFluxSNB` which calculates heat flux using the
  Shurtz-Nicolai-Busquet (SNB) model. Nonlocal (kinetic) corrections
  to the Spitzer-Harm collisional heat flux, which become important
  when the mean free path becomes a small (~1%) fraction of the
  temperature scale length
- The new `BOUT_SCOREP_REGION("name")` will automatically instrument a
  scope with [Score-P][scorep] if BOUT++ was compiled with Score-P
  support (see the [documentation][scorepdocs] for more information)
- `NYPE` may be given instead of `NXPE` in input files for decomposing
  the mesh in `(x, y)`
- Many fixes and improvements for Hypnotoad:
    - Use centred differencing to compute `dx` from `psi` instead of
      forward differencing
    - Fix for when the separatrix is exactly on a grid point
    - Fix for when the separatrix is very close to the X-point
    - Fix computation of `ShiftAngle`
    - Add a checkbox to the 'Output' tab, which if selected outputs
      metrics for orthogonal coordinates (i.e. using `ShiftedMetric`)
    - We now check what coordinate system was used to generate grid
      files in Hypnotoad when reading them in BOUT++
    - _Lots_ of fixes for non-orthogonal grids (see
      [\#1593](https://github.com/boutproject/BOUT-dev/pull/1593),
      [\#1596](https://github.com/boutproject/BOUT-dev/pull/1596), and
      [\#1636](https://github.com/boutproject/BOUT-dev/pull/1636))
    - Use double precision everywhere
    - Add option to write y-boundary guard cells
    - See [here for a complete
      list](https://github.com/boutproject/BOUT-dev/pulls?q=Hypnotoad+is%3Apr+is%3Amerged+created%3A2018-10-16..2019-10-24++base%3Anext)
- Many, many more tests! Unit test coverage since v4.2.0 has doubled
- We have begun to move parts of the codebase into a `bout::`
  namespace. This should help ensure we play nice with other
  libraries, as well as logically group related things across parts of
  the codebase
  
[mpark]: https://github.com/mpark/variant
[xlc]: https://bout-dev.readthedocs.io/en/latest/user_docs/advanced_install.html#issues
[scorep]: https://www.vi-hps.org/projects/score-p/
[scorepdocs]: https://bout-dev.readthedocs.io/en/latest/developer_docs/performance_profiling.html#scorep-scalasca-profiling

Deprecations:
- `invert_laplace`: create an instance of a `Laplacian` via
  `Laplacian::create` and use the `setCoef*` and `solve` methods
- Karniadakis time `Solver`: the current implementation is _very_ slow
  and will be removed in 5.0
- `MsgStack::setPoint`: use `MsgStack::push("")`
- The following `Mesh` methods were experimental, low-level
  communication routines that turned out to not be so useful:
    - `sendToProc`
    - `receiveFromProc`
    - `UpXSplitIndex`
    - `DownXSplitIndex`
    - `sendYOutIndest`
    - `sendYOutOutdest`
    - `sendYInIndest`
    - `sendYInOutdest`
    - `irecvYOutIndest`
    - `irecvYOutOutdest`
    - `irecvYInIndest`
    - `irecvYInOutdest`
- `Mesh::XGLOBAL` and `Mesh::YGLOBAL`: use `Mesh::getGlobalXIndex` and
  either `Mesh::getGlobalYIndexNoBoundaries` or
  `Mesh::getGlobalYIndex` instead. The former (`NoBoundaries`) is a
  direct replacement for `YGLOBAL`, whereas the latter includes the
  boundaries and so is consistent with `XGLOBAL` which does too
- `Laplacian::setFlags`: use `Laplacian::setGlobalFlags`,
  `Laplacian::setInnerBoundaryFlags` and
  `Laplacian::setOuterBoundaryFlags` instead
- The staggered parallel differential operators that end `CtoL` or
  `LtoC` (e.g. `Div_par_CtoL`, `Grad_par_LtoC`): the corresponding
  versions without the suffix now support staggering. For example,
  instead of `Div_par_CtoL(f)` use `Div_par(f, CELL_YLOW)` instead

Removed:
- The `serial` implementation of `parderiv`. The `cyclic` version
  works both serially and in parallel
- `comm_group`: not used internally and too low-level to be useful
- Support for the `scipy` and `scientific` netCDF libraries in
  `boututils` has been dropped. These were very slow and `scientific`
  is no longer available
- `Laplace3D`: use `Laplacian` instead


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
