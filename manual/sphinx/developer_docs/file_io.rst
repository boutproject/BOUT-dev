.. _sec-file-io:

File I/O
========

BOUT++ deals with essentially two types of binary files: "grid" files
which are input files read collectively and distributed over all
processes; and "dump" files which are used for input and/or output are
unique to each process. Dump files are accessed through the
`bout::OptionsNetCDF` interface, while grid files have another layer
over the top of this in order to handle the distribution over
processes.

`bout::OptionsNetCDF` provides an interface to read and write
`Options` to/from netCDF files. See :ref:`sec-options` for how to use
`Options`, and :ref:`sec-options-netcdf` for how to use
`bout::OptionsNetCDF`.

.. _sec-file-io-v5:

Changes in BOUT++ v5
--------------------

BOUT++ previously also supported the Portable Data Binary (PDB)
format, developed at LLNL, as well as HDF5 files directly, but support
for these formats were removed in BOUT++ v4 and v5, respectively.  The
whole file I/O system has been refactored in v5 to take advantage of
`bout::OptionsNetCDF` and drastically simplify the library code. This
has in turn changed how some I/O is handled. Previously, there was a
global ``bout::globals::dump`` object of type ``Datafile``, which
wrapped an internal ``Dataformat`` base class, which was an interface
to various concrete implementations of the different file formats
supported.

The ``dump`` instance stored non-``const`` pointers to variables, and
wrote them to disk periodically. This had some significant downsides:
writing ``const`` variables was painful, expressions required an
intermediate variable, and variables had to remain in scope until
``dump`` actually flushed to disk. This was easy to get wrong and
accidentally end up writing nonsense.

In v5, file IO is now done entirely through `Options` and
`bout::OptionsNetCDF`: data is stored in an ``Options`` instance and
then written to disk through ``bout::OptionsNetCDF``. The data is
copied into the ``Options`` instance, so there are no issues with
``const`` variables, expressions, or scope.

With the defunct ``dump`` being a global instance, any object anywhere
could write to file. Now that it has been removed, we have a new idiom
for objects that wish to write data: ``void outputVars(Options&)``. An
existing ``Options`` instance is passed in and the object sets whatever
data it likes. This allows the same ``Options`` instance to be passed
around and be added to before being written, as well as easily
enabling objects to write to separate files if required.

All of this has meant that some changes are necessary to physics
models when it comes to writing data. There are now two more
``virtual`` methods on `PhysicsModel`::

  /// Output additional variables other than the evolving variables
  virtual void outputVars(Options& options);
  /// Add additional variables other than the evolving variables to the restart files
  virtual void restartVars(Options& options);

`PhysicsModel::outputVars` is our new idiomatic method, while
`PhysicsModel::restartVars` is very similar just targeting the restart
file only. This allows the developer to add auxiliary variables to the
restart file. Instead of the global ``dump`` and member ``restart``
``Datafile`` instances, there are two ``Options`` members of
``PhysicsModel``: `PhysicsModel::output_file` and
`PhysicsModel::restart_file`. These are ``private`` and so can only be
accessed through the ``outputVars`` and ``restartVars`` methods.

As in previous versions, the `Solver` takes care of ensuring evolving
variables are added to the output and restart files through
`Solver::outputVars`, though now it takes a reference to an `Options`
rather than a ``Datafile``.

For non-evolving variables, the developer can ``override``
``PhysicsModel::outputVars`` to add further data, diagnostics, and
auxiliary variables to the output file.

The `PhysicsModel::PhysicsModelMonitor` is called by the solver every
``nout`` timesteps (see :ref:`sec-options-general` fore more details
on input options). This monitor calls `Solver::outputVars` to set the
evolving variables, and then `PhysicsModel::outputVars` to set any
model-specific data. The restart file is handled similarly.

To help upgrading to BOUT++ v5, there is a `bout::DataFileFacade`
member of `PhysicsModel` called ``dump``. This class provides a
similar interface to ``Datafile``, and just like ``Datafile`` it also
stores pointers to variables. This means that it still suffers from
all of the downsides of ``Datafile``, and developers are encouraged to
move to the ``outputVars`` approach. The `SAVE_ONCE`/`SAVE_REPEAT`
macros also work through ``DataFileFacade`` -- this means that they
cannot be used outside of ``PhysicsModel`` methods!


FieldPerp I/O
-------------

`FieldPerp` objects can be saved to output files and read from them. The `yindex` of a
`FieldPerp` is the local y-index on a certain processor, but is saved in output files as a
global y-index in the attribute `yindex_global`. The intention is that a `FieldPerp` being
saved should be a globally well-defined object, e.g. a set of values at one divertor
target boundary, that will only be saved from processors holding that global
y-index. The expectation is that the other processors would all save an invalid
`FieldPerp` variable, with a `yindex_global` that is more negative than the
lowest y-boundary guard cell [2]_. The reason for saving the invalid `FieldPerp` variables
is so that all variables are present in every dump file (even if they are not allocated or
used); in particular the Python `collect` routine assumes that any variable will be found
in the first output file, which `collect` uses to get its type and dimensions.

.. [2] Actually, the C++ I/O code should work fine even if a `FieldPerp` object is defined
       with different y-indices on different processors. This may be useful for diagnostic
       or debugging purposes. However, Python routines like `collect` and
       :py:`boutdata.restart.redistribute` will fail because they find inconsistent
       `yindex_global` values.
