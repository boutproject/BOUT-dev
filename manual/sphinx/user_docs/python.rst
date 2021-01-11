.. _sec-python-routines-list:

Python routines
===============

boututils
---------

-  ``class Datafile`` provides a convenient way to read and write NetCDF
   or HDF5 files. There are many different NetCDF libraries available
   for Python, so this class tries to provide a consistent interface to
   many of them, as well as to h5py.

-  ``deriv()``

-  ``determineNumberOfCPUs()``

-  ``file_import()`` reads the contents of a NetCDF file into a
   dictionary

-  ``integrate()``

-  ``launch()``

-  ``linear_regression()``

-  ``showdata()`` visualises and animates 2D data (time + 1 spatial dimension) or 3D data (time + 2 spatial dimensions). The animation object can be returned, or the animation can be saved to a file or displayed on screen.

-  ``boutwarnings`` contains functions to raise warning messages.
   ``alwayswarn()`` by default prints the warning every time it is called.
   ``defaultwarn()`` by default prints the warning only the first time an
   instance of it is called. This module is a wrapper for the Python
   ``warnings`` module, so printing the warnings can be controlled using
   ``warnings.simplefilter()`` or ``warnings.filterwarnings()``.

.. automodule:: boututils
   :members:
   :undoc-members:

boutdata
--------

-  ``collect()`` provides an interface to read BOUT++ data outputs,
   returning NumPy arrays of data. It deals with the processor layout,
   working out which file contains each part of the domain.

   .. code-block:: python

           from boutdata.collect import collect

           t = collect("t_array")  # Collect the time values


-  ``pol_slice()`` takes a 3 or 4-D data set for a toroidal equilibrium,
   and calculates a slice through it at fixed toroidal angle.

-  ``gen_surface()`` is a generator for iterating over flux surfaces

.. automodule:: boutdata
   :members:
   :undoc-members:

