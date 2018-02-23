.. _sec-python-routines-list:

Python routines
===============

.. _sec-python-requirements:

Requirements
------------

The Python tools provided with BOUT++ make heavy use of numpy_ and
scipy_, as well as matplotlib_ for the plotting routines. In order
to read BOUT++ output in Python, you will need either netcdf4_ or h5py_.

If you are developing BOUT++, you may also need Jinja2_ to edit some
of the generated code(see :ref:`sec-fieldops` for more information).

You can install most of the required Python modules by running

.. code-block:: console

   $ pip3 install --user --requirement requirements.txt

in the directory where you have unpacked BOUT++. This will install
supported versions of numpy, scipy, netcdf4, matplotlib and jinja2.

.. note:: If you have difficulties installing SciPy, please see their
          `installation instructions`_


.. _numpy: http://www.numpy.org/
.. _scipy: http://www.scipy.org/
.. _matplotlib: https://www.matplotlib.org
.. _netcdf4: http://unidata.github.io/netcdf4-python/
.. _h5py: http://www.h5py.org
.. _Jinja2: http://jinja.pocoo.org/
.. _installation instructions: https://www.scipy.org/install.html

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

.. _sec-bout_runners:

bout_runners
------------

``bout_runners`` contains classes which gives an alternative way of
running BOUT++ simulations either normally using the class
``basic_runner`` , or on a cluster through a generated Portable Batch
System (PBS) script using the child class ``PBS_runner`` . Examples can
be found in ``examples/bout_runners_example/``.

``bout_runners`` is especially useful if one needs to make several runs
with only small changes in the options (which is normally written in
``BOUT.inp`` or in the command-line), as is the case when performing a
parameter scan, or when performing a MMS test.

Instead of making several runs with several different input files with
only small changes in the option, one can with ``bout_runners`` specify
the changes as member data of an instance of the appropriate
``bout_runners`` class. One way to do this is to write a *driver* in the
same directory as the executable. The *driver* is just a python script
which imports ``bout_runners`` , creates an instance, specifies the
running option as member data of that instance and finally calls the
member function ``self.execute_runs()`` .

In addition, the ``bout_runners`` provides a way to run any python
post-processing script after finished simulations (as long as it accept
at least one parameter containing the folder name(s) of the run(s)). If
the simulations have been performed using the ``PBS_runner`` , the
post-processing function will be submitted to the cluster (although it
is possible to submit it to a different queue, using a different amount
of nodes etc.).

When the function ``self.execute_runs()`` is executed, a folder
structure like the one presented in :numref:`fig-folder-tree` is
created. ``BOUT.inp`` is copied to the folder of execution, where the
``BOUT.*.dmp`` files are stored. Secondly a list of combination of the
options specified in the driver is made. Eventually unset options are
obtained from ``BOUT.inp`` or given a default value if the option is
nowhere to be found.

.. _fig-folder-tree:
.. figure:: ../figs/folder_tree.*
   :alt: Longest possible folder tree

   Longest possible folder tree made by the ``self.execute_runs()``
   function.
