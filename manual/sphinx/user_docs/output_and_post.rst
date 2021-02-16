.. _sec-output:

Post-processing
===============

The recommended tool for analysing BOUT++ output is xBOUT, a Python
library that provides analysis, plotting and animation with
human-readable syntax (no magic numbers!) using `xarray
<http://xarray.pydata.org/en/stable/>`_. See the xBOUT documentation
`xbout.readthedocs.io <https://xbout.readthedocs.io/en/latest/>`_.

There is also older analysis and post-processing code, the majority
written in Python. Routines to read BOUT++ output data, usually called
"collect" because it collects data from multiple files, are also
available in IDL, Matlab, Mathematica and Octave. All these
post-processing routines are in the ``tools`` directory, with Python
modules in ``tools/pylib``. A summary of available routines is in
:ref:`sec-python-routines-list`; see below for how to install the
requirements.

.. _sec-pythonroutines:

Python routines
---------------

.. _sec-python-requirements:

Requirements
~~~~~~~~~~~~

The Python tools provided with BOUT++ make heavy use of numpy_ and
scipy_, as well as matplotlib_ for the plotting routines. In order
to read BOUT++ output in Python, you will need either netcdf4_ or h5py_.

While we try to ensure that the Python tools are compatible with both
Python 2 and 3, we officially only support Python 3.

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

Reading BOUT++ data
~~~~~~~~~~~~~~~~~~~

To read data from a BOUT++ simulation into Python, there is a ``collect`` routine.
This gathers together the data from multiple processors, taking care of the correct
layout.

.. code-block:: python

    from boutdata.collect import collect

    Ni = collect("Ni")  # Collect the variable "Ni"

The result is an up to 4D array, ``Ni`` in this case. The array is a BoutArray
object: BoutArray is a wrapper class for Numpy's ndarray which adds an
'attributes' member variable containing a dictionary of attributes.  The array
is ordered ``[t,x,y,z]``:

.. code-block:: python

    >>> Ni.shape
    [10,1,2,3]

so ``Ni`` would have 10 time slices, 1 point in x, 2 in y, and 3 in z.
This should correspond to the grid size used in the simulation.
Since the collected data is a NumPy array, all the useful routines
in NumPy, SciPy and Matplotlib can be used for further analysis.

The attributes of the data give:

- the ``bout_type`` of the variable

  - {``'Field3D_t'``, ``'Field2D_t'``, ``'scalar_t'``} for time-evolving variables

  - {``'Field3D'``, ``'Field2D'``, ``'scalar'``} for time-independent variables

- its location, one of {``'CELL_CENTRE'``, ``'CELL_XLOW'``, ``'CELL_YLOW'``, ``'CELL_ZLOW'``}. See :ref:`sec-staggergrids`.

.. code-block:: python

    >>> Ni.attributes("bout_type")
    'Field3D_t'
    >>> Ni.attributes("location")
    'CELL_CENTRE'

Attributes can also be read using the ``attributes`` routine:

.. code-block:: python

    from boutdata.collect import attributes

    attribs = attributes("Ni")

The result is a dictionary (map) of attribute name to attribute value.

If the data has less then 4 dimension, it can be checked with
``dimension`` what dimensions are available:

.. code-block:: python

    from boutdata.collect import dimension

    print(dimension("Ni"))
    print(dimension("dx"))

The first will print as expected ``[t, x, y, z]`` - while the second
will print ``[x, y]`` as dx is nether evolved in time, nor does it has
a ``z`` dependency.

To access both the input options (in the BOUT.inp file) and output data, there
is the ``BoutData`` class.

.. code-block:: pycon

    >>> from boutdata.data import BoutData
    >>> d = BoutData(path=".")

where the path is optional, and should point to the directory containing the BOUT.inp 
(input) and BOUT.dmp.* (output) files. This will return a dictionary with keys
"path" (the given path to the data), "options" (the input options) and "outputs" (the output data).
The tree of options can be printed:

.. code-block:: pycon

    >>> print d["options"]
      options
       |- timestep = 50
       |- myg = 0
       |- nout = 50
       |- mxg = 2
       |- all
       |   |- bndry_all = neumann
       |   |- scale = 0.0
       |- phisolver
       |   |- fourth_order = true        
       ...

and accessed as a tree of dictionaries:

.. code-block:: pycon

    >>> print d["options"]["phisolver"]["fourth_order"]
    true

Currently the values are either integers, floats, or strings, so in the above example "true" is a string,
not a Boolean.

In a similar way the outputs are available as dictionary keys:

.. code-block:: pycon

    >>> print d["outputs"]
    ZMAX
    rho_s
    zperiod
    BOUT_VERSION
    ...
    >>> d["outputs"]["rho_s"]
    0.00092165524660235405
    
There are several modules available for reading NetCDF files, so to
provide a consistent interface, file access is wrapped into a class
DataFile. This provides a simple interface for reading and writing files
from any of the following modules: ``netCDF4``;
``Scientific.IO.NetCDF``; and ``scipy.io.netcdf``. To open a file
using DataFile:

.. code-block:: python

    from boututils.datafile import DataFile

    f = DataFile("file.nc")  # Open the file
    var = f.read("variable") # Read a variable from the file
    f.close()                # Close the file

A more robust way to read from DataFiles is to use the context manager
syntax:

.. code-block:: python

    from boututils.datafile import DataFile

    with DataFile("file.nc") as f: # Open the file
        var = f.read("variable")     # Read a variable from the file

This way the DataFile is automatically closed at the end of the ``with``
block, even if there is an error in ``f.read``. To list the variables in
a file e.g.

.. code-block:: pycon

    >>> f = DataFile("test_io.grd.nc")
    >>> print(f.list())
    ['f3d', 'f2d', 'nx', 'ny', 'rvar', 'ivar']

and to list the names of the dimensions

.. code-block:: pycon

    >>> print(f.dimensions("f3d"))
    ('x', 'y', 'z')

or to get the sizes of the dimensions

.. code-block:: pycon

    >>> print(f.size("f3d"))
    [12, 12, 5]

or the dictionary of attributes

.. code-block:: pycon

    >>> print(f.attributes("f3d"))
    {}


To read in all variables in a file into a dictionary there is the
``file_import`` function

.. code-block:: python

    from boututils.file_import import file_import

    grid = file_import("grid.nc")

Python analysis routines
------------------------

The analysis and postprocessing routines are currently divided into two Python modules:
``boutdata``, which contains BOUT++ specific things like ``collect``, and ``boututils``
which contains more generic useful routines.

To plot data, a convenient wrapper around matplotlib is ``plotdata``

.. code-block:: python

    from boutdata import collect
    n = collect("n") # Read data as NumPy array [t,x,y,z]
    
    from boututils.plotdata import plotdata
    plotdata(n[-1,:,0,:])

If given a 2D array as in the above example, plotdata produces a contour plot
(using matplotlib pyplot.contourf) with colour bar. If given a 1D array then it will plot
a line plot (using pyplot.plot).

It is sometimes useful to see an animation of a simulation. To do this there is
``showdata``, which again is a wrapper around matplotlib:

.. code-block:: python

    from boutdata import collect
    n = collect("n") # Read data as NumPy array [t,x,y,z]
    
    from boututils.showdata import showdata
    showdata(n[:,:,0,:])

This always assumes that the first index is time and will be animated over. The above example
animates the variable ``n`` in time, at each time point plotting a contour plot in ``x`` and ``z`` dimensions.
The colour range is kept constant by default. If a 2D array is given to ``showdata`` then a line plot will be
drawn at each time, with the scale being kept constant.



Reading BOUT++ output into IDL
------------------------------

There are several routines provided for reading data from BOUT++
output into IDL. In the directory containing the BOUT++ output files
(usually ``data/``), you can list the variables available using

.. code-block:: idl

    IDL> print, file_list("BOUT.dmp.0.nc")
    Ajpar Apar BOUT_VERSION MXG MXSUB MYG MYSUB MZ NXPE NYPE Ni Ni0 Ni_x Te0 Te_x
    Ti0 Ti_x ZMAX ZMIN iteration jpar phi rho rho_s t_array wci

The ``file_list`` procedure just returns an array, listing all the
variables in a given file.

One thing new users can find confusing is that different simulations may
have very different outputs. This is because **BOUT++ is not a single
physics model**: the variables evolved and written to file are
determined by the model, and will be very different between (for
example) full MHD and reduced Braginskii models. There are however some
variables which all BOUT++ output files contain:

-  ``BOUT_VERSION``, which gives the version number of BOUT++ which
   produced the file. This is mainly to help output processing codes
   handle changes to the output file format. For example, BOUT++ version
   0.30 introduced 2D domain decomposition which needs to be handled
   when collecting data.

-  ``MXG``,\ ``MYG``. These are the sizes of the X and Y guard cells

-  ``MXSUB``, the number of X grid points in each processor. This does
   not include the guard cells, so the total X size of each field will
   be ``MXSUB + 2*MXG``.

-  ``MYSUB``, the number of Y grid points per processor (like MXSUB)

-  ``MZ``, the number of Z points

-  ``NXPE, NYPE``, the number of processors in the X and Y directions.
   ``NXPE * MXSUB + 2*MXG= NX``, ``NYPE * MYSUB = NY``

-  ``ZMIN``, ``ZMAX``, the range of Z in fractions of :math:`2\pi`.

-  ``iteration``, the last timestep in the file

-  ``t_array``, an array of times

Most of these - particularly those concerned with grid size and
processor layout - are used by post-processing routines such as
``collect``, and are seldom needed directly. To read a single variable
from a file, there is the ``file_read`` function:

.. code-block:: idl

    IDL> wci = file_read("BOUT.dmp.0.nc", "wci")
    IDL> print, wci
      9.58000e+06

To read in all the variables in a file into a structure, use the
``file_import`` function:

.. code-block:: idl

    IDL> d = file_import("BOUT.dmp.0.nc")
    IDL> print, d.wci
      9.58000e+06

This is often used to read in the entire grid file at once. Doing this
for output data files can take a long time and use a lot of memory.

Reading from individual files is fine for scalar quantities and time
arrays, but reading arrays which are spread across processors (i.e.
evolving variables) is tedious to do manually. Instead, there is the
``collect`` function to automate this:

.. code-block:: idl

    IDL> ni = collect(var="ni")
    Variable 'ni' not found
    -> Variables are case-sensitive: Using 'Ni'
    Reading from .//BOUT.dmp.0.nc: [0-35][2-6] -> [0-35][0-4]

This function takes care of the case, so that reading “ni” is
automatically corrected to “Ni”. The result is a 4D variable:

.. code-block:: idl

    IDL> help, ni
    NI              FLOAT     = Array[36, 5, 64, 400]

with the indices ``[X, Y, Z, T]``. Note that in the output files, these
variables are stored in ``[T, X, Y, Z]`` format instead but this is
changed by ``collect``. Sometimes you don’t want to read in the entire
array (which may be very large). To read in only a subset, there are
several optional keywords with ``[min,max]`` ranges:

.. code-block:: idl

    IDL> ni = collect(var="Ni", xind=[10,20], yind=[2,2], zind=[0,31],
    tind=[300,399])
    Reading from .//BOUT.dmp.0.nc: [10-20][4-4] -> [10-20][2-2]
    IDL> help, ni
    NI              FLOAT     = Array[11, 1, 32, 100]

Summary of IDL file routines
----------------------------

Functions file\_ can currently only read/write NetCDF files.

Open a NetCDF file:

.. code-block:: idl

    handle = file_open("filename", /write, /create)

Array of variable names:

.. code-block:: idl

    list = file_list(handle)
    list = file_list("filename")

Number of dimensions:

.. code-block:: idl

    nd = file_ndims(handle, "variable")
    nd = file_ndims("filename", "variable")

Read a variable from file. Inds = [xmin, xmax, ymin, ymax, ...]

.. code-block:: idl

    data = file_read(handle, "variable", inds=inds)
    data = file_read("filename", "variable", inds=inds)

Write a variable to file. For NetCDF it tries to match up dimensions,
and defines new dimensions when needed

.. code-block:: idl

    status = file_write(handle, "variable", data)

Close a file after use

.. code-block:: idl

    file_close, handle

To read in all the data in a file into a structure:

.. code-block:: idl

    data = file_import("filename")

and to write a structure to file:

.. code-block:: idl

    status = file_export("filename", data)

IDL analysis routines
---------------------

Now that the BOUT++ results have been read into IDL, all the usual
analysis and plotting routines can be used. In addition, there are many
useful routines included in the ``idllib`` subdirectory. There is a
``README`` file which describes what each of these routines, but some of
the most useful ones are listed here. All these examples assume there is
a variable ``P`` which has been read into IDL as a 4D [x,y,z,t]
variable:

-  ``fft_deriv`` and ``fft_integrate`` which differentiate and integrate
   periodic functions.

-  ``get_integer``, ``get_float``, and ``get_yesno`` request integers,
   floats and a yes/no answer from the user respectively.

-  ``showdata`` animates 1 or 2-dimensional variables. Useful for
   quickly displaying results in different ways. This is useful for
   taking a quick look at the data, but can also produce bitmap outputs
   for turning into a movie for presentation. To show an animated
   surface plot at a particular poloidal location (32 here):

   .. code-block:: idl

       IDL> showdata, p[*,32,*,*]

   To turn this into a contour plot,

   .. code-block:: idl

       IDL> showdata, p[*,32,*,*], /cont

   To show a slice through this at a particular toroidal location (0
   here):

   .. code-block:: idl

       IDL> showdata, p[*,32,0,*]

   There are a few other options, and ways to show data using this code;
   see the README file, or comments in ``showdata.pro``. Instead of
   plotting to screen, showdata can produce a series of numbered bitmap
   images by using the ``bmp`` option

   .. code-block:: idl

       IDL> showdata, p[*,32,*,*], /cont, bmp="result_"

   which will produce images called ``result_0000.bmp``,
   ``result_0001.bmp`` and so on. Note that the plotting should not be
   obscured or minimised, since this works by plotting to screen, then
   grabbing an image of the resulting plot.

-  ``moment_xyzt`` takes a 4D variable (such as those from ``collect``),
   and calculates RMS, DC and AC components in the Z direction.

-  ``safe_colors`` A general routine for IDL which arranges the color
   table so that colors are numbered 1 (black), 2 (red), 3 (green), 4
   (blue). Useful for plotting, and used by many other routines in this
   library.

There are many other useful routines in the ``idllib`` directory. See
the ``idllib/README`` file for a short description of each one.

Matlab routines
---------------

These are Matlab routines for collecting data, showing animation and
performing some basic analysis. To use these routines, either you may
copy these routines (from **tools/matlablib**) directly to your present
working directory or a path to **tools/matlablib** should be added
before analysis.

.. code-block:: matlab

    >> addpath <full_path_BOUT_directory>/tools/matlablib/

Now, the first routine to collect data and import it to Matlab for
further analysis is

.. code-block:: matlab

    >> var = import_dmp(path,var_name);

Here, *path* is the path where the output data in netcdf format has been
dumped. *var\_name* is the name of variable which user want to load for
further analysis. For example, to load “P” variable from present working
directory:

.. code-block:: matlab

    >> P = import_dmp('.','P');

Variable “P” can be any of [X,Y,Z,T]/[X,Y,Z]/[X,Y]/Constant formats. If
we are going to Import a large data set with [X,Y,Z,T] format. Normally
such data files are of very big size and Matlab goes out of memory/ or
may take too much time to load data for all time steps. To resolve this
limitation of above routine *import\_dmp*, another routine
*import\_data\_netcdf* is being provided. It serves all purposes the
routine *import\_dmp* does but also gives user freedom to import data at
only few/specific time steps.

.. code-block:: matlab

    >> var = import_data_netcdf(path,var_name,nt,ntsp);

Here, *path* and *var\_name* are same variables as described before.
*nt* is the number of time steps user wish to load data. *ntsp* is the
steps at which one wish to write data of of total simulation times the
data written.

.. code-block:: matlab

    >> P = import_data_netcdf('.','P',5,100);

Variable “P” has been imported from present working directory for 5 time
steps. As the original netcdf data contains time information of 500
steps (assume NT=500 in BOUT++ simulations), user will pick only 5 time
steps at steps of *ntsp* i.e. 100 here. Details of other Matlab routines
provided with BOUT++ package can be looked in to README.txt of
**tools/matlablib** directory. The Matlab users can develop their own
routines using ***ncread, ncinfo, ncwrite, ncdisp, netcdf etc.***
functions provided in Matlab package.

Mathematica routines
--------------------

A package to read BOUT++ output data into Mathematica is in
``tools/mathematicalib``. To read data into Mathematica, first add this
directory to Mathematica’s path by putting

.. code-block:: mathematica

       AppendTo[$Path,"/full/path/to/BOUT/tools/mathematicalib"]

in your Mathematica startup file (usually
``$HOME/.Mathematica/Kernel/init.m`` ). To use the package, call

.. code-block:: mathematica

       Import["BoutCollect.m"]

from inside Mathematica. Then you can use e.g.

.. code-block:: mathematica

       f=BoutCollect[variable,path->"data"]

or

.. code-block:: mathematica

       f=BoutCollect[variable,path->"data"]

’ ``bc``\ ’ is a shorthand for ’\ ``BoutCollect`` ’. All options
supported by the Python ``collect()`` function are included, though Info
does nothing yet.

Octave routines
---------------

There is minimal support for reading data into Octave, which has been
tested on Octave 3.2. It requires the ``octcdf`` library to access
NetCDF files.

.. code-block:: octave

    f = bcollect()  # optional path argument is "." by default

    f = bsetxrange(f, 1, 10) # Set ranges
    # Same for y, z, and t (NOTE: indexing from 1!)

    u = bread(f, "U")  # Finally read the variable

