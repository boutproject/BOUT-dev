.. _sec-zoidberg:

Zoidberg grid generator
=======================

The Zoidberg grid generator creates inputs for the Flux Coordinate Independent (FCI)
parallel transform (section :ref:`sec-parallel-transforms`). The domain is
divided into a set of 2D grids in the X-Z coordinates, and the magnetic field is followed 
along the Y coordinate from each 2D grid to where it either intersects the forward
and backward grid, or hits a boundary.

The simplest code which creates an output file is:

.. code:: python

   import zoidberg

   # Define the magnetic field
   field = zoidberg.field.Slab()
   # Define the grid points
   grid = zoidberg.grid.rectangular_grid(10,10,10)
   # Follow magnetic fields from each point
   maps = zoidberg.make_maps(grid, field)
   # Write everything to file
   zoidberg.write_maps(grid, field, maps, gridfile="grid.fci.nc")

As in the above code, creating an output file consists of the following steps:

1. Define a magnetic field
2. Define the grid points. This can be broken down into:
   
   a) Define 2D "poloidal" grids
   b) Form a 3D grid by putting 2D grids together along the Y direction

3. Create maps from each 2D grid to its neighbours
4. Save grids, fields and maps to file

Each of these stages can be customised to handle more complicated
magnetic fields, more complicated grids, and particular output formats. 
Details of the functionality available are described in sections below; 
here we will describe the key concepts.

An important input is the size of the domain in Y, and
whether the domain is periodic. By default ``rectangular_grid`` makes
a non-periodic rectangular box which is of length 10 in the Y direction.
This means that there are boundaries at :math:`y=0` and at :math:`y=10`.
``rectangular_grid`` puts the y slices at equally spaced intervals, and puts
the first and last points half an interval away from boundaries in y.
In this case with 10 points in y (second argument to ``rectangular_grid(nx,ny,nz)``)
the y locations are :math:`\left(0.5, 1.5, 2.5, \ldots, 9.5\right)`.

At each of these y locations ``rectangular_grid`` defines a rectangular 2D poloidal grid in
the X-Z coordinates, by default with a length of 1 in each direction and centred on :math:`x=0,z=0`. 
These 2D poloidal grids are then put together into a 3D ``Grid``. This process can be customised
by separating step 2 (the ``rectangular_grid`` call) into stages 2a) and 2b). 
For example, to create a periodic rectangular grid we could call:

.. code:: python
   
   import numpy as np

   # Create a 10x10 grid in X-Z with sides of length 1
   poloidal_grid = zoidberg.poloidal_grid.RectangularPoloidalGrid(10, 10, 1.0, 1.0)
   # Define the length of the domain in y
   ylength = 10.0
   # Define the y locations
   ycoords = np.linspace(0.0, ylength, 10, endpoint=False)
   # Create the 3D grid by putting together 2D poloidal grids
   grid = zoidberg.grid.Grid(poloidal_grid, ycoords, ylength, yperiodic=True)

In the above code the length of the domain in the y direction needs to be given to ``Grid``
so that it knows where to put boundaries (if not periodic), or where to wrap the domain
(if periodic). The array of y locations ycoords can be arbitrary, but note that finite
difference methods (like FCI) work best if grid point spacing varies smoothly.

In the last example only one poloidal grid was created (a ``RectangularPoloidalGrid``)
and then re-used for each y slice. We can instead define a different grid for each y
position. For example, to define a grid which expands along y (for some reason) we could do:

.. code:: python

   ylength = 10.0
   ycoords = np.linspace(0.0, ylength, 10, endpoint=False)
   # Create a list of poloidal grids, one for each y location
   poloidal_grids = [ RectangularPoloidalGrid(10, 10, 1.0 + y/10., 1.0 + y/10.)
                      for y in ycoords ]
   # Create the 3D grid by putting together 2D poloidal grids
   grid = zoidberg.grid.Grid(poloidal_grids, ycoords, ylength, yperiodic=True)

Note: Currently there is an assumption that the number of X and Z points is the
same on every poloidal grid. The shape of the grid can however be completely
different. The construction of a 3D ``Grid`` is the same in all cases, so for now
we will concentrate on producing different poloidal grids.

The FCI technique is not restricted to rectangular grids, and in particular
Zoidberg can handle structured grids in an annulus with quite complicated shapes.
The `StructuredPoloidalGrid` class handles more general geometries,
but still assumes that the grid is structured and logically rectangular.
Currently it also assumes that the z index is periodic.

One way to create this grid is to define the grid points manually e.g.:


.. code:: python
          
   import numpy as np
   import zoidberg
          
   r,theta = np.meshgrid(np.linspace(1,2,10),  # minor radius
                         np.linspace(0,2*np.pi, 10), # angle
                         indexing='ij')
   
   R = r * np.sin(theta)
   Z = r * np.cos(theta)

   poloidal_grid = zoidberg.poloidal_grid.StructuredPoloidalGrid(R,Z)

For more complicated shapes than circles, Zoidberg comes with an elliptic grid
generator which needs to be given only the inner and outer 
boundaries.

.. code:: python

   import zoidberg

   inner = zoidberg.rzline.shaped_line(R0=3.0, a=0.5,
                            elong=1.0, triang=0.0, indent=1.0,
                            n=50)
   
   outer = zoidberg.rzline.shaped_line(R0=2.8, a=1.5,
                            elong=1.0, triang=0.0, indent=0.2,
                            n=50)
   
   grid = zoidberg.poloidal_grid.grid_elliptic(inner, outer, 100, 100, show=True)


which should produce the figure below:

.. figure:: ../figs/zoidberg/elliptic_grid.png
   :name: elliptic
   :alt: 
   :scale: 50
   
   A grid produced by ``grid_elliptic`` from shaped inner and outer lines






There are several examples in the `examples/zoidberg` directory
   
Magnetic fields
---------------

The magnetic field is represented by a ``MagneticField`` class, in ``zoidberg.field``.

Slabs and curved slabs
~~~~~~~~~~~~~~~~~~~~~~

The simplest magnetic field is a straight slab geometry:

.. code:: python

   import zoidberg
   field = zoidberg.field.Slab()

By default this has a magnetic field :math:`\mathbf{B} = \left(0, 1, 0.1 + x\right)`.

A variant is a curved slab, which is defined in cylindrical coordinates
and has a given major radius (default 1):

.. code:: python

   import zoidberg
   field = zoidberg.field.CurvedSlab()

Note that this uses a large aspect-ratio approximation, so the major radius
is constant across the domain (independent of x). 
    
Straight stellarator
~~~~~~~~~~~~~~~~~~~~

This is generated by four coils with alternating currents arranged
on the edge of a circle, which spiral around the axis. 

This requires Sympy to generate the magnetic field, so if unavailable
an exception will be raised. 

.. code:: python
   
   import zoidberg
   field = zoidberg.StraightStellarator()



   

Plotting the magnetic field
---------------------------

Routines to plot the magnetic field are in ``zoidberg.plot``. 
To plot a Poincare plot, pass the ``MagneticField`` object,
start location(s) and periodicity information:

.. code:: python

   zoidberg.plot.plot_poincare(field, 0.1, 0.0, 1.0)

The inputs here are the starting location :math:`\left(x,z\right) = \left(0.1, 0.0\right)`,
and the periodicity in the y direction (1.0). By default this will
integrate from this given starting location 40 times (``revs`` option) around the y domain (0 to 1.0). 

An example generated by `this code <https://github.com/boutproject/BOUT-dev/blob/zoidberg-poloidal-grids/examples/zoidberg/poincare.py>`_ is shown in the poincare_ figure:

.. figure:: ../figs/zoidberg/poincare.png
   :name: poincare
   :alt: Points on four oval shaped flux surfaces in x-z at three locations along the y direction
   :scale: 50
   

   Poincare map of straight stellarator. Each colour corresponds to a different x-z plane
   in the y direction. Four flux surfaces are shown, each started at a point at :math:`y=0, z=0`.

   
         
Creating poloidal grids
-----------------------

The FCI technique is used for derivatives along the magnetic field
(in Y), and doesn't restrict the form of the grid in the X-Z
poloidal planes. A 3D grid created by Zoidberg is a collection of 2D planes
(poloidal grids), connected together by interpolations along
the magnetic field.To define a 3D grid we first need to define
the 2D poloidal grids.

Two types of poloidal grids can currently be created: Rectangular grids, and
curvilinear structured grids. All poloidal grids have the following
methods:

* `getCoordinate()` which returns the real space (R,Z) coordinates
  of a given (x,z) index, or derivatives thereof
* `findIndex()` which returns the (x,z) index of a given (R,Z) coordinate
  which in general is floating point
* `metric()` which returns the 2D metric tensor
* `plot()` which plots the grid

Rectangular grids
~~~~~~~~~~~~~~~~~

To create a rectangular grid, pass the number of points and lengths in the x and z directions
to ``RectangularPoloidalGrid``:

.. code:: python

   import zoidberg
   
   rect = zoidberg.poloidal_grid.RectangularPoloidalGrid( nx, nz, Lx, Lz )

By default the middle of the rectangle is at :math:`\left(R,Z\right) = \left(0,0\right)`
but this can be changed with the `Rcentre` and `Zcentre` options.



Curvilinear structured grids
~~~~~~~~~~~~~~~~~~~~~~~~~~~~


   
Here the ``shaped_line`` function creates RZline shapes with the following formula:

.. math::
   
   R = R_0 - b + \left(a + b \cos\left(\theta\right)\cos\left(\theta + \delta\sin\left(\theta\right)\right)\right)

   Z = \left(1 + \epsilon\right)a\sin\left(\theta\right)

where :math:`R_0` is the major radius, :math:`a` is the minor radius,
:math:`\epsilon` is the elongation (``elong``), :math:`\delta` the triangularity (``triang``), and :math:`b` the indentation (``indent``).

