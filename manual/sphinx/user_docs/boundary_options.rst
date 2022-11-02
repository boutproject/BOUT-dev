.. _sec-bndryopts:

Boundary conditions
===================

Like the variable initialisation, boundary conditions can be set for
each variable in individual sections, with default values in a section
``[All]``. Boundary conditions are specified for each variable, being
applied to variable itself during initialisation, and the
time-derivatives at each timestep. They are a combination of a basic
boundary condition, and optional modifiers.

When finding the boundary condition for a variable ``var`` on a boundary
region, the options are checked in order from most to least specific:

-  Section ``var``, ``bndry_`` + region name. Depending on the mesh
   file, regions of the grid are given labels. Currently these are
   ``core``, ``sol``, ``pf`` and ``target`` which are intended for
   tokamak edge simulations. Hence the variables checked are
   ``bndry_core``, ``bndry_pf`` etc.

-  Section ``var``, ``bndry_`` + boundary side. These names are ``xin``,
   ``xout``, ``yup`` and ``ydown``.

-  Section ``var``, variable ``bndry_all``

-  The same settings again except in section ``All``.

The default setting for everything is therefore ``bndry_all`` in the
``All`` section.

Boundary conditions are given names, with optional arguments in
brackets. Currently implemented boundary conditions are:

-  ``dirichlet`` - Set to zero

-  ``dirichlet(<number>)`` - Set to some number e.g. ``dirichlet(1)``
   sets the boundary to :math:`1.0`

-  ``neumann`` - Zero gradient

-  ``robin`` - A combination of zero-gradient and zero-value
   :math:`a f + b{{\frac{\partial f}{\partial x}}} = g` where the
   syntax is ``robin(a, b, g)``.

-  ``constgradient`` - Constant gradient across boundary

-  ``zerolaplace`` - Laplacian = 0, decaying solution (X boundaries
   only)

-  ``zerolaplace2`` - Laplacian = 0, using coefficients from the
   Laplacian inversion and Delp2 operator.

-  ``constlaplace`` - Laplacian = const, decaying solution (X boundaries
   only)

The zero- or constant-Laplacian boundary conditions works as follows:

.. math::

   \nabla_\perp^2 f &= 0 \\
   &\simeq g^{xx}\frac{\partial^2 f}{\partial x^2} + g^{zz}\frac{\partial^2 f}{\partial z^2}

which when Fourier transformed in :math:`z` becomes:

.. math::

   g^{xx}\frac{\partial^2 \hat{f}}{\partial x^2} - g^{zz}k_z^2 \hat{f} = 0

which has the solution

.. math::

   \hat{f} = Ae^{xk_z\sqrt{g^{zz}/g^{xx}}} + Be^{-xk_z\sqrt{g^{zz}/g^{xx}}}

Assuming that the solution should decay away from the domain, on the
inner :math:`x` boundary :math:`B = 0`, and on the outer boundary
:math:`A = 0`.

Boundary modifiers change the behaviour of boundary
conditions, and more than one modifier can be used. Currently the
following are available:

-  ``relax`` - Relaxing boundaries. Evolve the variable towards the
   given boundary condition at a given rate

-  ``shifted`` - Apply boundary conditions in orthogonal X-Z
   coordinates, rather than field-aligned

-  ``width`` - Modifies the width of the region over which the boundary
   condition is applied

These are described in later subsections.

Boundary conditions for non-orthogonal grids
--------------------------------------------

If non-orthogonal grids are used (meaning that the x- and y-directions are not orthogonal,
so ``g12 != 0.``), then corner cells may be required. The boundary conditions are applied
in corner cells[#disablecorners]_ by applying the y-boundary condition using x-boundary
values. This requires that x-boundary conditions are applied before y-boundary conditions.
The ordering is taken care of by the methods described in this section, but also needs to
be respected by any custom boundary conditions in user code (e.g. sheath boundary
conditions). Note that the iterators returned by the ``BoutMesh`` methods
``iterateBndryLowerY``, ``iterateBndryLowerInnerY``, ``iterateBndryLowerOuterY``,
``iterateBndryUpperY``, ``iterateBndryUpperInnerY``, and ``iterateBndryUpperOuterY``
do include the corner cells at the domain boundary corners.

.. [#disablecorners] although this may be disabled, reverting to the behaviour of BOUT++
                     up to v4, by setting the option ``mesh:include_corner_cells =
                     false``.

Relaxing boundaries
-------------------

All boundaries can be modified to be “relaxing” which are a combination
of zero-gradient time-derivative, and whatever boundary condition they
are applied to. The idea is that this prevents sharp discontinuities at
boundaries during transients, whilst maintaining the desired boundary
condition on longer time-scales. In some cases this can improve the
numerical stability and timestep.

For example, ``relax(dirichlet)`` will make a field :math:`f` at point
:math:`i` in the boundary follow a point :math:`i-1` in the domain:

.. math::

   .{{\frac{\partial f}{\partial t}}}|_i = .{{\frac{\partial f}{\partial t}}}|_{i-1}  - f_i / \tau

where :math:`\tau` is a time-scale for the boundary (currently set to
0.1, but will be a global option). When the time-derivatives are slow
close to the boundary, the boundary relaxes to the desired condition
(Dirichlet in this case), but when the time-derivatives are large then
the boundary approaches Neumann to reduce discontinuities.

By default, the relaxation rate is set to :math:`10` (i.e. a time-scale
of :math:`\tau=0.1`). To change this, give the rate as the second
argument e.g. ``relax(dirichlet, 2)`` would relax to a Dirichlet
boundary condition at a rate of :math:`2`.

Shifted boundaries
------------------

By default boundary conditions are applied in field-aligned coordinates,
where :math:`y` is along field-lines but :math:`x` has a discontinuity
at the twist-shift location. If radial derivatives are being done in
shifted coordinates where :math:`x` and :math:`z` are orthogonal, then
boundary conditions should also be applied in shifted coordinates. To do
this, the ``shifted`` boundary modifier applies a :math:`z` shift,
applies the boundary condition, then shifts back. For example::

    bndry_core = shifted( neumann )

would ensure that radial derivatives were zero in shifted coordinates on
the core boundary.

Changing the width of boundaries
--------------------------------

To change the width of a boundary region, the ``width`` modifier changes
the width of a boundary region before applying the boundary condition,
then changes the width back afterwards. To use, specify the boundary
condition and the width, for example

::

    bndry_core = width( neumann , 4 )

would apply a Neumann boundary condition on the innermost 4 cells in the
core, rather than the usual 2. When combining with other boundary
modifiers, this should be applied first e.g.

::

    bndry_sol = width( relax( dirichlet ), 3)

would relax the last 3 cells towards zero, whereas

::

    bndry_sol = relax( width( dirichlet, 3) )

would only apply to the usual 2, since relax didn’t use the updated
width.

Limitations:

#. Because it modifies then restores a globally-used BoundaryRegion,
   this code is not thread safe.

#. Boundary conditions can’t be applied across processors, and no checks
   are done that the width asked for fits within a single processor.

Examples
--------

This example is taken from the UEDGE benchmark test (in
``examples/uedge-benchmark``):

.. code-block:: cfg

    [All]
    bndry_all = neumann # Default for all variables, boundaries

    [Ni]
    bndry_target = neumann
    bndry_core = relax(dirichlet(1.))   # 1e13 cm^-3 on core boundary
    bndry_all  = relax(dirichlet(0.1))  # 1e12 cm^-3 on other boundaries

    [Vi]
    bndry_ydown = relax(dirichlet(-1.41648))   # -3.095e4/Vi_x
    bndry_yup   = relax(dirichlet( 1.41648))

The variable ``Ni`` (density) is set to a Neumann boundary condition on
the targets (yup and ydown), relaxes towards :math:`1` on the core
boundary, and relaxes to :math:`0.1` on all other boundaries. Note that
the ``bndry_target = neumann`` needs to be in the ``Ni`` section: If we
just had

.. code-block:: cfg

    [All]
    bndry_all = neumann # Default for all variables, boundaries

    [Ni]
    bndry_core = relax(dirichlet(1.))   # 1e13 cm^-3 on core boundary
    bndry_all  = relax(dirichlet(0.1))  # 1e12 cm^-3 on other boundaries

then the “target” boundary condition for ``Ni`` would first search in
the ``[Ni]`` section for ``bndry_target``, then for ``bndry_all`` in the
``[Ni]`` section. This is set to ``relax(dirichlet(0.1))``, not the
Neumann condition desired.

.. _sec-BoundaryRegion:

Boundary regions
----------------

The boundary condition code needs ways to loop over the boundary
regions, without needing to know the details of the mesh.

At the moment two mechanisms are provided: A RangeIterator over upper
and lower Y boundaries, and a vector of BoundaryRegion objects.

::

    // Boundary region iteration
    virtual const RangeIterator iterateBndryLowerY() const = 0;
    virtual const RangeIterator iterateBndryUpperY() const = 0;

    bool hasBndryLowerY();
    bool hasBndryUpperY();

    bool BoundaryOnCell; // NB: DOESN'T REALLY BELONG HERE

The `RangeIterator` class is an iterator which allows looping over a
set of indices. For example, in ``src/solver/solver.cxx`` to loop over
the upper Y boundary of a 2D variable ``var``::

    for(RangeIterator xi = mesh->iterateBndryUpperY(); !xi.isDone(); xi++) {
      ...
    }

The `BoundaryRegion` class is defined in
``include/boundary_region.hxx``

Boundary regions
----------------

Different regions of the boundary such as “core”, “sol” etc. are
labelled by the `Mesh` class (i.e. `BoutMesh`), which implements a
member function defined in ``mesh.hxx``::

      // Boundary regions
      virtual vector<BoundaryRegion*> getBoundaries() = 0;

This returns a vector of pointers to `BoundaryRegion` objects, each of
which describes a boundary region with a label, a ``BndryLoc``
location (i.e. inner x, outer x, lower y, upper y or all), and
iterator functions for looping over the points. This class is defined
in ``boundary_region.hxx``::

    /// Describes a region of the boundary, and a means of iterating over it
    class BoundaryRegion {
      public:
      BoundaryRegion();
      BoundaryRegion(const string &name, int xd, int yd);
      virtual ~BoundaryRegion();

      string label; // Label for this boundary region

      BndryLoc location; // Which side of the domain is it on?

      int x,y; // Indices of the point in the boundary
      int bx, by; // Direction of the boundary [x+dx][y+dy] is going outwards

      virtual void first() = 0;
      virtual void next() = 0; // Loop over every element from inside out (in X or
    Y first)
      virtual void nextX() = 0; // Just loop over X
      virtual void nextY() = 0; // Just loop over Y
      virtual bool isDone() = 0; // Returns true if outside domain. Can use this
    with nested nextX, nextY
    };

**Example:** To loop over all points in ``BoundaryRegion *bndry`` , use

::

      for(bndry->first(); !bndry->isDone(); bndry->next()) {
        ...
      }

Inside the loop, ``bndry->x`` and ``bndry->y`` are the indices of the
point, whilst ``bndry->bx`` and ``bndry->by`` are unit vectors out of
the domain. The loop is over all the points from the domain outwards
i.e. the point ``[bndry->x - bndry->bx][bndry->y - bndry->by]`` will
always be defined.

Sometimes it’s useful to be able to loop over just one direction along
the boundary. To do this, it is possible to use ``nextX()`` or
``nextY()`` rather than ``next()``. It is also possible to loop over
both dimensions using::

      for(bndry->first(); !bndry->isDone(); bndry->nextX())
        for(; !bndry->isDone(); bndry->nextY()) {
          ...
        }

Boundary operations
-------------------

On each boundary, conditions must be specified for each variable. The
different conditions are imposed by `BoundaryOp` objects. These set
the values in the boundary region such that they obey e.g. Dirichlet
or Neumann conditions. The `BoundaryOp` class is defined in
``boundary_op.hxx``::

    /// An operation on a boundary
    class BoundaryOp {
     public:
      BoundaryOp() {bndry = NULL;}
      BoundaryOp(BoundaryRegion *region)

      // Note: All methods must implement clone, except for modifiers (see below)
      virtual BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);

      /// Apply a boundary condition on field f
      virtual void apply(Field2D &f) = 0;
      virtual void apply(Field3D &f) = 0;

      virtual void apply(Vector2D &f);

      virtual void apply(Vector3D &f);

      /// Apply a boundary condition on ddt(f)
      virtual void apply_ddt(Field2D &f);
      virtual void apply_ddt(Field3D &f);
      virtual void apply_ddt(Vector2D &f);
      virtual void apply_ddt(Vector3D &f);

      BoundaryRegion *bndry;
    };

(where the implementations have been removed for clarity). Which has a
pointer to a `BoundaryRegion` object specifying which region this
boundary is operating on.

Boundary conditions need to be imposed on the initial conditions (after
`PhysicsModel::init`), and on the time-derivatives (after
`PhysicsModel::rhs`). The ``apply()`` functions are therefore called
during initialisation and given the evolving variables, whilst the
``apply_ddt`` functions are passed the time-derivatives.

To implement a boundary operation, as a minimum the ``apply(Field2D)``,
``apply(Field2D)`` and ``clone()`` need to be implemented: By default
the ``apply(Vector)`` will call the ``apply(Field)`` functions on each
component individually, and the ``apply_ddt()`` functions just call the
``apply()`` functions.

**Example**: Neumann boundary conditions are defined in
``boundary_standard.hxx``::

    /// Neumann (zero-gradient) boundary condition
    class BoundaryNeumann : public BoundaryOp {
     public:
      BoundaryNeumann() {}
     BoundaryNeumann(BoundaryRegion *region):BoundaryOp(region) { }
      BoundaryOp* clone(BoundaryRegion *region, const list<string> &args);
      void apply(Field2D &f);
      void apply(Field3D &f);
    };

and implemented in ``boundary_standard.cxx``

::

    void BoundaryNeumann::apply(Field2D &f) {
      // Loop over all elements and set equal to the next point in
      for(bndry->first(); !bndry->isDone(); bndry->next())
        f[bndry->x][bndry->y] = f[bndry->x - bndry->bx][bndry->y - bndry->by];
    }

    void BoundaryNeumann::apply(Field3D &f) {
      for(bndry->first(); !bndry->isDone(); bndry->next())
        for(int z=0;z<mesh->LocalNz;z++)
          f[bndry->x][bndry->y][z] = f[bndry->x - bndry->bx][bndry->y -
    bndry->by][z];
    }

This is all that’s needed in this case since there’s no difference
between applying Neumann conditions to a variable and to its
time-derivative, and Neumann conditions for vectors are just Neumann
conditions on each vector component.

To create a boundary condition, we need to give it a boundary region to
operate over::

    BoundaryRegion *bndry = ...
    BoundaryOp op = new BoundaryOp(bndry);

The ``clone`` function is used to create boundary operations given a
single object as a template in `BoundaryFactory`. This can take
additional arguments as a vector of strings - see explanation in
:ref:`sec-BoundaryFactory`.

Boundary modifiers
------------------

To create more complicated boundary conditions from simple ones (such
as Neumann conditions above), boundary operations can be modified by
wrapping them up in a `BoundaryModifier` object, defined in
``boundary_op.hxx``::

    class BoundaryModifier : public BoundaryOp {
     public:
      virtual BoundaryOp* clone(BoundaryOp *op, const list<string> &args) = 0;
     protected:
      BoundaryOp *op;
    };

Since `BoundaryModifier` inherits from `BoundaryOp`, modified boundary
operations are just a different boundary operation and can be treated
the same (Decorator pattern). Boundary modifiers could also be nested
inside each other to create even more complicated boundary
operations. Note that the ``clone`` function is different to the
`BoundaryOp` one: instead of a `BoundaryRegion` to operate on,
modifiers are passed a `BoundaryOp` to modify.

Currently the only modifier is `BoundaryRelax`, defined in
``boundary_standard.hxx``::

    /// Convert a boundary condition to a relaxing one
    class BoundaryRelax : public BoundaryModifier {
     public:
      BoundaryRelax(BoutReal rate) {r = fabs(rate);}
      BoundaryOp* clone(BoundaryOp *op, const list<string> &args);

      void apply(Field2D &f);
      void apply(Field3D &f);

      void apply_ddt(Field2D &f);
      void apply_ddt(Field3D &f);
     private:
      BoundaryRelax() {} // Must be initialised with a rate
      BoutReal r;
    };

.. _sec-BoundaryFactory:

Boundary factory
----------------

The boundary factory creates new boundary operations from input strings,
for example turning “relax(dirichlet)” into a relaxing Dirichlet
boundary operation on a given region. It is defined in
``boundary_factory.hxx`` as a Singleton, so to get a pointer to the
boundary factory use

::

      BoundaryFactory *bfact = BoundaryFactory::getInstance();

and to delete this singleton, free memory and clean-up at the end use::

      BoundaryFactory::cleanup();

Because users should be able to add new boundary conditions during
`PhysicsModel::init`, boundary conditions are not hard-wired into
`BoundaryFactory`. Instead, boundary conditions must be registered
with the factory, passing an instance which can later be cloned. This
is done in ``bout++.cxx`` for the standard boundary conditions::

      BoundaryFactory* bndry = BoundaryFactory::getInstance();
      bndry->add(new BoundaryDirichlet(), "dirichlet");
      ...
      bndry->addMod(new BoundaryRelax(10.), "relax");

where the ``add`` function adds BoundaryOp objects, whereas ``addMod``
adds `BoundaryModifier` objects. **Note**: The objects passed to
`BoundaryFactory` will be deleted when ``cleanup()`` is called.

When a boundary operation is added, it is given a name such as
“dirichlet”, and similarly for the modifiers (“relax” above). These
labels and object pointers are stored internally in `BoundaryFactory`
in maps defined in ``boundary_factory.hxx``::

      // Database of available boundary conditions and modifiers
      map<string, BoundaryOp*> opmap;
      map<string, BoundaryModifier*> modmap;

These are then used by `BoundaryFactory::create`::

      /// Create a boundary operation object
      BoundaryOp* create(const string &name, BoundaryRegion *region);
      BoundaryOp* create(const char* name, BoundaryRegion *region);

to turn a string such as “relax(dirichlet)” and a `BoundaryRegion`
pointer into a `BoundaryOp` object. These functions are implemented in
``boundary_factory.cxx``, starting around line 42. The parsing is done
recursively by matching the input string to one of:

-  ``modifier(<expression>, arg1, ...)``

-  ``modifier(<expression>)``

-  ``operation(arg1, ...)``

-  ``operation``

the ``<expression>`` variable is then resolved into a `BoundaryOp`
object by calling ``create(<expression>, region)``.

When an operator or modifier is found, it is created from the pointer
stored in the ``opmap`` or ``modmap`` maps using the ``clone`` method,
passing a ``list<string>`` reference containing any arguments. It’s up
to the operation implementation to ensure that the correct number of
arguments are passed, and to parse them into floats or other types.

**Example**: The Dirichlet boundary condition can take an optional
argument to change the value the boundary’s set to. In
``boundary_standard.cxx``::

    BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion *region, const list<string>
    &args) {
      if(!args.empty()) {
        // First argument should be a value
        stringstream ss;
        ss << args.front();

        BoutReal val;
        ss >> val;
        return new BoundaryDirichlet(region, val);
      }
      return new BoundaryDirichlet(region);
    }

If no arguments are passed i.e. the string was “dirichlet” or
“dirichlet()” then the ``args`` list is empty, and the default value
(0.0) is used. If one or more arguments is used then the first
argument is parsed into a `BoutReal` type and used to create a new
`BoundaryDirichlet` object. If more arguments are passed then these
are just ignored; probably a warning should be printed.

To set boundary conditions on a field, `FieldData` methods are defined
in ``field_data.hxx``::

    // Boundary conditions
      void setBoundary(const string &name); ///< Set the boundary conditions
      void setBoundary(const string &region, BoundaryOp *op); ///< Manually set
      virtual void applyBoundary() {}
      virtual void applyTDerivBoundary() {};
     protected:
      vector<BoundaryOp*> bndry_op; // Boundary conditions

The `FieldData::setBoundary` method is implemented in
``field_data.cxx``. It first gets a vector of pointers to
`BoundaryRegion`\ s from the mesh, then loops over these calling
`BoundaryFactory::createFromOptions` for each one and adding the
resulting boundary operator to the `FieldData::bndry_op` vector.
