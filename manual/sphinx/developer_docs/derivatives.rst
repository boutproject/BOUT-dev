.. _sec-derivatives:

Derivatives
===========

This is probably the part of the code most people will want to alter,
and is in ``bout++/src/sys/derivs.cxx``. The main task of this module is
to map functions on fields like ``DDX`` to direction-independent
differential methods on stencils such as :math:`4^{th}`-order central
differencing. This mapping depends on global settings in ``BOUT.inp``
and is illustrated in :numref:`fig-diffOverview`.

.. _fig-diffOverview:
.. figure:: ../figs/diffOverview.*
   :alt: Overview of ``derivs`` module

   Overview of ``derivs`` module, mapping derivative functions on fields
   to direction-independent differential methods

Four kinds of differencing methods are supported

#. | First derivative ``DDX``, ``DDY``, ``DDZ``
   | Central differencing type schemes for first-order derivatives

#. | Second derivatives ``D2DX2``, ``D2DZ2``, ``D2DZ2``
   | Central differencing second derivatives e.g. for :math:`\nabla^2`

#. | Upwinding ``VDDX``, ``VDDY``, ``VDDZ``
   | Terms like :math:`\mathbf{v}\cdot\nabla`

#. | Flux methods ``FDDX``, ``FDDY``, ``FDDZ``
   | Flux conserving, limiting methods for terms like
     :math:`\nabla\cdot\left(\mathbf{v}f\right)`

The differencing methods themselves are independent on direction, and
have types defined in :doc:`derivs.cxx<../_breathe_autogen/file/derivs_8cxx>`

::

    typedef BoutReal (*deriv_func)(stencil &); // f
    typedef BoutReal (*upwind_func)(stencil &, stencil &); // v, f

These operate on ``stencil`` objects. This class is in :doc:`stencils.hxx<../_breathe_autogen/file/stencils_8hxx>`

::

    class stencil {
      public:
        int jx, jy, jz;  // Central location
        BoutReal c, p, m, pp, mm; // stencil 2 each side of the centre
        Overloaded operators
          =,+,-,*,/
        Functions
          min, max, abs
    };

The main purpose of this class is to store a 5-element stencil. To
simplify some code this class also has a bunch of overloaded operators
on BoutReals and other stencil objects. There are also some functions to
calculate things like absolute, minimum, and maximum values.

Lookup tables
-------------

To convert between short variable names (“C2”), long descriptions
(“2nd order Central Differencing”), ``DIFF_METHOD`` enums used to
specify methods at runtime (DIFF\_C2, defined in
:doc:`bout_types.hxx<../_breathe_autogen/file/bout__types_8hxx>`), and
function pointers (``DDX_C2``), taking into account whether variables
are shifted or not, BOUT++ uses a set of lookup tables.

To find function pointers, tables of the following type are used:

::

    /// Translate between DIFF_METHOD codes, and functions
    struct DiffLookup {
      DIFF_METHOD method;
      deriv_func func;     // Single-argument differencing function
      upwind_func up_func; // Upwinding function
    };

Because the ``DiffLookup`` type contains a ``deriv_func`` and
``upwind_func`` pointer, it is used for all function lookup tables.
There is a separate table for each type of differencing method, so for
example the table of non-staggered upwinding methods is

::

    /// Upwinding functions lookup table
    static DiffLookup UpwindTable[] = { {DIFF_U1, NULL, VDDX_U1},
                        {DIFF_C2, NULL, VDDX_C2},
                        {DIFF_U4, NULL, VDDX_U4},
                        {DIFF_W3, NULL, VDDX_WENO3},
                        {DIFF_C4, NULL, VDDX_C4},
                        {DIFF_DEFAULT}};

The ``DIFF_DEFAULT`` at the end is used to terminate the array. These
tables are used by functions

::

    deriv_func lookupFunc(DiffLookup* table, DIFF_METHOD method);
    upwind_func lookupUpwindFunc(DiffLookup* table, DIFF_METHOD method);

which return the function pointer corresponding to the given method. If
the method isn’t in the table, then the first entry in the table is
used. These functions can be used at run-time to allow a user to specify
the method to use for specific operators.

When reading settings from the input file, they are specified as short
strings like “C2”, and a longer description of the method chosen should
be written to the output log. To do this, there is a name lookup table:

::

    /// Translate between short names, long names and DIFF_METHOD codes
    struct DiffNameLookup {
      DIFF_METHOD method;
      const char* label; // Short name
      const char* name;  // Long name
    };

    static DiffNameLookup DiffNameTable[] = { 
      {DIFF_U1, "U1", "First order upwinding"},
      {DIFF_C2, "C2", "Second order central"},
      {DIFF_W2, "W2", "Second order WENO"},
      {DIFF_W3, "W3", "Third order WENO"},
      {DIFF_C4, "C4", "Fourth order central"},
      {DIFF_U4, "U4", "Fourth order upwinding"},
      {DIFF_FFT, "FFT", "FFT"},
      {DIFF_DEFAULT}}; // Use to terminate the list

To search this table, there is the function

::

    DIFF_METHOD lookupFunc(DiffLookup *table, const string &label)

During initialisation, the lookup therefore works in two stages, shown
in :numref:`fig-diffLookup`. First the short description is turned into a
``DIFF_METHOD`` enum code, then this code is turned into a function
pointer.

.. _fig-diffLookup:
.. figure:: ../figs/diffLookup.*
   :alt: Lookup tables for differential method

   Lookup tables for mapping between differential method labels, codes,
   descriptions and function pointers

Staggered grids
---------------

By default, all quantities in BOUT++ are defined at cell centre, and all
derivative methods map cell-centred quantities to cell centres.
Switching on staggered grid support in BOUT.inp:

::

    StaggerGrids = true

allows quantities to be defined on cell boundaries. Functions such as
``DDX`` now have to handle all possible combinations of input and output
locations, in addition to the possible derivative methods.

Several things are not currently implemented, which probably should be:

-  Only 3D fields currently have a cell location attribute. The location
   (cell centre etc) of 2D fields is ignored at the moment. The
   rationale for this is that 2D fields are assumed to be slowly-varying
   equilibrium quantities for which it won’t matter so much. Still,
   needs to be improved in future

-  Twist-shift and X shifting still treat all quantities as
   cell-centred.

-  No boundary condition functions yet account for cell location.

Currently, BOUT++ does not support values at cell corners; values can
only be defined at cell centre, or at the lower X,Y, or Z boundaries.
This is

Once staggered grids are enabled, two types of stencil are needed: those
which map between the same cell location (e.g. cell-centred values to
cell-centred values), and those which map to different locations (e.g.
cell-centred to lower X).

.. figure:: ../figs/diffStencils.*
   :alt: Stencils with cell-centred and lower shifted values

   Stencils with cell-centred (solid) and lower shifted values (open).
   Processor boundaries marked by vertical dashed line

Central differencing using 4-point stencil:

.. math::

   \begin{aligned}
   y &=& \left(9y_{-1/2} + 9y_{1/2} - y_{-3/2} - y_{3/2}\right) / 16 \\
   {{\frac{\partial y}{\partial x}}} &=& \left( 27y_{1/2} - 27y_{-1/2} - y_{3/2} + y_{-3/2}\right) / 24\Delta x \\
   \frac{\partial^2 y}{\partial x^2} &=& \left(y_{3/2} + y_{-3/2} - y_{1/2} - y_{-1/2}\right) / 2\Delta x^2\end{aligned}

+----------+-------------------+----------------------------------------------------------------+
| Input    | Output            | Actions                                                        |
+==========+===================+================================================================+
|          | Central stencil   |                                                                |
+----------+-------------------+----------------------------------------------------------------+
| CENTRE   | XLOW              | Lower staggered stencil                                        |
+----------+-------------------+----------------------------------------------------------------+
| XLOW     | CENTRE            | Upper staggered stencil                                        |
+----------+-------------------+----------------------------------------------------------------+
| XLOW     | Any               | Staggered stencil to CENTRE, then interpolate                  |
+----------+-------------------+----------------------------------------------------------------+
| CENTRE   | Any               | Central stencil, then interpolate                              |
+----------+-------------------+----------------------------------------------------------------+
| Any      | Any               | Interpolate to centre, use central stencil, then interpolate   |
+----------+-------------------+----------------------------------------------------------------+

Table: DDX actions depending on input and output locations. Uses first
match.

