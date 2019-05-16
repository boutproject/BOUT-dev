.. _sec-laplacian:

Laplacian inversion
===================

A common problem in plasma models is to solve an equation of the form

.. math::
   :label: full_laplace_inv

   d\nabla^2_\perp x + \frac{1}{c_1}(\nabla_\perp c_2)\cdot\nabla_\perp x +
   a x = b

For example,

.. math::

   \nabla_\perp^2 x + a x = b

appears in reduced MHD for the vorticity inversion and :math:`j_{||}`.

Alternative formulations and ways to invert equation
:eq:`full_laplace_inv` can be found in section :ref:`sec-LaplaceXY` and
:ref:`sec-LaplaceXZ`

Several implementations of the Laplacian solver are available, which
are selected by changing the "type" setting.The currently available
implementations are listed in table :numref:`tab-laplacetypes`.

.. _tab-laplacetypes:
.. table:: Laplacian implementation types

   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | Name                   | Description                                                  | Requirements                             |
   +========================+==============================================================+==========================================+
   | cyclic                 | Serial/parallel. Gathers boundary rows onto one processor.   |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | `petsc                 | Serial/parallel. Lots of methods, no Boussinesq              | PETSc (section :ref:`sec-PETSc-install`) |
   | <sec-petsc-laplace_>`__|                                                              |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | multigrid              | Serial/parallel. Geometric multigrid, no Boussinesq          |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | `naulin                | Serial/parallel. Iterative treatment of non-Boussinesq terms |                                          |
   | <sec-naulin_>`__       |                                                              |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | `serial_tri            | Serial only. Thomas algorithm for tridiagonal system.        | Lapack (section :ref:`sec-lapack`)       |
   | <sec-tri_>`__          |                                                              |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | `serial_band           | Serial only. Enables 4th-order accuracy                      | Lapack (section :ref:`sec-lapack`)       |
   | <sec-band_>`__         |                                                              |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | `spt                   | Parallel only (NXPE>1). Thomas algorithm.                    |                                          |
   | <sec-spt_>`__          |                                                              |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | mumps                  | Serial/parallel. Direct solver                               | MUMPS (section :ref:`sec-mumps`)         |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | `pdd                   | Parallel Diagnonally Dominant algorithm. Experimental        |                                          |
   | <sec-pdd_>`__          |                                                              |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+
   | shoot                  | Shooting method. Experimental                                |                                          |
   +------------------------+--------------------------------------------------------------+------------------------------------------+

Usage of the laplacian inversion
--------------------------------

In BOUT++, equation :eq:`full_laplace_inv` can be solved in two
ways. The first method Fourier transforms in the :math:`z`-direction,
whilst the other is solving the full two dimensional problem by matrix
inversion. The derivation of :math:`\nabla_\perp^2f` for a general
coordinate system can be found in the ``coordinates`` manual. What is
important, is to note that if :math:`g_{xy}` and :math:`g_{yz}` are
non-zero, BOUT++ is neglecting the :math:`y`-parallel derivatives when
using the solvers `Laplacian` and `LaplaceXZ`.

By neglecting the :math:`y`-derivatives (or if
:math:`g_{xy}=g_{yz}=0`), one can solve equation
:eq:`full_laplace_inv` :math:`y` plane by :math:`y` plane.

The first approach utilizes that it is possible Fourier transform the
equation in :math:`z` (using some assumptions described in section
:ref:`sec-num-laplace`), and solve a tridiagonal system for each
mode. These inversion problems are band-diagonal (tri-diagonal in the
case of 2nd-order differencing) and so inversions can be very
efficient: :math:`O(n_z \log n_z)` for the FFTs, :math:`O(n_x)` for
tridiagonal inversion using the Thomas algorithm, where :math:`n_x`
and :math:`n_z` are the number of grid-points in the :math:`x` and
:math:`z` directions respectively.


In the second approach, the full :math:`2`\ -D system is being solved.
This requires PETSc to be built with BOUT++.


The `Laplacian` class is defined in ``invert_laplace.hxx`` and solves
problems formulated like equation :eq:`full_laplace_inv` To use this
class, first create an instance of it::

    Laplacian *lap = Laplacian::create();

By default, this will use the options in a section called “laplace”, but
can be given a different section as an argument. By default
:math:`d = 1`, :math:`a = 0`, and the :math:`c=1`. To set the values of
these coefficients, there are the ``setCoefA()``, ``setCoefC()``, and
``setCoefD()`` methods::

    Field2D a = ...;
    lap->setCoefA(a);
    lap->setCoefC(0.5);

arguments can be `Field2D`, `Field3D`, or `BoutReal` values.

Settings for the inversion can be set in the input file under the
section ``laplace`` (default) or whichever settings section name was
specified when the `Laplacian` class was created. Commonly used
settings are listed in tables :numref:`tab-laplacesettings` to
:numref:`tab-laplaceflags`.

In particular boundary conditions on the :math:`x` boundaries can be
set using the and ``outer_boundary_flags`` variables, as detailed in
table :numref:`tab-laplaceBCflags`. Note that DC (‘direct-current’)
refers to :math:`k = 0` Fourier component, AC (‘alternating-current’)
refers to :math:`k \neq 0` Fourier components. Non-Fourier solvers use
AC options (and ignore DC ones). Multiple boundary conditions can be
selected by adding together the required boundary condition flag
values together. For example, ``inner_boundary_flags = 3`` will set a
Neumann boundary condition on both AC and DC components.

It is pertinent to note here that the boundary in BOUT++ is defined by
default to be located half way between the first guard point and first
point inside the domain. For example, when a Dirichlet boundary
condition is set, using ``inner_boundary_flags = 0`` , ``16``, or
``32``, then the first guard point, :math:`f_{-}` will be set to
:math:`f_{-} = 2v - f_+`, where :math:`f_+` is the first grid point
inside the domain, and :math:`v` is the value to which the boundary is
being set to.

The ``global_flags``, ``inner_boundary_flags``,
``outer_boundary_flags`` and ``flags`` values can also be set from
within the physics module using ``setGlobalFlags``,
``setInnerBoundaryFlags`` , ``setOuterBoundaryFlags`` and
``setFlags``.

::

    lap->setGlobalFlags(Global_Flags_Value);
    lap->setInnerBoundaryFlags(Inner_Flags_Value);
    lap->setOuterBoundaryFlags(Outer_Flags_Value);
    lap->setFlags(Flags_Value);


.. _tab-laplacesettings:
.. table:: Laplacian inversion options

   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | Name                     | Meaning                                                                 | Default value                                |
   +==========================+=========================================================================+==============================================+
   | ``type``                 | Which implementation to use. See table :numref:`tab-laplacetypes`       | ``cyclic``                                   |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``filter``               | Filter out modes above :math:`(1-`\ ``filter``\                         | 0                                            |
   |                          | :math:`)\times k_{max}`, if using Fourier solver                        |                                              |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``maxmode``              | Filter modes with :math:`n >`\ ``maxmode``                              | ``MZ``/2                                     |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``all_terms``            | Include first derivative terms                                          | ``true``                                     |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``nonuniform``           | Include                                                                 | Same as global ``non_uniform``.              |
   |                          | :ref:`corrections for non-uniform meshes <sec-diffmethod-nonuniform>`   | See :ref:`here <sec-diffmethod-nonuniform>`  |
   |                          | (dx not constant)                                                       |                                              |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``global_flags``         | Sets global inversion options See table                                 | ``0``                                        |
   |                          | :ref:`Laplace global flags<tab-laplaceglobalflags>`                     |                                              |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``inner_boundary_flags`` | Sets boundary conditions on inner boundary. See table                   | ``0``                                        |
   |                          | :ref:`Laplace boundary flags<tab-laplaceBCflags>`                       |                                              |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``outer_boundary_flags`` | Sets boundary conditions on outer boundary. See table                   | ``0``                                        |
   |                          | :ref:`Laplace boundary flags<tab-laplaceBCflags>`                       |                                              |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``flags``                | DEPRECATED. Sets global solver options and boundary                     | ``0``                                        |
   |                          | conditions. See :ref:`Laplace flags<tab-laplaceflags>` or               |                                              |
   |                          | :doc:`invert_laplace.cxx<../_breathe_autogen/file/invert__laplace_8cxx>`|                                              |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+
   | ``include_yguards``      | Perform inversion in :math:`y`\ -boundary guard cells                   | ``false``                                    |
   +--------------------------+-------------------------------------------------------------------------+----------------------------------------------+

|

.. _tab-laplaceglobalflags:
.. table:: Laplacian inversion ``global_flags`` values: add the required quantities together.

   +--------+--------------------------------------------------------------------------------+-----------------------------+
   | Flag   | Meaning                                                                        | Code variable               |
   +========+================================================================================+=============================+
   | 0      | No global option set                                                           | :math:`-`                   |
   +--------+--------------------------------------------------------------------------------+-----------------------------+
   | 1      | zero DC component (Fourier solvers)                                            | ``INVERT_ZERO_DC``          |
   +--------+--------------------------------------------------------------------------------+-----------------------------+
   | 2      | set initial guess to 0 (iterative solvers)                                     | ``INVERT_START_NEW``        |
   +--------+--------------------------------------------------------------------------------+-----------------------------+
   | 4      | equivalent to                                                                  | ``INVERT_BOTH_BNDRY_ONE``   |
   |        | ``outer_boundary_flags = 128``,                                                |                             |
   |        | ``inner_boundary_flags = 128``                                                 |                             |
   +--------+--------------------------------------------------------------------------------+-----------------------------+
   | 8      | Use 4th order differencing (Apparently not actually implemented anywhere!!!)   | ``INVERT_4TH_ORDER``        |
   +--------+--------------------------------------------------------------------------------+-----------------------------+
   | 16     | Set constant component (:math:`k_x = k_z = 0`) to zero                         | ``INVERT_KX_ZERO``          |
   +--------+--------------------------------------------------------------------------------+-----------------------------+

|

.. _tab-laplaceBCflags:
.. table:: Laplacian inversion ``outer_boundary_flags`` or ``inner_boundary_flags`` values: add the required quantities together.

   +--------+----------------------------------------------------------------------+----------------------------+
   | Flag   | Meaning                                                              | Code variable              |
   +========+======================================================================+============================+
   | 0      | Dirichlet (Set boundary to 0)                                        | :math:`-`                  |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 1      | Neumann on DC component (set gradient to 0)                          | ``INVERT_DC_GRAD``         |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 2      | Neumann on AC component (set gradient to 0)                          | ``INVERT_AC_GRAD``         |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 4      | Zero or decaying Laplacian on AC components (                        | ``INVERT_AC_LAP``          |
   |        | :math:`\frac{\partial^2}{\partial x^2}+k_z^2` vanishes/decays)       |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 8      | Use symmetry to enforce zero value or gradient (redundant for 2nd    | ``INVERT_SYM``             |
   |        | order now)                                                           |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 16     | Set boundary condition to values in boundary guard cells of second   | ``INVERT_SET``             |
   |        | argument, ``x0``, of `Laplacian::solve(const Field3D &b,             |                            |
   |        | const Field3D &x0) <Laplacian::solve>`. May be combined with any     |                            |
   |        | combination of 0, 1 and 2, i.e. a Dirichlet or Neumann boundary      |                            |
   |        | condition set to values which are :math:`\neq 0` or :math:`f(y)`     |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 32     | Set boundary condition to values in boundary guard cells of RHS,     | ``INVERT_RHS``             |
   |        | ``b`` in `Laplacian::solve(const Field3D &b, const Field3D           |                            |
   |        | &x0) <Laplacian::solve>`. May be combined with any combination of 0, |                            |
   |        | 1 and 2, i.e. a Dirichlet or Neumann boundary condition set to values|                            |
   |        | which are :math:`\neq 0` or :math:`f(y)`                             |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 64     | Zero or decaying Laplacian on DC components                          | ``INVERT_DC_LAP``          |
   |        | (:math:`\frac{\partial^2}{\partial x^2}` vanishes/decays)            |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 128    | Assert that there is only one guard cell in the :math:`x`-boundary   | ``INVERT_BNDRY_ONE``       |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 256    | DC value is set to parallel gradient, :math:`\nabla_\parallel f`     | ``INVERT_DC_GRADPAR``      |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 512    | DC value is set to inverse of parallel gradient                      | ``INVERT_DC_GRADPARINV``   |
   |        | :math:`1/\nabla_\parallel f`                                         |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 1024   | Boundary condition for inner ‘boundary’ of cylinder                  | ``INVERT_IN_CYLINDER``     |
   +--------+----------------------------------------------------------------------+----------------------------+

|

.. _tab-laplaceflags:
.. table:: Laplacian inversion ``flags`` values (DEPRECATED!): add the required quantities together.

   +--------+------------------------------------------------------------------------------------------+
   | Flag   | Meaning                                                                                  |
   +========+==========================================================================================+
   | 1      | Zero-gradient DC on inner (X) boundary. Default is zero-value                            |
   +--------+------------------------------------------------------------------------------------------+
   | 2      | Zero-gradient AC on inner boundary                                                       |
   +--------+------------------------------------------------------------------------------------------+
   | 4      | Zero-gradient DC on outer boundary                                                       |
   +--------+------------------------------------------------------------------------------------------+
   | 8      | Zero-gradient AC on outer boundary                                                       |
   +--------+------------------------------------------------------------------------------------------+
   | 16     | Zero DC component everywhere                                                             |
   +--------+------------------------------------------------------------------------------------------+
   | 32     | Not used currently                                                                       |
   +--------+------------------------------------------------------------------------------------------+
   | 64     | Set width of boundary to 1 (default is ``MXG``)                                          |
   +--------+------------------------------------------------------------------------------------------+
   | 128    | Use 4\ :math:`^{th}`-order band solver (default is 2\ :math:`^{nd}` order tridiagonal)   |
   +--------+------------------------------------------------------------------------------------------+
   | 256    | Attempt to set zero laplacian AC component on inner boundary by combining                |
   |        | 2nd and 4th-order differencing at the boundary.                                          |
   |        | Ignored if tridiagonal solver used.                                                      |
   +--------+------------------------------------------------------------------------------------------+
   | 512    | Zero laplacian AC on outer boundary                                                      |
   +--------+------------------------------------------------------------------------------------------+
   | 1024   | Symmetric boundary condition on inner boundary                                           |
   +--------+------------------------------------------------------------------------------------------+
   | 2048   | Symmetric outer boundary condition                                                       |
   +--------+------------------------------------------------------------------------------------------+

To perform the inversion, there’s the ``solve`` method

::

    x = lap->solve(b);

If you prefer, there are functions compatible with older versions of the
BOUT++ code::

    Field2D a, c, d;
    invert_laplace(b, x, flags, &a, &c, &d);

and

::

    x = invert_laplace(b, flags, &a, &c, &d);

The input ``b`` and output ``x`` are 3D fields, and the coefficients
``a``, ``c``, and ``d`` are pointers to 2D fields. To omit any of the
three coefficients, set them to NULL.

.. _sec-num-laplace:

Numerical implementation
------------------------

We will here go through the implementation of the laplacian inversion
algorithm, as it is performed in BOUT++. We would like to solve the
following equation for :math:`f`

.. math::
   :label: to_invert

   d\nabla_\perp^2f + \frac{1}{c_1}(\nabla_\perp c_2)\cdot\nabla_\perp f + af = b

BOUT++ is neglecting the :math:`y`-parallel derivatives if
:math:`g_{xy}` and :math:`g_{yz}` are no-zero when using the solvers
`Laplacian` and `LaplaceXZ`. For these two
solvers, equation :eq:`to_invert` becomes (see ``coordinates`` manual
for derivation)

.. math::
   :label: invert_expanded

   \, &d (g^{xx} \partial_x^2 + G^x \partial_x + g^{zz} \partial_z^2 +
   G^z \partial_z + 2g^{xz} \partial_x \partial_z ) f \\
   +& \frac{1}{c_1}( {{\boldsymbol{e}}}^x \partial_x +
   {\boldsymbol{e}}^z \partial_z ) c_2 \cdot ({\boldsymbol{e}}^x
   \partial_x + {\boldsymbol{e}}^z \partial_z ) f \\ +& af = b


Using tridiagonal solvers
~~~~~~~~~~~~~~~~~~~~~~~~~

When using the tridiagonal solvers, :math:`c_1 = c_2` in equation
:eq:`to_invert`, hence, it is rather solving

.. math::
   :label: to_invert_tri

   d\nabla_\perp^2f + \frac{1}{c}(\nabla_\perp c)\cdot\nabla_\perp f + af = b

Since there are no parallel :math:`y`-derivatives if
:math:`g_{xy}=g_{yz}=0` (or if they are neglected), equation
:eq:`to_invert_tri` will only contain derivatives of :math:`x` and
:math:`z` for the dependent variable. The hope is that the modes in the
periodic :math:`z` direction will decouple, so that we in the end only
have to invert for the :math:`x` coordinate.

If the modes decouples when Fourier transforming equation
:eq:`invert_expanded`, we can use a tridiagonal solver to solve the
equation for each Fourier mode.

Using the discrete Fourier transform

.. math::

   F(x,y)_{k} = \frac{1}{N}\sum_{Z=0}^{N-1}f(x,y)_{Z}\exp(\frac{-2\pi i k
   Z}{N})

we see that the modes will not decouple if a term consist of a product
of two terms which depends on :math:`z`, as this would give terms like

.. math::

   \frac{1}{N}\sum_{Z=0}^{N-1} a(x,y)_Z f(x,y)_Z \exp(\frac{-2\pi i k
   Z}{N})

Thus, in order to use a tridiagonal solver, :math:`a`, :math:`c` and
:math:`d` cannot be functions of :math:`z`. Because of this, the
:math:`{{\boldsymbol{e}}}^z \partial_z c` term in equation
:eq:`invert_expanded` is zero. In principle the modes would still
decouple if the :math:`{{\boldsymbol{e}}}^z \partial_z f`
part of equation :eq:`invert_expanded` was kept, but currently this
part is also neglected in solvers using a tridiagonal matrix. Thus the
tridiagonal solvers are solving equations on the form

.. math::

   \, &d(x,y) ( g^{xx}(x,y) \partial_x^2 + G^x(x,y) \partial_x +
       g^{zz}(x,y) \partial_z^2 + G^z(x,y) \partial_z + 2g^{xz}(x,y)
       \partial_x \partial_z ) f(x,y,z) \\
     +& \frac{1}{c(x,y)}({{\boldsymbol{e}}}^x \partial_x ) c(x,y) \cdot (
       {{\boldsymbol{e}}}^x \partial_x ) f(x,y,z) \\
     +& a(x,y)f(x,y,z) = b(x,y,z)

after using the discrete Fourier transform (see section
:ref:`sec-derivatives-of-fft`), we get

.. math::

   \, &d (    g^{xx} \partial_x^2F_z + G^x \partial_xF_z + g^{zz} [i k]^2F_z
        + G^z [i k]F_z + 2g^{xz} \partial_x[i k]F_z ) \\
     +& \frac{1}{c}( {{\boldsymbol{e}}}^x \partial_x ) c \cdot ( {{\boldsymbol{e}}}^x
        \partial_xF_z ) \\
     +& aF_z = B_z

which gives

.. math::
   :label: FT_laplace_inversion

   \, &d ( g^{xx} \partial_x^2 + G^x \partial_x - k^2 g^{zz} + i
   kG^z + i k2g^{xz} \partial_x )F_z \\
   +& \frac{g^{xx}}{c} (\partial_x c ) \partial_xF_z \\
   +& aF_z = B_z

As nothing in equation :eq:`FT_laplace_inversion` couples points in
:math:`y` together (since we neglected the :math:`y`-derivatives if
:math:`g_{xy}` and :math:`g_{yz}` were non-zero). Also, as the modes are
decoupled, we may solve equation :eq:`FT_laplace_inversion` :math:`k`
mode by :math:`k` mode in addition to :math:`y`\ -plane by
:math:`y`\ -plane.

The second order centred approximation of the first and second
derivatives in :math:`x` reads

.. math::

   \partial_x f &\simeq \frac{-f_{n-1} + f_{n+1}}{2\text{d}x} \\
   \partial_x^2 f &\simeq \frac{f_{n-1} - f_{n} + f_{n+1}}{\text{d}x^2}

This gives

.. math::

   \, &d (    g^{xx} \frac{F_{z,n-1} - 2F_{z,n} + F_{z, n+1}}{\text{d}x^2} +
       G^x \frac{-F_{z,n-1} + F_{z,n+1}}{2\text{d}x} - k^2 g^{zz}F_{z,n} .\\
       &\quad.  + i kG^zF_{z,n} + i k2g^{xz} \frac{-F_{z,n-1} +
   F_{z,n+1}}{2\text{d}x} ) \\
       +& \frac{g^{xx}}{c} ( \frac{-c_{n-1} + c_{n+1}}{2\text{d}x} )
   \frac{-F_{z,n-1} + F_{z,n+1}}{2\text{d}x} \\
       +& aF_{z,n} = B_{z,n}

collecting point by point

.. math::
   :label: discretized_laplace

       &( \frac{dg^{xx}}{\text{d}x^2} - \frac{dG^x}{2\text{d}x} -
       \frac{g^{xx}}{c_{n}} \frac{-c_{n-1} + c_{n+1}}{4\text{d}x^2} - i\frac{d
       k2g^{xz}}{2\text{d}x} ) F_{z,n-1} \\
           +&( - \frac{ dg^{xx} }{\text{d}x^2} - dk^2 g^{zz} + a + idkG^z )
       F_{z,n} \\
           +&( \frac{dg^{xx}}{\text{d}x^2} + \frac{dG^x}{2\text{d}x} +
       \frac{g^{xx}}{c_{n}} \frac{-c_{n-1} + c_{n+1}}{4\text{d}x^2} +
       i\frac{dk2g^{xz}}{2\text{d}x} ) F_{z, n+1} \\
        =& B_{z,n}

We now introduce

.. math::

   c_1 = \frac{dg^{xx}}{\text{d}x^2}

   c_2 = dg^{zz}

   c_3 = \frac{2dg^{xz}}{2\text{d}x}

   c_4 = \frac{dG^x + g^{xx}\frac{-c_{n-1} + c_{n+1}}{2c_n\text{d}x}}{2\text{d}x}

   c_5 = dG^z

which inserted in equation :eq:`discretized_laplace` gives

.. math::

       ( c_1 - c_4 -ikc_3 ) F_{z,n-1} + ( -2c_1 - k^2c_2 +ikc_5 + a ) F_{z,n} + ( c_1 + c_4 + ikc_3 ) F_{z, n+1} = B_{z,n}

This can be formulated as the matrix equation

.. math::

   AF_z=B_z

where the matrix :math:`A` is tridiagonal. The boundary conditions are
set by setting the first and last rows in :math:`A` and :math:`B_z`.

.. _sec-petsc-laplace:

Using PETSc solvers
~~~~~~~~~~~~~~~~~~~

When using PETSc, all terms of equation :eq:`invert_expanded` is being
used when inverting to find :math:`f`. Note that when using PETSc, we
are not Fourier decomposing in the :math:`z`-direction, so it may take
substantially longer time to find the solution. As with the tridiagonal
solver, the fields are being sliced in the :math:`y`-direction, and a
solution is being found for one :math:`y` plane at the time.

Before solving, equation :eq:`invert_expanded` is rewritten to the
form
:math:`A{{\boldsymbol{x}}} ={{\boldsymbol{b}}}`
(however, the full :math:`A` is not expanded in memory). To do this, a
row :math:`i` in the matrix :math:`A` is indexed from bottom left of the
two dimensional field :math:`= (0,0) = 0` to top right
:math:`= (\texttt{meshx}-1,
\texttt{meshz}-1) = \texttt{meshx}\cdot\texttt{meshz}-1` of the two
dimensional field. This is done in such a way so that a row :math:`i` in
:math:`A` increments by :math:`1` for an increase of :math:`1` in the
:math:`z-`\ direction, and by :math:`\texttt{meshz}` for an increase of
:math:`1` in the :math:`x-`\ direction, where the variables
:math:`\texttt{meshx}` and :math:`\texttt{meshz}` represents the highest
value of the field in the given direction.

Similarly to equation :eq:`discretized_laplace`, the discretised
version of equation :eq:`invert_expanded` can be written. Doing the
same for the full two dimensional case yields:

Second order approximation

.. math::

       \; & c_{i,j} f_{i,j} \\
           &+ c_{i-1,j-1} f_{i-1,j-1} + c_{i-1,j} f_{i-1,j} \\
           &+ c_{i-1,j+1} f_{i-1,j+1} + c_{i,j-1} f_{i,j-1} \\
           &+ c_{i,j+1} f_{i,j+1} + c_{i+1,j-1} f_{i+1,j-1} \\
           &+ c_{i+1,j} f_{i+1,j} + c_{i+1,j+1} f_{i+1,j+1} \\
       =& b_{i,j}

Fourth order approximation

.. math::

       \; & c_{i,j} f_{i,j} \\
           &+ c_{i-2,j-2} f_{i-2,j-2} + c_{i-2,j-1} f_{i-2,j-1} \\
           &+ c_{i-2,j} f_{i-2,j} + c_{i-2,j+1} f_{i-2,j+1} \\
           &+ c_{i-2,j+2} f_{i-2,j+2} + c_{i-1,j-2} f_{i-1,j-2} \\
           &+ c_{i-1,j-1} f_{i-1,j-1} + c_{i-1,j} f_{i-1,j} \\
           &+ c_{i-1,j-1} f_{i-1,j-1} + c_{i-1,j} f_{i-1,j} \\
           &+ c_{i-1,j+1} f_{i-1,j+1} + c_{i-1,j+2} f_{i-1,j+2} \\
           &+ c_{i,j-2} f_{i,j-2} + c_{i,j-1} f_{i,j-1} \\
           &+ c_{i,j+1} f_{i,j+1} + c_{i,j+2} f_{i,j+2} \\
           &+ c_{i+1,j-2} f_{i+1,j-2} + c_{i+1,j-1} f_{i+1,j-1} \\
           &+ c_{i+1,j} f_{i+1,j} + c_{i+1,j+1} f_{i+1,j+1} \\
           &+ c_{i+1,j+2} f_{i+1,j+2} + c_{i+2,j-2} f_{i+2,j-2} \\
           &+ c_{i+2,j-1} f_{i+2,j-1} + c_{i+2,j} f_{i+2,j} \\
           &+ c_{i+2,j+1} f_{i+2,j+1} + c_{i+2,j+2} f_{i+2,j+2} \\
       =& b_{i,j}


To determine the coefficient for each node point, it is convenient to
introduce some quantities

.. math::
   :nowrap:

      \begin{align}
       &A_0 = a(x,y_{\text{current}},z)& &A_1 = dg^{xx}&\\
       &A_2 = dg^{zz}& &A_3 = 2dg^{xz}&
      \end{align}

In addition, we have:

Second order approximation (5-point stencil)

.. math::

       \texttt{ddx\_c} = \frac{\texttt{c2}_{x+1} - \texttt{c2}_{x-1} }{2\texttt{c1}\text{d}x}

       \texttt{ddz\_c} = \frac{\texttt{c2}_{z+1} - \texttt{c2}_{z-1} }{2\texttt{c1}\text{d}z}

Fourth order approximation (9-point stencil)

.. math::

       \texttt{ddx\_c} = \frac{-\texttt{c2}_{x+2} + 8\texttt{c2}_{x+1} -
       8\texttt{c2}_{x-1} + \texttt{c2}_{x-1} }{ 12\texttt{c1}\text{d}x} \\
       \texttt{ddz\_c} = \frac{-\texttt{c2}_{z+2} + 8\texttt{c2}_{z+1} -
       8\texttt{c2}_{z-1} + \texttt{c2}_{z-1} }{ 12\texttt{c1}\text{d}z}


This gives

.. math::
   A_4 = dG^x + g^{xx}\texttt{ddx\_c} + g^{xz}\texttt{ddz\_c}
   A_5 = dG^z + g^{xz}\texttt{ddx\_c} + g^{xx}\texttt{ddz\_c}

The coefficients :math:`c_{i+m,j+n}` are finally being set according
to the appropriate order of discretisation. The coefficients can be
found in the file ``petsc_laplace.cxx``.

Example: The 5-point stencil
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us now consider the 5-point stencil for a mesh with :math:`3` inner
points in the :math:`x`-direction, and :math:`3` inner points in the
:math:`z`-direction. The :math:`z` direction will be periodic, and the
:math:`x` direction will have the boundaries half between the grid-point
and the first ghost point (see :numref:`fig-lapl-inv-mesh`).

.. _fig-lapl-inv-mesh:
.. figure:: ../figs/5PointStencilMesh.*
   :alt: The mesh

   The mesh: The inner boundary points in :math:`x` are coloured in
   orange, whilst the outer boundary points in :math:`z` are coloured
   gray. Inner points are coloured blue.

Applying the :math:`5`-point stencil to point :math:`f_{22}` this mesh
will result in :numref:`fig-lapl-inv-mesh-w-stencil`.

.. _fig-lapl-inv-mesh-w-stencil:
.. figure:: ../figs/5PointStencilMeshWithStencil.*
   :alt: The 5-point stencil for the Laplacian

   The mesh with a stencil in point :math:`f_{22}`: The point under
   consideration is coloured blue. The point located :math:`+1` in
   :math:`z` direction (``zp``) is coloured yellow and the point
   located :math:`-1` in :math:`z` direction (``zm``) is coloured
   green. The point located :math:`+1` in :math:`x` direction (``xp``)
   is coloured purple and the point located :math:`-1` in :math:`x`
   direction (``xm``) is coloured red.

We want to solve a problem on the form
:math:`A{{\mathbf{x}}}={{\mathbf{b}}}`. We will order
:math:`{{\mathbf{x}}}` in a row-major order (so that :math:`z` is
varying faster than :math:`x`). Further, we put the inner :math:`x`
boundary points first in :math:`{{\mathbf{x}}}`, and the outer
:math:`x` boundary points last in :math:`{{\mathbf{x}}}`. The matrix
problem for our mesh can then be written like in
:numref:`fig-lapl-inv-matrix`.

.. _fig-lapl-inv-matrix:
.. figure:: ../figs/5PointStencilMatrix.*
   :alt: The matrix problem for the Laplacian inversion

   Matrix problem for our :math:`3\times3` mesh: The colors follow
   that of figure :numref:`fig-lapl-inv-mesh` and
   :numref:`fig-lapl-inv-mesh-w-stencil`.  The first index of the
   elements refers to the :math:`x`-position in figure
   :numref:`fig-lapl-inv-mesh`, and the last index of the elements
   refers to the :math:`z`-position in figure
   :numref:`fig-lapl-inv-mesh`. ``ig`` refers to "inner ghost point",
   ``og`` refers to "outer ghost point", and ``c`` refers to the point
   of consideration. Notice the "wrap-around" in :math:`z`-direction
   when the point of consideration neighbours the first/last
   :math:`z`-index.

As we are using a row-major implementation, the global indices of the
matrix will be as in :numref:`fig-lapl-inv-global`

.. _fig-lapl-inv-global:
.. figure:: ../figs/5PointStencilGlobalIndices.*
   :alt: Global indices of the matrix in figure
         :numref:`fig-lapl-inv-matrix`

   Global indices of the matrix in figure :numref:`fig-lapl-inv-matrix`


Implementation internals
------------------------

The Laplacian inversion code solves the equation:

.. math:: d\nabla^2_\perp x + \frac{1}{c}\nabla_\perp c\cdot\nabla_\perp x + a x = b

where :math:`x` and :math:`b` are 3D variables, whilst :math:`a`,
:math:`c` and :math:`d` are 2D variables. Several different algorithms
are implemented for Laplacian inversion, and they differ between
serial and parallel versions. Serial inversion can currently either be
done using a tridiagonal solver (Thomas algorithm), or a band-solver
(allowing :math:`4^{th}`-order differencing).

To support multiple implementations, a base class `Laplacian` is
defined in ``include/invert_laplace.hxx``. This defines a set of
functions which all implementations must provide::

    class Laplacian {
     public:
      virtual void setCoefA(const Field2D &val) = 0;
      virtual void setCoefC(const Field2D &val) = 0;
      virtual void setCoefD(const Field2D &val) = 0;

      virtual const FieldPerp solve(const FieldPerp &b) = 0;
    }

At minimum, all implementations must provide a way to set coefficients,
and a solve function which operates on a single FieldPerp (X-Y) object
at once. Several other functions are also virtual, so default code
exists but can be overridden by an implementation.

For convenience, the `Laplacian` base class also defines a function to
calculate coefficients in a Tridiagonal matrix::

      void tridagCoefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b,
                       dcomplex &c, const Field2D *ccoef = NULL,
                       const Field2D *d=NULL);

For the user of the class, some static functions are defined::

      static Laplacian* create(Options *opt = NULL);
      static Laplacian* defaultInstance();

The create function allows new Laplacian implementations to be created,
based on options. To use the options in the ``[laplace]`` section, just
use the default::

      Laplacian* lap = Laplacian::create();

The code for the `Laplacian` base class is in
``src/invert/laplace/invert_laplace.cxx``. The actual creation of new
Laplacian implementations is done in the `LaplaceFactory` class,
defined in ``src/invert/laplace/laplacefactory.cxx``. This file
includes all the headers for the implementations, and chooses which
one to create based on the “type” setting in the input options. This
factory therefore provides a single point of access to the underlying
Laplacian inversion implementations.

Each of the implementations is in a subdirectory of
``src/invert/laplace/impls`` and is discussed below.

.. _sec-tri:

Serial tridiagonal solver
~~~~~~~~~~~~~~~~~~~~~~~~~

This is the simplest implementation, and is in
``src/invert/laplace/impls/serial_tri/``

.. _sec-band:

Serial band solver
~~~~~~~~~~~~~~~~~~

This is band-solver which performs a :math:`4^{th}`-order inversion.
Currently this is only available when ``NXPE=1``; when more than one
processor is used in :math:`x`, the Laplacian algorithm currently
reverts to :math:`3^{rd}`-order.

.. _sec-spt:

SPT parallel tridiagonal
~~~~~~~~~~~~~~~~~~~~~~~~

This is a reference code which performs the same operations as the
serial code. To invert a single XZ slice (`FieldPerp` object), data
must pass from the innermost processor (``mesh->PE_XIND = 0``) to the
outermost ``mesh->PE_XIND = mesh->NXPE-1`` and back again.

Some parallelism is achieved by running several inversions
simultaneously, so while processor 1 is inverting Y=0, processor 0 is
starting on Y=1. This works ok as long as the number of slices to be
inverted is greater than the number of X processors
(``MYSUB > mesh->NXPE``). If ``MYSUB < mesh->NXPE`` then not all
processors can be busy at once, and so efficiency will fall sharply.
:numref:`fig-par-laplace` shows the useage of 4 processors inverting a
set of 3 poloidal slices (i.e. MYSUB=3)

.. _fig-par-laplace:
.. figure:: ../figs/par_laplace.*
   :alt: Parallel Laplacian inversion

   Parallel Laplacian inversion with MYSUB=3 on 4 processors. Red
   periods are where a processor is idle - in this case about 40% of the
   time

.. _sec-pdd:

PDD algorithm
~~~~~~~~~~~~~

This is the Parallel Diagonally Dominant (PDD) algorithm. It’s very
fast, but achieves this by neglecting some cross-processor terms. For
ELM simulations, it has been found that these terms are important, so
this method is not usually used.

.. _sec-naulin:

Naulin solver
~~~~~~~~~~~~~

This scheme was introduced for BOUT++ by Michael Løiten in the `CELMA code
<https://github.com/CELMA-project/CELMA>`_ and the iterative algoritm is detailed in
his thesis [Løiten2017]_.

The iteration can be under-relaxed (see naulin_laplace.cxx for more details of the
implementation). A factor 0<underrelax_factor<=1 is used, with a value of 1 corresponding
to no under-relaxation. If the iteration starts to diverge (the error increases on any
step) the underrelax_factor is reduced by a factor of 0.9, and the iteration is restarted
from the initial guess. The initial value of underrelax_factor, which underrelax_factor is
set to at the beginning of each call to ``solve`` can be set by the option
``initial_underrelax_factor`` (default is 1.0) in the appropriate section of the input
file (``[laplace]`` by default). Reducing the value of ``initial_underrelax_factor`` may
speed up convergence in some cases. Some statistics from the solver are written to the
output files to help in choosing this value. With ``<i>`` being the number of the
LaplaceNaulin solver, counting in the order they are created in the physics model:
- ``naulinsolver<i>_mean_underrelax_counts`` gives the mean number of times
  ``underrelax_factor`` had to be reduced to get the iteration to converge. If this is
  much above 0, it is probably worth reducing ``initial_underrelax_factor``.
- ``naulinsolver<i>_mean_its`` is the mean number of iterations taken to converge.  Try to
  minimise when adjusting ``initial_underrelax_factor``.

.. [Løiten2017] Michael Løiten, "Global numerical modeling of magnetized plasma
   in a linear device", 2017, https://celma-project.github.io/.

.. _sec-LaplaceXY:

LaplaceXY
---------

Perpendicular Laplacian solver in X-Y.

.. math::
   :label: nabl_perp_f

   \nabla_\perp f =& \nabla f - \mathbf{b}\left(\mathbf{b}\cdot\nabla\right)
       \nonumber \\ =& \left(\frac{\partial f}{\partial x} -
   \frac{g_{xy}}{g_{yy}}\frac{\partial f}{\partial y}\right)\nabla x +
   \left(\frac{\partial f}{\partial z} - \frac{g_{yz}}{g_{yy}}\frac{\partial
   f}{\partial y}\right)\nabla z

In 2D (X-Y), the :math:`g_{xy}` component can be dropped since this
depends on integrated shear :math:`I` which will cancel with the
:math:`g_{xz}` component. The :math:`z` derivative is zero and so this
simplifies to

.. math::

   \nabla_\perp f = \frac{\partial f}{\partial x}\nabla x -
   \frac{g_{yz}}{g_{yy}}\frac{\partial f}{\partial y}\nabla z

The divergence operator in conservative form is

.. math::

   \nabla\cdot\mathbf{A} = \frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}A_j\right)

and so the perpendicular Laplacian in X-Y is

.. math::

   \nabla_\perp^2f = \frac{1}{J}\frac{\partial}{\partial
   x}\left(Jg^{xx}\frac{\partial f}{\partial x}\right) -
   \frac{1}{J}\frac{\partial}{\partial
   y}\left(Jg^{yz}\frac{g_{yz}}{g_{yy}}\frac{\partial f}{\partial y}\right)

In field-aligned coordinates, the metrics in the :math:`y` derivative
term become:

.. math::

   g^{yz}\frac{g_{yz}}{g_{yy}} = \frac{B_{tor}^2}{B^2}\frac{1}{h_\theta^2}

In the LaplaceXY operator this is implemented in terms of fluxes at
cell faces.

.. math::

   \frac{1}{J}\frac{\partial}{\partial x}\left(Jg^{xx}\frac{\partial f}{\partial
   x}\right) \rightarrow
           \frac{1}{J_i\mathrm{dx_i}}\left[J_{i+1/2}g^{xx}_{i+1/2}\left(\frac{f_{i+1}
               - f_{i}}{\mathrm{dx}_{i+1/2}}\right) -
               J_{i-1/2}g^{xx}_{i-1/2}\left(\frac{f_{i} -
           f_{i-1}}{\mathrm{dx}_{i-1/2}}\right)\right]

Notes:

-  The ShiftXderivs option must be true for this to work, since it
   assumes that :math:`g^{xz} = 0`

.. _sec-LaplaceXZ:

LaplaceXZ
---------

This is a Laplacian inversion code in X-Z, similar to the
`Laplacian` solver described in :ref:`sec-laplacian`. The
difference is in the form of the Laplacian equation solved, and the
approach used to derive the finite difference formulae. The equation
solved is:

.. math::

     \nabla\cdot\left( A \nabla_\perp f \right) + Bf = b

where :math:`A` and :math:`B` are coefficients, :math:`b` is the known
RHS vector (e.g. vorticity), and :math:`f` is the unknown quantity to
be calculated (e.g. potential), and :math:`\nabla_\perp f` is the same
as equation (:eq:`nabl_perp_f`), but with negligible
:math:`y`-parallel derivatives if :math:`g_{xy}`, :math:`g_{yz}` and
:math:`g_{xz}` is non-vanishing. The Laplacian is written in
conservative form like the `LaplaceXY` solver, and
discretised in terms of fluxes through cell faces.

.. math::

     \frac{1}{J}\frac{\partial}{\partial x}\left(J A g^{xx}\frac{\partial
     f}{\partial x}\right) + \frac{1}{J}\frac{\partial}{\partial z}\left(J A
     g^{zz}\frac{\partial f}{\partial z}\right) + B f = b

The header file is ``include/bout/invert/laplacexz.hxx``. The solver is
constructed by using the `LaplaceXZ::create` function::

      LaplaceXZ *lap = LaplaceXZ::create(mesh);

Note that a pointer to a `Mesh` object must be given, which
for now is the global variable `mesh`. By default the
options section ``laplacexz`` is used, so to set the type of solver
created, set in the options

.. code-block:: cfg

      [laplacexz]
      type = petsc  # Set LaplaceXZ type

or on the command-line ``laplacexz:type=petsc`` .

The coefficients must be set using ``setCoefs`` . All coefficients must
be set at the same time::

      lap->setCoefs(1.0, 0.0);

Constants, `Field2D` or `Field3D` values can be passed. If the
implementation doesn’t support `Field3D` values then the average over
:math:`z` will be used as a `Field2D` value.

To perform the inversion, call the ``solve`` function::

      Field3D vort = ...;

      Field3D phi = lap->solve(vort, 0.0);

The second input to ``solve`` is an initial guess for the solution,
which can be used by iterative schemes e.g. using PETSc.

Implementations
~~~~~~~~~~~~~~~

The currently available implementations are:

- ``cyclic``: This implementation assumes coefficients are constant in
  :math:`Z`, and uses FFTs in :math:`z` and a complex tridiagonal solver
  in :math:`x` for each :math:`z` mode (the `CyclicReduction`
  solver). Code in ``src/invert/laplacexz/impls/cyclic/``.

- ``petsc``: This uses the PETSc KSP interface to solve a matrix with
  coefficients varying in both :math:`x` and :math:`z`. To improve
  efficiency of direct solves, a different matrix is used for
  preconditioning. When the coefficients are updated the
  preconditioner matrix is not usually updated. This means that LU
  factorisations of the preconditioner can be re-used. Since this
  factorisation is a large part of the cost of direct solves, this
  should greatly reduce the run-time.

Test case
~~~~~~~~~

The code in ``examples/test-laplacexz`` is a simple test case for
`LaplaceXZ` . First it creates a `LaplaceXZ`
object::

      LaplaceXZ *inv = LaplaceXZ::create(mesh);

For this test the ``petsc`` implementation is the default:

.. code-block:: cfg

      [laplacexz]
      type = petsc
      ksptype = gmres # Iterative method
      pctype  = lu  # Preconditioner

By default the LU preconditioner is used. PETSc’s built-in factorisation
only works in serial, so for parallel solves a different package is
needed. This is set using::

      factor_package = superlu_dist

This setting can be “petsc” for the built-in (serial) code, or one of
“superlu”, “superlu\_dist”, “mumps”, or “cusparse”.

Then we set the coefficients::

      inv->setCoefs(Field3D(1.0),Field3D(0.0));

Note that the scalars need to be cast to fields (Field2D or Field3D)
otherwise the call is ambiguous. Using the PETSc command-line flag
``-mat_view ::ascii_info`` information on the assembled matrix is
printed:

.. code-block:: bash

      $ mpirun -np 2 ./test-laplacexz -mat_view ::ascii_info
      ...
      Matrix Object: 2 MPI processes
      type: mpiaij
      rows=1088, cols=1088
      total: nonzeros=5248, allocated nonzeros=5248
      total number of mallocs used during MatSetValues calls =0
        not using I-node (on process 0) routines
      ...

which confirms that the matrix element pre-allocation is setting the
correct number of non-zero elements, since no additional memory
allocation was needed.

A field to invert is created using FieldFactory::

      Field3D rhs = FieldFactory::get()->create3D("rhs",
                                                  Options::getRoot(),
                                                  mesh);

which is currently set to a simple function in the options::

      rhs = sin(x - z)

and then the system is solved::

      Field3D x = inv->solve(rhs, 0.0);

Using the PETSc command-line flags ``-ksp_monitor`` to monitor the
iterative solve, and ``-mat_superlu_dist_statprint`` to monitor
SuperLU\_dist we get:

.. code-block:: bash

            Nonzeros in L       19984
            Nonzeros in U       19984
            nonzeros in L+U     38880
            nonzeros in LSUB    11900
            NUMfact space (MB) sum(procs):  L\U     0.45    all     0.61
            Total highmark (MB):  All       0.62    Avg     0.31    Max     0.36
            Mat conversion(PETSc->SuperLU_DIST) time (max/min/avg):
                                  4.69685e-05 / 4.69685e-05 / 4.69685e-05
            EQUIL time             0.00
            ROWPERM time           0.00
            COLPERM time           0.00
            SYMBFACT time          0.00
            DISTRIBUTE time        0.00
            FACTOR time            0.00
            Factor flops    1.073774e+06    Mflops    222.08
            SOLVE time             0.00
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     28.67
      0 KSP Residual norm 5.169560044060e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     60.50
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     49.86
      1 KSP Residual norm 1.359142853145e-12

So after the initial setup and factorisation, the system is solved in
one iteration using the LU direct solve.

As a test of re-using the preconditioner, the coefficients are then
modified::

      inv->setCoefs(Field3D(2.0),Field3D(0.1));

and solved again::

            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     84.15
      0 KSP Residual norm 5.169560044060e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     90.42
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     98.51
      1 KSP Residual norm 2.813291076609e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     94.88
      2 KSP Residual norm 1.688683980433e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     87.27
      3 KSP Residual norm 7.436784980024e+01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     88.77
      4 KSP Residual norm 1.835640800835e+01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     89.55
      5 KSP Residual norm 2.431147365563e+00
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     88.00
      6 KSP Residual norm 5.386963293959e-01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     93.50
      7 KSP Residual norm 2.093714782067e-01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     91.91
      8 KSP Residual norm 1.306701698197e-02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     89.44
      9 KSP Residual norm 5.838501185134e-04
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     81.47

Note that this time there is no factorisation step, but the direct solve
is still very effective.

Blob2d comparison
~~~~~~~~~~~~~~~~~

.. Use bash as the default language for syntax highlighting in this section
.. highlight:: console

The example ``examples/blob2d-laplacexz`` is the same as
``examples/blob2d`` but with ``LaplaceXZ`` rather than
`Laplacian`.

Tests on one processor: Using Boussinesq approximation, so that the
matrix elements are not changed, the cyclic solver produces output::

    1.000e+02        125       8.28e-01    71.8    8.2    0.4    0.6   18.9
    2.000e+02         44       3.00e-01    69.4    8.1    0.4    2.1   20.0

whilst the PETSc solver with LU preconditioner outputs::

    1.000e+02        146       1.15e+00    61.9   20.5    0.5    0.9   16.2
    2.000e+02         42       3.30e-01    58.2   20.2    0.4    3.7   17.5

so the PETSc direct solver seems to take only slightly longer than the
cyclic solver. For comparison, GMRES with Jacobi preconditioning gives::

    1.000e+02        130       2.66e+00    24.1   68.3    0.2    0.8    6.6
    2.000e+02         78       1.16e+00    33.8   54.9    0.3    1.1    9.9

and with SOR preconditioner::

    1.000e+02        124       1.54e+00    38.6   50.2    0.3    0.4   10.5
    2.000e+02         45       4.51e-01    46.8   37.8    0.3    1.7   13.4

When the Boussinesq approximation is not used, the PETSc solver with LU
preconditioning, re-setting the preconditioner every 100 solves gives::

    1.000e+02        142       3.06e+00    23.0   70.7    0.2    0.2    6.0
    2.000e+02         41       9.47e-01    21.0   72.1    0.3    0.6    6.1

i.e. around three times slower than the Boussinesq case. When using
jacobi preconditioner::

    1.000e+02        128       2.59e+00    22.9   70.8    0.2    0.2    5.9
    2.000e+02         68       1.18e+00    26.5   64.6    0.2    0.6    8.1

For comparison, the `Laplacian` solver using the
tridiagonal solver as preconditioner gives::

    1.000e+02        222       5.70e+00    17.4   77.9    0.1    0.1    4.5
    2.000e+02        172       3.84e+00    20.2   74.2    0.2    0.2    5.2

or with Jacobi preconditioner::

    1.000e+02        107       3.13e+00    15.8   79.5    0.1    0.2    4.3
    2.000e+02        110       2.14e+00    23.5   69.2    0.2    0.3    6.7

The `LaplaceXZ` solver does not appear to be dramatically faster **in
serial** than the `Laplacian` solver when the matrix coefficients are
modified every solve. When matrix elements are not modified then the
solve time is competitive with the tridiagonal solver.

As a test, timing only the ``setCoefs`` call for the non-Boussinesq case
gives::

    1.000e+02        142       1.86e+00    83.3    9.5    0.2    0.3    6.7
    2.000e+02         41       5.04e-01    83.1    8.0    0.3    1.2    7.3

so around 9% of the run-time is in setting the coefficients, and the
remaining :math:`\sim 60`\ % in the solve itself.
