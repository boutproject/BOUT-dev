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

Usage of the laplacian inversion
--------------------------------

In BOUT++, equation :eq:`full_laplace_inv` can be solved in two
ways. The first method Fourier transforms in the :math:`z`-direction,
whilst the other is solving the full two dimensional problem by matrix
inversion. The derivation of :math:`\nabla_\perp^2f` for a general
coordinate system can be found in the ``coordinates`` manual. What is
important, is to note that if :math:`g_{xy}` and :math:`g_{yz}` are
non-zero, BOUT++ is neglecting the :math:`y`-parallel derivatives when
using the solvers ``Laplacian`` and ``LaplaceXZ``.

By neglecting the :math:`y`-derivatives (or if
:math:`g_{xy}=g_{yz}=0`), one can solve equation
:eq:`full_laplace_inv` :math:`y` plane by :math:`y` plane.

The first approach utilizes that it is possible Fourier transform the
equation in :math:`z` (using some assumptions described in section
:ref:`sec-num-laplace`), and solve a tridiagonal system for each
mode. These inversion problems are band-diagonal (tri-diagonal in the
case of 2nd-order differencing) and so inversions can be very
efficient: :math:`O(n_z \log n_z)` for the FFTs,
:math:`O(n_x)` for tridiagonal inversion using the Thomas
algorithm [1]_, where :math:`n_x` and :math:`n_z` are the number of
grid-points in the :math:`x` and :math:`z` directions respectively.

.. [1] Numerical recipes in C. The art of scientific computing, Press, W H and Teukolsky, S A and Vetterling, W T and Flannery, B P

In the second approach, the full :math:`2`\ -D system is being solved.
This requires PETSc to be built with BOUT++.

The ``Laplacian`` class is defined in ``invert_laplace.hxx`` and solves
problems formulated like equation :eq:`full_laplace_inv` To use
this class, first create an instance of it::

    Laplacian *lap = Laplacian::create();

By default, this will use the options in a section called “laplace”, but
can be given a different section as an argument. By default
:math:`d = 1`, :math:`a = 0`, and the :math:`c=1`. To set the values of
these coefficients, there are the ``setCoefA()``, ``setCoefC()``, and
``setCoefD()`` methods:

::

    Field2D a = ...;
    lap->setCoefA(a);
    lap->setCoefC(0.5);

arguments can be ``Field2D``, ``Field3D``, or real values.

Settings for the inversion can be set in the input file under the
section ``laplace`` (default) or whichever settings section name was
specified when the ``Laplacian`` class was created. Commonly used
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
   | ``type``                 | Which implementation to use                                             | ``tri`` (serial), ``spt`` (parallel)         |
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
   |        | argument, ``x0``, of ``Laplacian::solve(const Field3D &b, const      |                            |
   |        | Field3D &x0)`` . May be combined with any combination of 0, 1 and 2, |                            |
   |        | i.e. a Dirichlet or Neumann boundary condition set to values which   |                            |
   |        | are :math:`\neq 0` or :math:`f(y)`                                   |                            |
   +--------+----------------------------------------------------------------------+----------------------------+
   | 32     | Set boundary condition to values in boundary guard cells of RHS,     | ``INVERT_RHS``             |
   |        | ``b`` in ``Laplacian::solve(const Field3D &b, const Field3D &x0)``   |                            |
   |        | . May be combined with any combination of 0, 1 and 2, i.e. a         |                            |
   |        | Dirichlet or Neumann boundary condition set to values which are      |                            |
   |        | :math:`\neq 0` or :math:`f(y)`                                       |                            |
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
BOUT++ code:

::

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
``Laplacian`` and ``LaplaceXZ``. For these two solvers, equation
:eq:`to_invert` becomes (see ``coordinates`` manual for derivation)

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
