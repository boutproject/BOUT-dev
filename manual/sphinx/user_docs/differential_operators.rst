.. _sec-diffops:

Differential operators
======================

There are a huge number of possible ways to perform differencing in
computational fluid dynamics, and BOUT++ is intended to be able to
implement a large number of them. This means that the way differentials
are handled internally is quite involved; see the developer’s manual for
full gory details. Much of the time this detail is not all that
important, and certainly not while learning to use BOUT++. Default
options are therefore set which work most of the time, so you can start
using the code without getting bogged down in these details.

In order to handle many different differencing methods and operations,
many layers are used, each of which handles just part of the problem.
The main division is between differencing methods (such as 4th-order
central differencing), and differential operators (such as
:math:`\nabla_{||}`).

.. _sec-diffmethod:

Differencing methods
--------------------

Methods are implemented on 5-point stencils, and are divided into three
categories:

-  Central-differencing methods, for diffusion operators
   :math:`\frac{df}{dx}`, :math:`\frac{d^2f}{dx^2}`. Each method has a
   short code, and currently include

   -  ``C2``: 2\ :math:`^{nd}` order :math:`f_{-1} - 2f_0 + f_1`

   -  ``C4``: 4\ :math:`^{th}` order
      :math:`(-f_{-2} + 16f_{-1} - 30f_0 + 16f_1 - f_2)/12`

   -  ``W2``: 2\ :math:`^{nd}` order CWENO

   -  ``W3``: 3\ :math:`^{rd}` order CWENO

   -  ``FFT``: Fourier Transform method in Z (axisymmetric) direction
      only

-  Upwinding methods for advection operators :math:`v_x\frac{df}{dx}`

   -  ``U1``: 1\ :math:`^{st}` order upwinding

   -  ``U4``: 4\ :math:`^{th}` order upwinding

   -  ``W3``: 3\ :math:`^{rd}` order Weighted Essentially
      Non-Oscillatory (WENO):raw-latex:`\cite{jiang-1997}`

-  Flux conserving and limiting methods for terms of the form
   :math:`\frac{d}{dx}(v_x f)`

   -  ``SPLIT``: split into upwind and central terms
      :math:`\frac{d}{dx}(v_x f) = v_x\frac{df}{dx} + f\frac{dv_x}{dx}`

   -  ``NND``: Non-oscillatory, containing No free parameters and
      Dissipative (NND) scheme:raw-latex:`\cite{nnd-2010}`

Both of these methods avoid overshoots (Gibbs phenomena) at sharp
gradients such as shocks, but the simple 1st-order method has very large
artificial diffusion. WENO schemes are a development of the ENO
reconstruction schemes which combine good handling of sharp-gradient
regions with high accuracy in smooth regions.

To use these differencing operators directly, add the following to the
top of your physics module

::

    #include <derivs.hxx>

+--------------+-----------------------------------------------+
| Function     | Formula                                       |
+==============+===============================================+
| DDX(f)       | :math:`\partial f / \partial x`               |
+--------------+-----------------------------------------------+
| DDY(f)       | :math:`\partial f / \partial y`               |
+--------------+-----------------------------------------------+
| DDZ(f)       | :math:`\partial f / \partial z`               |
+--------------+-----------------------------------------------+
| D2DX2(f)     | :math:`\partial^2 f / \partial x^2`           |
+--------------+-----------------------------------------------+
| D2DY2(f)     | :math:`\partial^2 f / \partial y^2`           |
+--------------+-----------------------------------------------+
| D2DZ2(f)     | :math:`\partial^2 f / \partial z^2`           |
+--------------+-----------------------------------------------+
| D2DX4(f)     | :math:`\partial^4 f / \partial x^4`           |
+--------------+-----------------------------------------------+
| D2DY4(f)     | :math:`\partial^4 f / \partial y^4`           |
+--------------+-----------------------------------------------+
| D2DZ4(f)     | :math:`\partial^4 f / \partial z^4`           |
+--------------+-----------------------------------------------+
| D2DXDZ(f)    | :math:`\partial^2 f / \partial x\partial z`   |
+--------------+-----------------------------------------------+
| D2DYDZ(f)    | :math:`\partial^2 f / \partial y\partial z`   |
+--------------+-----------------------------------------------+
| VDDX(f, g)   | :math:`f \partial g / \partial x`             |
+--------------+-----------------------------------------------+
| VDDY(f, g)   | :math:`f \partial g / \partial y`             |
+--------------+-----------------------------------------------+
| VDDZ(f, g)   | :math:`f \partial g / \partial z`             |
+--------------+-----------------------------------------------+
| FDDX(f, g)   | :math:`\partial/\partial x( f * g )`          |
+--------------+-----------------------------------------------+
| FDDY(f, g)   | :math:`\partial/\partial x( f * g )`          |
+--------------+-----------------------------------------------+
| FDDZ(f, g)   | :math:`\partial/\partial x( f * g )`          |
+--------------+-----------------------------------------------+

Table: Coordinate derivatives

By default the method used will be the one specified in the options
input file (see :ref:`sec-diffmethodoptions`), but most of these
methods can take an optional ``DIFF\_METHOD`` argument, specifying
exactly which method to use.

Non-uniform meshes
------------------

**examples/test-nonuniform seems to not work?** Setting
``non_uniform = true`` in the BOUT.inp options file enables corrections
to second derivatives in :math:`X` and :math:`Y`. This correction is
given by writing derivatives as:

.. math::

   {{\frac{\partial f}{\partial x}}} \simeq \frac{1}{\Delta x} {{\frac{\partial f}{\partial i}}}

where :math:`i` is the cell index number. The second derivative is
therefore given by

.. math::

   \frac{\partial^2 f}{\partial x^2} \simeq \frac{1}{\Delta x^2}\frac{\partial^2
   f}{\partial i^2} + \frac{1}{\Delta x}{{\frac{\partial f}{\partial x}}} \cdot
   {{\frac{\partial }{\partial i}}}(\frac{1}{\Delta x})

The correction factor :math:`\partial/\partial i(1/\Delta x)` can
be calculated automatically, but you can also specify ``d2x`` in the
grid file which is

.. math::

   \texttt{d2x} = {{\frac{\partial \Delta x}{\partial i}}} = \frac{\partial^2 x}{\partial i^2}

The correction factor is then calculated from ``d2x`` using

.. math::

   {{\frac{\partial }{\partial i}}}(\frac{1}{\Delta x}) = -\frac{1}{\Delta x^2} {{\frac{\partial \Delta x}{\partial i}}}

General operators
-----------------

These are differential operators which are for a general coordinate
system.

.. math::

   \begin{array}{rclrcl}
   \mathbf{v} =& \nabla f &\qquad {\texttt{Vector}} =& {\texttt{Grad(Field)}} \\
   f =& \nabla\cdot\mathbf{a} &\qquad {\texttt{Field}} =& {\texttt{Div(Vector)}} \\
   \mathbf{v} =& \nabla\times\mathbf{a} &\qquad {\texttt{Vector}} =&
   {\texttt{Curl(Vector)}} \\
   f =& \mathbf{v}\cdot\nabla g &\qquad {\texttt{Field}} =& {\texttt{V\_dot\_Grad(Vector,
   Field)}} \\
   \mathbf{v} =& \mathbf{a}\cdot\nabla\mathbf{c} &\qquad {\texttt{Vector}} =&
   {\texttt{V\_dot\_Grad(Vector, Vector)}} \\
   f =& \nabla^2 f &\qquad {\texttt{Field}} =& {\texttt{Laplace(Field)}}
   \end{array}

.. math::

   \nabla\phi =& {{\frac{\partial \phi}{\partial u^i}}}\nabla u^i arrow (\nabla\phi)_i =
       {{\frac{\partial \phi}{\partial u^i}}} \\ \nabla\cdot A =& =
       \frac{1}{J}{{\frac{\partial }{\partial u^i}}}(Jg^{ij}A_j) \\ \nabla^2\phi =&
       G^j{{\frac{\partial \phi}{\partial u^i}}} + g^{ij}\frac{\partial^2\phi}{\partial u^i\partial
       u^j}

where we have defined

.. math::

   G^j =& \frac{1}{J}{{\frac{\partial }{\partial u^i}}}(Jg^{ij})

**not** to be confused with the Christoffel symbol of the second kind
(see the coordinates manual for more details).

Clebsch operators
-----------------

Another set of operators assume that the equilibrium magnetic field is
written in Clebsch form as

.. math::

   \mathbf{B}_0 = \nabla z\times\nabla x \qquad B_0 = \frac{\sqrt{g_{yy}}}{J}

 where

.. math::

   \mathbf{B}_0 = |\mathbf{B}_0|\mathbf{b}_0 = B_0 \mathbf{b}_0

 is the background *equilibrium* magnetic field.

| l c Function & Formula
| ``Grad_par`` &
  :math:`\displaystyle\partial^0_{||} = \mathbf{b}_0\cdot\nabla =
  \frac{1}{\sqrt{g_{yy}}}{{\frac{\partial }{\partial y}}}`
| ``Div_par`` & :math:`\displaystyle \nabla^0_{||}f =
  B_0\partial^0_{||}(\frac{f}{B_0})`
| ``Grad2_par2`` & :math:`\displaystyle \partial^2_{||}\phi =
  \partial^0_{||}(\partial^0_{||}\phi) =
  \frac{1}{\sqrt{g_{yy}}}{{\frac{\partial }{\partial y}}}(\frac{1}{\sqrt{g_{yy}}}){{\frac{\partial 
  \phi}{\partial y}}} + \frac{1}{g_{yy}}\frac{\partial^2\phi}{\partial y^2}`
| ``Laplace_par`` & :math:`\displaystyle \nabla_{||}^2\phi =
  \nabla\cdot\mathbf{b}_0\mathbf{b}_0\cdot\nabla\phi =
  \frac{1}{J}{{\frac{\partial }{\partial y}}}(\frac{J}{g_{yy}}{{\frac{\partial \phi}{\partial y}}})`
| ``Laplace_perp`` &
  :math:`\displaystyle \nabla_\perp^2 = \nabla^2 - \nabla_{||}^2`
| ``Delp2`` & Perpendicular Laplacian, neglecting all :math:`y`
  derivatives
| & The ``Laplacian`` solver performs the inverse operation
| ``brackets`` & Poisson brackets
| & The Arakawa option, neglects the parallel :math:`y` derivatives if
  :math:`g_{xy}` and :math:`g_{yz}` are non-zero

We have that

.. math::

   \mathbf{b}_0\cdot\nabla\phi\times\nabla A =&
       \frac{1}{J\sqrt{g_{yy}}}[(g_{yy}{{\frac{\partial \phi}{\partial z}}} -
       g_{yz}{{\frac{\partial \phi}{\partial y}}}){{\frac{\partial A}{\partial x}}} + (g_{yz}{{\frac{\partial \phi}{\partial x}}} -
   g_{xy}{{\frac{\partial \phi}{\partial z}}}){{\frac{\partial A}{\partial y}}} + (g_{xy}{{\frac{\partial \phi}{\partial y}}} -
   g_{yy}{{\frac{\partial \phi}{\partial x}}}){{\frac{\partial A}{\partial z}}}]

.. math::

   \nabla_\perp \equiv \nabla - {{\boldsymbol{b}}}({{\boldsymbol{b}}}\cdot\nabla) \qquad
   {{\boldsymbol{b}}}\cdot\nabla = \frac{1}{JB}\frac{\partial}{\partial y}

.. math::

   {{\boldsymbol{b}}} = \frac{1}{JB}{{\boldsymbol{e}}}_y = \frac{1}{JB}[g_{xy}\nabla x + g_{yy}\nabla y
   + g_{yz}\nabla z]

In a Clebsch coordinate system
:math:`{{\boldsymbol{B}}} = \nabla z \times \nabla x = \frac{1}{J}{{\boldsymbol{e}}}_y`,
:math:`g_{yy} = {{\boldsymbol{e}}}_y\cdot{{\boldsymbol{e}}}_y = J^2B^2`,
and so the :math:`\nabla y` term cancels out:

.. math::

   \nabla_\perp =& \nabla x({{\frac{\partial }{\partial x}}} -
       \frac{g_{xy}}{(JB)^2}{{\frac{\partial }{\partial y}}}) + \nabla z({{\frac{\partial }{\partial z}}} -
       \frac{g_{yz}}{(JB)^2}{{\frac{\partial }{\partial y}}})

The bracket operators
---------------------

| The bracket operator ``brackets(phi, f, method)`` aims to
  differentiate equations on the form

  .. math::

         -\frac{\nabla\phi\times{{\boldsymbol{b}}}}{B}\cdot\nabla f

| Notice that when we use the Arakawa scheme, :math:`y`-derivatives are
  neglected if :math:`g_{xy}` and :math:`g_{yz}` are non-zero. An
  example of usage of the brackets can be found in for example
  ``examples/MMS/advection`` or ``examples/blob2d``.

Setting differencing method
---------------------------

