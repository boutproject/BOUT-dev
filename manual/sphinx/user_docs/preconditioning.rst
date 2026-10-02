.. default-role:: math

.. _sec-preconditioning:

=====================================
BOUT++ preconditioning and Jacobians
=====================================

Implicit solvers repeatedly solve linearised systems of the form
`(I - \gamma J) x = b`. BOUT++ supports several ways to help those solves:

- a user-supplied preconditioner registered with ``solver->setPrecon(...)``
- a user-supplied Jacobian-vector product registered with
  ``solver->setJacobian(...)``
- PETSc finite-difference Jacobians, usually accelerated with coloring
- extra Jacobian sparsity supplied from ``PetscCellOperator`` stencils with
  ``solver->addJacobianPattern(...)``

This page starts with the user-facing workflow and then moves into the current
implementation details and limitations.


Choosing an approach
--------------------

The most useful question is usually not "which feature exists?" but "which
approximation is easiest to provide for this model?".

- Use matrix-free methods when assembling a Jacobian is too expensive or the
  coupling is not well represented by a sparse local stencil.
- Add a user preconditioner when you know an approximate inverse for the stiff
  part of the physics.
- Use PETSc finite-difference Jacobians when the coupling is sparse and local
  enough that a structured sparsity pattern is a good approximation.
- Add ``PetscCellOperator`` sparsity when the solver's default stencil is too
  small, but the coupling is still naturally represented by a sparse
  cell-to-cell operator.

The solver options that choose between these modes are described in
:ref:`sec-time-integration`. PETSc option prefixes are described in
:ref:`sec-petsc`.


User-supplied preconditioners
-----------------------------

A user preconditioner is a callback that approximately solves
`(I - \gamma J) x = b` for the solver's current linearisation. The callback
signature is:

.. code-block:: C++

   int precon(BoutReal t, BoutReal gamma, BoutReal delta);

Register it during ``PhysicsModel::init``:

.. code-block:: C++

   int init(bool restarting) override {
     solver->setPrecon(precon);
     return 0;
   }

Inside the callback:

- the current state is stored in the usual evolving variables
- the vector to be preconditioned is stored in the corresponding
  ``ddt(variable)`` fields
- the callback should overwrite those ``ddt(...)`` fields with the
  preconditioned result

For CVODE, enable a user preconditioner with for example:

.. code-block:: cfg

   [solver]
   type = cvode
   cvode_precon_method = user

The ``examples/test-precon`` case shows the standard pattern. For the simple
wave system

.. math::

   \frac{\partial u}{\partial t} = \partial_{||} v
   \qquad
   \frac{\partial v}{\partial t} = \partial_{||} u

one useful approximation factors `(I - \gamma J)^{-1}` into inexpensive pieces:

.. math::

   \left(\begin{array}{cc}
   1 & -\gamma \partial_{||} \\
   -\gamma \partial_{||} & 1
   \end{array}\right)^{-1}
   =
   \left(\begin{array}{cc}
   1 & \gamma \partial_{||} \\
   0 & 1
   \end{array}\right)
   \left(\begin{array}{cc}
   1 & 0 \\
   0 & (1 - \gamma^2 \partial_{||}^2)^{-1}
   \end{array}\right)
   \left(\begin{array}{cc}
   1 & 0 \\
   \gamma \partial_{||} & 1
   \end{array}\right)

In BOUT++ that becomes a sequence of field operations and an ``InvertPar`` solve:

.. code-block:: C++

   int precon(BoutReal t, BoutReal gamma, BoutReal delta) {
     mesh->communicate(ddt(u));
     ddt(v) = gamma * Grad_par(ddt(u)) + ddt(v);

     inv->setCoefB(-SQ(gamma));
     ddt(v) = inv->solve(ddt(v));

     mesh->communicate(ddt(v));
     ddt(u) = ddt(u) + gamma * Grad_par(ddt(v));

     ddt(u).applyBoundary("dirichlet");
     ddt(v).applyBoundary("dirichlet");
     return 0;
   }

This is the main design pattern for analytic preconditioners: keep only the
stiff physics, approximate the rest, and return something that is cheap but
captures the dominant fast scales.


Jacobian-vector products
------------------------

If the model can apply the Jacobian to a vector directly, register that
operation with ``solver->setJacobian(...)``. This is primarily useful for
matrix-free methods, where the solver needs Jacobian-vector products but does
not assemble a full sparse Jacobian matrix.

For problems where a sparse matrix is more useful than a matrix-free operator,
PETSc-backed solvers can instead assemble the Jacobian by finite differences.


PETSc finite-difference Jacobians
---------------------------------

When a PETSc-backed implicit solver runs with ``matrix_free = false``, BOUT++
can assemble an explicit sparse Jacobian by finite differences. In most cases
``use_coloring = true`` is the right choice:

.. code-block:: cfg

   [solver]
   matrix_free = false
   use_coloring = true
   lag_jacobian = 5

   stencil:taxi = 2
   stencil:square = 0
   stencil:cross = 0

The default pattern comes from a solver-side stencil and is then passed to PETSc
coloring, which groups independent columns so the Jacobian can be built with far
fewer RHS evaluations than a brute-force finite-difference calculation.

The most important current assumption is that the Jacobian is sparse and local
in the solver ordering. If the RHS includes long-range couplings such as Fourier
transforms, matrix inversions, or other nonlocal operators, the default stencil
can miss real dependencies and the implicit solve may fail to converge. In those
cases the usual alternatives are:

- stay matrix-free
- reformulate the problem so the long-range solve becomes a constraint
- supply a better sparsity pattern when one exists

The main tuning options are:

- ``use_coloring = false`` to fall back to a brute-force finite-difference
  Jacobian
- ``stencil:taxi``, ``stencil:square``, and ``stencil:cross`` to widen the
  default solver stencil
- ``force_symmetric_coloring = true`` to make the coloring pattern symmetric
- ``lag_jacobian`` to reuse the assembled Jacobian across successive nonlinear
  iterations


Augmenting Jacobian structure with PetscOperator
------------------------------------------------

If the coupling you want to expose to the solver is already available as a
``PetscCellOperator``, you can register its sparsity pattern directly with the
solver before ``solver->init()``.

For example:

.. code-block:: C++

   #if BOUT_HAS_PETSC
   PetscOperators ops(mesh);
   auto parallel = ops.getParallel();

   auto n = solver->getVarRef("n");
   auto T = solver->getVarRef("T");

   // Insert the operator stencil into every Jacobian block
   solver->addJacobianPattern(parallel.Div_par_Grad_par);

   // Insert the operator stencil only into dF_n / dT
   solver->addJacobianPattern(parallel.Grad_par, n, T);
   #endif

The one-argument form is shorthand for:

.. code-block:: C++

   solver->addJacobianPattern(op, Solver::VarRef::All(), Solver::VarRef::All());

``Solver::VarRef::All()`` expands when the solver Jacobian is created:

- ``(All, All)`` inserts the stencil into every variable block
- ``(out, All)`` inserts it into one output-variable row block against all inputs
- ``(All, in)`` inserts it into all output-variable row blocks against one input
- ``(out, in)`` inserts it into a single block

This interface is guarded by ``BOUT_HAS_PETSC``. Even in PETSc-enabled builds,
``addJacobianPattern(...)`` may return ``false`` if the chosen solver does not
use the PETSc preconditioner/Jacobian path.


How operator-supplied sparsity is applied
-----------------------------------------

The operator itself lives on the full PETSc cell space, including boundary and
virtual-cell entries used by the operator construction. The solver Jacobian does
not: it contains only the evolving degrees of freedom.

To bridge those two layouts, BOUT++ uses the operator's stored cell mapping to
extract the evolving-cell submatrix first. That restricted submatrix is then
inserted into one or more solver Jacobian blocks using the helper documented in
``bout/petsc_jacobian.hxx``.

The important consequence is that ``addJacobianPattern(...)`` contributes only
the nonzero structure. The Jacobian entries themselves are still computed later
by the solver's finite-difference machinery.

At the moment the sequence is:

1. Model or component code calls ``solver->addJacobianPattern(...)`` during
   setup.
2. The solver stores the restricted operator submatrix and the requested output
   and input variable references.
3. During ``solver->init()``, a supported PETSc-backed solver creates its
   default Jacobian pattern.
4. The queued operator patterns are replayed into that Jacobian.
5. PETSc coloring and finite-difference Jacobian assembly use the augmented
   sparsity pattern.


Current limitations
-------------------

The helper used to insert operator sparsity into the solver Jacobian currently
assumes a uniform per-cell variable ordering of the form
``global_row = cell * nvars + var`` over evolving cells only. That means this
path is currently intended for:

- evolving ``Field3D`` variables
- no evolving boundary cells
- solver layouts that use a uniform per-cell interleaving of variables

It is not currently valid for mixed ``Field2D`` and ``Field3D`` solver layouts,
for layouts that include evolving boundary variables, or for any other
non-uniform cell-to-row mapping.


See also
--------

- :ref:`sec-time-integration`
- :ref:`sec-parallel-operators-petsc-fci`
- :ref:`sec-petsc`
