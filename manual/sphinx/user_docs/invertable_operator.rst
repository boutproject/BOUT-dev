.. _sec-invertable:

Invertable operators
====================

A common problem in physics models is solve a matrix equation of the
form

.. math::

   \underline{\underline{A}} \cdot \underline{x} = \underline{b}

for the unknown :math:`\underline{x}`. Here
:math:`\underline{\underline{A}}` represents some differential
operator subject to boundary conditions. A specific example is the set
of Laplacian operators described in :ref:`sec-laplacian`.

Whilst specific tools are provided to deal with Laplacian and parallel
Helmholtz like equations these do not capture all possible systems and
are typically implemented (at least partially) independently of the
finite difference representation of the forward operators provided by
the rest of BOUT++. To address this a class `InvertableOperator` has
been implemented that allows the user to define a generic differential
operator and provides a simple (for the user) method to invert the
operator to find :math:`\underline{x}`. This class currently relies on
PETSc to provide the inversion functionality and hence is not
available when configuring without PETSc support. It is available in
the namespace ``bout::inversion``.

There is an example in `examples/invertable_operator` that uses the
class to solve a simple Laplacian operator and compares to the
specific Laplacian inversion solvers.

The `InvertableOperator` class is templated on the field type of the
operator (essentially defining the domain over which the problem
exists).  To define the operator that the `InvertableOperator`
instances represents one should use the
`InvertableOperator::setOperatorFunction` method. This takes a
function of signature ``std::function<T(const T&)>``. This can be a
``std::function``, compatible function pointer, lambda or a
functor. The last of these allows more complicated functions that use
a local context. For example the following code snippet demonstrates a
functor that stores several auxilliary ``Field3D`` variables used in
the ``operator()`` call::

  struct myLaplacian {
    Field3D D = 1.0, C = 1.0, A = 0.0;

    // Drop C term for now
    Field3D operator()(const Field3D &input) {
      TRACE("myLaplacian::operator()");
      Timer timer("invertable_operator_operate");
      Field3D result = A * input + D * Delp2(input);

      // Ensure boundary points are set appropriately as given by the input field.
      result.setBoundaryTo(input);

      return result;
    };
  };



A more complete example is ::

  struct myLaplacian {
    Field3D D = 1.0, C = 1.0, A = 0.0;

    // Drop C term for now
    Field3D operator()(const Field3D &input) {
      TRACE("myLaplacian::operator()");
      Timer timer("invertable_operator_operate");
      Field3D result = A * input + D * Delp2(input);

      // Ensure boundary points are set appropriately as given by the input field.
      result.setBoundaryTo(input);

      return result;
    };
  };

  bout::inversion::InveratbleOperator<Field3D> solver;
  myLaplacian laplacianOperator;
  laplacianOperator.A = 1.0;
  laplacianOperator.B = 2.0;

  // Set the function defining the operator
  solver.setOperatorFunction(laplacianOperator);

  // Perform initial setup
  solver.setup();

  // Now invert the operator for a given right hand side.
  Field3D rhs = 3.0*x;
  auto solution = solver.invert(rhs);


The PETSc backend solver is an iterative solver. It can be controlled
through the usual PETSc command line options. Note we define the
options prefix here to be `-invertable`, so instead of `-ksp_type` one
would use `-invertable_ksp_type` for example.

By default the solver caches the result to use as the initial guess
for the next call to ``invert``. There is an overload of ``invert``
that takes a second field, which is used to set the initial guess to
use in that call.

The routine ``setOperatorFunction`` takes the function by value, and
hence subsequent changes to the functor will not be reflected in the
operator without a further call to ``setOperatorFunction``. For
example::

  using bout::inversion;
  InvertableOperator<Field3D> solver;
  myLaplacian laplacianOperator;
  laplacianOperator.A = 1.0;
  laplacianOperator.B = 2.0;

  // Set the function defining the operator
  solver.setOperatorFunction(laplacianOperator);

  // Perform initial setup
  solver.setup();

  // This does not change the operator represented by
  // solver yet.
  laplacianOperator.B = -1.0;

  // This call updates the function used by solver
  // and hence the operator is update to reflect the state
  // of laplacianOperator.
  solver.setOperatorFunction(laplacianOperator);

The class provides a ``reportTime`` method that reports the time spent
in various parts of the class. Note that by including ``Timer
timer("invertable_operator_operate");`` in the function representing
the operator ``reportTime`` will include the time spent actually
applying the operator.

The class provides both ``apply`` and ``operator()`` methods that can
be used to apply the operator to a field. For example the following
should be equivalent to no operation::

  // Here result should == input, at least in the main simulation domain
  auto result = solver(solver.invert(input));


The class provides a ``verify`` method that checks that applying the
operator to the calculated inverse returns the input field within some
tolerance.

It's also possible to register a function to use as a
preconditioner. By default this is the same as the full operator
function.
