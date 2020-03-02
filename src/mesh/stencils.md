Notes concerning the generation of stencils
================

We want to create a Tayler function
$f(x-x_0)=\sum_i=0^n \frac{1}{i!}f_i(x-x_0)^i$ where $n$
is the order of the function, $x_0$ is the point in the boundary
where we want to calculate the function. $f_i$ are some coefficients
that we need to determine. To be precise, only $f_0$ needs to be
determined.
We know that the function has at some points certain values. If the
value at some distance `spacing.f0` is a given value `val` then we
can build a linear system of equations using the above formula.
If rather the derivative is given, the above equations needs to be
differentiated once.

stencils_sympy.py calculates the coefficients of the above matrix
which represents our system of equations. The derivative is simply
one the factor of the next smaller term (or zero if the there is no
smaller one). This is what is calculated by `tayler`, `dirichlet`
and `neumann`, the respective matrix coefficients.

sympy does all the heavy lifting on analytically inverting the
matrix.

With the analytic inversion we can put in the numerical offsets
`spacing.f?` in C++ and get a fast expression for the respective
coefficients. As mentioned before, we do not need the full inverse,
just the first row, as we only care about the value, not about it's
derivative.
