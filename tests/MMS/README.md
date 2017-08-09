# BOUT++ Method of Manufactured Solutions tests

These tests are designed to check the numerical error scalings of
various BOUT++ operators match the expected theoretical scalings given
their discretisation schemes, etc.

Some of these tests are very expensive and aren't expected to be run
very often (they may require a cluster or supercomputer to get the
results in a reasonable amount of time). The following should be run
for at least every pull request:

* **MMS/diffusion** checks convergence for a diffusion equation in 1D
* **MMS/wave-1d** is an MMS test of a wave equation in X, including
  Neumann and Dirichlet boundaries
* **MMS/wave-1d-y** is an MMS test of a wave equation in Y, including
  Neumann and Dirichlet boundaries
