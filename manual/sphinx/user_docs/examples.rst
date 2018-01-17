
.. _sec-examples:

Examples
========

The code and input files in the ``examples/`` subdirectory are for
research, demonstrating BOUT++, and to check for broken functionality.
Some proper unit tests have been implemented, but this is something
which needs improving. The examples which were published in
[Dudson2009]_ were ``drift-instability``, ``interchange-instability``
and ``orszag-tang``.

.. [Dudson2009] https://www.sciencedirect.com/science/article/B6TJ5-4VTCM95-3/2/ed200cd23916d02f86fda4ce6887d798


advect1d
--------

The model in ``gas_compress.cxx`` solves the compressible gas dynamics
equations for the density :math:`n`, velocity :math:`\mathbf{V}`, and
pressure :math:`P`:

drift-instability
-----------------

The physics code ``2fluid.cxx`` implements a set of reduced Braginskii
2-fluid equations, similar to those solved by the original BOUT code.
This evolves 6 variables: Density, electron and ion temperatures,
parallel ion velocity, parallel current density and vorticity.

Input grid files are the same as the original BOUT code, but the output
format is different.

em-drift
--------

gyro-gem
--------

interchange-instability
-----------------------

.. figure:: ../figs/interchange_inst_test.*
   :alt: Interchange instability test

   Interchange instability test. Solid lines are from analytic theory,
   symbols from BOUT++ simulations, and the RMS density is averaged over
   :math:`z`. Vertical dashed line marks the reference point, where
   analytic and simulation results are set equal

jorek-compare
-------------

lapd-drift
----------

orszag-tang
-----------

The file ``mhd.cxx`` solves the full MHD equations for the full values
(perturbation + initial), whilst the file ``mhd_perturb.cxx`` solves for
a perturbation about the equilibrium.

shear-alfven-wave
-----------------

sod-shock
---------

.. figure:: ../figs/sod_result.*
   :alt: Sod shock-tube problem for testing shock-handling methods
   :width: 48.0%

   Sod shock-tube problem for testing shock-handling methods

uedge-benchmark
---------------
