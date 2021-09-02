3-field ELM simulation
======================

This version uses outer-loop operators, enabling higher performance
and use of RAJA to compile for GPUs. At the time of writing this is
experimental, and the canonical version is still the `elm-pb` example.

Compile-time settings
---------------------

To obtain high performance, it is necessary to move many run-time
options to compile time. The input options are the same as `elm-pb`, but
are now checked against the compile time settings used. This is to ensure
that the code is compiled to be consistent with the input settings.

To modify the settings used, use the constexpr settings near the top of the
`elm_pb_outerloop.cxx` file.

Inputs
------

`data` is the default input. It uses a low resolution input grid, default
solvers (i.e. whatever is available in the compiled BOUT++ library),
and runs a linear calculation to find the growth rate.

