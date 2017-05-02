test-staggered
==============

Test of staggered grids for wave equation

Runs the test case twice, once with staggered grids enabled, and once with
staggered grids disabled.

In each case, a top-hat density perturbation is used to give a mixture of
starting modes. By taking the FFT in both time and space, the frequency of each
mode is calculated.
