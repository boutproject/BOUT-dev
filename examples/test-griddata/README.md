test-griddata
=============

Test if variable references in the input file work correctly, along with
non-default output filenames.

This test calculates various geometry quantities in the input file that
reference other variables -- notably, `dx = dr * Bpxy * Rxy`. These are saved in
a non-default output file (i.e. not `BOUT.dmp.0.nc`). Some of these variables
are then read in and their consistency is checked with a tolerance of 1e-7.

This test should be split into two separate tests: one for checking variable
references in the input file, and one for checking non-default output filenames.
