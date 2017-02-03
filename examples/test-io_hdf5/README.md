test-io_hdf5
============

Test if variables can be read and written correctly for different number of
processes, using the HDF5 file format.

The test reads variables from gridfile, then writes them to the output file. The
output file is compared to existing benchmark files with a tolerance of 1e-10.
The test is run using 1, 2 and 4 MPI processes.

The following types are read in then immediately written to the output file:

- int
- BoutReal
- Field2D
- Field3D

Additionally, the following types are evolved:

- int
- BoutReal
- Vector2D
- Vector3D
