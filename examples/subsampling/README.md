Subsampling Example
==================

This example shows how the output monitors can be used for sampling a
field at a higher frequency.

In this example the output frequency for the fields is set to 2.
Therefore the fields, which are solved for, and any field that is
added with `SAVE_REPEAT()` will be written out every second timeunit.

In the code another monitor is added with an output timestep of 0.01.
This Monitor uses a separate output file to write the data. The
options for this output file are read from the `[probes]` section
in the input.

The `Monitor1dDump` just saves a single point of the 3D field.
As mentioned, it is set to an output timestep of 0.01.
This sets the internal timestep for the solver therefore to 0.01, as
this is smaller then the standard output timestep of 2.
Therefore 200 timesteps are written to the 1D datafile, for every
output to the `BOUT.dmp` datafile.

The timestep of the `Monitor` is by default independent of the output
timestep set in `BOUT.inp`. If this is not desired, it is easy to read
this value, and set an apropriate timestep.  If the monitors are used
for physics, e.g. by coupling to an external library for some
calulation or using it as an PID controller, it is probably desirable
to to not change the physics, if the output timestep is changed.
