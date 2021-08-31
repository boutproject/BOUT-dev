3-field ELM simulation
======================

Inputs
------

There are several subdirectories, each with a different input file:

`data` is the default. It uses a low resolution input grid, default
solvers (i.e. whatever is available in the compiled BOUT++ library),
and runs a linear calculation to find the growth rate.

The result for linear growth rate is shown in plots:
plot Elm_n15_pb_lineargrowthrate_diam_S8hyperRm4.ps and
Elm_n15_pb_lineargrowthrate_diam_S8hyperRm4.pdf.

    IDL>  p = collect(path="data", var="P")
    IDL> moment_xyzt, p, rms=rms
    IDL> plot, deriv(alog(rms[317,32,*]))


For low resolution mesh nx=68, data/cbm18_dens8.grid_nx68ny64.nc, 
the growth rate is 0.275002 *Alfven time.

For high resolution mesh nx=516, data/cbm18_8_y064_x516_090309.nc
the growth rate is 0.186655 *Alfven time.

The diference is 47%.

`data-hypre` is similar to `data`, but specifies that the Hypre3D
solver should be used to calculate the potential. If BOUT++ has been
compiled with Hypre then this can be used.

`data-nonlinear` is a nonlinear simulation of an ELM crash. It uses a high resolution
input grid, and turns on the nonlinear terms.

Running
-------

When running elm-pb with default options, on hopper.nersc.gov, 
for 30 min debug job, by using the following commend 
>qsub bout_hopper_debug.pbs
