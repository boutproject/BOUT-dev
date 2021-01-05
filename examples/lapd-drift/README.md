2-fluid Turbulence in a Linear device
=====================================

Benchmark cases for comparisons of BOUT 06 vs BOUT++ results

Electrostatic drift wave
LAPD plasma parameters/geometry


Analytic solution: independent solver for the eigenvalues for the
linear system of drift wave equations

BOUT: running BOUT with only linear terms included.  
Growth rate is calculated using RMS of Phi (averaged over volume). Phi
is then divided by the exponential (gamma t) and the frequency is
calculated in each point in space via matching the time series with a
sin(omega t) wave.  (fit_time program)

Everything except the first azimuthal harmonic is filtered out (zwindow = 0.) 
These solutions are for axial (y) harmonic n=1
Dominating radial mode in the solution: fundamental (~nx=0.5)

Running the BOUT++ cases
========================

1. Compile the lapd_drift executable (just run make)

2. Copy one of the BOUT.inp* input files into data/BOUT.inp
   (see below for cases) e.g. "$ cp BOUT.inp_nn data/BOUT.inp"

3. Run BOUT++ (e.g. mpirun -np 16 ./lapd_drift)

4. run the analysis script ($ idl runidl.pro)

Case 1
======
LAPD config, Ni0 profile from experiment, no phi0, no neutrals

BOUT-06 input: BOUT06.inp.nn
BOUT++ input : BOUT.inp.nn

Analytic solution: omega/OmCI  =( 0.035565860, 0.010584464)
BOUT-06:           omega/OmCI  =( 0.0348090  , 0.010362237)
BOUT++:            omega/OmCI  =( 0.0349127  , 0.010270213)

Case 2
======

LAPD config, Ni profile from experiment, no phi0, no neutrals, me=0
Same as 1), but with zero electrom mass

BOUT-06 input: BOUT06.inp.nn_zem
BOUT++ input : BOUT.inp.nn_zem

Analytic solution: omega/OmCI  =( 0.035565860, 0.010584464)
BOUT-06:           omega/OmCI  =( 0.0364727  , 0.0095402879)
BOUT++:            omega/OmCI  =( 0.0350513  , 0.010112489)

Case 3
======

LAPD config with neutrals (nu_in = 7.e-3 OmCI),
Ni profile from experiment, no phi0

BOUT-06 input: BOUT06.inp
BOUT++ input : BOUT.inp

Analytic solution: omega/OmCI  =( 0.035850606, 0.0067657669)
BOUT-06:           omega/OmCI  =( 0.0350248  , 0.0065545603)
BOUT++:            omega/OmCI  =( 0.0350987  , 0.0064631935)

