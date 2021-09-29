Tokamak edge turbulence
=======================

Equilibrium from DIII-D tokamak, discharge 129131

Running the case
----------------

To set up the case, run the following in this directory

    make

Then run the 2fluid executable on >= 16 processors

    mpirun -np 16 ./2fluid
