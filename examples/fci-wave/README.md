Wave equation in tokamak geometry using FCI scheme
==================================================

The input file for this test can be generated using Zoidberg,
by running `tokamak.py` in the `examples/zoidberg` directory. 
It is a MAST double-null equilibrium.

Two equations are evolved, the density and the momentum:

    dn/dt = -Div(n*v)

    d/dt(n*v) = -Div(n*v*v) - Grad_par(n) + Grad2_par2(n*v)
    
There are switches to control

* Whether density `n` is evolved, or the logarithm of the density `logn`
* The method used to calculate the divergence of the flux in the
  density equation.
  
Some alternative mod els are implemented in the `fci-wave-logn`
example, which always evolves log density, but is formulated as
a velocity equation rather than a momentum equation.
  
Zero-flux boundaries
--------------------

These cases set the velocity to zero at the boundaries, and so should
conserve total mass.

To run the tests run the python script, which will launch the simulations:

    $ python compare-density.py

If you want to run the tests individually, then they are as follows.

Evolve density, single interpolation calculation of divergence:

    mpirun -np 2 ./fci-wave -d div
    
Evolve density, area integration calculation of divergence:

    mpirun -np 2 ./fci-wave -d div-integration
    
Evolve log density, area integration calculation of divergence:

    mpirun -np 2 ./fci-wave -d logn
    
There is an analysis script which will plot the total mass as a
function of time:

    
Outgoing flow boundaries
------------------------

These set the parallel flow speed to +/- 1 on the boundaries. 

Evolve density, area integration:

    mpirun -np 2 ./fci-wave -d boundary
    
    
Evolve log density, area integration:

    mpirun -np 2 ./fci-wave -d boundary-logn

