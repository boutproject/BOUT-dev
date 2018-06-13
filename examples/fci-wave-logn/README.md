Wave equation in tokamak geometry using FCI scheme
==================================================

The input file for this test can be generated using Zoidberg,
by running `tokamak.py` in the `examples/zoidberg` directory. 
It is a MAST double-null equilibrium.

Two equations are evolved, the density and the momentum:

    dn/dt = -Div(n*v)

    d/dt(v) = -v * Grad_par(v) - Grad_par(n)/n + Grad2_par2(v)
    
The logarithm of the density is evolved, allowing the `Grad_par(n)/n`
term to be written as `Grad_par(logn)`. 

There are switches to control:

* Whether the density equation is solved by integrating the flux
  or by expanding the expression and using `logn` directly.
   
  
Some alternative models are implemented in the `fci-wave`
example, which is formulated as a momentum equation rather than a
velocity equation.
  
Zero-flux boundaries
--------------------

These cases set the velocity to zero at the boundaries, and so should
conserve total mass.
  
Expanded form of density equation *This is expected to fail*:


    mpirun -np 2 ./fci-wave -d expanded
    

Area integration calculation of divergence:

    mpirun -np 2 ./fci-wave -d div-integrate
    
There is an analysis script in the `fci-wave` directory 
which will plot the total mass as a function of time.
    
Outgoing flow boundaries
------------------------

These set the parallel flow speed to +/- 1 on the boundaries. 

Area integration:

    mpirun -np 2 ./fci-wave -d boundary
    
