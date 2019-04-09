Strang split Runge Kutta
========================

A second order Strang splitting scheme:

 - 2nd order Runge-Kutta-Legendre method for the diffusion (parabolic) part
   https://doi.org/10.1016/j.jcp.2013.08.021
   
 - 3rd order SSP-RK3 scheme for the advection (hyperbolic) part
   http://www.cscamm.umd.edu/tadmor/pub/linear-stability/Gottlieb-Shu-Tadmor.SIREV-01.pdf

Each timestep consists of

 - A half timestep of the diffusion part
 - A full timestep of the advection part
 - A half timestep of the diffusion part

