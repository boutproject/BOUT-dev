Advection-Reaction equation
===========================

Split into advective and reaction parts. Can be simulated using unsplit methods
(the two parts are just combined), but intended for testing split schemes.

Currently one of the RHS functions has to be called `physics_run`, so here
`physics_run` contains advection term only.

Grid file simple_xz.nc contains:
- `nx = 68`
- `ny = 5`
- `dx = 1. / 64`   so X domain has length 1

In `data/BOUT.inp`:
- Domain is set to periodic in X
- The Z domain is set to size 1 (`1 / 2*pi`-th of a torus)
