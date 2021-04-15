Advection-Reaction equation
===========================

Split into advective and reaction parts. Can be simulated using unsplit methods
(the two parts are just combined), but intended for testing split schemes.

`Split_operator::convective` contains the advective piece, while
`Split_operator::diffusive` contains the reaction part.

Grid file simple_xz.nc contains:
- `nx = 68`
- `ny = 5`
- `dx = 1. / 64`   so X domain has length 1

In `data/BOUT.inp`:
- Domain is set to periodic in X
- The Z domain is set to size 1 (`1 / 2*pi`-th of a torus)
