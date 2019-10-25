# Advection-Diffusion in 2D

Doubly periodic 2D domain in X-Z, evolving a single variable `U`.
An advection velocity `Vx` in the X direction, and diffusion in Z.

The default settings use the `splitrk` method, which uses Strang splitting
to combine RK3-SSP (advection) with Runge-Kutta-Legendre (diffusion). 

```
$ ./imex
1.000e-02         31              181       9.23e-02    60.7    0.0    0.7    8.6   30.1
2.000e-02         31              181       9.00e-02    61.6    0.0    0.7    7.2   30.6
...
```

This test case can be used to try different solvers, for example the PVODE solver
(adaptive implicit BDF scheme):

```
$ ./imex solver:type=pvode
1.000e-02        338              338       3.68e-01    62.8    0.0    0.5    2.0   34.7
2.000e-02        156              156       1.78e-01    59.8    0.0    0.5    3.5   36.1
...
```

In this case the `splitrk` solver is faster, but note that it has no error
control. Increasing the tolerances for `CVODE` also speeds up the calculation:
```
$ ./imex solver:type=pvode solver:atol=1e-4 solver:rtol=1e-2
1.000e-02        116              116       1.35e-01    58.6    0.0    0.5    5.5   35.4
2.000e-02         56               56       7.09e-02    54.5    0.0    0.5    8.4   36.6
3.000e-02         34               34       4.67e-02    49.5    0.0    0.5   13.3   36.7
...
```
