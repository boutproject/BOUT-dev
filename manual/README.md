BOUT++ documentation
====================
The documentation has been moved to sphinx.
It can be build localy with `make`.
To read the documentation run e.g. `firefox html/index.html`.
A online version can be found at https://bout-dev.readthedocs.io/en/latest/


To build the (deprecated) LaTeX files into a PDF document, run "make old":

```bash
$ make old
```


## User manual
This describes how to use BOUT++, starting with installation
and simple examples.


## Developer manual
This describes in more detail how parts of the code work internally.
The developer manual is written assuming some level of C++ knowledge,
and that you have read the user manual.


## Coordinates
This manual describes the field-aligned coordinate system usually used
in BOUT++ simulations of tokamak plasmas. Commonly used quantities
and identities are derived.


## BOUT Gradperp op
Memo written for understanding approximations when the field-aligned coordinates are
used for the conventional turbulence ordering (k_par and k_perp).


## Preconditioning
Describes how to use preconditioning in BOUT++
