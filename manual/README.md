#BOUT++ documentation

The most up to date documentation consists of reStructuredText in the "sphinx"
subdirectory. This is used to create the online manual at https://bout-dev.readthedocs.io/en/latest/

These documents can be built into a PDF using "sphinx-build":

```bash
$ make
```

To use e.g. "sphinx-build-3" instead of "sphinx-build", run
```bash
$ make sphinx-build=sphinx-build-3
```

This should create a file "BOUT.pdf" in the "manual" directory.

To get a local html version, run

```bash
$ make html
```

This should create a file "index.html" in the "manual/html" directory.

## LaTeX documents

The LaTeX documents under "tex_files" are not updated much, but include things like
a detailed derivation of the coordinate system, which are not (yet) in the RST docs.
To build the LaTeX files into a PDF document, run "make":

```bash
$ make old
```

This should build:

### User manual
This describes how to use BOUT++, starting with installation
and simple examples.


### Developer manual
This describes in more detail how parts of the code work internally.
The developer manual is written assuming some level of C++ knowledge,
and that you have read the user manual.


### Coordinates
This manual describes the field-aligned coordinate system usually used
in BOUT++ simulations of tokamak plasmas. Commonly used quantities
and identities are derived.


### BOUT Gradperp op
Memo written for understanding approximations when the field-aligned coordinates are
used for the conventional turbulence ordering (k_par and k_perp).


### Preconditioning
Describes how to use preconditioning in BOUT++
