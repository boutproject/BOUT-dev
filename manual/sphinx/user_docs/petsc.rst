.. default-role:: math

.. _sec-petsc:


PETSc solvers
=============

Options for PETSc solvers can be passed in the input file (or on the command line).
Global options are set in the ``[petsc]`` section. To set options specific to a
particular PETSc-based solver, the options can be set in a ``petsc`` subsection of the
solver's options, e.g. for a LaplaceXY solver (using the default options section) use the
``[laplacexy:petsc]`` section. Note that the global options, including any
passed on the command line [*]_, will be ignored for that solver if the subsection
is created.

Any options that can be passed on the command line to PETSc can be set, with no preceding
hyphen. Flags passed with no value can be passed as options with no value. So
for example, if the command line options would be::

    -ksp_monitor -ksp_type gmres

in the input file you would put::

    [petsc]
    ksp_monitor
    ksp_type = gmres


.. [*] The object-specific options are passed to PETSc by creating an object-specific
       prefix ``boutpetsclib#_``, where ``#`` is replaced with an integer counter,
       counting the number of PetscLib instances. So an option could in principle be
       passed to a particular solver if you work out what the counter is for that solver,
       e.g.::

            -boutpetsclib1_ksp_type gmres

       The PETSc arguments ``-options_view`` and ``options_left`` might be helpful for
       this - they will show what options have been set, so will show the prefixes used.
