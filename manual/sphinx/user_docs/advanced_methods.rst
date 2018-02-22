Advanced methods
================

This section describes the more advanced methods which can be used to
speed up simulations using implicit time stepping schemes. At the time
of writing (Dec ’12), they can be used with either the SUNDIALS CVODE or
PETSc solvers.

Global field gather / scatter
-----------------------------

In BOUT++ each processor performs calculations on a sub-set of the mesh,
and communicates with other processors primarily through exchange of
guard cells (the ``mesh->commmunicate`` function). If you need to gather
data from the entire mesh onto a single processor, then this can be done
using either 2D or 3D ``GlobalFields`` .

First include the header file

::

    #include <bout/globalfield.hxx>

which defines both ``GlobalField2D`` and ``GlobalField3D`` . To create a
3D global field, pass it the mesh pointer:

::

      GlobalField3D g3d(mesh);

By default all data will be gathered onto processor 0. To change this,
specify which processor the data should go to as the second input

::

      GlobalField3D g3d(mesh, processor);

Gather and scatter methods are defined:

::

      Field3D localData;
      // Set local data to some value

      g3d.gather(localData);  // Gathers all data onto one processor

      localData = g3d.scatter(); // Scatter data back

**Note:** Boundary guard cells are **not** handled by the scatter step,
as this would mean handling branch-cuts etc. To obtain valid data in the
guard and Y boundary cells, you will need to communicate and set Y
boundaries.

**Note:** Gather and Scatter are global operations, so all processors
must call these functions.

Once data has been gathered, it can be used on one processor. To check
if the data is available, call the method ``dataIsLocal()``, which will
return ``true`` only on one processor

::

      if(g3d.dataIsLocal()) {
        // Data is available on this processor

      }

The sizes of the global array are available through ``xSize()``,
``ySize()`` and ``zSize()`` methods. The data itself can be accessed
indirectly using ``(x,y,z)`` operators:

::

      for(int x=0; x<g3d.xSize(); x++)
        for(int y=0; y<g3d.ySize(); y++)
          for(int z=0; z<g3d.zSize(); z++)
            output.write("Value at (%d,%d,%d) is %e\n",
            x,y,z,
            g3d(x,y,z) );

or by getting a pointer to the underlying data, which is stored as a 1D
array:

::

      BoutReal *data = g3d.getData();
      nx = g3d.xSize();
      ny = g3d.ySize();
      nz = g3d.zSize();

      data[x*ny*nz + y*nz + z]; // Value at g3d(x,y,z)

See the example ``examples/test-globalfield`` for more examples.

.. _sec-LaplaceXY:

LaplaceXY
---------

Perpendicular Laplacian solver in X-Y.

.. math::
   :label: nabl_perp_f

   \nabla_\perp f =& \nabla f - \mathbf{b}\left(\mathbf{b}\cdot\nabla\right)
       \nonumber \\ =& \left(\frac{\partial f}{\partial x} -
   \frac{g_{xy}}{g_{yy}}\frac{\partial f}{\partial y}\right)\nabla x +
   \left(\frac{\partial f}{\partial z} - \frac{g_{yz}}{g_{yy}}\frac{\partial
   f}{\partial y}\right)\nabla z

In 2D (X-Y), the :math:`g_{xy}` component can be dropped since this
depends on integrated shear :math:`I` which will cancel with the
:math:`g_{xz}` component. The :math:`z` derivative is zero and so this
simplifies to

.. math::

   \nabla_\perp f = \frac{\partial f}{\partial x}\nabla x -
   \frac{g_{yz}}{g_{yy}}\frac{\partial f}{\partial y}\nabla z

The divergence operator in conservative form is

.. math::

   \nabla\cdot\mathbf{A} = \frac{1}{J}\frac{\partial}{\partial
   u^i}\left(Jg^{ij}A_j\right)

and so the perpendicular Laplacian in X-Y is

.. math::

   \nabla_\perp^2f = \frac{1}{J}\frac{\partial}{\partial
   x}\left(Jg^{xx}\frac{\partial f}{\partial x}\right) -
   \frac{1}{J}\frac{\partial}{\partial
   y}\left(Jg^{yz}\frac{g_{yz}}{g_{yy}}\frac{\partial f}{\partial y}\right)

In field-aligned coordinates, the metrics in the :math:`y` derivative
term become:

.. math::

   g^{yz}\frac{g_{yz}}{g_{yy}} = \frac{B_{tor}^2}{B^2}\frac{1}{h_\theta^2}

In the LaplaceXY operator this is implemented in terms of fluxes at
cell faces.

.. math::

   \frac{1}{J}\frac{\partial}{\partial x}\left(Jg^{xx}\frac{\partial f}{\partial
   x}\right) \rightarrow
           \frac{1}{J_i\mathrm{dx_i}}\left[J_{i+1/2}g^{xx}_{i+1/2}\left(\frac{f_{i+1}
               - f_{i}}{\mathrm{dx}_{i+1/2}}\right) -
               J_{i-1/2}g^{xx}_{i-1/2}\left(\frac{f_{i} -
           f_{i-1}}{\mathrm{dx}_{i-1/2}}\right)\right]

Notes:

-  The ShiftXderivs option must be true for this to work, since it
   assumes that :math:`g^{xz} = 0`

.. _sec-LaplaceXZ:

LaplaceXZ
---------

This is a Laplacian inversion code in X-Z, similar to the ``Laplacian``
solver described in :ref:`sec-laplacian`. The difference is in the
form of the Laplacian equation solved, and the approach used to derive
the finite difference formulae. The equation solved is:

.. math::

     \nabla\cdot\left( A \nabla_\perp f \right) + Bf = b

where :math:`A` and :math:`B` are coefficients, :math:`b` is the known
RHS vector (e.g. vorticity), and :math:`f` is the unknown quantity to be
calculated (e.g. potential), and :math:`\nabla_\perp f` is the same as
equation (:eq:`nabl_perp_f`), but with negligible :math:`y`-parallel
derivatives if :math:`g_{xy}`, :math:`g_{yz}` and :math:`g_{xz}` is
non-vanishing. The Laplacian is written in conservative form like the
``LaplaceXY`` solver, and discretised in terms of fluxes through cell
faces.

.. math::

     \frac{1}{J}\frac{\partial}{\partial x}\left(J A g^{xx}\frac{\partial
     f}{\partial x}\right) + \frac{1}{J}\frac{\partial}{\partial z}\left(J A
     g^{zz}\frac{\partial f}{\partial z}\right) + B f = b

The header file is ``include/bout/invert/laplacexz.hxx``. The solver is
constructed by using the ``LaplaceXZ::create`` function:

::

      LaplaceXZ *lap = LaplaceXZ::create(mesh);

Note that a pointer to a ``Mesh`` object must be given, which for now is
the global variable ``mesh`` . By default the options section
``laplacexz`` is used, so to set the type of solver created, set in the
options

.. code-block:: cfg

      [laplacexz]
      type = petsc  # Set LaplaceXZ type

or on the command-line ``laplacexz:type=petsc`` .

The coefficients must be set using ``setCoefs`` . All coefficients must
be set at the same time:

::

      lap->setCoefs(1.0, 0.0);

Constants, ``Field2D`` or ``Field3D`` values can be passed. If the
implementation doesn’t support ``Field3D`` values then the average over
:math:`z` will be used as a ``Field2D`` value.

To perform the inversion, call the ``solve`` function:

::

      Field3D vort = ...;

      Field3D phi = lap->solve(vort, 0.0);

The second input to ``solve`` is an initial guess for the solution,
which can be used by iterative schemes e.g. using PETSc.

Implementations
~~~~~~~~~~~~~~~

The currently available implementations are:

-  ``cyclic``: This implementation assumes coefficients are constant in
   :math:`Z`, and uses FFTs in :math:`z` and a complex tridiagonal
   solver in :math:`x` for each :math:`z` mode (the ``CyclicReduction``
   solver). Code in ``src/invert/laplacexz/impls/cyclic/``.

-  ``petsc``: This uses the PETSc KSP interface to solve a matrix with
   coefficients varying in both :math:`x` and :math:`z`. To improve
   efficiency of direct solves, a different matrix is used for
   preconditioning. When the coefficients are updated the preconditioner
   matrix is not usually updated. This means that LU factorisations of
   the preconditioner can be re-used. Since this factorisation is a
   large part of the cost of direct solves, this should greatly reduce
   the run-time.

Test case
~~~~~~~~~

The code in ``examples/test-laplacexz`` is a simple test case for
``LaplaceXZ`` . First it creates a ``LaplaceXZ`` object:

::

      LaplaceXZ *inv = LaplaceXZ::create(mesh);

For this test the ``petsc`` implementation is the default:

.. code-block:: cfg

      [laplacexz]
      type = petsc
      ksptype = gmres # Iterative method
      pctype  = lu  # Preconditioner

By default the LU preconditioner is used. PETSc’s built-in factorisation
only works in serial, so for parallel solves a different package is
needed. This is set using:

::

      factor_package = superlu_dist

This setting can be “petsc” for the built-in (serial) code, or one of
“superlu”, “superlu\_dist”, “mumps”, or “cusparse”.

Then we set the coefficients:

::

      inv->setCoefs(Field3D(1.0),Field3D(0.0));

Note that the scalars need to be cast to fields (Field2D or Field3D)
otherwise the call is ambiguous. Using the PETSc command-line flag
``-mat_view ::ascii_info`` information on the assembled matrix is
printed:

.. code-block:: bash

      $ mpirun -np 2 ./test-laplacexz -mat_view ::ascii_info
      ...
      Matrix Object: 2 MPI processes
      type: mpiaij
      rows=1088, cols=1088
      total: nonzeros=5248, allocated nonzeros=5248
      total number of mallocs used during MatSetValues calls =0
        not using I-node (on process 0) routines
      ...

which confirms that the matrix element pre-allocation is setting the
correct number of non-zero elements, since no additional memory
allocation was needed.

A field to invert is created using FieldFactory:

::

      Field3D rhs = FieldFactory::get()->create3D("rhs",
                                                  Options::getRoot(),
                                                  mesh);

which is currently set to a simple function in the options:

::

      rhs = sin(x - z)

and then the system is solved:

::

      Field3D x = inv->solve(rhs, 0.0);

Using the PETSc command-line flags ``-ksp_monitor`` to monitor the
iterative solve, and ``-mat_superlu_dist_statprint`` to monitor
SuperLU\_dist we get:

.. code-block:: bash

            Nonzeros in L       19984
            Nonzeros in U       19984
            nonzeros in L+U     38880
            nonzeros in LSUB    11900
            NUMfact space (MB) sum(procs):  L\U     0.45    all     0.61
            Total highmark (MB):  All       0.62    Avg     0.31    Max     0.36
            Mat conversion(PETSc->SuperLU_DIST) time (max/min/avg):
                                  4.69685e-05 / 4.69685e-05 / 4.69685e-05
            EQUIL time             0.00
            ROWPERM time           0.00
            COLPERM time           0.00
            SYMBFACT time          0.00
            DISTRIBUTE time        0.00
            FACTOR time            0.00
            Factor flops    1.073774e+06    Mflops    222.08
            SOLVE time             0.00
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     28.67
      0 KSP Residual norm 5.169560044060e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     60.50
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     49.86
      1 KSP Residual norm 1.359142853145e-12

So after the initial setup and factorisation, the system is solved in
one iteration using the LU direct solve.

As a test of re-using the preconditioner, the coefficients are then
modified:

::

      inv->setCoefs(Field3D(2.0),Field3D(0.1));

and solved again:

::

            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     84.15
      0 KSP Residual norm 5.169560044060e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     90.42
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     98.51
      1 KSP Residual norm 2.813291076609e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     94.88
      2 KSP Residual norm 1.688683980433e+02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     87.27
      3 KSP Residual norm 7.436784980024e+01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     88.77
      4 KSP Residual norm 1.835640800835e+01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     89.55
      5 KSP Residual norm 2.431147365563e+00
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     88.00
      6 KSP Residual norm 5.386963293959e-01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     93.50
      7 KSP Residual norm 2.093714782067e-01
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     91.91
      8 KSP Residual norm 1.306701698197e-02
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     89.44
      9 KSP Residual norm 5.838501185134e-04
            SOLVE time             0.00
            Solve flops     8.245800e+04    Mflops     81.47

Note that this time there is no factorisation step, but the direct solve
is still very effective.

Blob2d comparison
~~~~~~~~~~~~~~~~~

The example ``examples/blob2d-laplacexz`` is the same as
``examples/blob2d`` but with ``LaplaceXZ`` rather than ``Laplacian``.

Tests on one processor: Using Boussinesq approximation, so that the
matrix elements are not changed, the cyclic solver produces output

::

    1.000e+02        125       8.28e-01    71.8    8.2    0.4    0.6   18.9
    2.000e+02         44       3.00e-01    69.4    8.1    0.4    2.1   20.0

whilst the PETSc solver with LU preconditioner outputs

::

    1.000e+02        146       1.15e+00    61.9   20.5    0.5    0.9   16.2
    2.000e+02         42       3.30e-01    58.2   20.2    0.4    3.7   17.5

so the PETSc direct solver seems to take only slightly longer than the
cyclic solver. For comparison, GMRES with Jacobi preconditioning gives:

::

    1.000e+02        130       2.66e+00    24.1   68.3    0.2    0.8    6.6
    2.000e+02         78       1.16e+00    33.8   54.9    0.3    1.1    9.9

and with SOR preconditioner

::

    1.000e+02        124       1.54e+00    38.6   50.2    0.3    0.4   10.5
    2.000e+02         45       4.51e-01    46.8   37.8    0.3    1.7   13.4

When the Boussinesq approximation is not used, the PETSc solver with LU
preconditioning, re-setting the preconditioner every 100 solves gives:

::

    1.000e+02        142       3.06e+00    23.0   70.7    0.2    0.2    6.0
    2.000e+02         41       9.47e-01    21.0   72.1    0.3    0.6    6.1

i.e. around three times slower than the Boussinesq case. When using
jacobi preconditioner:

::

    1.000e+02        128       2.59e+00    22.9   70.8    0.2    0.2    5.9
    2.000e+02         68       1.18e+00    26.5   64.6    0.2    0.6    8.1

For comparison, the ``Laplacian`` solver using the tridiagonal solver as
preconditioner gives:

::

    1.000e+02        222       5.70e+00    17.4   77.9    0.1    0.1    4.5
    2.000e+02        172       3.84e+00    20.2   74.2    0.2    0.2    5.2

or with Jacobi preconditioner:

::

    1.000e+02        107       3.13e+00    15.8   79.5    0.1    0.2    4.3
    2.000e+02        110       2.14e+00    23.5   69.2    0.2    0.3    6.7

The ``LaplaceXZ`` solver does not appear to be dramatically faster **in
serial** than the ``Laplacian`` solver when the matrix coefficients are
modified every solve. When matrix elements are not modified then the
solve time is competitive with the tridiagonal solver.

As a test, timing only the ``setCoefs`` call for the non-Boussinesq case
gives

::

    1.000e+02        142       1.86e+00    83.3    9.5    0.2    0.3    6.7
    2.000e+02         41       5.04e-01    83.1    8.0    0.3    1.2    7.3

so around 9% of the run-time is in setting the coefficients, and the
remaining :math:`\sim 60`\ % in the solve itself.

