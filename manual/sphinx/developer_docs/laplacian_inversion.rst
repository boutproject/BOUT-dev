Laplacian inversion
===================

The Laplacian inversion code solves the equation:

.. math:: d\nabla^2_\perp x + \frac{1}{c}\nabla_\perp c\cdot\nabla_\perp x + a x = b

where :math:`x` and :math:`b` are 3D variables, whilst :math:`a`,
:math:`c` and :math:`d` are 2D variables. Several different algorithms
are implemented for Laplacian inversion, and they differ between
serial and parallel versions. Serial inversion can currently either be
done using a tridiagonal solver (Thomas algorithm), or a band-solver
(allowing :math:`4^{th}`-order differencing).

To support multiple implementations, a base class ``Laplacian`` is
defined in ``include/invert_laplace.hxx``. This defines a set of
functions which all implementations must provide::

    class Laplacian {
     public:
      virtual void setCoefA(const Field2D &val) = 0;
      virtual void setCoefC(const Field2D &val) = 0;
      virtual void setCoefD(const Field2D &val) = 0;
     
      virtual const FieldPerp solve(const FieldPerp &b) = 0;
    }

At minimum, all implementations must provide a way to set coefficients,
and a solve function which operates on a single FieldPerp (X-Y) object
at once. Several other functions are also virtual, so default code
exists but can be overridden by an implementation.

For convenience, the ``Laplacian`` base class also defines a function to
calculate coefficients in a Tridiagonal matrix::

      void tridagCoefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b,
                       dcomplex &c, const Field2D *ccoef = NULL,
                       const Field2D *d=NULL);

For the user of the class, some static functions are defined::

      static Laplacian* create(Options *opt = NULL);
      static Laplacian* defaultInstance();

The create function allows new Laplacian implementations to be created,
based on options. To use the options in the ``[laplace]`` section, just
use the default::

      Laplacian* lap = Laplacian::create();

The code for the ``Laplacian`` base class is in
``src/invert/laplace/invert_laplace.cxx``. The actual creation of new
Laplacian implementations is done in the ``LaplaceFactory`` class,
defined in ``src/invert/laplace/laplacefactory.cxx``. This file includes
all the headers for the implementations, and chooses which one to create
based on the “type” setting in the input options. This factory therefore
provides a single point of access to the underlying Laplacian inversion
implementations.

Each of the implementations is in a subdirectory of
``src/invert/laplace/impls`` and is discussed below.

Serial tridiagonal solver
-------------------------

This is the simplest implementation, and is in
``src/invert/laplace/impls/serial_tri/``

Serial band solver
------------------

This is band-solver which performs a :math:`4^{th}`-order inversion.
Currently this is only available when ``NXPE=1``; when more than one
processor is used in :math:`x`, the Laplacian algorithm currently
reverts to :math:`3^{rd}`-order.

SPT parallel tridiagonal
------------------------

This is a reference code which performs the same operations as the
serial code. To invert a single XZ slice (``FieldPerp`` object), data
must pass from the innermost processor (``mesh->PE_XIND = 0``) to the
outermost ``mesh->PE_XIND = mesh->NXPE-1`` and back again.

Some parallelism is achieved by running several inversions
simultaneously, so while processor 1 is inverting Y=0, processor 0 is
starting on Y=1. This works ok as long as the number of slices to be
inverted is greater than the number of X processors
(``MYSUB > mesh->NXPE``). If ``MYSUB < mesh->NXPE`` then not all
processors can be busy at once, and so efficiency will fall sharply.
:numref:`fig-par-laplace` shows the useage of 4 processors inverting a
set of 3 poloidal slices (i.e. MYSUB=3)

.. _fig-par-laplace:
.. figure:: ../figs/par_laplace.*
   :alt: Parallel Laplacian inversion

   Parallel Laplacian inversion with MYSUB=3 on 4 processors. Red
   periods are where a processor is idle - in this case about 40% of the
   time

PDD algorithm
-------------

This is the Parallel Diagonally Dominant (PDD) algorithm. It’s very
fast, but achieves this by neglecting some cross-processor terms. For
ELM simulations, it has been found that these terms are important, so
this method is not usually used.

