The python boutcore module
==========================

Installing
----------

Installing boutcore can be tricky.
Ideally it should be just

.. code-block:: bash

   ./configure --enable-shared
   make -j 4
   make python


but getting all the
dependencies can be difficult.
``make python`` creates the python3 module.

If problems arise, it might be worth checking a copy of the bout
module out, to reduce the risk of causing issues with the old bout
installation. This is especially true if you are trying to run
boutcore not on compute nodes of a super computer but rather on
post-processing/login/... nodes.

To use boutcore on the login node, a self compiled version of mpi may be
required, as the provided one may be only for the compute nodes.
Further, numpy header files are required, therfore numpy needs to be
compiled as well.
Further, the header files need to be exposed to the boutcore cython
compilation, e.g. by adding them to ``_boutcore_build/boutcore.pyx.in``.
It seems both ``NUMPY/numpy/core/include`` and
``NUMPY/build/src.linux-x86_64-2.7/numpy/core/include/numpy`` need to be
added, where ``NUMPY`` is the path of the numpy directory.
For running boutcore on the postprocessing nodes, fftw3 needs to be
compiled as well, if certain fftw routines are used. Note, fftw needs
to be configured with ``--enable-shared``.

After installing mpi e.g. in ``~/local/mpich``, bout needs to be
configured with something like:
``./configure CXXFLAGS="-fPIC" MPICC=~/local/mpich/bin/mpicc MPICXX=~/local/mpich/bin/mpicxx --with-fftw=~/local/fftw/``

``-fPIC`` is required, so that pvode etc is compiles as position
independend code. This can be achieved by configuring bout with
``--enable-shared``.

To only build the python2 module, run ``make python2``.

To build both, run ``make python-all``

Purpose
-------

The boutcore module exposes (part) of the BOUT++ C++ library to python.
It allows to calculate e.g. BOUT++ derivatives in python.


State
-----
It is still incomplete.
Field3D is mostly working. The other fields are not exposed.
Field3D can be accessed directly using the [] operators, and give a list of slice objects.
The get all data, ``f3d.getAll()`` is equivalent to ``f3d[:,:,]`` and returns a numpy array.
This array can be adressed with
e.g. ``[]`` operators, and then the field can be set again with
``f3d.setAll(numpyarray)``.
It is also possible to set a part of an Field3D with the ``[]`` operators.
Addition, multiplication etc should all be available.
Most of the derivatives are available, if something is missing open a bug.

Functions
---------
To be added
 - only in boutcore

 - from BOUT++

Examples
--------
To be added
