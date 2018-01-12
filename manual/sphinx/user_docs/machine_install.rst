.. _sec-machine-specific:

Machine-specific installation
=============================

Archer
------

As of 30th April 2014, the following configuration should work

.. code-block:: bash

    $ module swap PrgEnv-cray PrgEnv-gnu/5.1.29
    $ module load fftw
    $ module load netcdf/4.1.3

KNL
~~~

To use the KNL system, configure BOUT++ as follows:

.. code-block:: bash

    ./configure MPICXX=CC --host=knl --with-netcdf --with-pnetcdf=no --with-hypre=no CXXFLAGS="-xMIC-AVX512 -D_GLIBCXX_USE_CXX11_ABI=0"
