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

KNL @ Archer
~~~~~~~~~~~~

To use the KNL system, configure BOUT++ as follows:

.. code-block:: bash

    ./configure MPICXX=CC --host=knl --with-netcdf --with-pnetcdf=no --with-hypre=no CXXFLAGS="-xMIC-AVX512 -D_GLIBCXX_USE_CXX11_ABI=0"

Atlas
-----

.. code-block:: bash

   ./configure --with-netcdf=/usr/local/tools/hdf5-gnu-serial-1.8.1/lib --with-fftw=/usr/local --with-pdb=/usr/gapps/pact/new/lnx-2.5-ib/gnu


Cab
---

.. code-block:: bash

   ./configure --with-netcdf=/usr/local/tools/hdf5-gnu-serial-1.8.1/lib --with-fftw=/usr/local/tools/fftw3-3.2 --with-pdb=/usr/gapps/pact/new/lnx-2.5-ib/gnu



Edison
------

.. code-block:: bash

   module swap PrgEnv-intel PrgEnv-gnu
   module load fftw
   ./configure MPICC=cc MPICXX=CC --with-netcdf=/global/u2/c/chma/PUBLIC/netcdf_edison/netcdf --with-fftw=/opt/fftw/3.3.0.1/x86_64


Franklin
--------

.. code-block:: bash

   module load netcdf
   module load fftw
   ./configure MPICC=cc MPICXX=CC --with-pdb=no --with-fftw=/opt/fftw/3.2.2.1 --with-netcdf=/opt/cray/netcdf/4.1.1.0/netcdf-gnu

Hoffman2
--------

.. code-block:: bash

   ./configure --with-netcdf=/u/local/apps/netcdf/current --with-fftw=/u/local/apps/fftw3/current --with-cvode=/u/local/apps/sundials/2.4.0 --with-lapack=/u/local/apps/lapack/current

Hopper
------

.. code-block:: bash

    module swap PrgEnv-pgi PrgEnv-gnu
    module load netcdf
    module swap netcdf netcdf/4.1.3
    module swap gcc gcc/4.6.3
    ./configure MPICC=cc MPICXX=CC --with-fftw=/opt/fftw/3.2.2.1 --with-pdb=/global/homes/u/umansky/PUBLIC/PACT_HOPP2/pact

Hyperion
--------

With the bash shell use

.. code-block:: bash

   export PETSC_DIR=~farley9/projects/petsc/petsc-3.2-p1
   export PETSC_ARCH=arch-c
   ./configure --with-netcdf=/usr/local/tools/netcdf-gnu-4.1 --with-fftw=/usr/local MPICXX=mpiCC EXTRA_LIBS=-lcurl --with-petsc --with-cvode=~farley9/local --with-ida=~farley9/local

With the tcsh shell use

.. code-block:: tcsh

   setenv PETSC_DIR ~farley9/projects/petsc/petsc-3.2-p1
   setenv PETSC_ARCH arch-c
   ./configure --with-netcdf=/usr/local/tools/netcdf-gnu-4.1 --with-fftw=/usr/local MPICXX=mpiCC EXTRA_LIBS=-lcurl --with-petsc --with-cvode=~farley9/local --with-ida=~farley9/local

Ubgl
----

.. code-block:: bash

   ./configure --with-netcdf CXXFLAGS=-DMPICH_IGNORE_CXX_SEEK CFLAGS=-DMPICH_IGNORE_CXX_SEEK --with-pdb=/usr/gapps/pact/new_s/lnx-2.5-ib --with-netcdf=/usr/local/tools/netcdf/netcdf-4.1_c++



