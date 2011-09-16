setenv PETSC_DIR ~farley9/projects/petsc/petsc-3.2-p1
setenv PETSC_ARCH arch-c
./configure --with-netcdf=/usr/local/tools/netcdf-gnu-4.1 --with-fftw=/usr/local MPICXX=mpiCC EXTRA_LIBS=-lcurl --with-petsc --with-cvode=~farley9/local --with-ida=~farley9/local
