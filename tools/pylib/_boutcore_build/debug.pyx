# distutils: language=c++
# distutils: include_dirs = /usr/lib64/python3.6/site-packages/numpy/core/include  /usr/include /usr/include
# distutils: libraries =  m fftw3 netcdf_c++4 netcdf m dl z hdf5_hl hdf5 hdf5_hl blas lapack
# distutils: library_dirs =  /home/dave/soft/bout-dev/cython/lib /usr/lib64 /usr/lib64
# distutils: extra_compile_args =  

cimport debug as c
from libc.stdlib cimport malloc, free

cdef bla():
    args=["/usr/bin/python3","fubar"]
    cdef char **argv = <char **>malloc(len(args) * sizeof(char*))
    fu=[]
    cdef char * tmp
    for i in range(len(args)):
        t2=str.encode(args[i])
        tmp=t2
        fu.append(tmp)
        argv[i]=<char*>fu[i]
    cdef int argc=len(args)
    c.MPI_Init(&argc,&argv);

def blas():
    bla()

cdef finis():
    c.MPI_Finalize()

def finish():
    finis()

cdef ran():
    return c.get_rank()

def rank():
    return ran()