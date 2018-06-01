# distutils: language=c++

from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "mpi.h":
    int MPI_Init(int *argc, char ***argv)
    int MPI_Finalize()

cdef extern from "debug_helper.hxx":
    int get_rank()