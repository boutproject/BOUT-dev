# distutils: language = c++
# distutils: include_dirs = ../../include
# distutils: libraries = m bout++ fftw3 netcdf_c++ netcdf m dl z hdf5_hl hdf5 hdf5_hl blas lapack pvode pvpre
# distutils: library_dirs = ../../lib


##### distutils: sources = Rectangle.cpp

cdef extern from "field3d.hxx":
    cppclass Field3D:
        Field3D();
        double & operator()(int,int,int)
    
print "hello World"

import numpy as np


cdef class Field3D_:
    """A wrapper for the Field3D_"""
    cdef Field3D * cobj
    def __init__(self,data):
        dims=np.shape(data)
        if len(dims) != 3:
            raise MemoryError("wrong number of dimensions")
        self.cobj = new Field3D()
        if self.cobj == NULL:
            raise MemoryError('Not enough memory, allocation failed.')
        
    
    def __del__(self):
        del self.cobj

cdef extern from "bout/mesh.hxx":
    cppclass Mesh:
        Mesh()

cdef extern from "meshfactory.hxx":
    cppclass MeshFactory:
        @staticmethod
        Mesh * createMesh()
cdef class Mesh_:
    cdef Mesh * cobj;
    def __init__(self):
        self.cobj = MeshFactory.createMesh();
        if self.cobj == NULL:
            raise MemoryError('Not enough memory, allocation failed.')
