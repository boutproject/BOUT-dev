# distutils: language = c++
# distutils: include_dirs = ../../include
# distutils: libraries = m bout++ fftw3 netcdf_c++ netcdf m dl z hdf5_hl hdf5 hdf5_hl blas lapack pvode pvpre
# distutils: library_dirs = ../../lib
# distutils: sources = helper.cxx
# distutils: extra_compile_args = -std=c++14

cdef extern from "field3d.hxx":
    cppclass Field3D:
        Field3D(Mesh * mesh);
        double & operator()(int,int,int)
    
print "hello World"
import numpy as np
cimport numpy as np

cdef extern from "helper.h":
     void c_set_f3d_all(Field3D * f3d, double * data)

cdef class Field3D_:
    """A wrapper for the Field3D_"""
    cdef Field3D * cobj
    def __init__(self,mesh):
        self.cobj = new Field3D((<Mesh_?>mesh).cobj)
        if self.cobj == NULL:
            raise MemoryError('Not enough memory, allocation failed.')
        
    def setAll(self,data):
        cdef np.ndarray[double, mode="c", ndim=3] data_ = np.ascontiguousarray(data)
        c_set_f3d_all(self.cobj,&data_[0,0,0]);
        #self.cobj.setData(0,0,0,data[0,0,0])
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
