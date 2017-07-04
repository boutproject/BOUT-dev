# distutils: language = c++
# distutils: include_dirs = ../../include
# distutils: libraries = m bout++ fftw3 netcdf_c++ netcdf m dl z hdf5_hl hdf5 hdf5_hl blas lapack pvode pvpre
# distutils: library_dirs = ../../lib
# distutils: sources = helper.cxx
# distutils: extra_compile_args = -std=c++14

cimport boutcpp as c
print "hello World"
import numpy as np
cimport numpy as np
import atexit

from libc.stdlib cimport malloc, free

cdef extern from "helper.h":
     void c_set_f3d_all(c.Field3D * f3d, double * data)
     void c_get_f3d_all(c.Field3D * f3d, double * data)
     c.Field3D * fadd( c.Field3D*,c.Field3D*)
     c.Field3D * fmul( c.Field3D*,c.Field3D*)
     c.Mesh * c_get_global_mesh()


     
cdef class Field3D:
    """A wrapper for the Field3D"""
    cdef c.Field3D * cobj
    @classmethod
    def fromMesh(cls,mesh):
        fu=Field3D()
        fu.cobj=new c.Field3D((<Mesh?>mesh).cobj)
        return fu
    
    def __cinit__(self,Field3D obj=None):
        self.cobj=NULL
        if obj:
            self.cobj=obj.cobj
        #self.cobj = (<c.Field3D * ?> cobj_)
        #if self.cobj == NULL:
        #    raise MemoryError('Not enough memory, allocation failed.')

    def setAll(self,data):
        cdef np.ndarray[double, mode="c", ndim=3] data_ = np.ascontiguousarray(data)
        c_set_f3d_all(self.cobj,&data_[0,0,0]);
        #self.cobj.setData(0,0,0,data[0,0,0])

    def getAll(self):
        nx=c.getNx(self.cobj)
        ny=c.getNy(self.cobj)
        nz=c.getNz(self.cobj)
        print(nx,ny,nz)
        cdef np.ndarray[double, mode="c", ndim=3] data_ = np.ascontiguousarray(np.zeros((nx,ny,nz)))
        c_get_f3d_all(self.cobj,&data_[0,0,0]);
        return data_
    def __add__(self,other):
        fu=Field3D()
        fu.cobj=fadd((<Field3D?>self).cobj , (<Field3D?>other).cobj)
        return fu
    def __mul__(self,other):
        fu=Field3D()
        fu.cobj=fmul((<Field3D?>self).cobj , (<Field3D?>other).cobj)
        return fu

    def __del__(self):
        del self.cobj
        
cdef class Mesh:
    cdef c.Mesh * cobj;
    def __init__(self, create=True):
        if create:
            self.cobj = c.MeshFactory.getInstance().createMesh();
            if self.cobj == NULL:
                raise MemoryError('Not enough memory, allocation failed.')
    @classmethod
    def getGlobal(cls):
        msh = Mesh(create=False);
        msh.cobj = c_get_global_mesh();
        return msh

cdef extern from "bout.hxx":
    int BoutInitialise(int&, char **&)
    void BoutFinalise()
    
isInit=False
def init(args):
    cdef char **string_buf = <char **>malloc(len(args) * sizeof(char*))
    fu=[]
    cdef char * tmp
    for i in range(len(args)):
        tmp=args[i]
        fu.append(tmp)
        string_buf[i]=<char*>fu[i]
    cdef int fuu=len(args)
    ret=BoutInitialise(fuu,string_buf)
    free(string_buf)
    if ret:
        BoutFinalise()
    else:
        isInit=True
    atexit.register(finalise)

def finalise():
    BoutFinalise()
    isInit=False

