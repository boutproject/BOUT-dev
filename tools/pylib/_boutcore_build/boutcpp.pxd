from libcpp cimport bool
from libcpp.string cimport string

cdef extern from "field3d.hxx":
    cppclass Field3D:
        Field3D(Mesh * mesh);
        Field3D(const Field3D &)
        double & operator()(int,int,int)
        int getNx()
        int getNy()
        int getNz()
        bool isAllocated()
    Field3D sqrt(Field3D)
    Field3D exp(Field3D)


cdef extern from "bout/mesh.hxx":
    cppclass Mesh:
        Mesh()
        @staticmethod
        Mesh * create(Options * option)
        void load()
        void setParallelTransform()

cdef extern from "invert_laplace.hxx":
    cppclass Laplacian:
        @staticmethod
        Laplacian * create()
        Field3D solve(Field3D,Field3D)

cimport resolve_enum as benum

cdef extern from "difops.hxx":
    Field3D Div_par(Field3D, benum.CELL_LOC, benum.DIFF_METHOD)
    Field3D Grad_par(Field3D, benum.CELL_LOC, benum.DIFF_METHOD)
    Field3D Vpar_Grad_par(Field3D, Field3D, benum.CELL_LOC, benum.DIFF_METHOD)
    Field3D bracket(Field3D,Field3D, benum.BRACKET_METHOD, benum.CELL_LOC)
    Field3D Delp2(Field3D,double)

cdef extern from "options.hxx":
    cppclass Options:
        Options()
        @staticmethod
        Options * getRoot()
        Options * getSection(string fu)
        void set(string ,string,string)
    
cdef extern from "derivs.hxx":
    Field3D DDZ(Field3D, benum.CELL_LOC, benum.DIFF_METHOD,bool)

cdef extern from "interpolation.hxx":
    Field3D interp_to(Field3D, benum.CELL_LOC)


cdef extern from "field_factory.hxx":
    cppclass FieldFactory:
        FieldFactory(Mesh*,Options*)
        Field3D create3D(string bla, Options * o, Mesh * m,benum.CELL_LOC loc, double t)

