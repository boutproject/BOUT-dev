
cdef extern from "field3d.hxx":
    cppclass Field3D:
        Field3D(Mesh * mesh);
        Field3D(const Field3D &)
        double & operator()(int,int,int)
        int getNx()
        int getNy()
        int getNz()


cdef extern from "bout/mesh.hxx":
    cppclass Mesh:
        Mesh()
        @staticmethod
        Mesh * create()
        void load()
        void setParallelTransform()

cdef extern from "invert_laplace.hxx":
    cppclass Laplacian:
        @staticmethod
        Laplacian * create()
        void solve(Field3D,Field3D)

cimport resolve_enum as benum

cdef extern from "difops.hxx":
    Field3D Div_par(Field3D, benum.CELL_LOC)
    Field3D Vpar_Grad_par(Field3D, Field3D, benum.CELL_LOC, benum.DIFF_METHOD)
    Field3D bracket(Field3D,Field3D, benum.BRACKET_METHOD, benum.CELL_LOC)

cdef extern from "derivs.hxx":
    Field3D DDZ(Field3D, benum.CELL_LOC, benum.DIFF_METHOD,bool)

cdef extern from "interpolation.hxx":
    Field3D interp_to(Field3D, benum.CELL_LOC)