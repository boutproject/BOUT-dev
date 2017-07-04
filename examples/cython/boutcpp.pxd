
cdef extern from "field3d.hxx":
    cppclass Field3D:
        Field3D(Mesh * mesh);
        double & operator()(int,int,int)
    Field3D operator+(Field3D,Field3D)
    int getNx(Field3D *)
    int getNy(Field3D *)
    int getNz(Field3D *)


cdef extern from "bout/mesh.hxx":
    cppclass Mesh:
        Mesh()
        @staticmethod
        Mesh * create()
        void load()
        void setParallelTransform()

