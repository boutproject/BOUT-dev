# distutils: language=c++

from libcpp cimport bool
from libcpp.string cimport string


cdef extern from "boutexception_helper.hxx":
     cdef void raise_bout_py_error()


cdef extern from "options_netcdf.hxx" namespace "bout":
    cdef void writeDefaultOutputFile();
    cdef void writeDefaultOutputFile(Options& options);
    cppclass OptionsNetCDF:
        enum FileMode:
             replace
             append
        OptionsNetCDF() except +raise_bout_py_error
        OptionsNetCDF(string filename) except +raise_bout_py_error
        OptionsNetCDF(string filename, FileMode mode) except +raise_bout_py_error
        OptionsNetCDF(const OptionsNetCDF&);
        OptionsNetCDF(OptionsNetCDF&&);
        OptionsNetCDF& operator=(const OptionsNetCDF&);
        OptionsNetCDF& operator=(OptionsNetCDF&&);
        Options read();
        void write(const Options& options);
        void write(const Options& options, string time_dim);
        void verifyTimesteps() const;


cdef extern from "options.hxx":
    cppclass Options:
        Options()
        @staticmethod
        Options* getRoot()
        @staticmethod
        Options& root()
        Options& operator[](string)
        T operator=[T](T)
        T as[T]() const
        void assign[T](T val)
        void assign[T](T val, string souce)
        void force[T](T val)
        void force[T](T val, string souce)
        void assignRepeat[T](T val)
        void assignRepeat[T](T val, string time_dimension, bool save_repeat, string source)
        bool isSet()
        T withDefault[T](T default)
        Options* getSection(string section)
        void set(string, string, string, bool)
        void get(string, string&, string)
        void get(string, double&, double)
        void get(string, bool&, bool)
        void cleanCache()


cdef extern from "optionsreader.hxx":
    cppclass OptionsReader:
        @staticmethod
        OptionsReader *getInstance()
        void read(Options *options, const char *file, ...);
