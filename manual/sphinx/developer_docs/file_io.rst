File I/O
========

BOUT++ needs to deal with binary format files to read the grid; read
and write restart restart files; and write dump files. The two parts
of the code which need to read and write data are therefore the grid
routines (:doc:`grid.hxx<../_breathe_autogen/file/griddata_8hxx>`),
and the `Datafile` class
(:doc:`datafile.hxx<../_breathe_autogen/file/datafile_8hxx>` and
:doc:`datafile.cxx<../_breathe_autogen/file/datafile_8cxx>`). All
other parts which need to read or write data go through these methods.

Several different file formats are commonly used, such as HDF, HDF5,
and netCDF. For historical reasons (inherited from BOUT), BOUT++
originally used the Portable Data Binary (PDB) format developed at
LLNL [1]_. To separate the basic file format functions from the higher
level grid and Datafile classes, these use an abstract class
`DataFormat`. Any class which implements the functions listed in
:doc:`dataformat.hxx<../_breathe_autogen/file/dataformat_8hxx>` can
therefore be passed to grid or datafile. This makes implementing a new
file format, and switching between formats at run-time, relatively
straightforward.

Access to data in files is provided using a Bridge pattern: The
`Datafile` class provides an interface to the rest of the code to read
and write variables, whilst file formats implement the `Dataformat`
interface.

::

    class Datafile {
     public:
      Datafile();
      Datafile(DataFormat *format);
      ~Datafile();
      
      /// Set the file format by passing an interface class
      void setFormat(DataFormat *format);

      void setLowPrecision(); ///< Only output floats

      void add(var, const char *name, int grow = 0);

      int read(const char *filename, ...);
      int write(const char *filename, ...);
      int append(const char *filename, ...);
      bool write(const string &filename, bool append=false);

      /// Set this to false to switch off all data writing
      static bool enabled;
    };

The important bits of the DataFormat interface are::

    class DataFormat {
     public:
      bool openr(const char *name);
      bool openw(const char *name, bool append=false);
      
      bool is_valid();
      
      void close();
      
      const char* filename();

      const vector<int> getSize(const char *var);
      const vector<int> getSize(const string &var);

      // Set the origin for all subsequent calls
      bool setOrigin(int x = 0, int y = 0, int z = 0); 
      bool setRecord(int t); // negative -> latest
      
      // Read / Write simple variables up to 3D

      bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
      bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);

      bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
      bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);

      // Read / Write record-based variables

      bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0);
      bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0);

      bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0);
      bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0);

      // Optional functions
      
      void setLowPrecision();
    };

.. [1] Support for PDB files was removed in BOUT++ 4.0.0
