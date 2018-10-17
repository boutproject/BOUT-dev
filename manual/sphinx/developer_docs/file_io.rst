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
      Datafile(Options *opt = NULL);
      Datafile(Datafile &&other);
      ~Datafile();
      
      Datafile& operator=(Datafile &&rhs);
      Datafile& operator=(const Datafile &rhs) = delete;

      bool openr(const char *filename, ...);
      bool openw(const char *filename, ...); // Overwrites existing file
      bool opena(const char *filename, ...); // Appends if exists

      bool isValid();  // Checks if the data source is valid

      void close();

      void setLowPrecision(); ///< Only output floats
      template <typename t>
      void addRepeat(t &value, std::string name){
        add(value,name.c_str(),true);
      }
      template <typename t>
      void addOnce(t &value, std::string name){
        add(value,name.c_str(),false);
      }
      void add(int &i, const char *name, bool save_repeat = false);
      void add(BoutReal &r, const char *name, bool save_repeat = false);
      void add(Field2D &f, const char *name, bool save_repeat = false);
      void add(Field3D &f, const char *name, bool save_repeat = false);
      void add(Vector2D &f, const char *name, bool save_repeat = false);
      void add(Vector3D &f, const char *name, bool save_repeat = false);

      bool read();  ///< Read data into added variables
      bool write(); ///< Write added variables

      bool write(const char *filename, ...) const; ///< Opens, writes, closes file

      // Write a variable to the file now
      DEPRECATED(bool writeVar(const int &i, const char *name));
      DEPRECATED(bool writeVar(BoutReal r, const char *name));

      void setAttribute(const string &varname, const string &attrname, const string &text) {
        attrib_string[varname][attrname] = text;
      }
      void setAttribute(const string &varname, const string &attrname, int value) {
        attrib_int[varname][attrname] = value;
      }
    }

The important bits of the DataFormat interface are::

    class DataFormat {
     public:
      virtual ~DataFormat() { }
      // File opening routines
      virtual bool openr(const char *name) = 0;
      virtual bool openr(const string &name) {
        return openr(name.c_str());
      }
      virtual bool openr(const string &base, int mype);
      virtual bool openw(const char *name, bool append=false) = 0;
      virtual bool openw(const string &name, bool append=false) {
        return openw(name.c_str(), append);
      }
      virtual bool openw(const string &base, int mype, bool append=false);
      
      virtual bool is_valid() = 0;
      
      virtual void close() = 0;

      virtual void flush() = 0;

      virtual const vector<int> getSize(const char *var) = 0;
      virtual const vector<int> getSize(const string &var) = 0;

      // Set the origin for all subsequent calls
      virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) = 0;
      virtual bool setLocalOrigin(int x = 0, int y = 0, int z = 0, int offset_x = 0, int offset_y = 0, int offset_z = 0);
      virtual bool setRecord(int t) = 0; // negative -> latest
      
      // Read / Write simple variables up to 3D

      virtual bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;

      virtual bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;

      // Read / Write record-based variables

      virtual bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec(int *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec(BoutReal *var, const string &name, int lx = 1, int ly = 0, int lz = 0) = 0;

      virtual bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec(int *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec(BoutReal *var, const string &name, int lx = 0, int ly = 0, int lz = 0) = 0;

      // Optional functions
      
      virtual void setLowPrecision() { }  // By default doesn't do anything

      // Attributes

      /// Sets a string attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist
      /// @param[in] attrname    Attribute name
      /// @param[in] text        A string attribute to attach to the variable
      virtual void setAttribute(const string &UNUSED(varname), const string &UNUSED(attrname),
                                const string &UNUSED(text)) {}

      /// Sets an integer attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist
      /// @param[in] attrname    Attribute name
      /// @param[in] value       A string attribute to attach to the variable
      virtual void setAttribute(const string &UNUSED(varname), const string &UNUSED(attrname),
                                int UNUSED(value)) {}
    };

.. [1] Support for PDB files was removed in BOUT++ 4.0.0
