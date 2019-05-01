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
      Datafile(Options *opt = nullptr, Mesh* mesh_in = nullptr);
      Datafile(Datafile &&other) noexcept;
      ~Datafile(); // need to delete filename

      Datafile& operator=(Datafile &&rhs) noexcept;
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
      void add(FieldPerp &f, const char *name, bool save_repeat = false);
      void add(Vector2D &f, const char *name, bool save_repeat = false);
      void add(Vector3D &f, const char *name, bool save_repeat = false);

      bool read();  ///< Read data into added variables
      bool write(); ///< Write added variables

      /// Opens, writes, closes file
      bool write(const char* filename, ...) const;

      void setAttribute(const std::string &varname, const std::string &attrname, const std::string &text);
      void setAttribute(const std::string &varname, const std::string &attrname, int value);
      void setAttribute(const std::string &varname, const std::string &attrname, BoutReal value);
    };

The important bits of the DataFormat interface are::

    class DataFormat {
     public:
      DataFormat(Mesh* mesh_in = nullptr);
      virtual ~DataFormat() { }
      // File opening routines
      virtual bool openr(const char *name) = 0;
      virtual bool openr(const std::string &name) {
        return openr(name.c_str());
      }
      virtual bool openr(const std::string &base, int mype);
      virtual bool openw(const char *name, bool append=false) = 0;
      virtual bool openw(const std::string &name, bool append=false) {
        return openw(name.c_str(), append);
      }
      virtual bool openw(const std::string &base, int mype, bool append=false);

      virtual bool is_valid() = 0;

      virtual void close() = 0;

      virtual void flush() = 0;

      virtual const std::vector<int> getSize(const char *var) = 0;
      virtual const std::vector<int> getSize(const std::string &var) = 0;

      // Set the origin for all subsequent calls
      virtual bool setGlobalOrigin(int x = 0, int y = 0, int z = 0) = 0;
      virtual bool setLocalOrigin(int x = 0, int y = 0, int z = 0, int offset_x = 0, int offset_y = 0, int offset_z = 0);
      virtual bool setRecord(int t) = 0; // negative -> latest

      // Add a variable to the file
      virtual bool addVarInt(const std::string &name, bool repeat) = 0;
      virtual bool addVarBoutReal(const std::string &name, bool repeat) = 0;
      virtual bool addVarField2D(const std::string &name, bool repeat) = 0;
      virtual bool addVarField3D(const std::string &name, bool repeat) = 0;
      virtual bool addVarFieldPerp(const std::string &name, bool repeat) = 0;

      // Read / Write simple variables up to 3D

      virtual bool read(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_perp(BoutReal *var, const std::string &name, int lx = 1, int lz = 0) = 0;

      virtual bool write(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_perp(BoutReal *var, const std::string &name, int lx = 0, int lz = 0) = 0;

      // Read / Write record-based variables

      virtual bool read_rec(int *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec(int *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec(BoutReal *var, const char *name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec(BoutReal *var, const std::string &name, int lx = 1, int ly = 0, int lz = 0) = 0;
      virtual bool read_rec_perp(BoutReal *var, const std::string &name, int lx = 1, int lz = 0) = 0;

      virtual bool write_rec(int *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec(int *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec(BoutReal *var, const char *name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec(BoutReal *var, const std::string &name, int lx = 0, int ly = 0, int lz = 0) = 0;
      virtual bool write_rec_perp(BoutReal *var, const std::string &name, int lx = 0, int lz = 0) = 0;

      // Optional functions

      virtual void setLowPrecision() { }  // By default doesn't do anything

      // Attributes

      /// Sets a string attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist. If
      ///                        varname is the empty string "" then the attribute
      ///                        will be added to the file instead of to a
      ///                        variable.
      /// @param[in] attrname    Attribute name
      /// @param[in] text        A string attribute to attach to the variable
      virtual void setAttribute(const std::string &varname, const std::string &attrname,
                                const std::string &text) = 0;

      /// Sets an integer attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist. If
      ///                        varname is the empty string "" then the attribute
      ///                        will be added to the file instead of to a
      ///                        variable.
      /// @param[in] attrname    Attribute name
      /// @param[in] value       An int attribute to attach to the variable
      virtual void setAttribute(const std::string &varname, const std::string &attrname,
                                int value) = 0;

      /// Sets a BoutReal attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist. If
      ///                        varname is the empty string "" then the attribute
      ///                        will be added to the file instead of to a
      ///                        variable.
      /// @param[in] attrname    Attribute name
      /// @param[in] value       A BoutReal attribute to attach to the variable
      virtual void setAttribute(const std::string &varname, const std::string &attrname,
                                BoutReal value) = 0;

      /// Gets a string attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist. If
      ///                        varname is the empty string "" then get the
      ///                        attribute from the top-level of the file instead
      ///                        of from a variable.
      /// @param[in] attrname    Attribute name
      ///
      /// Returns
      /// -------
      /// text                   A string attribute of the variable
      virtual bool getAttribute(const std::string &varname, const std::string &attrname, std::string &text) = 0;

      /// Gets an integer attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist. If
      ///                        varname is the empty string "" then get the
      ///                        attribute from the top-level of the file instead
      ///                        of from a variable.
      /// @param[in] attrname    Attribute name
      ///
      /// Returns
      /// -------
      /// value                  An int attribute of the variable
      virtual bool getAttribute(const std::string &varname, const std::string &attrname, int &value) = 0;

      /// Gets a BoutReal attribute
      ///
      /// Inputs
      /// ------
      ///
      /// @param[in] varname     Variable name. The variable must already exist. If
      ///                        varname is the empty string "" then get the
      ///                        attribute from the top-level of the file instead
      ///                        of from a variable.
      /// @param[in] attrname    Attribute name
      ///
      /// Returns
      /// -------
      /// value                  A BoutReal attribute of the variable
      virtual bool getAttribute(const std::string &varname, const std::string &attrname, BoutReal &value) = 0;

      /// Write out the meta-data of a field as attributes of the variable
      void writeFieldAttributes(const std::string& name, const Field& f);
      /// Overload for FieldPerp so we can also write 'yindex'
      void writeFieldAttributes(const std::string& name, const FieldPerp& f);

      /// Read the attributes of a field
      void readFieldAttributes(const std::string& name, Field& f);
      /// Overload for FieldPerp so we can also read 'yindex'
      void readFieldAttributes(const std::string& name, FieldPerp& f);
    };

.. [1] Support for PDB files was removed in BOUT++ 4.0.0

FieldPerp I/O
-------------

`FieldPerp` objects can be saved to output files and read from them. The `yindex` of a
`FieldPerp` is the local y-index on a certain processor, but is saved in output files as a
global y-index in the attribute `yindex_global`. The intention is that a `FieldPerp` being
saved should be a globally well-defined object, e.g. a set of values at one divertor
target boundary, that will only be saved from processors holding that global y-index.
Actually, the C++ I/O code should work fine even if a `FieldPerp` object is defined with
different y-indices on different processors, but Python routines like `collect` and
`restart.redistribute` will fail because they find inconsistent `yindex_global` values.

`FieldPerp` objects should be added to the output files on every processor even if they
are not allocated or used, otherwise `collect` cannot find the variable in the first
output file to get its dimensions.
