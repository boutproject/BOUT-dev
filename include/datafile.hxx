/*!
 * \file datafile.hxx
 *
 * \brief Data file handling object definition
 *
 * \author B.Dudson
 * \date   April 2009
 *
 * 26th Sep 2009: Modified to use varargs
 */

class Datafile;

#ifndef __DATAFILE_H__
#define __DATAFILE_H__

#include "bout_types.hxx"
#include "bout/macro_for_each.hxx"

#include "dataformat.hxx"
#include "bout/format.hxx"

#include <cstdarg>
#include <cstdio>
class Mesh;
class Field;
class Field2D;
class Field3D;
class Options;
class Vector2D;
class Vector3D;

#include <vector>
#include <string>
#include <memory>

/*!
  Uses a generic interface to file formats (DataFormat)
  and provides an interface for reading/writing simulation data.
*/
class Datafile {
 public:
  Datafile(Options *opt = nullptr, Mesh* mesh_in = nullptr);
  Datafile(Datafile &&other) noexcept;
  ~Datafile(); // need to delete filename
  
  Datafile& operator=(Datafile &&rhs) noexcept;
  Datafile& operator=(const Datafile &rhs) = delete;

  bool openr(const char *filename, ...)
    BOUT_FORMAT_ARGS( 2, 3);
  bool openw(const char *filename, ...)
    BOUT_FORMAT_ARGS( 2, 3); // Overwrites existing file
  bool opena(const char *filename, ...)
    BOUT_FORMAT_ARGS( 2, 3); // Appends if exists
  
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

  /// Opens, writes, closes file
  bool write(const char* filename, ...) const BOUT_FORMAT_ARGS(2, 3);

  void setAttribute(const std::string &varname, const std::string &attrname, const std::string &text);
  void setAttribute(const std::string &varname, const std::string &attrname, int value);
  void setAttribute(const std::string &varname, const std::string &attrname, BoutReal value);

 private:
  Mesh* mesh;
  bool parallel; // Use parallel formats?
  bool flush;    // Flush after every write?
  bool guards;   // Write guard cells?
  bool floats;   // Low precision?
  bool openclose; // Open and close file for each write
  int Lx,Ly,Lz; // The sizes in the x-, y- and z-directions of the arrays to be written
  bool enabled;  // Enable / Disable writing
  bool init_missing; // Initialise missing variables?
  bool shiftOutput; // Do we want to write out in shifted space?
  bool shiftInput;  // Read in shifted space?
  int flushFrequencyCounter; //Counter used in determining when next openclose required
  int flushFrequency; //How many write calls do we want between openclose

  std::unique_ptr<DataFormat> file;
  size_t filenamelen;
  static const size_t FILENAMELEN=512;
  char *filename;
  bool writable; // is file open for writing?
  bool appending;
  bool first_time; // is this the first time the data will be written?

  /// Shallow copy, not including dataformat, therefore private
  Datafile(const Datafile& other);

  /// A structure to hold a pointer to a class, and associated name and flags
  template <class T>
  struct VarStr {
    T *ptr;             ///< Pointer to the data.
                        ///< Note that this may be a user object, not a copy, so must not be destroyed
    std::string name;        ///< Name as it appears in the output file
    bool save_repeat;   ///< If true, has a time dimension and is saved every time step
    bool covar;         ///< For vectors, true if a covariant vector, false if contravariant
  };

  // one set per variable type
  std::vector<VarStr<int>> int_arr;
  std::vector<VarStr<BoutReal>> BoutReal_arr;
  std::vector<VarStr<Field2D>> f2d_arr;
  std::vector<VarStr<Field3D>> f3d_arr;
  std::vector<VarStr<Vector2D>> v2d_arr;
  std::vector<VarStr<Vector3D>> v3d_arr;

  bool read_f2d(const std::string &name, Field2D *f, bool save_repeat);
  bool read_f3d(const std::string &name, Field3D *f, bool save_repeat);

  bool write_int(const std::string &name, int *f, bool save_repeat);
  bool write_real(const std::string &name, BoutReal *f, bool save_repeat);
  bool write_f2d(const std::string &name, Field2D *f, bool save_repeat);
  bool write_f3d(const std::string &name, Field3D *f, bool save_repeat);

  /// Write out the meta-data of a field as attributes of the variable in
  /// 'file'.
  void writeFieldAttributes(const std::string& name, const Field& f);

  /// Check if a variable has already been added
  bool varAdded(const std::string &name);

  /// Get the pointer to the variable, nullptr if not added
  /// This is used to check if the same variable is being added
  void* varPtr(const std::string &name);
};

/// Write this variable once to the grid file
#define SAVE_ONCE1(var) bout::globals::dump.add(var, #var, 0);
#define SAVE_ONCE2(var1, var2) { \
  bout::globals::dump.add(var1, #var1, 0); \
  bout::globals::dump.add(var2, #var2, 0);}
#define SAVE_ONCE3(var1, var2, var3) {\
  bout::globals::dump.add(var1, #var1, 0); \
  bout::globals::dump.add(var2, #var2, 0); \
  bout::globals::dump.add(var3, #var3, 0);}
#define SAVE_ONCE4(var1, var2, var3, var4) { \
  bout::globals::dump.add(var1, #var1, 0); \
  bout::globals::dump.add(var2, #var2, 0); \
  bout::globals::dump.add(var3, #var3, 0); \
  bout::globals::dump.add(var4, #var4, 0);}
#define SAVE_ONCE5(var1, var2, var3, var4, var5) {\
  bout::globals::dump.add(var1, #var1, 0); \
  bout::globals::dump.add(var2, #var2, 0); \
  bout::globals::dump.add(var3, #var3, 0); \
  bout::globals::dump.add(var4, #var4, 0); \
  bout::globals::dump.add(var5, #var5, 0);}
#define SAVE_ONCE6(var1, var2, var3, var4, var5, var6) {\
  bout::globals::dump.add(var1, #var1, 0); \
  bout::globals::dump.add(var2, #var2, 0); \
  bout::globals::dump.add(var3, #var3, 0); \
  bout::globals::dump.add(var4, #var4, 0); \
  bout::globals::dump.add(var5, #var5, 0); \
  bout::globals::dump.add(var6, #var6, 0);}

#define SAVE_ONCE(...)                          \
  { MACRO_FOR_EACH(SAVE_ONCE1, __VA_ARGS__) }

/// Write this variable every timestep
#define SAVE_REPEAT1(var) bout::globals::dump.add(var, #var, 1);
#define SAVE_REPEAT2(var1, var2) { \
  bout::globals::dump.add(var1, #var1, 1); \
  bout::globals::dump.add(var2, #var2, 1);}
#define SAVE_REPEAT3(var1, var2, var3) {\
  bout::globals::dump.add(var1, #var1, 1); \
  bout::globals::dump.add(var2, #var2, 1); \
  bout::globals::dump.add(var3, #var3, 1);}
#define SAVE_REPEAT4(var1, var2, var3, var4) { \
  bout::globals::dump.add(var1, #var1, 1); \
  bout::globals::dump.add(var2, #var2, 1); \
  bout::globals::dump.add(var3, #var3, 1); \
  bout::globals::dump.add(var4, #var4, 1);}
#define SAVE_REPEAT5(var1, var2, var3, var4, var5) {\
  bout::globals::dump.add(var1, #var1, 1); \
  bout::globals::dump.add(var2, #var2, 1); \
  bout::globals::dump.add(var3, #var3, 1); \
  bout::globals::dump.add(var4, #var4, 1); \
  bout::globals::dump.add(var5, #var5, 1);}
#define SAVE_REPEAT6(var1, var2, var3, var4, var5, var6) {\
  bout::globals::dump.add(var1, #var1, 1); \
  bout::globals::dump.add(var2, #var2, 1); \
  bout::globals::dump.add(var3, #var3, 1); \
  bout::globals::dump.add(var4, #var4, 1); \
  bout::globals::dump.add(var5, #var5, 1); \
  bout::globals::dump.add(var6, #var6, 1);}

#define SAVE_REPEAT(...)                        \
  { MACRO_FOR_EACH(SAVE_REPEAT1, __VA_ARGS__) }

#endif // __DATAFILE_H__
