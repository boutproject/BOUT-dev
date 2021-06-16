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

#include <fmt/core.h>

#include <cstdarg>
#include <cstdio>
class Mesh;
class Field;
class Field2D;
class Field3D;
class FieldPerp;
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
  ~Datafile() = default;
  
  Datafile& operator=(Datafile &&rhs) noexcept;
  Datafile& operator=(const Datafile &rhs) = delete;

  /// Open read-only
  bool openr(const std::string& filename);

  template <class S, class... Args>
  bool openr(const S& format, const Args&... args) {
    return openr(fmt::format(format, args...));
  }

  /// Overwrites existing file
  bool openw(const std::string& filename);

  template <class S, class... Args>
  bool openw(const S& format, const Args&... args) {
    return openw(fmt::format(format, args...));
  }

  /// Appends if exists
  bool opena(const std::string& filename);

  template <class S, class... Args>
  bool opena(const S& format, const Args&... args) {
    return opena(fmt::format(format, args...));
  }

  
  bool isValid();  // Checks if the data source is valid

  void close();

  void setLowPrecision(); ///< Only output floats
  template <typename T>
  void addRepeat(T& value, std::string name) {
    add(value, name.c_str(), true);
  }
  template <typename T>
  void addOnce(T& value, std::string name) {
    add(value, name.c_str(), false);
  }
  void add(int &i, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(std::vector<int> &ivec, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(std::string &s, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(BoutReal &r, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(bool &b, const char* name, bool save_repeat = false,
           const std::string &description = "");
  void add(Field2D &f, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(Field3D &f, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(FieldPerp &f, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(Vector2D &f, const char *name, bool save_repeat = false,
           const std::string &description = "");
  void add(Vector3D &f, const char *name, bool save_repeat = false,
           const std::string &description = "");
  
  bool read();  ///< Read data into added variables 
  bool write(); ///< Write added variables

  /// Opens, writes, closes file
  bool write(const std::string& filename) const;
  template <typename S, typename... Args>
  bool write(const S& format, const Args&... args) const {
    return write(fmt::format(format, args...));
  }

  void setAttribute(const std::string &varname, const std::string &attrname, const std::string &text);
  void setAttribute(const std::string &varname, const std::string &attrname, int value);
  void setAttribute(const std::string &varname, const std::string &attrname, BoutReal value);

 private:
  Mesh* mesh;
  bool parallel{false}; // Use parallel formats?
  bool flush{true};     // Flush after every write?
  bool guards{true};    // Write guard cells?
  bool floats{false};   // Low precision?
  bool openclose{true}; // Open and close file for each write
  int Lx,Ly,Lz; // The sizes in the x-, y- and z-directions of the arrays to be written
  bool enabled{true}; // Enable / Disable writing
  bool init_missing; // Initialise missing variables?
  bool shiftoutput{false}; // Do we want to write out in shifted space?
  bool shiftinput{false};  // Read in shifted space?
  // Counter used in determining when next openclose required
  int flushFrequencyCounter{0};
  int flushfrequency{1}; // How many write calls do we want between openclose

  std::unique_ptr<DataFormat> file;
  std::string filename;
  bool writable{false}; // is file open for writing?
  bool appending{false};
  bool first_time{true}; // is this the first time the data will be written?

  /// Shallow copy, not including dataformat, therefore private
  Datafile(const Datafile& other);

  /// A structure to hold a pointer to a class, and associated name and flags
  template <class T>
  struct VarStr {
    T *ptr;                       ///< Pointer to the data.
                                  ///< Note that this may be a user object, not a copy, so must not be destroyed
    std::string name;             ///< Name as it appears in the output file
    bool save_repeat;             ///< If true, has a time dimension and is saved every time step
    bool covar;                   ///< For vectors, true if a covariant vector, false if contravariant
    size_t size;                  ///< Size of a stored vector or string, to check it does not change after being added
    std::string description{""};  ///< Documentation of what the variable is
  };

  // one set per variable type
  std::vector<VarStr<int>> int_arr;
  std::vector<VarStr<std::vector<int>>> int_vec_arr;
  std::vector<VarStr<std::string>> string_arr;
  std::vector<VarStr<BoutReal>> BoutReal_arr;
  std::vector<VarStr<bool>> bool_arr;
  std::vector<VarStr<Field2D>> f2d_arr;
  std::vector<VarStr<Field3D>> f3d_arr;
  std::vector<VarStr<FieldPerp>> fperp_arr;
  std::vector<VarStr<Vector2D>> v2d_arr;
  std::vector<VarStr<Vector3D>> v3d_arr;

  bool read_f2d(const std::string &name, Field2D *f, bool save_repeat);
  bool read_f3d(const std::string &name, Field3D *f, bool save_repeat);
  bool read_fperp(const std::string &name, FieldPerp *f, bool save_repeat);

  bool write_int(const std::string &name, int *f, bool save_repeat);
  bool write_int_vec(const std::string &name, std::vector<int> *f, bool save_repeat);
  bool write_string(const std::string &name, std::string *f, bool save_repeat);
  bool write_real(const std::string &name, BoutReal *f, bool save_repeat);
  bool write_f2d(const std::string &name, Field2D *f, bool save_repeat);
  bool write_f3d(const std::string &name, Field3D *f, bool save_repeat);
  bool write_fperp(const std::string &name, FieldPerp *f, bool save_repeat);

  /// Check if a variable has already been added
  bool varAdded(const std::string &name);

  /// Get the pointer to the variable, nullptr if not added
  /// This is used to check if the same variable is being added
  void* varPtr(const std::string &name);
};

#endif // __DATAFILE_H__
