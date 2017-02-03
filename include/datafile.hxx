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
#include "field2d.hxx"
#include "field3d.hxx"
#include "vector2d.hxx"
#include "vector3d.hxx"
#include "options.hxx"

#include "dataformat.hxx"

#include <stdarg.h>
#include <stdio.h>

#include <vector>
#include <string>

/*!
  Uses a generic interface to file formats (DataFormat)
  and provides an interface for reading/writing simulation data.
*/
class Datafile {
 public:
  Datafile(Options *opt = NULL);
  Datafile(const Datafile &other);
  ~Datafile(); // need to delete filename
  
  Datafile& operator=(const Datafile &rhs);

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
  
 private:
  bool parallel; // Use parallel formats?
  bool flush;    // Flush after every write?
  bool guards;   // Write guard cells?
  bool floats;   // Low precision?
  bool openclose; // Open and close file for each write
  int Lx,Ly,Lz; // The sizes in the x-, y- and z-directions of the arrays to be written
  bool enabled;  // Enable / Disable writing
  bool init_missing; // Initialise missing variables?
  bool shiftOutput; //Do we want to write out in shifted space?

  DataFormat *file;
  int filenamelen;
  static const int FILENAMELEN=512;
  char *filename;
  bool appending;

  /// A structure to hold a pointer to a class, and associated name and flags
  template <class T>
    struct VarStr {
      T *ptr;
      string name;
      bool save_repeat;
      bool covar;
    };

  // one set per variable type
  vector< VarStr<int> >      int_arr;
  vector< VarStr<BoutReal> > BoutReal_arr;
  vector< VarStr<Field2D> >  f2d_arr;
  vector< VarStr<Field3D> >  f3d_arr;
  vector< VarStr<Vector2D> > v2d_arr;
  vector< VarStr<Vector3D> > v3d_arr;

  bool read_f2d(const string &name, Field2D *f, bool save_repeat);
  bool read_f3d(const string &name, Field3D *f, bool save_repeat);

  bool write_int(const string &name, int *f, bool save_repeat);
  bool write_real(const string &name, BoutReal *f, bool save_repeat);
  bool write_f2d(const string &name, Field2D *f, bool save_repeat);
  bool write_f3d(const string &name, Field3D *f, bool save_repeat);

  bool varAdded(const string &name); // Check if a variable has already been added
};

/// Write this variable once to the grid file
#define SAVE_ONCE(var) dump.add(var, #var, 0)
#define SAVE_ONCE2(var1, var2) { \
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0);}
#define SAVE_ONCE3(var1, var2, var3) {\
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0);}
#define SAVE_ONCE4(var1, var2, var3, var4) { \
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0); \
    dump.add(var4, #var4, 0);}
#define SAVE_ONCE5(var1, var2, var3, var4, var5) {\
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0); \
    dump.add(var4, #var4, 0); \
    dump.add(var5, #var5, 0);}
#define SAVE_ONCE6(var1, var2, var3, var4, var5, var6) {\
    dump.add(var1, #var1, 0); \
    dump.add(var2, #var2, 0); \
    dump.add(var3, #var3, 0); \
    dump.add(var4, #var4, 0); \
    dump.add(var5, #var5, 0); \
    dump.add(var6, #var6, 0);}

/// Write this variable every timestep
#define SAVE_REPEAT(var) dump.add(var, #var, 1)
#define SAVE_REPEAT2(var1, var2) { \
    dump.add(var1, #var1, 1); \
    dump.add(var2, #var2, 1);}
#define SAVE_REPEAT3(var1, var2, var3) {\
    dump.add(var1, #var1, 1); \
    dump.add(var2, #var2, 1); \
    dump.add(var3, #var3, 1);}
#define SAVE_REPEAT4(var1, var2, var3, var4) { \
    dump.add(var1, #var1, 1); \
    dump.add(var2, #var2, 1); \
    dump.add(var3, #var3, 1); \
    dump.add(var4, #var4, 1);}
#define SAVE_REPEAT5(var1, var2, var3, var4, var5) {\
    dump.add(var1, #var1, 1); \
    dump.add(var2, #var2, 1); \
    dump.add(var3, #var3, 1); \
    dump.add(var4, #var4, 1); \
    dump.add(var5, #var5, 1);}
#define SAVE_REPEAT6(var1, var2, var3, var4, var5, var6) {\
    dump.add(var1, #var1, 1); \
    dump.add(var2, #var2, 1); \
    dump.add(var3, #var3, 1); \
    dump.add(var4, #var4, 1); \
    dump.add(var5, #var5, 1); \
    dump.add(var6, #var6, 1);}

#endif // __DATAFILE_H__
