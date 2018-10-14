/*******************************************************************************
 * Functions to read input variables on a grid/mesh
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 *******************************************************************************/

#ifndef __GRIDDATA_DIRECTION__
#define __GRIDDATA_DIRECTION__
/// enum to use to specify Direction
enum Direction { X = 1, Y = 2, Z = 3 };
#endif //__GRIDDATA_DIRECTION__

class GridDataSource;

#ifndef __GRIDDATA_H__
#define __GRIDDATA_H__

#include "options.hxx"

#include "bout_types.hxx"
#include "dataformat.hxx"

#include "mesh.hxx"

#include <field2d.hxx>
#include <field3d.hxx>

#include <list>

/// Interface class to serve grid data
/*!
 * Provides a generic interface for sources of
 * equilibrium data.
 * Could be used to simplify interfacing between BOUT++ and other codes
 */
class GridDataSource {
public:
  virtual ~GridDataSource() {}

  virtual bool hasVar(const string &name) = 0; ///< Test if source can supply a variable

  virtual bool get(Mesh *m, int &ival, const string &name) = 0; ///< Get an integer
  virtual bool get(Mesh *m, BoutReal &rval,
                   const string &name) = 0; ///< Get a BoutReal number
  virtual bool get(Mesh *m, Field2D &var, const string &name, BoutReal def = 0.0) = 0;
  virtual bool get(Mesh *m, Field3D &var, const string &name, BoutReal def = 0.0) = 0;

  virtual bool get(Mesh *m, vector<int> &var, const string &name, int len, int offset = 0,
                   Direction dir = Direction::X) = 0;
  virtual bool get(Mesh *m, vector<BoutReal> &var, const string &name, int len,
                   int offset = 0, Direction dir = Direction::X) = 0;

  /// Test if grid data source includes y-boundary guard cells.
  /// Older grid files may not: this then requires some special handling.
  virtual bool hasYGuards() { return true; }
};

/// Interface to grid data in a file
/*!
 * This is a thin wrapper around a DataFormat object. Only needs to implement
 * reading routines.
 */
class GridFile : public GridDataSource {
public:
  GridFile() = delete;
  GridFile(std::unique_ptr<DataFormat> format, string gridfilename);
  ~GridFile() override;

  bool hasVar(const string &name) override;

  bool get(Mesh *m, int &ival, const string &name) override; ///< Get an integer
  bool get(Mesh *m, BoutReal &rval,
           const string &name) override; ///< Get a BoutReal number
  bool get(Mesh *m, Field2D &var, const string &name, BoutReal def = 0.0) override;
  bool get(Mesh *m, Field3D &var, const string &name, BoutReal def = 0.0) override;

  bool get(Mesh *m, vector<int> &var, const string &name, int len, int offset = 0,
           Direction dir = Direction::X) override;
  bool get(Mesh *m, vector<BoutReal> &var, const string &name, int len, int offset = 0,
           Direction dir = Direction::X) override;

  // Grid files from hypnotoad don't have y-guard cells at all. But even if
  // they did, they don't have guard cells for the y-boundaries of the core
  // region because the grid structure means the grid cells adjacent to the
  // first/last core cells in the input arrays in the file store data for the
  // PF region. So until a GridFile supports input with different domains, it
  // cannot contain *all* y-guard cells
  bool hasYGuards() override { return false; }

private:
  std::unique_ptr<DataFormat> file;
  string filename;

  bool readgrid_3dvar_fft(Mesh *m, const string &name, int yread, int ydest, int ysize,
                          int xge, int xlt, Field3D &var);

  bool readgrid_3dvar_real(Mesh *m, const string &name, int yread, int ydest, int ysize,
                           int xge, int xlt, Field3D &var);
};

/*!
 * Provides a way to create variables from options, which can
 * be set in the input file or on the command line. This is done
 * using FieldFactory to convert string expressions into fields.
 */
class GridFromOptions : public GridDataSource {
public:
  /*!
   * Constructor, passing optional Options object
   *
   * @param[in] opt  Options section to use as input. By default the "mesh" section under
   * root will be used.
   */
  GridFromOptions(Options *opt = nullptr) : options(opt) {}

  /*!
   * Checks if the options has a given variable
   */
  bool hasVar(const string &name) override;

  /*!
   * Reads integers from options. Uses Options::get to handle
   * expressions
   *
   * @param[in] mesh   Not used
   * @param[in] name   Name of variable
   * @param[out] ival  Always given a value, defaults to 0
   *
   * @return True if option is set, false if ival is default (0)
   */
  bool get(Mesh *mesh, int &ival, const string &name) override;

  /*!
   * Reads BoutReal from options. Uses Options::get to handle
   * expressions
   *
   * @param[in] mesh   Not used
   * @param[in] name   Name of variable
   * @param[out] rval  Always given a value, defaults to 0
   *
   * @return True if option is set, false if ival is default (0)
   */
  bool get(Mesh *mesh, BoutReal &rval, const string &name) override;

  /*!
   * Get a Field2D object by finding the option with the given name,
   * and passing the string to FieldFactory
   *
   * @param[in] mesh  The Mesh object over which the field is defined
   * @param[out] var  The variable which will be set
   * @param[in] name  The name in the options. Not case sensitive
   * @param[in] def   Default value to use if option not found
   */
  bool get(Mesh *mesh, Field2D &var, const string &name, BoutReal def = 0.0) override;

  /*!
   * Get a Field3D object by finding the option with the given name,
   * and passing the string to FieldFactory
   *
   * @param[in] mesh  The Mesh object over which the field is defined
   * @param[out] var  The variable which will be set
   * @param[in] name  The name in the options. Not case sensitive
   * @param[in] def   Default value to use if option not found
   */
  bool get(Mesh *mesh, Field3D &var, const string &name, BoutReal def = 0.0) override;

  /*!
   * Get an array of integers. Currently reads a single
   * integer, and sets the whole array to the same value
   *
   * @param[in] mesh  Mesh object
   * @param[out] var  A vector which will be resized to length len
   * @param[in] name  The name of the option
   * @param[in] len   The length of the vector
   * @param[in] offset Not currently used
   * @param[in] dir  The direction (X,Y,Z) of the array
   */
  bool get(Mesh *mesh, vector<int> &var, const string &name, int len, int offset = 0,
           Direction dir = Direction::X) override;

  /*!
   * Get an array of BoutReals. Uses FieldFactory to generate
   * an expression, then evaluates it at indices depending
   * on the direction (dir) and length (len)
   *
   * @param[in] mesh  Mesh object
   * @param[out] var  A vector which will be resized to length len
   * @param[in] name  The name of the option
   * @param[in] len   The length of the vector
   * @param[in] offset The index where this vector starts i.e. var[0] is at x=offset if
   *                   dir is X.
   * @param[in] dir  The direction (X,Y,Z) of the array
   */
  bool get(Mesh *mesh, vector<BoutReal> &var, const string &name, int len, int offset = 0,
           Direction dir = Direction::X) override;

private:
  /// The options section to use. Could be nullptr
  Options *options;
};

#endif // __GRIDDATA_H__
