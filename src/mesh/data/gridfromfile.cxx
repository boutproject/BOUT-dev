
#include "bout/traits.hxx"
#include <bout/griddata.hxx>

#include <bout/constants.hxx>
#include <bout/sys/timer.hxx>
#include <boutexception.hxx>
#include <fft.hxx>
#include <msg_stack.hxx>
#include <options_netcdf.hxx>
#include <output.hxx>
#include <unused.hxx>
#include <utility>
#include <utils.hxx>

GridFile::GridFile(std::string gridfilename)
  : GridDataSource(true), data(bout::OptionsNetCDF(gridfilename).read()), filename(std::move(gridfilename)) {
  TRACE("GridFile constructor");

  // Get number of y-boundary guard cells saved in the grid file
  grid_yguards = data["y_boundary_guards"].withDefault<int>(0);
  // Get number ny_inner from the grid file.
  // Is already read in BoutMesh, but this way we don't have to the Mesh API to
  // get it from there.
  ny_inner = data["ny_inner"].withDefault<int>(0);
}

/*!
 * Tests whether a variable exists in the file
 *
 */
bool GridFile::hasVar(const std::string& name) { return data.isSet(name); }

/*!
 * Read a string from file. If the string is not
 * found, then string is set to "" and false is returned.
 *
 * Inputs
 * ------
 *
 *   m     Pointer to mesh, not used
 *   name  String containing name of variable
 *
 * Outputs
 * -------
 *
 *   sval   Reference to string
 *
 * Returns
 * -------
 *
 *   Boolean. True on success.
 *
 */
bool GridFile::get(Mesh* UNUSED(m), std::string& sval, const std::string& name,
                   const std::string& def) {
  Timer timer("io");
  TRACE("GridFile::get(std::string)");
  const bool success = data.isSet(name);
  if (not success) {
    // Override any previously set defaults
    data[name].force(def, data.getDefaultSource());
  }
  sval = data[name];
  return success;
}

/*!
 * Read a single integer from file. If the integer is not
 * found, then ival is set to zero and false is returned.
 *
 * Inputs
 * ------
 *
 *   m     Pointer to mesh, not used
 *   name  String containing name of variable
 *
 * Outputs
 * -------
 *
 *   ival   Reference to integer
 *
 * Returns
 * -------
 *
 *   Boolean. True on success.
 *
 */
bool GridFile::get(Mesh* UNUSED(m), int& ival, const std::string& name, int def) {
  Timer timer("io");
  TRACE("GridFile::get(int)");
  const bool success = data.isSet(name);
  if (not success) {
    // Override any previously set defaults
    data[name].force(def);
  }
  ival = data[name];
  return success;
}

/*!
 *
 *
 */
bool GridFile::get(Mesh *UNUSED(m), BoutReal &rval, const std::string &name, BoutReal def) {
  Timer timer("io");
  TRACE("GridFile::get(BoutReal)");
  const bool success = data.isSet(name);
  if (not success) {
    // Override any previously set defaults
    data[name].force(def);
  }
  rval = data[name];
  return success;
}

/*!
 * Reads a 2D, 3D or FieldPerp field variable from a file
 * 
 * Successfully reads Field2D or FieldPerp if the variable in the file is 0-D or 2-D.
 * Successfully reads Field3D if the variable in the file is 0-D, 2-D or 3-D.
 */
bool GridFile::get(Mesh *m, Field2D &var, const std::string &name, BoutReal def, CELL_LOC location) {
  return getField(m, var, name, def, location);
}
bool GridFile::get(Mesh *m, Field3D &var, const std::string &name, BoutReal def, CELL_LOC location) {
  return getField(m, var, name, def, location);
}

namespace {
/// Visitor that returns the shape of its argument
struct GetDimensions {
  std::vector<int> operator()(MAYBE_UNUSED(bool value)) { return {1}; }
  std::vector<int> operator()(MAYBE_UNUSED(int value)) { return {1}; }
  std::vector<int> operator()(MAYBE_UNUSED(BoutReal value)) { return {1}; }
  std::vector<int> operator()(MAYBE_UNUSED(const std::string& value)) { return {1}; }
  std::vector<int> operator()(const Array<BoutReal>& array) { return {array.size()}; }
  std::vector<int> operator()(const Matrix<BoutReal>& array) {
    const auto shape = array.shape();
    return {std::get<0>(shape), std::get<1>(shape)};
  }
  std::vector<int> operator()(const Tensor<BoutReal>& array) {
    const auto shape = array.shape();
    return {std::get<0>(shape), std::get<1>(shape), std::get<2>(shape)};
  }
  std::vector<int> operator()(const Field& array) {
    return {array.getNx(), array.getNy(), array.getNz()};
  }
};
} // namespace

template <typename T>
bool GridFile::getField(Mesh* m, T& var, const std::string& name, BoutReal def,
                        CELL_LOC location) {
  static_assert(bout::utils::is_Field<T>::value,
                "templated GridFile::get only works for Field2D, Field3D or FieldPerp");

  Timer timer("io");
  AUTO_TRACE();

  if (not data.isSet(name)) {
    // Variable not found
    output_warn.write("\tWARNING: Could not read '{:s}' from grid. Setting to {:e}\n", name, def);
    var = def;
    var.setLocation(location);
    return false;
  }

  Options option = data[name];

  // Global (x, y, z) dimensions of field
  const std::vector<int> size = bout::utils::visit(GetDimensions{}, option.value);

  switch(size.size()) {
  case 1: {
    // 0 or 1 dimension
    if (size[0] != 1) {
      throw BoutException(
          "Expecting a 2D variable, but '{:s}' is 1D with {:d} elements\n", name,
          size[0]);
    }
    var = option.as<BoutReal>();
    var.setLocation(location);
    return true;
  }
  case 2: {
    // Check size
    break;
  }
  case 3: {
    // Check size if getting Field3D
    if (bout::utils::is_Field2D<T>::value or bout::utils::is_FieldPerp<T>::value) {
      output_warn.write("WARNING: Variable '{:s}' should be 2D, but has {:d} dimensions. Ignored\n", name, size.size());
      var = def;
      var.setLocation(location);
      return false;
    }
    break;
  }
  default: {
    output_warn.write("WARNING: Variable '{:s}' should be 2D or 3D, but has {:d} dimensions. Ignored\n", name, size.size());
    var = def;
    var.setLocation(location);
    return false;
  }
  };

  ///Ghost region widths.
  const int mxg = (m->LocalNx - (m->xend - m->xstart + 1)) / 2;
  const int myg = (m->LocalNy - (m->yend - m->ystart + 1)) / 2;
  ///Check that ghost region widths are in fact integers
  ASSERT1((m->LocalNx - (m->xend - m->xstart + 1)) % 2 == 0);
  ASSERT1((m->LocalNy - (m->yend - m->ystart + 1)) % 2 == 0);

  // Index offsets into source array
  int xs = m->OffsetX;
  // Need to increase offset by 2*(# boundary guards) for each target position
  // we pass
  int ys = m->OffsetY;

  // Total number of y-boundary cells in grid file, used for check later.
  // Value depends on if we are double-null or not.
  int total_grid_yguards = 2*grid_yguards;
  if (m->numberOfXPoints > 1) {
    ASSERT1(m->numberOfXPoints == 2);
    // Need to check if we are before or after the target in the middle of the
    // y-domain, and increase ys for the extra boundary guard cells at that
    // target if we are after it.
    if (m->OffsetY >= ny_inner) {
      // Note: neither ny_inner nor OffsetY include guard cells
      ys += 2*grid_yguards;
    }

    // Add y-boundary guard cells at upper target
    total_grid_yguards += 2*grid_yguards;
  }

  // Index offsets into destination
  int xd = -1;
  int yd = -1;

  // Number of points to read.
  int nx_to_read = -1;
  int ny_to_read = -1;

  ///Check if field dimensions are correct. x-direction
  int grid_xguards = (size[0] - (m->GlobalNx - 2*mxg)) / 2;
  // Check there is no rounding in calculation of grid_xguards
  ASSERT1( (size[0] - (m->GlobalNx - 2*mxg)) % 2 == 0 );
  if (grid_xguards >= 0) { ///including ghostpoints
    nx_to_read = m->LocalNx;
    xd = grid_xguards - mxg;
    ASSERT1(xd >= 0);
  } else if (grid_xguards == 0) { ///excluding ghostpoints
    nx_to_read = m->LocalNx - 2*mxg;
    xd = mxg;
  } else {
    throw BoutException(
        "Could not read '{:s}' from file: number of x-boundary guard cells "
        "in the grid file grid_xguards={:d} neither matches grid_xguards >= mxg={:d} "
        "nor grid_xguards = 0",
        name, grid_xguards, mxg);
  }

  if (not bout::utils::is_FieldPerp<T>::value) {
    // Check if field dimensions are correct. y-direction
    if (grid_yguards > 0) {
      // including ghostpoints
      ASSERT1(size[1] == m->GlobalNy);
      ny_to_read = m->LocalNy;
      yd = grid_yguards - myg;
      ASSERT1(yd >= 0);
    } else if (grid_yguards == 0) {
      // excluding ghostpoints
      ASSERT1(size[1] == m->GlobalNy - m->numberOfYBoundaries()*2*myg);
      ny_to_read = m->LocalNy - 2*myg;
      yd = myg;
    } else {
      throw BoutException(
          "Could not read '{:s}' from file: number of y-boundary guard cells "
          "in the grid file grid_yguards={:d} neither matches grid_yguards >= myg={:d} "
          "nor grid_yguards = 0",
          name, grid_yguards, myg);
    }
  }

  // Now read data from file
  readField(m, name, ys, yd, ny_to_read, xs, xd, nx_to_read, size, var);

  if (location != CELL_DEFAULT and var.getLocation() != location) {
    throw BoutException("Incorrect location of field {:s} in grid file, expected {:s}, "
                        "got {:s}.", name, toString(location),
                        toString(var.getLocation()));
  }

  if (var.isAllocated()) {
    // FieldPerps might not be allocated if they are not read on this processor

    // If field does not include ghost points in x-direction ->
    // Upper and lower X boundaries copied from nearest point
    if (size[0] == m->GlobalNx - 2*mxg ) {
      for (int x=0; x<m->xstart; x++) {
        for (int y=0; y<m->LocalNy; y++) {
          for (int z=0; z<var.getNz(); z++) {
            var(x, y, z) = var(m->xstart, y, z);
          }
        }
      }
      for (int x=m->xend+1;x<m->LocalNx;x++) {
        for (int y=0; y<m->LocalNy; y++) {
          for (int z=0; z<var.getNz(); z++) {
            var(x, y, z) = var(m->xend, y, z);
          }
        }
      }
    }

    if (not bout::utils::is_FieldPerp<T>::value) {
      ///If field does not include ghost points in y-direction ->
      ///Upper and lower Y boundaries copied from nearest point
      if (grid_yguards == 0) {
        for(int x=0; x<m->LocalNx; x++) {
          for(int y=0; y<m->ystart; y++) {
            for (int z=0; z<var.getNz(); z++) {
              var(x, y, z) = var(x, m->ystart, z);
            }
          }
          for(int y=m->yend+1; y<m->LocalNy; y++) {
            for (int z=0; z<var.getNz(); z++) {
              var(x, y, z) = var(x, m->yend, z);
            }
          }
        }
      }
    }
  }

  return true;
}

namespace {
void readFieldAttributes(Options& field_options, Field& f) {
  if (field_options.hasAttribute("cell_location")) {
    const std::string location_string = field_options.attributes.at("cell_location");
    f.setLocation(CELL_LOCFromString(location_string));
  }

  if (field_options.hasAttribute("direction_y")) {
    const std::string direction_y_string = field_options.attributes.at("direction_y");
    f.setDirectionY(YDirectionTypeFromString(direction_y_string));
  }

  if (field_options.hasAttribute("direction_z")) {
    const std::string direction_z_string = field_options.attributes.at("direction_z");
    f.setDirectionZ(ZDirectionTypeFromString(direction_z_string));
  }
}
}

void GridFile::readField(Mesh* UNUSED(m), const std::string& name, int ys, int yd,
                         int ny_to_read, int xs, int xd, int nx_to_read,
                         const std::vector<int>& UNUSED(size), Field2D& var) {

  readFieldAttributes(data[name], var);

  var.allocate();

  const auto full_var = data[name].as<Matrix<BoutReal>>();

  for (int x = xs; x < xs + nx_to_read; ++x) {
    for (int y = ys; y < ys + ny_to_read; ++y) {
      var(x - xs + xd, y - ys + yd) = full_var(x, y);
    }
  }
}

void GridFile::readField(Mesh* m, const std::string& name, int ys, int yd,
    int ny_to_read, int xs, int xd, int nx_to_read, const std::vector<int>& size,
    Field3D& var) {

  readFieldAttributes(data[name], var);

  var.allocate();

  // Check whether "nz" is defined
  if (hasVar("nz")) {
    // Check the array is the right size
    if (size[2] != m->LocalNz) {
      throw BoutException("3D variable '{:s}' has incorrect size {:d} (expecting {:d})",
                          name, size[2], m->LocalNz);
    }

    if (!readgrid_3dvar_real(name,
          ys,// Start reading at global y-index
          yd,// Insert data starting from y=yd
          ny_to_read,// Length of data in Y
          xs,// Start reading at global x-index
          xd,// Insert data starting from x=xd
          nx_to_read, // Length of data in X
          var) ) {
      throw BoutException("\tWARNING: Could not read '{:s}' from grid. Setting to zero\n",
                          name);
    }
  } else {
    // No Z size specified in file. Assume FFT format
    if (!readgrid_3dvar_fft(m, name,
          ys,// Start reading at global y-index
          yd,// Insert data starting from y=yd
          ny_to_read,// Length of data in Y
          xs,// Start reading at global x-index
          xd,// Insert data starting from x=xd
          nx_to_read, // Length of data in X
          var) ) {
      throw BoutException("\tWARNING: Could not read '{:s}' from grid. Setting to zero\n",
                          name);
    }
  }
}

void GridFile::readField(Mesh* m, const std::string& name, int UNUSED(ys), int UNUSED(yd),
    int UNUSED(ny_to_read), int xs, int xd, int nx_to_read, const std::vector<int>& size,
    FieldPerp& var) {

  readFieldAttributes(data[name], var);

  int yindex = var.getIndex();

  if (yindex >= 0 and yindex <= m->LocalNy) {
    // Only read if yindex is on this processor

    var.allocate();

    // Check whether "nz" is defined
    if (hasVar("nz")) {
      // Check the array is the right size
      if (size[2] != m->LocalNz) {
        throw BoutException(
            "FieldPerp variable '{:s}' has incorrect size {:d} (expecting {:d})", name,
            size[2], m->LocalNz);
      }

      if (!readgrid_perpvar_real(name,
            xs,// Start reading at global x-index
            xd,// Insert data starting from x=xd
            nx_to_read, // Length of data in X
            var) ) {
        throw BoutException(
            "\tWARNING: Could not read '{:s}' from grid. Setting to zero\n", name);
      }
    } else {
      // No Z size specified in file. Assume FFT format
      if (!readgrid_perpvar_fft(m, name,
            xs,// Start reading at global x-index
            xd,// Insert data starting from x=xd
            nx_to_read, // Length of data in X
            var) ) {
        throw BoutException(
            "\tWARNING: Could not read '{:s}' from grid. Setting to zero\n", name);
      }
    }
  }
}

bool GridFile::get(MAYBE_UNUSED(Mesh* m), MAYBE_UNUSED(std::vector<int>& var),
                   MAYBE_UNUSED(const std::string& name), MAYBE_UNUSED(int len),
                   MAYBE_UNUSED(int offset),
                   MAYBE_UNUSED(GridDataSource::Direction dir)) {
  TRACE("GridFile::get(vector<int>)");

  return false;
}

bool GridFile::get(Mesh *UNUSED(m), std::vector<BoutReal> &var, const std::string &name,
                   int len, int offset, GridDataSource::Direction UNUSED(dir)) {
  TRACE("GridFile::get(vector<BoutReal>)");
  
  if (not data.isSet(name)) { return false;}

  const auto full_var = data[name].as<Array<BoutReal>>();
  const auto* it = std::begin(full_var);
  std::advance(it, offset);
  std::copy_n(it, len, std::begin(var));

  return true;
}

bool GridFile::hasXBoundaryGuards(Mesh* m) {
  // Global (x,y) dimensions of some field
  // a grid file should always contain "dx"
  Options option = data["dx"];
  const std::vector<int> size = bout::utils::visit(GetDimensions{}, option.value);

  if (size.empty()) {
    // handle case where "dx" is not present - non-standard grid file
    // - e.g. for tests
    return false;
  }

  return size[0] > m->GlobalNx - 2*m->xstart;
}

/////////////////////////////////////////////////////////////
// Private routines

/// Reads in a portion of the X-Y domain
/*
  Data stored as toroidal FFTs in BoutReal space at each X-Y point.
  In toroidal direction, array must have an odd number of points.
  Format is:

  DC, r1,i1, r2,i2, ... , rn,in

  with the BoutReal and imaginary parts of each (positive) frequency
  up to the nyquist frequency.
 */
bool GridFile::readgrid_3dvar_fft(Mesh *m, const std::string &name, 
				 int yread, int ydest, int ysize, 
				 int xread, int xdest, int xsize, Field3D &var) {
  /// Check the arguments make sense
  if ((yread < 0) || (ydest < 0) || (ysize < 0) || (xread < 0) || (xdest < 0)
      || (xsize < 0)) {
    return false;
  }
  
  /// Check the size of the data
  Options option = data[name];
  const std::vector<int> size = bout::utils::visit(GetDimensions{}, option.value);
  
  if (size.size() != 3) {
    output_warn.write("\tWARNING: Number of dimensions of {:s} incorrect\n", name);
    return false;
  }

  int maxmode = (size[2] - 1)/2; ///< Maximum mode-number n

  int ncz = m->LocalNz;

  /// we should be able to replace the following with
  /// var.getCoordinates()->zlength();
  /// but don't do it yet as we don't assert that m == var.getMesh()
  /// Expect the assertion to be true, in which case we probably don't
  /// need to pass m as can just use var.getMesh()
  BoutReal zlength = m->getCoordinates(var.getLocation())->zlength();
  
  int zperiod = ROUND(TWOPI / zlength); /// Number of periods in 2pi

  // Print out which modes are going to be read in
  if (zperiod > maxmode) {
    // Domain is too small: Only DC
    output_warn.write("zperiod ({:d}) > maxmode ({:d}) => Only reading n = 0 component\n", zperiod, maxmode);
  } else {
    // Get maximum mode in the input which is a multiple of zperiod
    int mm = (maxmode / zperiod) * zperiod;
    if ( (ncz/2)*zperiod < mm )
      mm = (ncz/2)*zperiod; // Limited by Z resolution
    
    if (mm == zperiod) {
      output_info.write(" => Reading n = 0, {:d}\n", zperiod);
    } else {
      output_info.write(" => Reading n = 0, {:d} ... {:d}\n", zperiod, mm);
    }
  }

  /// Data for FFT. Only positive frequencies
  Array<dcomplex> fdata(ncz / 2 + 1);
  Array<BoutReal> zdata(size[2]);

  const auto full_var = data[name].as<Tensor<BoutReal>>();

  for (int jx = xread; jx < xread+xsize; jx++) {
    // jx is global x-index to start from
    for (int jy = yread; jy < yread+ysize; jy++) {
      // jy is global y-index to start from
      for (int jz = 0; jz < size[2]; ++jz) {
        zdata[jz] = full_var(jx, jy, jz);
      }

      /// Load into dcomplex array
      fdata[0] = zdata[0]; // DC component

      for(int i=1;i<=ncz/2;i++) {
        int modenr = i*zperiod; // Z mode number

        if (modenr <= maxmode) {
          // Have data for this mode
          fdata[i] = dcomplex(zdata[modenr*2 - 1], zdata[modenr*2]);
        } else {
          fdata[i] = 0.0;
        }
      }
      irfft(std::begin(fdata), ncz, &var(jx-xread+xdest, jy-yread+ydest, 0));
    }
  }

  return true;
}

/*!
 * Reads a 3D variable directly from the file, without
 * any processing
 */
bool GridFile::readgrid_3dvar_real(const std::string& name, int yread, int ydest,
                                   int ysize, int xread, int xdest, int xsize,
                                   Field3D& var) {
  /// Check the arguments make sense
  if ((yread < 0) || (ydest < 0) || (ysize < 0) || (xread < 0) || (xdest < 0)
      || (xsize < 0)) {
    return false;
  }

  Options option = data[name];

  /// Check the size of the data
  const std::vector<int> size = bout::utils::visit(GetDimensions{}, option.value);

  if (size.size() != 3) {
    output_warn.write("\tWARNING: Number of dimensions of {:s} incorrect\n", name);
    return false;
  }

  const auto full_var = option.as<Tensor<BoutReal>>();

  for (int jx = xread; jx < xread + xsize; jx++) {
    // jx is global x-index to start from
    for (int jy = yread; jy < yread + ysize; jy++) {
      // jy is global y-index to start from
      for (int jz = 0; jz < size[2]; ++jz) {
        var(jx - xread + xdest, jy - yread + ydest, jz) = full_var(jx, jy, jz);
      }
    }
  }

  return true;
}

/*
  Data stored as toroidal FFTs in BoutReal space at each X point.
  In toroidal direction, array must have an odd number of points.
  Format is:

  DC, r1,i1, r2,i2, ... , rn,in

  with the BoutReal and imaginary parts of each (positive) frequency
  up to the nyquist frequency.
 */
bool GridFile::readgrid_perpvar_fft(Mesh *m, const std::string &name,
    int xread, int xdest, int xsize, FieldPerp &var) {
  /// Check the arguments make sense
  if ((xread < 0) || (xdest < 0) || (xsize < 0)) {
    return false;
  }

  /// Check the size of the data
  Options option = data[name];
  const std::vector<int> size = bout::utils::visit(GetDimensions{}, option.value);

  if (size.size() != 2) {
    output_warn.write("\tWARNING: Number of dimensions of {:s} incorrect\n", name);
    return false;
  }

  int maxmode = (size[1] - 1)/2; ///< Maximum mode-number n

  int ncz = m->LocalNz;

  /// we should be able to replace the following with
  /// var.getCoordinates()->zlength();
  /// but don't do it yet as we don't assert that m == var.getMesh()
  /// Expect the assertion to be true, in which case we probably don't
  /// need to pass m as can just use var.getMesh()
  BoutReal zlength = m->getCoordinates(var.getLocation())->zlength();

  int zperiod = ROUND(TWOPI / zlength); /// Number of periods in 2pi

  // Print out which modes are going to be read in
  if (zperiod > maxmode) {
    // Domain is too small: Only DC
    output_warn.write("zperiod ({:d}) > maxmode ({:d}) => Only reading n = 0 component\n", zperiod, maxmode);
  } else {
    // Get maximum mode in the input which is a multiple of zperiod
    int mm = (maxmode / zperiod) * zperiod;
    if ( (ncz/2)*zperiod < mm )
      mm = (ncz/2)*zperiod; // Limited by Z resolution

    if (mm == zperiod) {
      output_info.write(" => Reading n = 0, {:d}\n", zperiod);
    } else {
      output_info.write(" => Reading n = 0, {:d} ... {:d}\n", zperiod, mm);
    }
  }

  /// Data for FFT. Only positive frequencies
  Array<dcomplex> fdata(ncz / 2 + 1);
  Array<BoutReal> zdata(size[1]);

  const auto full_var = option.as<Matrix<BoutReal>>();

  for (int jx = xread; jx < xread+xsize; jx++) {
    // jx is global x-index to start from
    for (int jz = 0; jz < size[1]; ++jz) {
      zdata[jz] = full_var(jx, jz);
    }

    /// Load into dcomplex array
    fdata[0] = zdata[0]; // DC component

    for(int i=1;i<=ncz/2;i++) {
      int modenr = i*zperiod; // Z mode number

      if (modenr <= maxmode) {
        // Have data for this mode
        fdata[i] = dcomplex(zdata[modenr*2 - 1], zdata[modenr*2]);
      } else {
        fdata[i] = 0.0;
      }
    }
    irfft(std::begin(fdata), ncz, &var(jx-xread+xdest, 0));
  }

  return true;
}

/*!
 * Reads a FieldPerp variable directly from the file, without
 * any processing
 */
bool GridFile::readgrid_perpvar_real(const std::string &name,
    int xread, int xdest, int xsize, FieldPerp &var) {
  /// Check the arguments make sense
  if ((xread < 0) || (xdest < 0) || (xsize < 0)) {
    return false;
  }

  /// Check the size of the data
  Options option = data[name];
  const std::vector<int> size = bout::utils::visit(GetDimensions{}, option.value);

  if (size.size() != 2) {
    output_warn.write("\tWARNING: Number of dimensions of {:s} incorrect\n", name);
    return false;
  }

  const auto full_var = option.as<Matrix<BoutReal>>();

  for (int jx = xread; jx < xread+xsize; jx++) {
    // jx is global x-index to start from
    for (int jz = 0; jz < size[1]; ++jz) {
      var(jx - xread + xdest, jz) = full_var(jx, jz);
    }
  }

  return true;
}

