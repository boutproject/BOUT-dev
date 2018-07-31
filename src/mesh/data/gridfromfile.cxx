
#include <bout/griddata.hxx>

#include <msg_stack.hxx>
#include <bout/sys/timer.hxx>

#include <boutexception.hxx>

#include <output.hxx>

#include <bout/constants.hxx>

#include <utils.hxx>  // for ROUND function

#include <fft.hxx>

#include <unused.hxx>

#include <utility>

/*!
 * Creates a GridFile object
 * 
 * format     Pointer to DataFormat. This will be deleted in
 *            destructor
 */
GridFile::GridFile(std::unique_ptr<DataFormat> format, string gridfilename)
    : file(std::move(format)), filename(std::move(gridfilename)) {
  TRACE("GridFile constructor");

  if (! file->openr(filename) ) {
    throw BoutException("Could not open file '%s'", filename.c_str());
  }

  file->setGlobalOrigin(); // Set default global origin

}

GridFile::~GridFile() {
  file->close();
}

/*!
 * Tests whether a variable exists in the file
 * 
 * Currently this is done by getting the variable's size,
 * and testing for zero size.
 */
bool GridFile::hasVar(const string &name) {
  if (!file->is_valid()) {
    return false;
  }
  
  /// Get the size of the variable
  vector<int> s = file->getSize(name);
  
  /// Test if the variable has zero size
  return s.size() != 0;
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
bool GridFile::get(Mesh *UNUSED(m), int &ival,      const string &name) {
  Timer timer("io");
  TRACE("GridFile::get(int)");
  
  if (!file->is_valid()) {
    throw BoutException("File cannot be read");
  }
  
  bool success = file->read(&ival, name);
  if (success) {
    output_info << "\tOption " << name  << " = " << ival << " (" << filename <<")" << endl;
  }

  return success;

}

/*!
 *
 *
 */
bool GridFile::get(Mesh *UNUSED(m), BoutReal &rval, const string &name) {
  Timer timer("io");
  TRACE("GridFile::get(BoutReal)");
  
  if (!file->is_valid()) {
    throw BoutException("File cannot be read");
  }
  bool success = file->read(&rval, name);
  if (success) {
    output_info << "\tOption " << name  << " = " << rval << " (" << filename <<")" << endl;
  }

  return success;

}

/*!
 * Reads a 2D field variable from a file
 * 
 * Succeeds if the variable in the file is 0-D or 2-D
 */
bool GridFile::get(Mesh *m, Field2D &var,   const string &name, BoutReal def) {
  Timer timer("io");
  TRACE("GridFile::get(Field2D)");

  if (!file->is_valid()) {
    throw BoutException("Could not read '%s' from file: File cannot be read", name.c_str());
  }
  vector<int> size = file->getSize(name);
  
  switch(size.size()) {
  case 0: {
    // Variable not found
    output_warn.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name.c_str(), def);
    var = def;
    return false;
  }
  case 1: {
    // 0 or 1 dimension
    if (size[0] != 1) {
      throw BoutException("Expecting a 2D variable, but '%s' is 1D with %d elements\n", name.c_str(), size[0]);
    }
    BoutReal rval;
    if (!file->read(&rval, name)) {
      throw BoutException("Couldn't read 0D variable '%s'\n", name.c_str());
    }
    var = rval;
    return true;
  }
  case 2: {
    // Check size
    break;
  }
  default: {
    output_warn.write("WARNING: Variable '%s' should be 2D, but has %d dimensions. Ignored\n", 
                      name.c_str(), size.size());
    var = def;
    return false;
  }
  };

  var.allocate(); // Make sure data allocated

  // Index offsets into source array
  int xs = m->OffsetX;
  int ys = m->OffsetY;

  // Index offsets into destination
  int xd = -1;
  int yd = -1;

  ///Ghost region widths.
  int mxg = (m->LocalNx - (m->xend - m->xstart + 1)) / 2;
  int myg = (m->LocalNy - (m->yend - m->ystart + 1)) / 2;
  ///Check that ghost region widths are in fact integers
  ASSERT1((m->LocalNx - (m->xend - m->xstart + 1)) % 2 == 0);
  ASSERT1((m->LocalNy - (m->yend - m->ystart + 1)) % 2 == 0);

  ///Global (x,y) dimensions of field
  const vector<int> field_dimensions = file->getSize(name);

  // Number of points to read.
  int nx_to_read = -1;
  int ny_to_read = -1;

  ///Check if field dimensions are correct. x-direction
  if (field_dimensions[0] == m->GlobalNx) { ///including ghostpoints
    nx_to_read = m->LocalNx;
    xd = 0;
  } else if ( field_dimensions[0] == m->GlobalNx - 2*mxg ) {///including ghostpoints
    nx_to_read = m->LocalNx - 2*mxg;
    xd = mxg;
  } else {
    throw BoutException("Could not read '%s' from file: x-dimension = %i do neither match nx = %i"
                "nor nx-2*mxg = %i ", name.c_str(), field_dimensions[0], m->GlobalNx, m->GlobalNx-2*mxg);
  }

  ///Check if field dimensions are correct. y-direction
  if (field_dimensions[1] == m->GlobalNy) { ///including ghostpoints
    ny_to_read = m->LocalNy;
    yd = 0;
  } else if ( field_dimensions[1] == m->GlobalNy - 2*myg ) {///including ghostpoints
    ny_to_read = m->LocalNy - 2*myg;
    yd = myg;
  } else {
    throw BoutException("Could not read '%s' from file: y-dimension = %i do neither match ny = %i"
                "nor ny-2*myg = %i ", name.c_str(), field_dimensions[1], m->GlobalNy, m->GlobalNy-2*myg);
  }

  ///Now read data from file
  for(int x=xs;x < xs+nx_to_read; x++) {
    file->setGlobalOrigin(x,ys,0);
    if (!file->read(&var(x-xs+xd, yd), name, 1, ny_to_read) ) {
      throw BoutException("Could not fetch data for '%s'", name.c_str());
    }
  }

  ///If field does not include ghost points in y-direction ->
  ///Upper and lower Y boundaries copied from nearest point
  if (field_dimensions[1] == m->GlobalNy - 2*myg ) {
    for(int x=0;x<m->LocalNx;x++) {
      for(int y=0;y<m->ystart;y++)
        var(x, y) = var(x, m->ystart);
      for(int y=m->yend+1;y<m->LocalNy;y++)
        var(x, y) = var(x, m->yend);
    }
  }
  file->setGlobalOrigin();

  return true;
}

/*!
 * Reads a 3D variable from a file
 * 
 * 
 */
bool GridFile::get(Mesh *m, Field3D &var,   const string &name, BoutReal def) {
  Timer timer("io");
  TRACE("GridFile::get(Field3D)");

  // Check that the file can be read
  
  if (!file->is_valid()) {
    throw BoutException("Could not read '%s' from file: File cannot be read", name.c_str());
  }

  // Check the size of the variable in the file
  
  vector<int> size = file->getSize(name);
  switch(size.size()) {
  case 0: {
    // Variable not found
    output_warn.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name.c_str(), def);
    var = def;
    return false;
  }
  case 1: {
    // 0 or 1 dimension
    if (size[0] != 1) {
      throw BoutException("Expecting a 3D variable, but '%s' is 1D with %d elements\n", name.c_str(), size[0]);
    }
    BoutReal rval;
    if (!file->read(&rval, name)) {
      throw BoutException("Couldn't read 0D variable '%s'\n", name.c_str());
    }
    var = rval;
    return true;
  }
  case 2: {
    // Read as 2D

    Field2D var2d(m);
    if (!get(m, var2d, name, def)) {
      throw BoutException("Couldn't read 2D variable '%s'\n", name.c_str());
    }
    var = var2d;
    return true;
  }
  case 3: {
    // Check whether "nz" is defined
    if (hasVar("nz")) {
      // Read directly into arrays
      
      // Check the array is the right size
      
      if (size[2] != m->LocalNz) {
        throw BoutException("3D variable '%s' has incorrect size %d (expecting %d)", name.c_str(), size[2], m->LocalNz);
      }
      
      if (! readgrid_3dvar_real(m, name,
			       m->OffsetY,// Start reading at global index
			       m->ystart,// Insert data starting from y=ystart
			       m->yend-m->ystart+1, // Length of data in Y
			       0, m->LocalNx, // All x indices (local indices)
			       var) ) {
	      throw BoutException("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name.c_str());
      }
    } else {
      // No Z size specified in file. Assume FFT format
      if (! readgrid_3dvar_fft(m, name,
			      m->OffsetY,// Start reading at global index
			      m->ystart,// Insert data starting from y=ystart
			      m->yend-m->ystart+1, // Length of data in Y
			      0, m->LocalNx, // All x indices (local indices)
			      var) ) {
	      throw BoutException("\tWARNING: Could not read '%s' from grid. Setting to zero\n", name.c_str());
      }
    }
    
    break;
  }
  default: {
    throw BoutException("Error: Variable '%s' should be 3D, but has %d dimensions\n", 
			name.c_str(), size.size());
  }
  };

  // Upper and lower Y boundaries copied from nearest point
  for(int x=0;x<m->LocalNx;x++) {
    for(int y=0;y<m->ystart;y++)
      for(int z=0;z<m->LocalNz;z++)
	var(x, y, z) = var(x, m->ystart, z);
    for(int y=m->yend+1;y<m->LocalNy;y++)
      for(int z=0;z<m->LocalNz;z++)
	var(x, y, z) = var(x, m->yend, z);
  }
  
  return true;
}

bool GridFile::get(Mesh *UNUSED(m), vector<int> &var, const string &name,
                   int len, int offset, GridDataSource::Direction UNUSED(dir)) {
  TRACE("GridFile::get(vector<int>)");
  
  if (!file->is_valid()) {
    return false;
  }

  file->setGlobalOrigin(offset);

  if (!file->read(&var[0], name, len)){
    return false;
  }
  
  file->setGlobalOrigin();
  return true;
}

bool GridFile::get(Mesh *UNUSED(m), vector<BoutReal> &var, const string &name,
                   int len, int offset, GridDataSource::Direction UNUSED(dir)) {
  TRACE("GridFile::get(vector<BoutReal>)");
  
  if (!file->is_valid()){
    return false;
  }

  file->setGlobalOrigin(offset);

  if (!file->read(&var[0], name, len)) {
    return false;
  }
  
  file->setGlobalOrigin();
  return true;
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
bool GridFile::readgrid_3dvar_fft(Mesh *m, const string &name, 
				 int yread, int ydest, int ysize, 
				 int xge, int xlt, Field3D &var) {
  /// Check the arguments make sense
  if ((yread < 0) || (ydest < 0) || (ysize < 0) || (xge < 0) || (xlt < 0)) {
    return false;
  }
  
  /// Check the size of the data
  vector<int> size = file->getSize(name);
  
  if (size.size() != 3) {
    output_warn.write("\tWARNING: Number of dimensions of %s incorrect\n", name.c_str());
    return false;
  }

  int maxmode = (size[2] - 1)/2; ///< Maximum mode-number n

  int ncz = m->LocalNz;

  BoutReal zlength = m->coordinates(var.getLocation())->zlength();
  
  int zperiod = ROUND(TWOPI / zlength); /// Number of periods in 2pi

  // Print out which modes are going to be read in
  if (zperiod > maxmode) {
    // Domain is too small: Only DC
    output_warn.write("zperiod (%d) > maxmode (%d) => Only reading n = 0 component\n", zperiod, maxmode);
  } else {
    // Get maximum mode in the input which is a multiple of zperiod
    int mm = (maxmode / zperiod) * zperiod;
    if ( (ncz/2)*zperiod < mm )
      mm = (ncz/2)*zperiod; // Limited by Z resolution
    
    if (mm == zperiod) {
      output_info.write(" => Reading n = 0, %d\n", zperiod);
    } else {
      output_info.write(" => Reading n = 0, %d ... %d\n", zperiod, mm);
    }
  }

  /// Data for FFT. Only positive frequencies
  Array<dcomplex> fdata(ncz / 2 + 1);
  Array<BoutReal> zdata(size[2]);

  for(int jx=xge;jx<xlt;jx++) {
    // Set the global X index

    for(int jy=0; jy < ysize; jy++) {
      /// Read data

      int yind = yread + jy; // Global location to read from

      file->setGlobalOrigin(jx + m->OffsetX, yind);
      if (!file->read(std::begin(zdata), name, 1, 1, size[2])) {
        return true;
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
      irfft(std::begin(fdata), ncz, &var(jx, ydest + jy, 0));
    }
  }

  file->setGlobalOrigin();
  
  return true;
}

/*!
 * Reads a 3D variable directly from the file, without 
 * any processing
 */ 
bool GridFile::readgrid_3dvar_real(Mesh *m, const string &name, 
				   int yread, int ydest, int ysize, 
				   int xge, int xlt, Field3D &var) {
  /// Check the arguments make sense
  if ((yread < 0) || (ydest < 0) || (ysize < 0) || (xge < 0) || (xlt < 0)) {
    return false;
  }
  
  /// Check the size of the data
  vector<int> size = file->getSize(name);
  
  if (size.size() != 3) {
    output_warn.write("\tWARNING: Number of dimensions of %s incorrect\n", name.c_str());
    return false;
  }
  
  for(int jx=xge;jx<xlt;jx++) {
    // Set the global X index
    
    for(int jy=0; jy < ysize; jy++) {
      /// Read data
      
      int yind = yread + jy; // Global location to read from
      
      file->setGlobalOrigin(jx + m->OffsetX, yind);
      if (!file->read(&var(jx,ydest+jy,0), name, 1, 1, size[2])) {
        return false;
      }
    }
  }
  file->setGlobalOrigin();
  
  return true;
}

