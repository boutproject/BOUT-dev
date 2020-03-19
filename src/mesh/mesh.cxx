#include <globals.hxx>
#include <bout/mesh.hxx>
#include <bout/coordinates.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include <boutcomm.hxx>
#include <output.hxx>

#include "impls/bout/boutmesh.hxx"

MeshFactory::ReturnType MeshFactory::create(Options* options, GridDataSource* source) {
  if (source != nullptr) {
    return Factory::create(getType(options), source, options);
  }

  if (options == nullptr) {
    options = Options::getRoot()->getSection(section_name);
  }

  if (options->isSet("file") or Options::root().isSet("grid")) {
    // Specified mesh file
    const auto grid_name =
        (*options)["file"].withDefault(Options::root()["grid"].withDefault(""));
    output << "\nGetting grid data from file " << grid_name << "\n";

    // Create a grid file, using specified format if given
    const auto grid_ext =
        (*options)["format"].withDefault(Options::root()["format"].withDefault(""));

    // Create a grid file
    source = static_cast<GridDataSource*>(new GridFile(
        data_format((grid_ext.empty()) ? grid_name.c_str() : grid_ext.c_str()),
        grid_name));
  } else {
    output << "\nGetting grid data from options\n";
    source = static_cast<GridDataSource*>(new GridFromOptions(options));
  }

  return Factory::create(getType(options), source, options);
}

Mesh* Mesh::create(GridDataSource *s, Options *opt) {
  return MeshFactory::getInstance().create(opt, s).release();
}

Mesh *Mesh::create(Options *opt) { return create(nullptr, opt); }

Mesh::Mesh(GridDataSource *s, Options* opt)
  : source(s), options(opt == nullptr ? Options::getRoot()->getSection("mesh") : opt),
    include_corner_cells((*options)["include_corner_cells"]
                         .doc("Communicate corner guard and boundary cells. Can be set "
                               "to false if you are sure that you will not need these "
                               "cells, for mixed derivatives D2DXDY (or anything else), "
                               "for example if your grid has orthogonal x- and "
                               "y-directions.  This might slightly reduce communication "
                               "time.")
                         .withDefault(true)) {
  if(s == nullptr)
    throw BoutException("GridDataSource passed to Mesh::Mesh() is NULL");
  
  /// Get mesh options
  OPTION(options, StaggerGrids,   false); // Stagger grids
  OPTION(options, maxregionblocksize, MAXREGIONBLOCKSIZE);
  OPTION(options, calcParallelSlices_on_communicate, true);
  // Initialise derivatives
  derivs_init(options);  // in index_derivs.cxx for now
}

Mesh::~Mesh() { delete source; }

/**************************************************************************
 * Functions for reading data from external sources
 *
 * These functions are delegated to a GridDataSource object,
 * which may then read from a file, options, or other sources.
 **************************************************************************/

namespace {
// Wrapper for writing nicely to the screen
template <class T>
void warn_default_used(const T& value, const std::string& name) {
  output_warn << "\tWARNING: Mesh has no source. Setting '" << name << "' = " << value
              << std::endl;
}
} // namespace

int Mesh::get(std::string& sval, const std::string& name, const std::string& def) {
  TRACE("Mesh::get(sval, {:s})", name);

  if (source == nullptr) {
    warn_default_used(def, name);
    sval = def;
    return true;
  }

  return !source->get(this, sval, name, def);
}

int Mesh::get(int &ival, const std::string &name, int def) {
  TRACE("Mesh::get(ival, {:s})", name);

  if (source == nullptr) {
    warn_default_used(def, name);
    ival = def;
    return true;
  }

  return !source->get(this, ival, name, def);
}

int Mesh::get(BoutReal& rval, const std::string& name, BoutReal def) {
  TRACE("Mesh::get(rval, {:s})", name);

  if (source == nullptr) {
    warn_default_used(def, name);
    rval = def;
    return true;
  }

  return !source->get(this, rval, name, def);
}

int Mesh::get(bool &bval, const std::string &name, bool def) {
  TRACE("Mesh::get(bval, {:s})", name);

  if (source == nullptr) {
    warn_default_used(def, name);
    bval = def;
    return true;
  }

  int bval_as_int = 0;
  bool success = source->get(this, bval_as_int, name, def);
  bval = bool(bval_as_int);
  return !success;
}

int Mesh::get(Field2D& var, const std::string& name, BoutReal def, bool communicate) {
  TRACE("Loading 2D field: Mesh::get(Field2D, {:s})", name);

  if (source == nullptr or !source->get(this, var, name, def)) {
    // set val to default in source==nullptr too:
    var = def;
    return 1;
  }

  // Communicate to get guard cell data
  if (communicate) {
    Mesh::communicate(var);
  }

  // Check that the data is valid
  checkData(var);

  return 0;
}

int Mesh::get(Field3D &var, const std::string &name, BoutReal def, bool communicate) {
  TRACE("Loading 3D field: Mesh::get(Field3D, {:s})", name);

  if (source == nullptr or !source->get(this, var, name, def)) {
    // set val to default in source==nullptr too:
    var = def;
    return 1;
  }

  // Communicate to get guard cell data
  if(communicate) {
    Mesh::communicate(var);
  }

  // Check that the data is valid
  checkData(var);

  return 0;
}

int Mesh::get(FieldPerp &var, const std::string &name, BoutReal def,
    bool UNUSED(communicate)) {
  TRACE("Loading FieldPerp: Mesh::get(FieldPerp, {:s})", name);

  if (source == nullptr or !source->get(this, var, name, def)) {
    // set val to default in source==nullptr too:
    var = def;
    return 1;
  }

  int yindex = var.getIndex();
  if (yindex >= 0 and yindex < var.getMesh()->LocalNy) {
    // Communicate to get guard cell data
    Mesh::communicate(var);

    // Check that the data is valid
    checkData(var);
  }

  return 0;
}

/**************************************************************************
 * Data get routines
 **************************************************************************/

int Mesh::get(Vector2D& var, const std::string& name, BoutReal def, bool communicate) {
  TRACE("Loading 2D vector: Mesh::get(Vector2D, {:s})", name);

  if(var.covariant) {
    output << _("\tReading covariant vector ") << name << endl;

    get(var.x, name + "_x", def, communicate);
    get(var.y, name + "_y", def, communicate);
    get(var.z, name + "_z", def, communicate);

  }else {
    output << _("\tReading contravariant vector ") << name << endl;

    get(var.x, name + "x", def, communicate);
    get(var.y, name + "y", def, communicate);
    get(var.z, name + "z", def, communicate);
  }

  return 0;
}

int Mesh::get(Vector3D& var, const std::string& name, BoutReal def, bool communicate) {
  TRACE("Loading 3D vector: Mesh::get(Vector3D, {:s})", name);

  if(var.covariant) {
    output << _("\tReading covariant vector ") << name << endl;

    get(var.x, name + "_x", def, communicate);
    get(var.y, name + "_y", def, communicate);
    get(var.z, name + "_z", def, communicate);

  }else {
    output << ("\tReading contravariant vector ") << name << endl;

    get(var.x, name + "x", def, communicate);
    get(var.y, name + "y", def, communicate);
    get(var.z, name + "z", def, communicate);
  }

  return 0;
}

bool Mesh::isDataSourceGridFile() const {
  return source != nullptr and source->is_file;
}

bool Mesh::sourceHasVar(const std::string &name) {
  TRACE("Mesh::sourceHasVar({:s})", name);
  if (source == nullptr)
    return false;
  return source->hasVar(name);
}

/// Wrapper for GridDataSource::hasXBoundaryGuards
bool Mesh::sourceHasXBoundaryGuards() {
  return source->hasXBoundaryGuards(this);
}

/// Wrapper for GridDataSource::hasYBoundaryGuards
bool Mesh::sourceHasYBoundaryGuards() {
  return source->hasYBoundaryGuards();
}

/**************************************************************************
 * Communications
 **************************************************************************/

void Mesh::communicateXZ(FieldGroup &g) {
  TRACE("Mesh::communicate(FieldGroup&)");

  // Send data
  comm_handle h = sendX(g);

  // Wait for data from other processors
  wait(h);
}

void Mesh::communicateYZ(FieldGroup &g) {
  TRACE("Mesh::communicate(FieldGroup&)");

  // Send data
  comm_handle h = sendY(g);

  // Wait for data from other processors
  wait(h);

  // Calculate yup and ydown fields for 3D fields
  if (calcParallelSlices_on_communicate) {
    for(const auto& fptr : g.field3d()) {
      fptr->calcParallelSlices();
    }
  }
}

void Mesh::communicate(FieldGroup &g) {
  TRACE("Mesh::communicate(FieldGroup&)");

  if (include_corner_cells) {
    // Send data in y-direction
    comm_handle h = sendY(g);

    // Wait for data from other processors
    wait(h);

    // Send data in x-direction
    h = sendX(g);

    // Wait for data from other processors
    wait(h);
  } else {
    // Send data
    comm_handle h = send(g);

    // Wait for data from other processors
    wait(h);
  }

  // Calculate yup and ydown fields for 3D fields
  if (calcParallelSlices_on_communicate) {
    for(const auto& fptr : g.field3d()) {
      fptr->calcParallelSlices();
    }
  }
}

/// This is a bit of a hack for now to get FieldPerp communications
/// The FieldData class needs to be changed to accomodate FieldPerp objects
void Mesh::communicate(FieldPerp &f) {
  comm_handle recv[2];
  
  int nin = xstart; // Number of x points in inner guard cell
  int nout = LocalNx-xend-1; // Number of x points in outer guard cell

  // Post receives for guard cell regions

  recv[0] = irecvXIn(f[0],       nin*LocalNz, 0);
  recv[1] = irecvXOut(f[xend+1], nout*LocalNz, 1);
  
  // Send data
  sendXIn(f[xstart], nin*LocalNz, 1);
  sendXOut(f[xend-nout+1], nout*LocalNz, 0);
 
  // Wait for receive
  wait(recv[0]);
  wait(recv[1]);
}

int Mesh::msg_len(const std::vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt) {
  int len = 0;

  /// Loop over variables
  for(const auto& var : var_list) {
    if(var->is3D()) {
      len += (xlt - xge) * (ylt - yge) * LocalNz * var->BoutRealSize();
    } else {
      len += (xlt - xge) * (ylt - yge) * var->BoutRealSize();
    }
  }

  return len;
}

bool Mesh::periodicY(int jx) const {
  BoutReal ts; return periodicY(jx, ts);
}

int Mesh::ySize(int jx) const {
  // Get the size of a surface in Y using MPI communicator
  MPI_Comm comm = getYcomm(jx);

  int local = yend - ystart + 1;
  int all;
  mpi->MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_SUM, comm);
  return all;
}

bool Mesh::hasBndryLowerY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = static_cast<int>(!(iterateBndryLowerY().isDone()));
  int allbndry;
  mpi->MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm(yend));
  answer = static_cast<bool>(allbndry);
  calc = true;
  return answer;
}

bool Mesh::hasBndryUpperY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = static_cast<int>(!(iterateBndryUpperY().isDone()));
  int allbndry;
  mpi->MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm(ystart));
  answer = static_cast<bool>(allbndry);
  calc = true;
  return answer;
}

const std::vector<int> Mesh::readInts(const std::string &name, int n) {
  TRACE("Mesh::readInts({:s})", name);

  if (source == nullptr) {
    throw BoutException("Can't read integer array {:s} as 'Mesh::source' is nullptr\n",
                        name);
  }

  std::vector<int> result;

  if(source->hasVar(name)) {
    if(!source->get(this, result, name, n, 0)) {
      // Error reading
      throw BoutException(_("Could not read integer array '{:s}'\n"), name.c_str());
    }
  }else {
    // Not found
    throw BoutException(_("Missing integer array {:s}\n"), name.c_str());
  }

  return result;
}

std::shared_ptr<Coordinates> Mesh::createDefaultCoordinates(const CELL_LOC location,
    bool force_interpolate_from_centre) {

  if (location == CELL_CENTRE || location == CELL_DEFAULT) {
    // Initialize coordinates from input
    return std::make_shared<Coordinates>(this, options);
  } else {
    // Interpolate coordinates from CELL_CENTRE version
    return std::make_shared<Coordinates>(this, options, location,
        getCoordinates(CELL_CENTRE), force_interpolate_from_centre);
  }
}

const Region<>& Mesh::getRegion3D(const std::string& region_name) const {
  const auto found = regionMap3D.find(region_name);
  if (found == end(regionMap3D)) {
    throw BoutException(_("Couldn't find region {:s} in regionMap3D"), region_name);
  }
  return found->second;
}

const Region<Ind2D>& Mesh::getRegion2D(const std::string& region_name) const {
  const auto found = regionMap2D.find(region_name);
  if (found == end(regionMap2D)) {
    throw BoutException(_("Couldn't find region {:s} in regionMap2D"), region_name);
  }
  return found->second;
}

const Region<IndPerp>& Mesh::getRegionPerp(const std::string& region_name) const {
  const auto found = regionMapPerp.find(region_name);
  if (found == end(regionMapPerp)) {
    throw BoutException(_("Couldn't find region {:s} in regionMapPerp"), region_name);
  }
  return found->second;
}

bool Mesh::hasRegion3D(const std::string& region_name) const {
  return regionMap3D.find(region_name) != std::end(regionMap3D);
}

bool Mesh::hasRegion2D(const std::string& region_name) const {
  return regionMap2D.find(region_name) != std::end(regionMap2D);
}

bool Mesh::hasRegionPerp(const std::string& region_name) const {
  return regionMapPerp.find(region_name) != std::end(regionMapPerp);
}

void Mesh::addRegion3D(const std::string &region_name, const Region<> &region) {
  if (regionMap3D.count(region_name)) {
    throw BoutException(_("Trying to add an already existing region {:s} to regionMap3D"),
                        region_name);
  }
  regionMap3D[region_name] = region;
  output_verbose.write(_("Registered region 3D {:s}"),region_name);
  output_verbose << "\n:\t" << region.getStats() << "\n";
}

void Mesh::addRegion2D(const std::string &region_name, const Region<Ind2D> &region) {
  if (regionMap2D.count(region_name)) {
    throw BoutException(_("Trying to add an already existing region {:s} to regionMap2D"),
                        region_name);
  }
  regionMap2D[region_name] = region;
  output_verbose.write(_("Registered region 2D {:s}"),region_name);
  output_verbose << "\n:\t" << region.getStats() << "\n";
}

void Mesh::addRegionPerp(const std::string &region_name, const Region<IndPerp> &region) {
  if (regionMapPerp.count(region_name)) {
    throw BoutException(
        _("Trying to add an already existing region {:s} to regionMapPerp"), region_name);
  }
  regionMapPerp[region_name] = region;
  output_verbose.write(_("Registered region Perp {:s}"),region_name);
  output_verbose << "\n:\t" << region.getStats() << "\n";
}

void Mesh::createDefaultRegions(){
  //3D regions
  addRegion3D("RGN_ALL", Region<Ind3D>(0, LocalNx - 1, 0, LocalNy - 1, 0, LocalNz - 1,
                                       LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOBNDRY", Region<Ind3D>(xstart, xend, ystart, yend, zstart, zend,
                                           LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOX", Region<Ind3D>(xstart, xend, 0, LocalNy - 1, 0, LocalNz - 1,
                                       LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOY", Region<Ind3D>(0, LocalNx - 1, ystart, yend, 0, LocalNz - 1,
                                       LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOZ", Region<Ind3D>(0, LocalNx - 1, 0, LocalNy - 1, zstart, zend,
                                       LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_GUARDS", mask(getRegion3D("RGN_ALL"), getRegion3D("RGN_NOBNDRY")));
  addRegion3D("RGN_XGUARDS", Region<Ind3D>(0, xstart - 1, ystart, yend, zstart, zend,
          LocalNy, LocalNz, maxregionblocksize)
      + Region<Ind3D>(xend + 1, LocalNx - 1, ystart, yend, zstart, zend,
          LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_YGUARDS", Region<Ind3D>(xstart, xend, 0, ystart - 1, zstart, zend,
          LocalNy, LocalNz, maxregionblocksize)
      + Region<Ind3D>(xstart, xend, yend + 1, LocalNy - 1, zstart, zend,
          LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_ZGUARDS", Region<Ind3D>(xstart, xend, ystart, yend, 0, zstart - 1,
          LocalNy, LocalNz, maxregionblocksize)
      + Region<Ind3D>(xstart, xend, ystart, yend, zend + 1, LocalNz - 1,
          LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOCORNERS",
      (getRegion3D("RGN_NOBNDRY") + getRegion3D("RGN_XGUARDS") +
        getRegion3D("RGN_YGUARDS") + getRegion3D("RGN_ZGUARDS")).unique());

  addRegion3D("RGN_UPPER_Y_THIN",	
              Region<Ind3D>(xstart, xend, yend + 1, yend + 1, 0, LocalNz - 1, LocalNy, LocalNz,	
                            maxregionblocksize));
  addRegion3D("RGN_LOWER_Y_THIN",	
              Region<Ind3D>(xstart, xend, ystart + 1, ystart + 1, 0, LocalNz - 1, LocalNy, LocalNz,	
                            maxregionblocksize));
  
  
  //2D regions
  addRegion2D("RGN_ALL", Region<Ind2D>(0, LocalNx - 1, 0, LocalNy - 1, 0, 0, LocalNy, 1,
                                       maxregionblocksize));
  addRegion2D("RGN_NOBNDRY", Region<Ind2D>(xstart, xend, ystart, yend, 0, 0, LocalNy, 1,
                                           maxregionblocksize));
  addRegion2D("RGN_NOX", Region<Ind2D>(xstart, xend, 0, LocalNy - 1, 0, 0, LocalNy, 1,
                                       maxregionblocksize));
  addRegion2D("RGN_NOY", Region<Ind2D>(0, LocalNx - 1, ystart, yend, 0, 0, LocalNy, 1,
                                       maxregionblocksize));
  addRegion2D("RGN_NOZ", Region<Ind2D>(0, LocalNx - 1, 0, LocalNy - 1, 0, 0, LocalNy, 1,
                                       maxregionblocksize));
  addRegion2D("RGN_GUARDS", mask(getRegion2D("RGN_ALL"), getRegion2D("RGN_NOBNDRY")));
  addRegion2D("RGN_XGUARDS", Region<Ind2D>(0, xstart - 1, ystart, yend, 0, 0, LocalNy, 1,
          maxregionblocksize)
      + Region<Ind2D>(xend + 1, LocalNx - 1, ystart, yend, 0, 0, LocalNy, 1,
          maxregionblocksize));
  addRegion2D("RGN_YGUARDS", Region<Ind2D>(xstart, xend, 0, ystart - 1, 0, 0, LocalNy, 1,
          maxregionblocksize)
      + Region<Ind2D>(xstart, xend, yend + 1, LocalNy - 1, 0, 0, LocalNy, 1,
          maxregionblocksize));
  addRegion2D("RGN_ZGUARDS", Region<Ind2D>(xstart, xend, ystart, yend, 0, -1, LocalNy, 1,
          maxregionblocksize)
      + Region<Ind2D>(xstart, xend, ystart, yend, 0, -1, LocalNy, 1,
          maxregionblocksize));
  addRegion2D("RGN_NOCORNERS",
      (getRegion2D("RGN_NOBNDRY") + getRegion2D("RGN_XGUARDS") +
        getRegion2D("RGN_YGUARDS") + getRegion2D("RGN_ZGUARDS")).unique());

  addRegion2D("RGN_UPPER_Y_THIN",	
              Region<Ind2D>(xstart, xend, yend + 1, yend + 1, 0, -1, LocalNy, 1,
			    maxregionblocksize));
  addRegion2D("RGN_LOWER_Y_THIN",	
              Region<Ind2D>(xstart, xend, ystart + 1, ystart + 1, 0, -1, LocalNy, 1,	
                            maxregionblocksize));

  // Perp regions
  addRegionPerp("RGN_ALL", Region<IndPerp>(0, LocalNx - 1, 0, 0, 0, LocalNz - 1, 1,
                                           LocalNz, maxregionblocksize));
  addRegionPerp("RGN_NOBNDRY", Region<IndPerp>(xstart, xend, 0, 0, zstart, zend, 1,
                                               LocalNz, maxregionblocksize));
  addRegionPerp("RGN_NOX", Region<IndPerp>(xstart, xend, 0, 0, 0, LocalNz - 1, 1, LocalNz,
                                           maxregionblocksize)); // Same as NOBNDRY
  addRegionPerp("RGN_NOY", Region<IndPerp>(0, LocalNx - 1, 0, 0, 0, LocalNz - 1, 1,
                                           LocalNz, maxregionblocksize));

  addRegionPerp("RGN_NOZ", Region<IndPerp>(0, LocalNx - 1, 0, 0, zstart, zend, 1, LocalNz,
                                           maxregionblocksize));
  addRegionPerp("RGN_GUARDS", mask(getRegionPerp("RGN_ALL"), getRegionPerp("RGN_NOBNDRY")));
  addRegionPerp("RGN_XGUARDS", Region<IndPerp>(0, xstart - 1, 0, 0, zstart, zend, 1,
          LocalNz, maxregionblocksize)
      + Region<IndPerp>(xend + 1, LocalNx - 1, 0, 0, zstart, zend, 1,
          LocalNz, maxregionblocksize));
  addRegionPerp("RGN_YGUARDS", Region<IndPerp>(xstart, xend, 0, -1, zstart, zend, 1,
          LocalNz, maxregionblocksize)
      + Region<IndPerp>(xstart, xend, 0, -1, zstart, zend, 1,
          LocalNz, maxregionblocksize));
  addRegionPerp("RGN_ZGUARDS", Region<IndPerp>(xstart, xend, 0, 0, 0, zstart - 1, 1,
          LocalNz, maxregionblocksize)
      + Region<IndPerp>(xstart, xend, 0, 0, zend + 1, LocalNz - 1, 1,
          LocalNz, maxregionblocksize));
  addRegionPerp("RGN_NOCORNERS",
      (getRegionPerp("RGN_NOBNDRY") + getRegionPerp("RGN_XGUARDS") +
        getRegionPerp("RGN_YGUARDS") + getRegionPerp("RGN_ZGUARDS")).unique());

  // Construct index lookup for 3D-->2D
  indexLookup3Dto2D = Array<int>(LocalNx*LocalNy*LocalNz);
  BOUT_FOR(ind3D, getRegion3D("RGN_ALL")) {
    indexLookup3Dto2D[ind3D.ind] = ind3Dto2D(ind3D).ind;
  }
}

void Mesh::recalculateStaggeredCoordinates() {
  for (auto &i : coords_map) {
    CELL_LOC location = i.first;

    if (location == CELL_CENTRE) {
      // Only reset staggered locations
      continue;
    }

    *coords_map[location] = std::move(*createDefaultCoordinates(location, true));
  }
}

constexpr decltype(MeshFactory::type_name) MeshFactory::type_name;
constexpr decltype(MeshFactory::section_name) MeshFactory::section_name;
constexpr decltype(MeshFactory::option_name) MeshFactory::option_name;
constexpr decltype(MeshFactory::default_type) MeshFactory::default_type;
