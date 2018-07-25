
#include <globals.hxx>
#include <bout/mesh.hxx>
#include <bout/coordinates.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include "meshfactory.hxx"

#include <output.hxx>

#include "parallel/fci.hxx"

Mesh* Mesh::create(GridDataSource *s, Options *opt) {
  return MeshFactory::getInstance()->createMesh(s, opt);
}

Mesh *Mesh::create(Options *opt) { return create(nullptr, opt); }

Mesh::Mesh(GridDataSource *s, Options* opt) : source(s), coords(nullptr), options(opt) {
  if(s == nullptr)
    throw BoutException("GridDataSource passed to Mesh::Mesh() is NULL");
  
  if (options == nullptr) {
    options = Options::getRoot()->getSection("mesh");
  }

  /// Get mesh options
  OPTION(options, StaggerGrids,   false); // Stagger grids
  OPTION(options, maxregionblocksize, MAXREGIONBLOCKSIZE);
  // Initialise derivatives
  derivs_init(options);  // in index_derivs.cxx for now
}

Mesh::~Mesh() {
  if (source) {
    delete source;
  }

  if (coords) {
    delete coords;
  }
}

/**************************************************************************
 * Functions for reading data from external sources
 *
 * These functions are delegated to a GridDataSource object,
 * which may then read from a file, options, or other sources.
 **************************************************************************/

/// Get an integer
int Mesh::get(int &ival, const string &name) {
  TRACE("Mesh::get(ival)");

  if(!source->get(this, ival, name))
    return 1;

  return 0;
}

/// A BoutReal number
int Mesh::get(BoutReal &rval, const string &name) {
  TRACE("Mesh::get(rval)");

  if(!source->get(this, rval, name))
    return 1;

  return 0;
}

int Mesh::get(Field2D &var, const string &name, BoutReal def) {
  TRACE("Loading 2D field: Mesh::get(Field2D)");

  // Ensure data allocated
  var.allocate();

  if(!source->get(this, var, name, def))
    return 1;

  // Communicate to get guard cell data
  Mesh::communicate(var);

  // Check that the data is valid
  checkData(var);

  return 0;
}

int Mesh::get(Field3D &var, const string &name, BoutReal def, bool communicate) {
  TRACE("Loading 3D field: Mesh::get(Field3D)");

  // Ensure data allocated
  var.allocate();

  if(!source->get(this, var, name, def))
    return 1;

  // Communicate to get guard cell data
  if(communicate) {
    Mesh::communicate(var);
  }

  // Check that the data is valid
  checkData(var);

  return 0;
}

/**************************************************************************
 * Data get routines
 **************************************************************************/

int Mesh::get(Vector2D &var, const string &name) {
  TRACE("Loading 2D vector: Mesh::get(Vector2D, %s)", name.c_str());

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;

    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");

  }else {
    output << "\tReading contravariant vector " << name << endl;

    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }

  return 0;
}

int Mesh::get(Vector3D &var, const string &name) {
  TRACE("Loading 3D vector: Mesh::get(Vector3D, %s)", name.c_str());

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;

    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");

  }else {
    output << "\tReading contravariant vector " << name << endl;

    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }

  return 0;
}

bool Mesh::sourceHasVar(const string &name) {
  return source->hasVar(name);
}

/**************************************************************************
 * Communications
 **************************************************************************/

void Mesh::communicateXZ(FieldGroup &g) {
  TRACE("Mesh::communicate(FieldGroup&)");

  // Send data
  comm_handle h = send(g);

  // Wait for data from other processors
  wait(h);
}

void Mesh::communicate(FieldGroup &g) {
  TRACE("Mesh::communicate(FieldGroup&)");

  // Send data
  comm_handle h = send(g);

  // Wait for data from other processors
  wait(h);

  // Calculate yup and ydown fields for 3D fields
  for(const auto& fptr : g.field3d())
    getParallelTransform().calcYUpDown(*fptr);
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

int Mesh::msg_len(const vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt) {
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
  MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_SUM, comm);
  return all;
}

bool Mesh::hasBndryLowerY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = static_cast<int>(!(iterateBndryLowerY().isDone()));
  int allbndry;
  MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm(yend));
  answer = static_cast<bool>(allbndry);
  calc = true;
  return answer;
}

bool Mesh::hasBndryUpperY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = static_cast<int>(!(iterateBndryUpperY().isDone()));
  int allbndry;
  MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm(ystart));
  answer = static_cast<bool>(allbndry);
  calc = true;
  return answer;
}

const vector<int> Mesh::readInts(const string &name, int n) {
  vector<int> result;

  if(source->hasVar(name)) {
    if(!source->get(this, result, name, n, 0)) {
      // Error reading
      throw BoutException("Could not read integer array '%s'\n", name.c_str());
    }
  }else {
    // Not found
    throw BoutException("Missing integer array %s\n", name.c_str());
  }

  return result;
}

void Mesh::setParallelTransform() {

  string ptstr;
  options->get("paralleltransform", ptstr, "identity");

  // Convert to lower case for comparison
  ptstr = lowercase(ptstr);
    
  if(ptstr == "identity") {
    // Identity method i.e. no transform needed
    transform = std::unique_ptr<ParallelTransform>(new ParallelTransformIdentity());
      
  }else if(ptstr == "shifted") {
    // Shifted metric method
    transform = std::unique_ptr<ParallelTransform>(new ShiftedMetric(*this));
      
  }else if(ptstr == "fci") {

    Options *fci_options = Options::getRoot()->getSection("fci");
    // Flux Coordinate Independent method
    bool fci_zperiodic;
    fci_options->get("z_periodic", fci_zperiodic, true);
    transform = std::unique_ptr<ParallelTransform>(new FCITransform(*this, fci_zperiodic));
      
  }else {
    throw BoutException("Unrecognised paralleltransform option.\n"
                        "Valid choices are 'identity', 'shifted', 'fci'");
  }
}

ParallelTransform& Mesh::getParallelTransform() {
  if(!transform) {
    // No ParallelTransform object yet. Set from options
    setParallelTransform();
  }
  
  // Return a reference to the ParallelTransform object
  return *transform;
}

Coordinates *Mesh::createDefaultCoordinates() {
  return new Coordinates(this);
}


Region<> & Mesh::getRegion3D(const std::string &region_name){
   auto found = regionMap3D.find(region_name);
   if (found == end(regionMap3D)) {
     throw BoutException("Couldn't find region %s in regionMap3D", region_name.c_str());
   }
   return found->second;
}

Region<Ind2D> & Mesh::getRegion2D(const std::string &region_name){
   auto found = regionMap2D.find(region_name);
   if (found == end(regionMap2D)) {
     throw BoutException("Couldn't find region %s in regionMap2D", region_name.c_str());
   }
   return found->second;
}

Region<IndPerp> &Mesh::getRegionPerp(const std::string &region_name) {
  auto found = regionMapPerp.find(region_name);
  if (found == end(regionMapPerp)) {
    throw BoutException("Couldn't find region %s in regionMapPerp", region_name.c_str());
  }
  return found->second;
}

void Mesh::addRegion3D(const std::string &region_name, const Region<> &region) {
  if (regionMap3D.count(region_name)) {
    throw BoutException("Trying to add an already existing region %s to regionMap3D");
  }
  regionMap3D[region_name] = region;
}

void Mesh::addRegion2D(const std::string &region_name, const Region<Ind2D> &region) {
  if (regionMap2D.count(region_name)) {
    throw BoutException("Trying to add an already existing region %s to regionMap2D");
  }
  regionMap2D[region_name] = region;
}

void Mesh::addRegionPerp(const std::string &region_name, const Region<IndPerp> &region) {
  if (regionMapPerp.count(region_name)) {
    throw BoutException("Trying to add an already existing region %s to regionMapPerp");
  }
  regionMapPerp[region_name] = region;
}

void Mesh::createDefaultRegions() {
  //3D regions
  addRegion3D("RGN_ALL", Region<Ind3D>(0, LocalNx - 1, 0, LocalNy - 1, 0, LocalNz - 1,
                                       LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOBNDRY", Region<Ind3D>(xstart, xend, ystart, yend, 0, LocalNz - 1,
                                           LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOX", Region<Ind3D>(xstart, xend, 0, LocalNy - 1, 0, LocalNz - 1,
                                       LocalNy, LocalNz, maxregionblocksize));
  addRegion3D("RGN_NOY", Region<Ind3D>(0, LocalNx - 1, ystart, yend, 0, LocalNz - 1,
                                       LocalNy, LocalNz, maxregionblocksize));

  //2D regions
  addRegion2D("RGN_ALL", Region<Ind2D>(0, LocalNx - 1, 0, LocalNy - 1, 0, 0, LocalNy, 1,
                                       maxregionblocksize));
  addRegion2D("RGN_NOBNDRY", Region<Ind2D>(xstart, xend, ystart, yend, 0, 0, LocalNy, 1,
                                           maxregionblocksize));
  addRegion2D("RGN_NOX", Region<Ind2D>(xstart, xend, 0, LocalNy - 1, 0, 0, LocalNy, 1,
                                       maxregionblocksize));
  addRegion2D("RGN_NOY", Region<Ind2D>(0, LocalNx - 1, ystart, yend, 0, 0, LocalNy, 1,
                                       maxregionblocksize));

  //Perp regions
  addRegionPerp("RGN_ALL", Region<IndPerp>(0, LocalNx - 1, 0, 0, 0, LocalNz - 1, 1, LocalNz,
                                       maxregionblocksize));
  addRegionPerp("RGN_NOBNDRY", Region<IndPerp>(xstart, xend, 0, 0, 0, LocalNz - 1, 1, LocalNz,
                                           maxregionblocksize));
  addRegionPerp("RGN_NOX", Region<IndPerp>(xstart, xend, 0, 0, 0, LocalNz - 1, 1, LocalNz,
                                       maxregionblocksize));
  addRegionPerp("RGN_NOY", Region<IndPerp>(0, LocalNx - 1, 0, 0, 0, LocalNz - 1, 1, LocalNz,
                                       maxregionblocksize));

  // Construct index lookup for 3D-->2D
  indexLookup3Dto2D = Array<int>(LocalNx*LocalNy*LocalNz);
  for (const auto &ind3D: getRegion3D("RGN_ALL")){
    indexLookup3Dto2D[ind3D.ind] = ind3Dto2D(ind3D).ind;
  }

  indexLookupPerpto3D = Matrix<int>(LocalNx*LocalNz,LocalNy);
  for (const auto &ind3D: getRegion3D("RGN_ALL")){
    int jperp = ind3DtoPerp(ind3D).ind;
    int jy = (ind3D.ind / LocalNz) % LocalNy;
    indexLookupPerpto3D(jperp,jy) = ind3D.ind;
  }
}
