
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

Mesh* Mesh::create(Options *opt) {
  return create(NULL, opt);
}

Mesh::Mesh(GridDataSource *s, Options* opt) : source(s), coords(0), options(opt) {
  if(s == NULL)
    throw BoutException("GridDataSource passed to Mesh::Mesh() is NULL");
  
  /// Get mesh options
  OPTION(options, StaggerGrids,   false); // Stagger grids

  // Will be set to true if any variable has a free boundary condition applied to the corresponding boundary
  freeboundary_xin = false;
  freeboundary_xout = false;
  freeboundary_ydown = false;
  freeboundary_yup = false;
  
  // Initialise derivatives
  derivs_init(options);  // in index_derivs.cxx for now
}

Mesh::~Mesh() {
  delete source;

  if(coords) {
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
  msg_stack.push("Loading 2D vector: Mesh::get(Vector2D, %s)", name.c_str());

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

  msg_stack.pop();

  return 0;
}

int Mesh::get(Vector3D &var, const string &name) {
  msg_stack.push("Loading 3D vector: Mesh::get(Vector3D, %s)", name.c_str());

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

  msg_stack.pop();

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

  int mybndry = (int) !(iterateBndryLowerY().isDone());
  int allbndry;
  MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm(yend));
  answer = (bool) allbndry;
  calc = true;
  return answer;
}

Coordinates* Mesh::coordinates() {
  if(!coords) {
    // No coordinate system set. Create default
    coords = new Coordinates(this);
  }
  return coords;
}

bool Mesh::hasBndryUpperY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = (int) !(iterateBndryUpperY().isDone());
  int allbndry;
  MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm(ystart));
  answer = (bool) allbndry;
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
    bool fci_yperiodic;
    fci_options->get("y_periodic", fci_yperiodic, true);
    bool fci_zperiodic;
    fci_options->get("z_periodic", fci_zperiodic, true);
    transform = std::unique_ptr<ParallelTransform>(new FCITransform(*this, fci_yperiodic, fci_zperiodic));
      
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
