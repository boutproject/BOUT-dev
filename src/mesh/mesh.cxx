
#include <globals.hxx>
#include <bout/mesh.hxx>
#include <bout/coordinates.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include "meshfactory.hxx"

#include <output.hxx>

Mesh* Mesh::create(GridDataSource *s, Options *opt) {
  return MeshFactory::getInstance()->createMesh(s, opt);
}

Mesh* Mesh::create(Options *opt) {
  return create(NULL, opt);
}

Mesh::Mesh(GridDataSource *s, Options* options) : source(s), coords(0) {
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
  
  /*
  // Gaussj working arrays
  if(ilen > 0) {
    ivfree(indxc);
    ivfree(indxr);
    ivfree(ipiv);
  }
  */
}

/**************************************************************************
 * Default functions for getting scalars
 **************************************************************************/

/// Get an integer
int Mesh::get(int &ival, const string &name) {
  int msg_pos = msg_stack.push("Loading integer: Mesh::get(int, %s)", name.c_str());

  if(!source->hasVar(name)) {
    msg_stack.pop(msg_pos);
    output << "Source doesn't have var: " << name << endl;
    return 1;
  }
  
  source->open(name);
  bool success = source->fetch(&ival, name);
  source->close();
  
  msg_stack.pop(msg_pos);

  if(!success) {
    output << "Failed to read var: " << name << endl;
    return 2;
  }
  return 0;
}

/// A BoutReal number
int Mesh::get(BoutReal &rval, const string &name) {
  if(!source->hasVar(name))
    return 1;
  
  source->open(name);
  bool success = source->fetch(&rval, name);
  source->close();
  
  if(!success)
    return 2;
  return 0;
}

int Mesh::get(Field2D &var, const string &name, BoutReal def) {
  int msg_pos = msg_stack.push("Loading 2D field: BoutMesh::get(Field2D, %s)", name.c_str());
  
  output << "GET: " << name << " = " << def << endl;

  if(!source->hasVar(name)) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name.c_str(), def);
    var = def;
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 2;
  }
  
  var.allocate(); // Make sure data allocated
  
  BoutReal **data = var.getData(); // pointer for faster access
  
  // Send an open signal to the source
  source->open(name);
  
  // Get the size of the variable
  vector<int> size = source->getSize(name);
  switch(size.size()) {
  case 1: {
    // 0 or 1 dimension
    if(size[0] != 1) {
      output.write("Expecting a 2D variable, but '%s' is 1D with %d elements\n", name.c_str(), size[0]);
      source->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    BoutReal val;
    if(!source->fetch(&val, name)) {
      output.write("Couldn't read 0D variable '%s'\n", name.c_str());
      source->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    var = val;
    // Close source
    source->close();
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 0;
  }
  case 2: {
    // Check size? More complicated now...
    break;
  }
  default: {
    output.write("Error: Variable '%s' should be 2D, but has %d dimensions\n", 
                 name.c_str(), size.size());
    source->close();
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 1;
  }
  }
  
  // Read bulk of points
  read2Dvar(source, name, 
            OffsetX, OffsetY,  // Coordinates in grid file
            xstart, ystart,    // Coordinates in this processor
            xend-xstart+1, yend-ystart+1,      // Number of points to read
            data);
  
  // Close the data source
  source->close();
  
  // Communicate to get guard cell data
  Mesh::communicate(var);
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return 0;
}

int Mesh::get(Field3D &var, const string &name) {
  return 1;  // Not yet implemented
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

/**************************************************************************
 * Communications
 **************************************************************************/

int Mesh::communicate(FieldData &f) {
  FieldGroup group;
  group.add(f);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2) {
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3) {
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4) {
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  group.add(f4);
  return communicate(group);
}

comm_handle Mesh::send(FieldData &f) {
  FieldGroup group;
  group.add(f);
  return send(group);
}

/// This is a bit of a hack for now to get FieldPerp communications
/// The FieldData class needs to be changed to accomodate FieldPerp objects
int Mesh::communicate(FieldPerp &f) {
  comm_handle recv[2];
  
  BoutReal **fd = f.getData();
  
  int nin = xstart; // Number of x points in inner guard cell
  int nout = ngx-xend-1; // Number of x points in outer guard cell
  
  // Post receives for guard cell regions
  recv[0] = irecvXIn(fd[0],       nin*ngz, 0);
  recv[1] = irecvXOut(fd[xend+1], nout*ngz, 1);
  
  // Send data
  sendXIn(fd[xstart], nin*ngz, 1);
  sendXOut(fd[xend-nout+1], nout*ngz, 0);
 
  // Wait for receive
  wait(recv[0]);
  wait(recv[1]);

  return 0;
}

int Mesh::msg_len(const vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt) {
  int len = 0;

  /// Loop over variables
  for(std::vector<FieldData*>::const_iterator it = var_list.begin(); it != var_list.end(); it++) {
    if((*it)->is3D()) {
      len += (xlt - xge) * (ylt - yge) * (ngz-1) * (*it)->BoutRealSize();
    }else
      len += (xlt - xge) * (ylt - yge) * (*it)->BoutRealSize();
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
    source->open(name);
    source->setGlobalOrigin();
    result.resize(n);
    if(!source->fetch(&(result.front()), name, n)) {
      // Error reading
      source->close();
      throw BoutException("Could not read integer array '%s'\n", name.c_str());
    }
    source->close();
  }else {
    // Not found
    throw BoutException("Missing integer array %s\n", name.c_str());
  }
  
  return result;
}

void Mesh::read2Dvar(GridDataSource *s, const string &name, 
                     int xs, int ys,
                     int xd, int yd,
                     int nx, int ny,
                     BoutReal **data) {
  
  for(int x=xs;x < xs+nx; x++) {
    // Set source origin
    s->setGlobalOrigin(x, ys);
    // Read in block of data for this x value
    if(!s->fetch(&data[x-xs+xd][yd], name, 1, ny))
      throw BoutException("Could not fetch data for '%s'", name.c_str());
  }
  s->setGlobalOrigin();
}

