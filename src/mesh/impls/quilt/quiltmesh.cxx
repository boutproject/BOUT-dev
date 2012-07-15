
#include <globals.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <boutcomm.hxx>

#include <bout/sys/timer.hxx>

#include <algorithm>
#include <numeric>
using std::accumulate;
#include <cmath>

#include "quiltmesh.hxx"
#include "../partition.hxx"

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

QuiltMesh::~QuiltMesh() {
  
}

int QuiltMesh::load() {
  return load(BoutComm::get());
}

int QuiltMesh::load(MPI_Comm comm) {
  // Get number of processors
  int NPES, MYPE;
  MPI_Comm_size(comm, &NPES);
  MPI_Comm_rank(comm, &MYPE);
  
  // Get number of regions
  int nregions;
  if(Mesh::get(nregions, "nregions"))
    throw BoutException("Mesh file doesn't have nregions variable. Incorrect format?\n");
  
  if(nregions > NPES)
    throw BoutException("Can't divide %d regions between %d processors\n", nregions, NPES);
  
  // Get array of grid sizes
  vector<int> nx = readInts("nx", nregions);
  vector<int> ny = readInts("ny", nregions);
  
  // Create domains
  vector<QuiltDomain*> domains;
  
  for(int i=0;i<nregions;i++)
    domains.push_back(new QuiltDomain(nx[i], ny[i]));
  
  // Join domains together
  
  // Partition the mesh (all domains connected together)
  partitionAll(domains[0], NPES);
  
  // Assign domains to processors
  int p = 0;
  for(Domain::iterator it=domains[0]->begin(); it != domains[0]->end(); it++) {
    // For now very simple numbering. Should cluster nearby domains
    QuiltDomain *qd = static_cast<QuiltDomain*>( &(*it) );
    qd->proc = p; p++;
    if(p == MYPE)
      mydomain = qd; // Domain for this processor
  }
  
  // Iterate through boundaries to get guard cell regions
  for(QuiltDomain::bndry_iterator it = mydomain->bndry_begin(); it != mydomain->bndry_end(); it++) {
    Domain::Bndry* b = *it;
     
    QuiltDomain* to = (QuiltDomain*) b->getNeighbour(mydomain);
    
    // Check which boundary(s) this is on
    if(b->onSide(mydomain, Domain::xlow)) {
      // get range of indices
      int min = b->getMin(Domain::xlow);
      int max = b->getMax(Domain::xlow);
      
      // Create new GuardRange 
      GuardRange *g = new GuardRange();
      g->xmin = 0;
      g->xmax = MXG;
      g->ymin = min;
      g->ymax = max;
      g->proc = to->proc;
    }
  }
}

/****************************************************************
 * Getting variables
 ****************************************************************/
  
int QuiltMesh::get(Field2D &var, const char *name, BoutReal def) {
  if(name == NULL)
    return 1;
  
#ifdef CHECK
  int msg_pos = msg_stack.push("Loading 2D field: BoutMesh::get(Field2D, %s)", name);
#endif
  
  GridDataSource *s = findSource(name);
  if(s == NULL) {
    output.write("\tWARNING: Could not read '%s' from grid. Setting to %le\n", name, def);
    var = def;
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 2;
  }
  
  var.allocate(); // Make sure data allocated
  
  BoutReal **data = var.getData(); // pointer for faster access
  
  // Send an open signal to the source
  s->open(name);
  
  // Get the size of the variable
  vector<int> size = s->getSize(name);
  switch(size.size()) {
  case 1: {
    // 0 or 1 dimension
    if(size[0] != 1) {
      output.write("Expecting a 2D variable, but '%s' is 1D with %d elements\n", name, size[0]);
      s->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    BoutReal val;
    if(!s->fetch(&val, name)) {
      output.write("Couldn't read 0D variable '%s'\n", name);
      s->close();
#ifdef CHECK
      msg_stack.pop(msg_pos);
#endif
      return 1;
    }
    var = val;
    // Close source
    s->close();
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
                 name, size.size());
    s->close();
#ifdef CHECK
    msg_stack.pop(msg_pos);
#endif
    return 1;
  }
  }
  
  // Read bulk of points
  read2Dvar(s, name, 
            mydomain->xOrigin(), mydomain->yOrigin(),  // Coordinates in grid file
            MXG, MYG,                    // Coordinates in this processor
            mydomain->xSize(), mydomain->ySize(),  // Number of points to read
            data);
  
  // Close the data source
  s->close();
  
  // Communicate to get guard cell data
  Mesh::communicate(var);
  
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return 0;
}

int QuiltMesh::get(Field2D &var, const string &name, BoutReal def) {
  return get(var, name.c_str());
}

int QuiltMesh::get(Field3D &var, const char *name) {
  
}

int QuiltMesh::get(Field3D &var, const string &name) {
  return get(var, name.c_str());
}

/****************************************************************
 * Communications
 ****************************************************************/

int QuiltMesh::communicate(FieldGroup &g) {
  comm_handle c = send(g);
  return wait(c);
}

comm_handle QuiltMesh::send(FieldGroup &g) {
  /// Start timer
  Timer timer("comms");
  
  /// Get the list of variables to send
  vector<FieldData*> var_list = g.get();
  
  // Iterate over boundaries
  for(QuiltDomain::bndry_iterator it = mydomain->bndry_begin(); it != mydomain->bndry_end(); it++) {
    // Post recieves
    Domain::Bndry* b = *it;
     
    QuiltDomain* to = (QuiltDomain*) b->getNeighbour(mydomain);
    
    // Check which boundary(s) this is on
    for(int i = 0; i < 4; i++) {
      Domain::BndrySide s = static_cast<Domain::BndrySide>(i);
      if(b->onSide(mydomain, s)) {
        // get range of indices
        int min = b->getMin(s);
        int max = b->getMax(s);
      }
    }
  }
  
  for(QuiltDomain::bndry_iterator it = mydomain->bndry_begin(); it != mydomain->bndry_end(); it++) {
    // Send data
    
  }
 
  
}

int QuiltMesh::wait(comm_handle handle) {
  if(handle == NULL)
    return 1;
  
  QMCommHandle *ch = (QMCommHandle*) handle;
  
  if(!ch->inProgress)
    return 2;
  
  /*
  BoutReal t = MPI_Wtime(); // Starting time
  
  if(ch->var_list.empty()) {
    // just waiting for communications, not unpacking into variables
    
    MPI_Waitall(ch->request.size(),
		&ch->request[0],
		MPI_STATUSES_IGNORE);
    
    freeHandle(ch);
    
    // Add the time elapsed to the communications wall time
    wtime_comms += MPI_Wtime() - t;
    
    return 0;
  }
  
  // Waiting for data, and unpacking the result into variables
  
  do {
    int ind;
    MPI_Waitany(ch->request.size(),   // How many requests 
		&ch->request[0],      // Array of requests
		&ind,   // Index into the array: which one arrived?
		MPI_STATUSES_IGNORE); // Don't store statuses for now
    
    if(ind == MPI_UNDEFINED)
      break; // Finished
    
    // Unpack data into variables
    unpackData(ch->buffer[ind], ch->range[ind], ch->var_list);
    
    // Twist-shift condition
    if(TwistShift && (mesh->TwistOrder == 0) && ch->range[ind]->zshift) {
      GuardRange *r = ch->range[ind];
      for(vector<FieldData*>::iterator it = ch->var_list.begin(); it != ch->var_list.end(); it++)
	if((*it)->is3D()) {
	  // Only need to shift 3D variables
	  for(int jx=r->xmin;jx<=r->xmax;jx++)
	    for(int jy=r->ymin;jy <= r->ymax; jy++)
	      (*it)->shiftZ(jx, jy, r->shiftAngle[jx-r->xmin]);
	}
    }
    
    // Mark the request as NULL
    ch->request[ind] = MPI_REQUEST_NULL;
  }while(1);
  
  freeHandle(ch);

  wtime_comms += MPI_Wtime() - t;
  */
  return 0;
}


/****************************************************************
 * X communications
 ****************************************************************/

bool QuiltMesh::firstX() {
  
}

bool QuiltMesh::lastX() {
  
}

int QuiltMesh::sendXOut(BoutReal *buffer, int size, int tag) {
  if(lastX())
    return 1;

  return 0;
}

int QuiltMesh::sendXIn(BoutReal *buffer, int size, int tag) {
  if(firstX())
    return 1;
  
  return 0;
}

comm_handle QuiltMesh::irecvXOut(BoutReal *buffer, int size, int tag) {
  if(lastX())
    return NULL; // No processor to receive from
  
  BoutReal t = MPI_Wtime();

  QMCommHandle *ch = getHandle(1); // Just one request
  
  return (comm_handle) ch;
}

comm_handle QuiltMesh::irecvXIn(BoutReal *buffer, int size, int tag) {
  if(firstX())
    return NULL;
  
  BoutReal t = MPI_Wtime();
  
  QMCommHandle *ch = getHandle(1);
 
  return (comm_handle) ch;
}

MPI_Comm QuiltMesh::getXcomm() const {
  
}

SurfaceIter* QuiltMesh::iterateSurfaces() {

}

const Field2D QuiltMesh::averageY(const Field2D &f) {

}

const Field3D QuiltMesh::averageY(const Field3D &f) {

}

bool QuiltMesh::surfaceClosed(int jx, BoutReal &ts) {
  
}

const RangeIterator QuiltMesh::iterateBndryLowerY() const {

}

const RangeIterator QuiltMesh::iterateBndryUpperY() const {
  
}

vector<BoundaryRegion*> QuiltMesh::getBoundaries() {
  
}

BoutReal QuiltMesh::GlobalX(int jx) {

}

BoutReal QuiltMesh::GlobalY(int jy) {
  
}

void QuiltMesh::outputVars(Datafile &file) {
  
}

int QuiltMesh::XGLOBAL(int xloc) {
  
}

int QuiltMesh::YGLOBAL(int yloc) {
  
}

/****************************************************************
 * Private functions
 ****************************************************************/

void QuiltMesh::freeHandle(QMCommHandle *h) {
  h->inProgress = false;
  comm_list.push_back(h);
}

QuiltMesh::QMCommHandle* QuiltMesh::getHandle(int n) {
  QMCommHandle* h;
  h->inProgress = false;

  return h;
}

void QuiltMesh::packData(const vector<FieldData*> &vars, GuardRange* range, vector<BoutReal> &data) {
  
  /// Get size of buffer needed
  int len = msg_len(vars, range->xmin, range->xmax+1, range->ymin, range->ymax+1);
  
  data.resize(len);
  
  len = 0;
  // Loop over variables
  for(vector<FieldData*>::const_iterator it = vars.begin(); it != vars.end(); it++) {
    if((*it)->is3D()) {
      // 3D variable
      for(int x=range->xmin; x <= range->xmax; x++)
	for(int y=range->ymin; y <= range->ymax; y++)
	  for(int z=0;z<ngz-1;z++)
	    len += (*it)->getData(x,y,z,&data[len]);
    }else {
      // 2D variable
      for(int x=range->xmin; x <= range->xmax; x++)
	for(int y=range->ymin; y <= range->ymax; y++)
	  len += (*it)->getData(x,y,0,&data[len]);
    }
  }
}

void QuiltMesh::unpackData(vector<BoutReal> &data, GuardRange* range, vector<FieldData*> &vars) {
  int len = 0;
  // Loop over variables
  for(vector<FieldData*>::iterator it = vars.begin(); it != vars.end(); it++) {
    if((*it)->is3D()) {
      // 3D variable
      for(int x=range->xmin; x <= range->xmax; x++)
	for(int y=range->ymin; y <= range->ymax; y++)
	  for(int z=0;z<ngz-1;z++)
	    len += (*it)->setData(x,y,z,&data[len]);
    }else {
      // 2D variable
      for(int x=range->xmin; x <= range->xmax; x++)
	for(int y=range->ymin; y <= range->ymax; y++)
	  len += (*it)->setData(x,y,0,&data[len]);
    }
  }
}

void QuiltMesh::read2Dvar(GridDataSource *s, const char *name, 
                          int xs, int ys,
                          int xd, int yd,
                          int nx, int ny,
                          BoutReal **data) {
  
}
