
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
#include <utils.hxx>

#include "quiltmesh.hxx"
#include "../partition.hxx"

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

QuiltMesh::QuiltMesh(GridDataSource *s, Options *opt) : Mesh(s), options(opt) {
  if(options == NULL)
    options = Options::getRoot()->getSection("mesh");
}

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
  
  // Read some options
  options->get("MXG", MXG, 2);
  options->get("MYG", MYG, 2);
  if(options->isSet("MZ")) {
    options->get("MZ", ngz, 65);
  }else
    Options::getRoot()->get("MZ", ngz, 65);
  if(!is_pow2(ngz-1)) {
    if(is_pow2(ngz)) {
      ngz++;
      output.write("WARNING: Number of toroidal points increased to %d\n", ngz);
    }else {
      output.write("ERROR: Number of toroidal points must be 2^n + 1");
      return 1;
    }
  }

  // Create domains
  output << "\tCreating domains...";
  vector<QuiltDomain*> domains;
  
  for(int i=0;i<nregions;i++)
    domains.push_back(new QuiltDomain(nx[i], ny[i]));
  

  // Join domains together
  
  
  // Partition the mesh (all domains connected together)
  output << "done\n\tPartitioning...";
  partitionAll(domains[0], NPES);
  

  // Assign domains to processors
  output << "done\n\tAssigning...";
  int p = 0;
  for(Domain::iterator it=domains[0]->begin(); it != domains[0]->end(); it++) {
    // For now very simple numbering. Should cluster nearby domains
    QuiltDomain *qd = static_cast<QuiltDomain*>( &(*it) );
    
    // Check that domain is large enough
    if( (qd->xSize() < MXG) || (qd->ySize() < MYG) )
      throw BoutException("Domain size too small (< %dx%d). Try fewer processors", MXG, MYG);
    
    qd->proc = p; p++;
    if(p == MYPE)
      mydomain = qd; // Domain for this processor
  }
  
  // Get size of the mesh
  ngx = mydomain->xSize() + 2*MXG;
  ngy = mydomain->ySize() + 2*MYG;
  
  // Local ranges
  xstart = MXG;
  xend   = ngx-1-MXG;
  ystart = MYG;
  yend   = ngy-1-MYG;
  
  // Create ranges for boundary iteration
  ylow_range = yhigh_range = RangeIterator(0, ngx-1);
  xlow_range = xhigh_range = RangeIterator(0, ngy-1);
  
  // Iterate through boundaries to get guard cell regions
  for(QuiltDomain::bndry_iterator it = mydomain->bndry_begin(); it != mydomain->bndry_end(); it++) {
    Domain::Bndry* b = *it;
     
    QuiltDomain* to = (QuiltDomain*) b->getNeighbour(mydomain);
    
    // Check which boundary(s) this is on
    // (multiple boundaries in case of periodic domain)
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
      
      all_guards.push_back(g);
      xlow_guards.push_back(g);
      xlow_range -= RangeIterator(min, max);
    }
    if(b->onSide(mydomain, Domain::xhigh)) {
      // get range of indices
      int min = b->getMin(Domain::xhigh);
      int max = b->getMax(Domain::xhigh);
      
      // Create new GuardRange 
      GuardRange *g = new GuardRange();
      g->xmin = ngx-1-MXG;
      g->xmax = ngx-1;
      g->ymin = min;
      g->ymax = max;
      g->proc = to->proc;
      
      all_guards.push_back(g);
      xhigh_guards.push_back(g);
      xhigh_range -= RangeIterator(min, max);
    }
    if(b->onSide(mydomain, Domain::ylow)) {
      // get range of indices
      int min = b->getMin(Domain::ylow);
      int max = b->getMax(Domain::ylow);
      
      // Create new GuardRange 
      GuardRange *g = new GuardRange();
      g->ymin = 0;
      g->ymax = MYG;
      g->xmin = min;
      g->xmax = max;
      g->proc = to->proc;
      
      all_guards.push_back(g);
      
      // Remove range from boundary
      ylow_range -= RangeIterator(min, max);
    }
    if(b->onSide(mydomain, Domain::yhigh)) {
      // get range of indices
      int min = b->getMin(Domain::yhigh);
      int max = b->getMax(Domain::yhigh);
      
      // Create new GuardRange 
      GuardRange *g = new GuardRange();
      g->ymin = ngy-1-MYG;
      g->ymax = ngy-1;
      g->xmin = min;
      g->xmax = max;
      g->proc = to->proc;
      
      all_guards.push_back(g);
      
      // Remove range from boundary
      yhigh_range -= RangeIterator(min, max);
    }
  }
  
  // Turn boundary ranges into BoundaryRegion objects
  for(RangeIterator *r = &xlow_range; r != 0; r=r->nextRange())
    boundaries.push_back(new BoundaryRegionXIn("xin", r->min(), r->max()));
  
  for(RangeIterator *r = &xhigh_range; r != 0; r=r->nextRange())
    boundaries.push_back(new BoundaryRegionXOut("xout", r->min(), r->max()));
  
  for(RangeIterator *r = &ylow_range; r != 0; r=r->nextRange())
    boundaries.push_back(new BoundaryRegionYDown("ydown", r->min(), r->max()));
  
  for(RangeIterator *r = &yhigh_range; r != 0; r=r->nextRange())
    boundaries.push_back(new BoundaryRegionYUp("yup", r->min(), r->max()));
  
  // Get communicators for X and Y communication
  

  output << "done\n";
  
  return 0;
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
  
  // Get an empty communication handle to return
  QMCommHandle *handle = getCommHandle();
  
  /// Get the list of variables to send
  handle->var_list = g.get();
  
  // Iterate over guard cells
  for(vector<GuardRange*>::const_iterator it=all_guards.begin(); it != all_guards.end(); it++) {
    // Post recieves
    GuardRange* g = *it;
    
    // Calculate length of the message
    int len = msg_len(handle->var_list, g->xmin, g->xmax+1, g->ymin, g->ymax+1);
    
    // Get a request handle of this size
    QMRequest* rq = getRequestHandle(len);
    rq->guard = g;
    
    MPI_Irecv(&(rq->buffer[0]),
	      len,
	      PVEC_REAL_MPI_TYPE,
	      g->proc,
	      1234,
	      BoutComm::get(),
	      &(rq->request));
    
    // Add request to the CommHandle
    handle->request.push_back(rq);
    handle->mpi_rq.push_back(rq->request);
  }
  
  for(vector<GuardRange*>::const_iterator it=all_guards.begin(); it != all_guards.end(); it++) {
    // Send data
    GuardRange* g = *it;
    
    // Calculate length of the message
    int len = msg_len(handle->var_list, g->xmin, g->xmax+1, g->ymin, g->ymax+1);
    
    // Get a request handle of this size
    QMRequest* rq = getRequestHandle(len);
    rq->guard = NULL; // Signals that no data is to be copied on success
    
    // Copy data into the buffer
    packData(handle->var_list, g, rq->buffer);
    
    if(async_send) {
      // Asyncronous sending
      MPI_Isend(&(rq->buffer[0]), // Buffer to send
                len,              // Length of buffer in BoutReals
                PVEC_REAL_MPI_TYPE,  // Real variable type
                g->proc,          // Destination processor
                1234,             // Label (tag) for the message
                BoutComm::get(),  // Communicator
                &(rq->request));  // MPI Request handle
      
      // Add request to Comm Handle
      handle->request.push_back(rq);
      handle->mpi_rq.push_back(rq->request);
    }else {
      MPI_Send(&(rq->buffer[0]),
               len,
               PVEC_REAL_MPI_TYPE, 
               g->proc,
               1234,
               BoutComm::get());
      
      // Free request handle
      freeRequestHandle(rq);
    }
  }
  
  handle->inProgress = true;
  
  return static_cast<comm_handle>(handle);
}

int QuiltMesh::wait(comm_handle handle) {
  if(handle == NULL)
    return 1;

  QMCommHandle *ch = (QMCommHandle*) handle;
  if(!ch->inProgress)
    return 2;
  
  /// Start timer
  Timer timer("comms");

  do {
    int ind;
    MPI_Waitany(ch->mpi_rq.size(),   // How many requests 
		&ch->mpi_rq[0],      // Array of requests
		&ind,   // Index into the array: which one arrived?
		MPI_STATUSES_IGNORE); // Don't store statuses for now
    if(ind == MPI_UNDEFINED)
      break; // Finished
    
    // Check if we need to unpack data
    if(ch->request[ind]->guard != NULL) {
      unpackData(ch->request[ind]->buffer, ch->request[ind]->guard, ch->var_list);
      
      /*
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
      */
    }
    // Mark the request as NULL
    ch->mpi_rq[ind] = MPI_REQUEST_NULL;
  }while(1);
  
  // Free handle
  freeCommHandle(ch);
  
  return 0;
}

/***************************************************************
 *             Non-Local Communications
 ***************************************************************/

MPI_Request QuiltMesh::sendToProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) {
  throw BoutException("sendToProc() is not implemented in QuiltMesh");
}

comm_handle QuiltMesh::receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) {
  throw BoutException("receiveFromProc() is not implemented in QuiltMesh");
}

int QuiltMesh::getNXPE() {
  throw BoutException("getNXPE() makes no sense in QuiltMesh");
}

int QuiltMesh::getNYPE() {
  throw BoutException("getNYPE() makes no sense in QuiltMesh");
}

int QuiltMesh::getXProcIndex() {
  throw BoutException("getXProcIndex() makes no sense in QuiltMesh");
}

int QuiltMesh::getYProcIndex() {
  throw BoutException("getYProcIndex() makes no sense in QuiltMesh");
}

/****************************************************************
 * X communications
 ****************************************************************/

bool QuiltMesh::firstX() {
  return xlow_guards.size() == 0;
}

bool QuiltMesh::lastX() {
  return xhigh_guards.size() == 0;
}

int QuiltMesh::sendXOut(BoutReal *buffer, int size, int tag) {
  if(lastX())
    return 1;

  // Get an empty communication handle to return
  QMCommHandle *handle = getCommHandle();
  
  

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
  
  

  QMCommHandle *ch = getCommHandle(); // Just one request
  
  return (comm_handle) ch;
}

comm_handle QuiltMesh::irecvXIn(BoutReal *buffer, int size, int tag) {
  if(firstX())
    return NULL;
  
  
  
  QMCommHandle *ch = getCommHandle();
 
  return (comm_handle) ch;
}

MPI_Comm QuiltMesh::getXcomm() const {
  
}

MPI_Comm QuiltMesh::getYcomm(int jx) const {
  
}

const Field2D QuiltMesh::averageY(const Field2D &f) {
  
}

const Field3D QuiltMesh::averageY(const Field3D &f) {
  
}

bool QuiltMesh::periodicY(int jx, BoutReal &ts) const {
  return yperiodic.intersects(jx);
}

const RangeIterator QuiltMesh::iterateBndryLowerY() const {
  return ylow_range;
}

const RangeIterator QuiltMesh::iterateBndryUpperY() const {
  return yhigh_range;
}

vector<BoundaryRegion*> QuiltMesh::getBoundaries() {
  return boundaries;
}

BoutReal QuiltMesh::GlobalX(int jx) const {
  return ((BoutReal) XGLOBAL(jx));
}

BoutReal QuiltMesh::GlobalY(int jy) const {
  return ((BoutReal) YGLOBAL(jy));
}

void QuiltMesh::outputVars(Datafile &file) {
  
}

int QuiltMesh::XGLOBAL(int xloc) const {
  return mydomain->xOrigin() + xloc - MXG;
}

int QuiltMesh::YGLOBAL(int yloc) const {
  return mydomain->yOrigin() + yloc - MYG;
}

/****************************************************************
 * Private functions
 ****************************************************************/

QuiltMesh::QMCommHandle* QuiltMesh::getCommHandle() {
  return NULL;
}

void QuiltMesh::freeCommHandle(QMCommHandle *h) {
  h->inProgress = false;
  comm_list.push_back(h);
}

QuiltMesh::QMRequest* QuiltMesh::getRequestHandle(int n) {
  return NULL;
}

void QuiltMesh::freeRequestHandle(QMRequest* r) {
  
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


bool QuiltMesh::firstY() {
  throw BoutException("firstY() is not implemented in QuiltMesh");
}

bool QuiltMesh::lastY() {
  throw BoutException("lastY() is not implemented in QuiltMesh");
}

int QuiltMesh::UpXSplitIndex() {
  throw BoutException("UpXSplitIndex() is not implemented in QuiltMesh");
}

int QuiltMesh::DownXSplitIndex() {
  throw BoutException("DownXSplitIndex() is not implemented in QuiltMesh");
}

int QuiltMesh::sendYOutIndest(BoutReal *buffer, int size, int tag) {
  throw BoutException("sendYOutIndest() is not implemented in QuiltMesh");
}

int QuiltMesh::sendYOutOutdest(BoutReal *buffer, int size, int tag) {
  throw BoutException("sendYOutOutdest() is not implemented in QuiltMesh");
}

int QuiltMesh::sendYInIndest(BoutReal *buffer, int size, int tag) {
  throw BoutException("sendYInIndest() is not implemented in QuiltMesh");
}

int QuiltMesh::sendYInOutdest(BoutReal *buffer, int size, int tag) {
  throw BoutException("sendYInOutdest() is not implemented in QuiltMesh");
}

comm_handle QuiltMesh::irecvYOutIndest(BoutReal *buffer, int size, int tag) {
  throw BoutException("irecvYOutIndest() is not implemented in QuiltMesh");
}

comm_handle QuiltMesh::irecvYOutOutdest(BoutReal *buffer, int size, int tag) {
  throw BoutException("irecvYOutOutdest() is not implemented in QuiltMesh");
}

comm_handle QuiltMesh::irecvYInIndest(BoutReal *buffer, int size, int tag) {
  throw BoutException("irecvYInIndest() is not implemented in QuiltMesh");
}

comm_handle QuiltMesh::irecvYInOutdest(BoutReal *buffer, int size, int tag) {
  throw BoutException("irecvYInOutdest() is not implemented in QuiltMesh");
}
