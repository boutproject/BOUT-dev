
#include <globals.hxx>
#include <quiltmesh.hxx>
#include <boutexception.hxx>

#include <algorithm>
#include <numeric>
using std::accumulate;
#include <cmath>

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

QuiltMesh::~QuiltMesh()
{
  
}


int QuiltMesh::load()
{
  // Get number of processors
  int NPES, MYPE;
  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);
  
  // Get number of regions
  int nregions;
  if(Mesh::get(nregions, "nregions"))
    throw BoutException("Mesh file doesn't have nregions variable. Incorrect format?\n");
  
  if(nregions > NPES)
    throw BoutException("Can't divide %d regions between %d processors\n", nregions, NPES);
  
  // Get array of grid sizes
  vector<int> nx = readInts("nx", nregions);
  vector<int> ny = readInts("ny", nregions);
  
  // Partition the mesh
  domains = partition(nx, ny, NPES);
  
  
}

/****************************************************************
 * Getting variables
 ****************************************************************/
  
int QuiltMesh::get(Field2D &var, const char *name, BoutReal def) {

}

int QuiltMesh::get(Field2D &var, const string &name, BoutReal def) {

}

int QuiltMesh::get(Field3D &var, const char *name) {

}

int QuiltMesh::get(Field3D &var, const string &name) {

}

/****************************************************************
 * Communications
 ****************************************************************/

int QuiltMesh::communicate(FieldGroup &g) {
  comm_handle c = send(g);
  return wait(c);
}

comm_handle QuiltMesh::send(FieldGroup &g) {
  /// Record starting wall-time
  BoutReal t = MPI_Wtime();
  
  /// Get the list of variables to send
  vector<FieldData*> var_list = g.get();
  
  QMCommHandle *ch = getHandle(0);

  //////////////////////////////
  /// Post recieves
  
  if(mydomain->xin != NULL) {
    // Inner processor
    
    /*
    MPI_Irecv(ch.umsg_recvbuff,
	      msg_len(var_list, 0, MXG, 0, mydomain->ny),
              PVEC_REAL_MPI_TYPE,
              xin->destination->proc; // Destination processor
	      OUT_SENT_IN,
	      BoutComm::get(),
	      &ch.request[0]);
    */
  }
  if(mydomain->xout != NULL) {
    // Outer processor
    
  }
  
  vector<GuardRange*>::iterator it;
  // Iterate over the upper guard cell ranges
  for(it = mydomain->yup.begin(); it != mydomain->yup.end(); it++) {
    
  }
  
  // Iterate over lower guard cell ranges
  for(it = mydomain->ydown.begin(); it != mydomain->ydown.end(); it++) {
    
  }

  //////////////////////////////
  /// Send data
  
  if(mydomain->xin != NULL) {
    // Inner processor
    
  }
  if(mydomain->xout != NULL) {
    // Outer processor
    
  }
  
  // Iterate over the upper guard cell ranges
  for(it = mydomain->yup.begin(); it != mydomain->yup.end(); it++) {
    
  }
  
  // Iterate over lower guard cell ranges
  for(it = mydomain->ydown.begin(); it != mydomain->ydown.end(); it++) {
    
  }
}

int QuiltMesh::wait(comm_handle handle) {
  if(handle == NULL)
    return 1;
  
  QMCommHandle *ch = (QMCommHandle*) handle;
  
  if(!ch->inProgress)
    return 2;
  
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
  
  return 0;
}


/****************************************************************
 * X communications
 ****************************************************************/

bool QuiltMesh::firstX() {
  return mydomain->xin == NULL;
}

bool QuiltMesh::lastX() {
  return mydomain->xout == NULL;
}

int QuiltMesh::sendXOut(BoutReal *buffer, int size, int tag) {
  if(lastX())
    return 1;
  
  BoutReal t = MPI_Wtime();

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   mydomain->xout->proc,
	   tag,
	   BoutComm::get());
  
  wtime_comms += MPI_Wtime() - t;

  return 0;
}

int QuiltMesh::sendXIn(BoutReal *buffer, int size, int tag) {
  if(firstX())
    return 1;
  
  BoutReal t = MPI_Wtime();

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE,
	   mydomain->xin->proc,
	   tag,
	   BoutComm::get());
  
  wtime_comms += MPI_Wtime() - t;

  return 0;
}

comm_handle QuiltMesh::irecvXOut(BoutReal *buffer, int size, int tag) {
  if(lastX())
    return NULL; // No processor to receive from
  
  BoutReal t = MPI_Wtime();

  QMCommHandle *ch = getHandle(1); // Just one request
  
  ch->tag[0] = tag;
  ch->var_list.clear();

  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    mydomain->xout->proc, // Processor number
	    tag,
	    BoutComm::get(),
	    &ch->request[0]);
  
  ch->inProgress = true;
  
  wtime_comms += MPI_Wtime() - t;
  
  return (comm_handle) ch;
}

comm_handle QuiltMesh::irecvXIn(BoutReal *buffer, int size, int tag) {
  if(firstX())
    return NULL;
  
  BoutReal t = MPI_Wtime();
  
  QMCommHandle *ch = getHandle(1);
  
  ch->tag[0] = tag;
  ch->var_list.clear();
  
  MPI_Irecv(buffer,
	    size,
	    PVEC_REAL_MPI_TYPE,
	    mydomain->xin->proc, // Processor number
	    tag,
	    BoutComm::get(),
	    &ch->request[0]);
  
  ch->inProgress = true;
  
  wtime_comms += MPI_Wtime() - t;
  
  return (comm_handle) ch;
}

/****************************************************************
 * Private functions
 ****************************************************************/

const vector<int> QuiltMesh::readInts(const string &name, int n)
{
  vector<int> result;
  
  // First get a data source
  GridDataSource* s = findSource(name);
  if(s) {
    s->open(name);
    s->setOrigin();
    result.resize(n);
    if(!s->fetch(&(result.front()), name, n)) {
      // Error reading
      s->close();
      throw BoutException("Could not read integer array '%s'\n", name.c_str());
    }
    s->close();
  }else {
    // Not found
    throw BoutException("Missing integer array %s\n", name.c_str());
  }
  
  return result;
}


///////////////////////////////////////////////////////////
// Partition this grid between NPES processors.
// Try to make the number of points on each processor approximately
// equal, and make each domain as square as possible (minimise comms)
vector<QuiltMesh::MeshDomain*> QuiltMesh::partition(const vector<int> &nx, 
                                         const vector<int> &ny, 
                                         int NPES) {
  int nregions = nx.size();
  
  // Number of points in each region ntot = nx * ny
  vector<int> ntot(nregions);
  transform(nx.begin(), nx.end(), 
            ny.begin(),
            ntot.begin(),
            std::multiplies<int>());
  
  // Sum ntot to get total number of points in the grid
  int alln =  accumulate( ntot.begin(), ntot.end(), 0 );
  
  // Allocate processors to regions
  vector<int> nproc(nregions);
  for(int i=0;i<nregions;i++) {
    nproc[i] = 1; // Need at least one per region
    // Share the others out
    nproc[i] += (int) ((NPES-nregions) * ntot[i]) / alln;
  }
  // Find number of processors remaining
  int npremain = NPES - accumulate( nproc.begin(), nproc.end(), 0 );
  
  if(npremain != 0) {
    // Assign extra points to the regions with highest load
    vector<BoutReal> load(nregions);
    for(int i=0;i<nregions;i++)
      load[i] = ((BoutReal) ntot[i]) / ((BoutReal) nproc[i]);
    
    for(int p=0;p<npremain;p++) {
      // Get maximum load index
      int ind = 0;
      BoutReal maxload = load[0];
      for(int i=1;i<nregions;i++)
        if(load[i] > maxload) {
          ind = i;
          maxload = load[i];
        }else if(fabs(load[i] - maxload) < 1.e-5) {
          // Equal loads, so add to the one with the largest number of points
          // (minimise difference in loads after)
          if(ntot[i] > ntot[ind])
            ind = i;
        }
      // Add a processor to this region
      nproc[ind]++;
      load[ind] = ((BoutReal) ntot[ind]) / ((BoutReal) nproc[ind]);
    }
  }
  
  ///////////////////////////////////////////////////////////
  // Each region now has ntot points and nproc processors
  // Divide up each region to produce a set of MeshDomains
  
  vector<MeshDomain*> result;
  vector<int> nxproc(nregions); // Number of X processors in each region

  int proc = 0; // Processor number
  for(int r = 0; r < nregions; r++) { // Loop over the regions
    
    // Allocate all processors first to make it easier to link together
    for(int i=0;i<nproc[r];i++)
      result.push_back(new MeshDomain());
    
    // Calculate number of processors in X for a square domain (mxsub = mysub)
    int nxp_sq = (int) (((BoutReal) nx[r])*sqrt(((BoutReal) nproc[r]) / ((BoutReal) ntot[r])));
    
    int nxp;
    for(nxp=nxp_sq; nxp > 0; nxp--) {
      if(nproc[r] % nxp == 0) {
        break;
      }
    }
    nxproc[r] = nxp;
    
    int nyp = nproc[r] / nxp;

    // Divide up into nyp strips
    int nyall = (int) ny[r] / nyp;
    int nyextra = ny[r] % nyp;
    
    int nxall = (int) nx[r] / nxp;
    int nxextra = nx[r] % nxp;
    
    // Loop over X processors fastest so X processors are close
    // together (for inversions, tightly coupled in X).
    
    int y0 = 0;
    for(int yp = 0; yp < nyp; yp++) {
      int nysub = nyall; // ny on these processors
      if(yp < nyextra)
        nysub++; // Extra y point in this processor
      
      int x0 = 0;
      for(int xp = 0; xp < nxp; xp++) {
        int nxsub = nxall; // nx on these processors
        if(xp < nxextra)
          nxsub++;
        
        // Create a new domain
        MeshDomain *d = result[proc];
        d->proc = proc;
        d->nx = nxsub;
        d->ny = nysub;
        d->x0 = x0;
        d->y0 = y0;
        
        // Inside a region, fill in the communications.
        // Communication between regions is more complicated
        // and needs to be done after all regions have been created

        if(xp == 0) {
          d->xin = NULL; // No inner X processor
        }else
          d->xin = result[proc-1];
        if(xp == (nxp-1)) {
          d->xout = NULL;
        }else
          d->xout = result[proc+1];
        
        if(yp > 0) {
          // Communications going down
          GuardRange *gd = new GuardRange;
          gd->xmin = 0;
          gd->xmax = nxsub-1;
          gd->destination = result[proc-nxp];
          gd->xshift = 0;
          d->ydown.push_back(gd);
        }
        if(yp < (nyp-1)) {
          // Communications going up
          GuardRange *gu = new GuardRange;
          gu->xmin = 0;
          gu->xmax = nxsub-1;
          gu->destination = result[proc+nxp];
          gu->xshift = 0;
          d->yup.push_back(gu);
        }
        
        proc++;
        x0 += nxsub;
      }
      y0 += nysub;
    }
  }
  
  // Have all regions, connected internally.
  // Now connect regions together

  return result;
}

void QuiltMesh::freeHandle(QMCommHandle *h) {
  h->inProgress = false;
  comm_list.push_back(h);
}

QuiltMesh::QMCommHandle* QuiltMesh::getHandle(int n) {
  QMCommHandle* h;
  if(!comm_list.empty()) {
    h = comm_list.front();
    comm_list.pop_front();
  }else
    h = new QMCommHandle;
  
  h->buffer.resize(n);
  h->request.resize(n);
  h->tag.resize(n);
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
