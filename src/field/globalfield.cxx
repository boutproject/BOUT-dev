
#include <bout/globalfield.hxx>
#include <boutexception.hxx>
#include <boutcomm.hxx>

GlobalField::GlobalField(Mesh *m, int xsize, int ysize, int zsize) 
  : mesh(m), nx(xsize), ny(ysize), nz(zsize) {
  
  comm = BoutComm::get(); // This should come from Mesh
  
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);
  
}

void GlobalField::proc_origin(int proc, int *x, int *y, int *z) {
  // This is a hack, and should be improved with a better Mesh interface
  
  // Get the number of processors in X and Y
  int nxpe = mesh->getNXPE();
  int nype = mesh->getNYPE();
  
  // Get the X and Y indices
  int pex = proc % nxpe;
  int pey = proc / nxpe;
  
  // Get the size of the processor domain
  int nx = mesh->xend - mesh->xstart + 1;
  int ny = mesh->yend - mesh->ystart + 1;

  // Set the origin values
  *x = pex * nx;
  *y = pey * ny;
  if(z != NULL)
    *z = 0;
}

void GlobalField::proc_size(int proc, int *lx, int *ly, int *lz) {
  // Get the size of the processor domain. 
  // Assumes every processor has the same size domain
  
  *lx = mesh->xend - mesh->xstart + 1;
  *ly = mesh->yend - mesh->ystart + 1;
  if(lz != NULL)
    *lz = mesh->ngz-1;
}

///////////////////////////////////////////////////////////////////////////////////////////

GlobalField2D::GlobalField2D(Mesh *m) : GlobalField(m, mesh->GlobalNx, mesh->GlobalNy, 1){
  
  buffer = new BoutReal*[npes];
  for(int p=0;p<npes;p++)
    buffer[p] = new BoutReal[msg_len(p)];
}

GlobalField2D::~GlobalField2D() {
  for(int p=0;p<npes;p++)
    delete[] buffer[p];
  delete[] buffer;
}

void GlobalField2D::gather(const Field2D &f, int proc) {
  // Gather all data onto processor 'proc'
  
  if((proc < 0) || (proc >= npes))
    throw BoutException("Processor out of range");

  MPI_Request req[npes]; // Array of receive handles
  
  if(mype == proc) {
    // This processor will receive the data
    
    // Post receives
    for(int p = 0; p < npes; p++) {
      if( p != mype ) {
        // Check size of the array
        
        MPI_Irecv(buffer[p], msg_len(p), MPI_DOUBLE, p, 3141, comm, &req[p]);
      }
    }
    
    // Copy data from this processor
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) {
        (*this)(x+xorig,y+yorig) =  f(local_xorig+x, local_yorig+y);
      }
    
    if(npes > 1) {
      // Wait for receives, process as they arrive
      int pe;
      MPI_Status status;
      do {
        MPI_Waitany(npes-1, req, &pe, &status);
        
        // Unpack data from processor 'pe'
        int xorig, yorig;
        proc_origin(pe, &xorig, &yorig);
        int xsize, ysize;
        proc_size(pe, &xsize, &ysize);
        
        for(int x=0;x<xsize;x++)
          for(int y=0;y<ysize;y++) {
            (*this)(x+xorig,y+yorig) =  buffer[pe][x*ysize + y];
          }
        
        if(pe != MPI_UNDEFINED)
          req[pe] = MPI_REQUEST_NULL;
      }while(pe != MPI_UNDEFINED);
    }
  }else {
    // Sending data to proc
    BoutReal *d = f.getData()[0]; // Pointer to start of data
    
    MPI_Send(d, msg_len(mype), MPI_DOUBLE, proc, 3141, comm);
  }
  
  data_on_proc = proc;
}

const Field2D GlobalField2D::scatter() const {
  Field2D result;
  result.allocate();
  
  MPI_Status status;
  if(mype == data_on_proc) {
    // Data is on this processor. Send to other processors
    
    for(int p = 0; p < npes; p++) {
      if(p == mype) continue;
      
      int xorig, yorig;
      proc_origin(pe, &xorig, &yorig);
      int xsize, ysize;
      proc_size(pe, &xsize, &ysize);
      
      // Send to processor p
      for(int x=0;x<xsize;x++)
        for(int y=0;y<ysize;y++) {
          buffer[p][x*ysize + y] = (*this)(x+xorig,y+yorig);
        }
      
      MPI_Send(buffer[p], xsize*ysize, MPI_DOUBLE, p, 1413, comm);
    }

    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    
    // Copy to result
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) {
        result(local_xorig+x,local_yorig+y) = (*this)(x+xorig,y+yorig);
      }
  }else {
    // Receive data
    MPI_Recv(data, msg_len(mype), MPI_DOUBLE, data_on_proc, 1413, comm, &status);
    
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) {
        result(local_xorig+x,local_yorig+y) = data[x*ysize + y];
      }
  }
  return result;
}

int GlobalField2D::msg_len(int proc) const {
  return 1;
}
