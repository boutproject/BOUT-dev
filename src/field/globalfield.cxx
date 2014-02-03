
#include <bout/globalfield.hxx>
#include <boutexception.hxx>
#include <boutcomm.hxx>

GlobalField::GlobalField(Mesh *m, int proc, int xsize, int ysize, int zsize) 
  : mesh(m), data_on_proc(proc), nx(xsize), ny(ysize), nz(zsize) {
  
  comm = BoutComm::get(); // This should come from Mesh
  
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);
  
  if(nx*ny*nz <= 0)
    throw BoutException("GlobalField data must have non-zero size");

  if(mype == proc) {
    // Allocate memory
    data = new BoutReal[nx*ny*nz];
  }else {
    data = NULL;
  }
}

GlobalField::~GlobalField() {
  if(data)
    delete[] data;
}

void GlobalField::proc_local_origin(int proc, int *x, int *y, int *z) const {
  
  int nxpe = mesh->getNXPE();
  if(proc % nxpe == 0) {
    *x = 0;
  }else
    *x = mesh->xstart;
  
  *y = mesh->ystart;
  
  if(z != NULL)
    *z = 0;
}

void GlobalField::proc_origin(int proc, int *x, int *y, int *z) const {
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

  if(pex != 0)
    *x += mesh->xstart;
}

void GlobalField::proc_size(int proc, int *lx, int *ly, int *lz) const {
  // Get the size of the processor domain. 
  // Assumes every processor has the same size domain
  
  *lx = mesh->xend - mesh->xstart + 1;
  *ly = mesh->yend - mesh->ystart + 1;
  if(lz != NULL)
    *lz = mesh->ngz-1;
  
  int nxpe = mesh->getNXPE();
  int pex = proc % nxpe;
  if(pex == 0)
    *lx += mesh->xstart;
  if(pex == (nxpe-1))
    *lx += mesh->xstart;
}

///////////////////////////////////////////////////////////////////////////////////////////

GlobalField2D::GlobalField2D(Mesh *m, int proc) : GlobalField(m, proc, m->GlobalNx, m->GlobalNy-2*m->ystart, 1), 
  data_valid(false) {
  
  if((proc < 0) || (proc >= npes))
    throw BoutException("Processor out of range");
  
  if(mype == data_on_proc) {
    // Gathering onto this processor
    buffer = new BoutReal*[npes];
    for(int p=0;p<npes;p++)
      buffer[p] = new BoutReal[msg_len(p)];
  }else {
    buffer = new BoutReal*[1];
    buffer[0] = new BoutReal[msg_len(mype)];
  }
}

GlobalField2D::~GlobalField2D() {
  if(mype == data_on_proc) {
    for(int p=0;p<npes;p++)
      delete[] buffer[p];
  }else
    delete[] buffer[0];
  delete[] buffer;
}

void GlobalField2D::gather(const Field2D &f) {
  // Gather all data onto processor 'proc'
  
  if(mype == data_on_proc) {
    // This processor will receive the data
    
    MPI_Request req[npes]; // Array of receive handles

    // Post receives
    for(int p = 0; p < npes; p++) {
      if( p != mype ) {
        // Check size of the array
        MPI_Irecv(buffer[p], msg_len(p), MPI_DOUBLE, p, 3141, comm, &req[p]);
      }
    }
    
    req[mype] = MPI_REQUEST_NULL; // Mark as not used
    
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
        MPI_Waitany(npes, req, &pe, &status);

        if(pe != MPI_UNDEFINED) {
          // Unpack data from processor 'pe'
          int xorig, yorig;
          proc_origin(pe, &xorig, &yorig);
          int xsize, ysize;
          proc_size(pe, &xsize, &ysize);
          
          for(int x=0;x<xsize;x++)
            for(int y=0;y<ysize;y++) {
              (*this)(x+xorig,y+yorig) =  buffer[pe][x*ysize + y];
            }
          
          req[pe] = MPI_REQUEST_NULL;
        }
      }while(pe != MPI_UNDEFINED);
    }
  }else {
    // Sending data to proc
    
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) {
        buffer[0][x*ysize + y] = f(local_xorig+x, local_yorig+y);
      }

    MPI_Send(buffer[0], msg_len(mype), MPI_DOUBLE, data_on_proc, 3141, comm);
  }
  data_valid = true;
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
      proc_origin(p, &xorig, &yorig);
      int xsize, ysize;
      proc_size(p, &xsize, &ysize);
      
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
    MPI_Recv(buffer[0], msg_len(mype), MPI_DOUBLE, data_on_proc, 1413, comm, &status);
    
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) {
        result(local_xorig+x,local_yorig+y) = buffer[0][x*ysize + y];
      }
  }
  return result;
}

int GlobalField2D::msg_len(int proc) const {
  int xsize, ysize;
  proc_size(proc, &xsize, &ysize);
  return xsize * ysize;
}

///////////////////////////////////////////////////////////////////////////////////////////

GlobalField3D::GlobalField3D(Mesh *m, int proc) : GlobalField(m, proc, m->GlobalNx, m->GlobalNy-2*m->ystart, m->ngz-1), 
  data_valid(false) {
  
  if((proc < 0) || (proc >= npes))
    throw BoutException("Processor out of range");
  
  if(mype == data_on_proc) {
    // Gathering onto this processor
    buffer = new BoutReal*[npes];
    for(int p=0;p<npes;p++)
      buffer[p] = new BoutReal[msg_len(p)];
  }else {
    buffer = new BoutReal*[1];
    buffer[0] = new BoutReal[msg_len(mype)];
  }
}

GlobalField3D::~GlobalField3D() {
  if(mype == data_on_proc) {
    for(int p=0;p<npes;p++)
      delete[] buffer[p];
  }else
    delete[] buffer[0];
  delete[] buffer;
}

void GlobalField3D::gather(const Field3D &f) {
  // Gather all data onto processor 'proc'
  
  if(mype == data_on_proc) {
    // This processor will receive the data
    
    MPI_Request req[npes]; // Array of receive handles

    // Post receives
    for(int p = 0; p < npes; p++) {
      if( p != mype ) {
        // Check size of the array
        MPI_Irecv(buffer[p], msg_len(p), MPI_DOUBLE, p, 3141, comm, &req[p]);
      }
    }
    
    req[mype] = MPI_REQUEST_NULL; // Mark as not used
    
    // Copy data from this processor
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) 
        for(int z=0;z<mesh->ngz-1;z++) {
          (*this)(x+xorig,y+yorig,z) =  f(local_xorig+x, local_yorig+y, z);
        }
    
    if(npes > 1) {
      // Wait for receives, process as they arrive
      int pe;
      MPI_Status status;
      do {
        MPI_Waitany(npes, req, &pe, &status);

        if(pe != MPI_UNDEFINED) {
          // Unpack data from processor 'pe'
          int xorig, yorig;
          proc_origin(pe, &xorig, &yorig);
          int xsize, ysize;
          proc_size(pe, &xsize, &ysize);
          int zsize = mesh->ngz-1;
          
          for(int x=0;x<xsize;x++)
            for(int y=0;y<ysize;y++)
              for(int z=0;z<mesh->ngz-1;z++) {
                (*this)(x+xorig,y+yorig,z) =  buffer[pe][x*ysize*zsize + y*zsize + z];
              }
          
          req[pe] = MPI_REQUEST_NULL;
        }
      }while(pe != MPI_UNDEFINED);
    }
  }else {
    // Sending data to proc
    
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    int zsize = mesh->ngz-1;
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++)
        for(int z=0;z<mesh->ngz-1;z++) {
          buffer[0][x*ysize*zsize + y*zsize + z] = f(local_xorig+x, local_yorig+y, z);
        }
    
    MPI_Send(buffer[0], msg_len(mype), MPI_DOUBLE, data_on_proc, 3141, comm);
  }
  data_valid = true;
}

const Field3D GlobalField3D::scatter() const {
  Field3D result;
  result.allocate();
  
  MPI_Status status;
  if(mype == data_on_proc) {
    // Data is on this processor. Send to other processors
    
    for(int p = 0; p < npes; p++) {
      if(p == mype) continue;
      
      int xorig, yorig;
      proc_origin(p, &xorig, &yorig);
      int xsize, ysize;
      proc_size(p, &xsize, &ysize);
      int zsize = mesh->ngz-1;
      
      // Send to processor p
      for(int x=0;x<xsize;x++)
        for(int y=0;y<ysize;y++) 
          for(int z=0;z<zsize;z++) {
            buffer[p][x*ysize*zsize + y*zsize + z] = (*this)(x+xorig,y+yorig,z);
          }
      
      MPI_Send(buffer[p], xsize*ysize*zsize, MPI_DOUBLE, p, 1413, comm);
    }

    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    int zsize = mesh->ngz-1;
    
    // Copy to result
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++)
        for(int z=0;z<zsize;z++) {
          result(local_xorig+x,local_yorig+y,z) = (*this)(x+xorig,y+yorig,z);
        }
  }else {
    // Receive data
    MPI_Recv(buffer[0], msg_len(mype), MPI_DOUBLE, data_on_proc, 1413, comm, &status);
    
    int local_xorig, local_yorig;
    proc_local_origin(mype, &local_xorig, &local_yorig);
    int xorig, yorig;
    proc_origin(mype, &xorig, &yorig);
    int xsize, ysize;
    proc_size(mype, &xsize, &ysize);
    int zsize = mesh->ngz-1;
    
    for(int x=0;x<xsize;x++)
      for(int y=0;y<ysize;y++) 
        for(int z=0;z<zsize;z++) {
          result(local_xorig+x,local_yorig+y,z) = buffer[0][x*ysize*zsize + y*zsize + z];
        }
  }
  return result;
}

int GlobalField3D::msg_len(int proc) const {
  int xsize, ysize;
  proc_size(proc, &xsize, &ysize);
  int zsize = mesh->ngz-1;
  return xsize * ysize * zsize;
}
