
#include <bout/globalfield.hxx>


GlobalField::GlobalField(Mesh *m, int xsize, int ysize, int zsize) 
  : mesh(m), nx(xsize), ny(ysize), nz(zsize) {
  
  comm = BoutComm::get(); // This should come from Mesh
  
  MPI_Comm_size(comm, &npes);
  MPI_Comm_rank(comm, &mype);
  
  
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
        
        MPI_Irecv(d, len, MPI_DOUBLE, 3141, comm, &req[p]);
      }
    }
    
    // Copy data from this processor

    if(npes > 1) {
      // Wait for receives, process as they arrive
      do {
        int pe, status;
        MPI_Waitany(npes-1, req, &pe, status);
        
        // Unpack data from processor 'pe'
        
        

        if(pe != MPI_UNDEFINED)
          req[pe] = MPI_REQUEST_NULL;
      }while(pe != MPI_UNDEFINED);
    }
  }else {
    // Sending data to proc
    BoutReal *d = f.getData()[0]; // Pointer to start of data
    
    MPI_Send(d, len, MPI_DOUBLE, proc, 3141, comm);
  }
}

const Field2D GlobalField2D::scatter() const {
  Field2D result;
  
}

