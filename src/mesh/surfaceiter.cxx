
#include <bout/surfaceiter.hxx>

#include <boutexception.hxx>
#include <unused.hxx>

int SurfaceIter::ySize() {
  return m->ySize(xpos);
}

bool SurfaceIter::closed() {
  return m->periodicY(xpos);
}

bool SurfaceIter::closed(BoutReal &ts) {
  return m->periodicY(xpos, ts);
}

MPI_Comm SurfaceIter::communicator() {
  return m->getYcomm(xpos);
}

int SurfaceIter::yGlobal(int UNUSED(yloc)) {
  // Communicator for this surface
  MPI_Comm comm = communicator();
  
  // Get number of processors and processor rank
  int np;
  bout::globals::mpi->MPI_Comm_size(comm, &np);
  int myp;
  bout::globals::mpi->MPI_Comm_rank(comm, &myp);
  
  // Need a generic method

  throw BoutException("SurfaceIter::yGlobal not implemented");
  
  return 0;
}

bool SurfaceIter::firstY() { ///< Is this processor at the lower end?
  if(closed())
    return false;
  
  // Communicator for this surface
  MPI_Comm comm = communicator();
  
  // Get processor rank
  int myp;
  bout::globals::mpi->MPI_Comm_rank(comm, &myp);
  
  return myp == 0;
}

bool SurfaceIter::lastY() {
  if(closed())
    return false;
  
  // Communicator for this surface
  MPI_Comm comm = communicator();
  
  // Get number of processors and processor rank
  int np;
  bout::globals::mpi->MPI_Comm_size(comm, &np);
  int myp;
  bout::globals::mpi->MPI_Comm_rank(comm, &myp);
  
  return myp == np-1;
}

void SurfaceIter::first() {
  xpos = 0;
}

void SurfaceIter::next() {
  if(xpos < 0)
    return;
  
  xpos++;
  if(xpos >= m->LocalNx)
    xpos = -1;
}

bool SurfaceIter::isDone() {
  return xpos < 0;
}
