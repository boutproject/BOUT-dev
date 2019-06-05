#include <boutcomm.hxx>
#include <bout_types.hxx>

BoutComm* BoutComm::instance = nullptr;

BoutComm::BoutComm() : comm(MPI_COMM_NULL) {}

BoutComm::~BoutComm() {
  if(comm != MPI_COMM_NULL)
    MPI_Comm_free(&comm);
  
  if(!isSet()) {
    // If BoutComm was set, then assume that MPI_Finalize is called elsewhere
    // but might need to revisit if that isn't the case
    MPI_Finalize();
  }
}

void BoutComm::setComm(MPI_Comm c) {
  if(comm != MPI_COMM_NULL)
    MPI_Comm_free(&comm);
  MPI_Comm_dup(c, &comm);
  hasBeenSet = true;
}

MPI_Comm BoutComm::getComm() {
  if(comm == MPI_COMM_NULL) {
    // No communicator set. Initialise MPI
    MPI_Init(pargc,pargv);
    
    // Duplicate MPI_COMM_WORLD
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
  }
  return comm;
}

bool BoutComm::isSet() {
  return hasBeenSet;
}

// Static functions below. Must use getInstance()
MPI_Comm BoutComm::get() {
  return getInstance()->getComm();
}

void BoutComm::setArgs(int &c, char**&v) {
  getInstance()->pargc = &c;
  getInstance()->pargv = &v;
}

int BoutComm::rank() {
  int MYPE;
  MPI_Comm_rank(get(), &MYPE);
  return MYPE;
}

int BoutComm::size() {
  int NPES;
  MPI_Comm_size(get(), &NPES);
  return NPES;
}

BoutComm* BoutComm::getInstance() {
  if(instance == nullptr) {
    // Create the singleton object
    instance = new BoutComm();
  }
  return instance;
}

void BoutComm::cleanup() {
  delete instance;
  instance = nullptr;
}
