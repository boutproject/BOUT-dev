#include <boutcomm.hxx>
#include <bout_types.hxx>

BoutComm* BoutComm::instance = nullptr;

BoutComm::BoutComm()
    : pargc(nullptr), pargv(nullptr), hasBeenSet(false), comm(MPI_COMM_NULL) {}

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
  if(instance != nullptr) delete instance;
  instance = nullptr;
}

// Communicator interface

Communicator::Request BoutComm::irecvBytes(void *data, int length, int destination, int identifier) {
  MPI_Request request;
  
  MPI_Irecv(data,
            length,
            MPI_BYTE, // Just sending raw data, unknown type
            destination, // Destination processor
            identifier,
            get(),    // Communicator
            &request);

  return request;
}

void BoutComm::sendBytes(void *data, int length, int destination, int identifier) {
  MPI_Send(data,     // Data pointer
           length,   // Number
           MPI_BYTE, // Type
           destination, // Destination processor
           identifier,
           get());   // Communicator
}

int BoutComm::waitAny(const std::vector<Communicator::Request> &requests) {
  std::vector<MPI_Request> mpi_requests;
  std::vector<int> indices;

  for (std::vector<Communicator::Request>::size_type i = 0; i < requests.size(); i++) {
    const auto& req = requests[i];
    
    if (req != Communicator::REQUEST_NULL) {
      mpi_requests.push_back(bout::utils::get<MPI_Request>(req));
      indices.push_back(i); // Keep track of original index
    }
  }

  if (mpi_requests.empty()) {
    return Communicator::UNDEFINED;
  }
  
  int ind;
  MPI_Status stat;
  MPI_Waitany(mpi_requests.size(), mpi_requests.data(), &ind, &stat);

  if (ind == MPI_UNDEFINED) {
    return Communicator::UNDEFINED;
  }

  // Look up original index
  return indices[ind];
}
