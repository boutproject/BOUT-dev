#include "boutcomm.h"

BoutComm* BoutComm::instance = NULL;

BoutComm::BoutComm() : comm(MPI_COMM_WORLD) {
}

void BoutComm::setComm(MPI_Comm c) {
  this->comm = c;
}

MPI_Comm BoutComm::getComm() {
  return this->comm;
}

// Static functions below. Must use getInstance()
MPI_Comm BoutComm::get() {
  return getInstance()->getComm();
}

BoutComm* BoutComm::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new BoutComm();
  }
  return instance;
}
