#include "boutcomm.hxx"

BoutComm* BoutComm::instance = NULL;

BoutComm::BoutComm() : comm(MPI_COMM_WORLD), hasBeenSet(false) {
}

void BoutComm::setComm(MPI_Comm c) {
  hasBeenSet = true;
  this->comm = c;
}

MPI_Comm BoutComm::getComm() {
  return this->comm;
}

bool BoutComm::isSet() {
  return hasBeenSet;
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
