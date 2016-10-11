#include <mpi.h>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <iostream>
#include <stdarg.h>
#include <output.hxx>
#include <utils.hxx>

void BoutParallelThrowRhsFail(int &status, const char* message) {
  int allstatus;
  MPI_Allreduce(&status,&allstatus,1,MPI_INT,MPI_LOR,BoutComm::get());
  
  if (allstatus) throw BoutRhsFail(message);
}

BoutException::~BoutException() throw() {
  delete[] buffer;
}

BoutException::BoutException(const char* s, ...) {
  buflen=1024;
  buffer=new char[buflen];
  
  if(s == (const char*) NULL)
    return;

  myvsnprintf(buffer, buflen, s, ap);

  message.assign(buffer);
}

const char* BoutException::what() const throw() {
  #ifdef CHECK
    /// Print out the message stack to help debugging
    msg_stack.dump();
  #else
    output.write("Enable checking (-DCHECK flag) to get a trace\n");
  #endif
  return message.c_str();
}

BoutRhsFail::BoutRhsFail(const char* s, ...)  : BoutException::BoutException(s) {
}

BoutIterationFail::BoutIterationFail(const char* s, ...) : BoutException::BoutException(s) {
}
