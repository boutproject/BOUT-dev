#include <mpi.h>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <iostream>
#include <stdarg.h>
#include <output.hxx>

void BoutParallelThrowRhsFail(int &status, const char* message) {
  int allstatus;
  MPI_Allreduce(&status,&allstatus,1,MPI_INT,MPI_LOR,BoutComm::get());
  
  if (allstatus) throw BoutRhsFail(message);
}

BoutException::~BoutException() throw()
{
}

BoutException::BoutException(const char* s, ...)
{
  va_list ap;  // List of arguments

  if(s == (const char*) NULL)
    return;
  
  char buffer[1024];
  va_start(ap, s);
    vsprintf(buffer, s, ap);
  va_end(ap);
  
  message.assign(buffer);
}

const char* BoutException::what() const throw()
{
  #ifdef CHECK
    /// Print out the message stack to help debugging
    msg_stack.dump();
  #else
    output.write("Enable checking (-DCHECK flag) to get a trace\n");
  #endif
  return message.c_str();
}

BoutRhsFail::BoutRhsFail(const char* s, ...)  : BoutException::BoutException(s){
  va_list ap;  // List of arguments

  if(s == (const char*) NULL)
    return;
  
  char buffer[1024];
  va_start(ap, s);
    vsprintf(buffer, s, ap);
  va_end(ap);
  
  message.assign(buffer);
}

BoutIterationFail::BoutIterationFail(const char* s, ...) : BoutException::BoutException(s) {
  va_list ap;  // List of arguments

  if(s == (const char*) NULL)
    return;
  
  char buffer[1024];
  va_start(ap, s);
    vsprintf(buffer, s, ap);
  va_end(ap);
  
  message.assign(buffer);
}
