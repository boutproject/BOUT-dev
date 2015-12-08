#include <mpi.h>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <iostream>
#include <stdarg.h>
#include <output.hxx>
#define BACKTRACE
#ifdef BACKTRACE
#include <execinfo.h>
#endif

void BoutParallelThrowRhsFail(int &status, const char* message) {
  int allstatus;
  MPI_Allreduce(&status,&allstatus,1,MPI_INT,MPI_LOR,BoutComm::get());
  
  if (allstatus) throw BoutRhsFail(message);
}

BoutException::~BoutException() throw()
{
}

void BoutException::Backtrace(){
#ifdef BACKTRACE
  void *trace[64];
  char **messages = (char **)NULL;
  int i, trace_size = 0;
  
  trace_size = backtrace(trace, 64);
  messages = backtrace_symbols(trace, trace_size);
  /* skip first stack frame (points here) */
  output.write("\n[bt] Execution path:\n");
  for (i=1; i<trace_size; ++i)
    {
      output.write("[bt] #%d %s\n", i, messages[i]);
      
      /* find first occurence of '(' or ' ' in message[i] and assume
       * everything before that is the file name. (Don't go beyond 0 though
       * (string terminator)*/
      size_t p = 0;
      while(messages[i][p] != '(' && messages[i][p] != ' '
	    && messages[i][p] != 0)
	++p;
      
      char syscom[256];
      sprintf(syscom,"addr2line %p -Cfpie %.*s", trace[i], p, messages[i]);
      //last parameter is the file name of the symbol
      FILE * fp = popen(syscom, "r");
      if (fp!= NULL){
	char out[1024];
	fgets(out, sizeof(out)-1, fp);
	if (pclose(fp) == 0){
	  output.write(out);
	}
      } else {
	output.write(syscom);
      }
    }
#else
  output.write("Backtrace not enabled!\n");
#endif
}

BoutException::BoutException(const char* s, ...)
{
  this->Backtrace();
  
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
