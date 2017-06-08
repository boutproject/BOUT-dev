#include <mpi.h>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <iostream>
#include <stdarg.h>
#include <output.hxx>
#ifdef BACKTRACE
#include <execinfo.h>
#endif
#include <utils.hxx>

void BoutParallelThrowRhsFail(int &status, const std::string &message) {
  int allstatus;
  MPI_Allreduce(&status, &allstatus, 1, MPI_INT, MPI_LOR, BoutComm::get());

  if (allstatus) {
    throw BoutRhsFail(message);
  }
}

void BoutException::Backtrace() {
#if CHECK > 1
  /// Print out the message stack to help debugging
  std::string tmp=msg_stack.getDump();
  message+=tmp;
#else
  message+="Enable checking (configure with --enable-check or set flag -DCHECK > 1) to get a trace\n";
#endif
#ifdef BACKTRACE
  void *trace[64];
  char **messages = (char **)NULL;
  int i, trace_size = 0;

  trace_size = backtrace(trace, 64);
  messages = backtrace_symbols(trace, trace_size);
  /* skip first stack frame (points here) */
  //output.write("\n[bt] Execution path:\n");
  message+=("====== Exception path ======\n");
  char buf[1024];
  for (i=1; i<trace_size; ++i)
    {
      snprintf(buf,1023,"[bt] #%d %s\n", i, messages[i]);
      message+=buf;
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
          message+=out;
	}
      } else {
        message+=syscom;
      }
    }
#else
  message+="Stacktrace not enabled.\n";
#endif
}

BoutException::BoutException(const std::string &msg) {
  message = "====== Exception thrown ======\n" + msg + "\n";

  this->Backtrace();
}
