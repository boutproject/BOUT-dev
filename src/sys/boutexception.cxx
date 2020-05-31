#include <mpi.h>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <iostream>
#include <msg_stack.hxx>
#include <output.hxx>

#ifdef BACKTRACE
#include <execinfo.h>
#include <dlfcn.h>
#endif

#include <utils.hxx>

#include <fmt/format.h>

void BoutParallelThrowRhsFail(int status, const char *message) {
  int allstatus;
  MPI_Allreduce(&status, &allstatus, 1, MPI_INT, MPI_LOR, BoutComm::get());

  if (allstatus) {
    throw BoutRhsFail(message);
  }
}

BoutException::~BoutException() {
  // If an exception is thrown while a TRACE is active, we won't clear
  // up the msg_stack. We also won't know how many messages to pop, so
  // just clear everything
  msg_stack.clear();
#ifdef BACKTRACE
  free(messages);
#endif
}

std::string BoutException::getBacktrace() const {
  std::string backtrace_message;
#ifdef BACKTRACE
  backtrace_message = "====== Exception path ======\n";
  // skip first stack frame (points here)
  for (int i = trace_size - 1; i > 1; --i) {
    backtrace_message += fmt::format(FMT_STRING("[bt] #{:d} {:s}\n"), i - 1, messages[i]);
    // find first occurence of '(' or ' ' in message[i] and assume
    // everything before that is the file name. (Don't go beyond 0 though
    // (string terminator)
    int p = 0; // snprintf %.*s expects int
    while (messages[i][p] != '(' && messages[i][p] != ' ' && messages[i][p] != 0) {
      ++p;
    }

    // If we are compiled as PIE, need to get base pointer of .so and substract
    Dl_info info;
    void * ptr=trace[i];
    if (dladdr(trace[i],&info)){
      // Additionally, check whether this is the default offset for an executable
      if (info.dli_fbase != reinterpret_cast<void*>(0x400000))
        ptr=reinterpret_cast<void*>(reinterpret_cast<size_t>(trace[i])-reinterpret_cast<size_t>(info.dli_fbase)));
    }

    // Pipe stderr to /dev/null to avoid cluttering output
    // when addr2line fails or is not installed
    const auto syscom =
      fmt::format(FMT_STRING("addr2line {:p} -Cfpie {:.{}s} 2> /dev/null"), ptr, messages[i], p);
    // last parameter is the file name of the symbol
    FILE *fp = popen(syscom.c_str(), "r");
    if (fp != nullptr) {
      char out[1024];
      char *retstr;
      std::string buf;
      do {
        retstr = fgets(out, sizeof(out) - 1, fp);
        if (retstr != nullptr)
          buf+=retstr;
      } while (retstr != nullptr);
      int status = pclose(fp);
      if (status == 0) {
        backtrace_message += buf;
      }
    }
  }
#else
  backtrace_message = "Stacktrace not enabled.\n";
#endif

  return backtrace_message + msg_stack.getDump() + "\n" + header + message + "\n";
}

void BoutException::makeBacktrace() {
#ifdef BACKTRACE
  trace_size = backtrace(trace, TRACE_MAX);
  messages = backtrace_symbols(trace, trace_size);
#endif
}
