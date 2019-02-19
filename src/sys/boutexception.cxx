#include <mpi.h>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <iostream>
#include <msg_stack.hxx>
#include <output.hxx>
#include <stdarg.h>

#ifdef BACKTRACE
#include <execinfo.h>
#include <dlfcn.h>
#endif

#include <utils.hxx>

void BoutParallelThrowRhsFail(int status, const char *message) {
  int allstatus;
  MPI_Allreduce(&status, &allstatus, 1, MPI_INT, MPI_LOR, BoutComm::get());

  if (allstatus) {
    throw BoutRhsFail("%s",message);
  }
}

BoutException::~BoutException() {
  if (buffer != nullptr) {
    delete[] buffer;
    buffer = nullptr;
  }
}

void BoutException::Backtrace() {
#if CHECK > 1
  /// Print out the message stack to help debugging
  std::string tmp = msg_stack.getDump();
  message += tmp;
#else
  message += "Enable checking (configure with --enable-check or set flag -DCHECK > 1) to "
             "get a trace\n";
#endif
  
#ifdef BACKTRACE

  trace_size = backtrace(trace, TRACE_MAX);
  messages = backtrace_symbols(trace, trace_size);

#else // BACKTRACE
  message += "Stacktrace not enabled.\n";
#endif
}

std::string BoutException::BacktraceGenerate() const{
  std::string message;
#ifdef BACKTRACE
    // skip first stack frame (points here)
  message = ("====== Exception path ======\n");
  char buf[1024];
  for (int i = 1; i < trace_size; ++i) {
    snprintf(buf, sizeof(buf) - 1, "[bt] #%d %s\n", i, messages[i]);
    message += buf;
    // find first occurence of '(' or ' ' in message[i] and assume
    // everything before that is the file name. (Don't go beyond 0 though
    // (string terminator)
    int p = 0; // snprintf %.*s expects int
    while (messages[i][p] != '(' && messages[i][p] != ' ' && messages[i][p] != 0) {
      ++p;
    }

    char syscom[256];
    // If we are compiled as PIE, need to get base pointer of .so and substract
    Dl_info info;
    void * ptr=trace[i];
    if (dladdr(trace[i],&info)){
      // Additionally, check whether this is the default offset for an executable
      if (info.dli_fbase != (void*)0x400000)
        ptr=(void*) ((size_t)trace[i]-(size_t)info.dli_fbase);
    }

    // Pipe stderr to /dev/null to avoid cluttering output
    // when addr2line fails or is not installed
    snprintf(syscom, sizeof(syscom) - 1, "addr2line %p -Cfpie %.*s 2> /dev/null",
             ptr, p, messages[i]);
    // last parameter is the file name of the symbol
    FILE *fp = popen(syscom, "r");
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
        message += buf;
      }
    }
  }
#endif
  return message;
}

/// Common set up for exceptions
///
/// Formats the message s using C-style printf formatting
#define INIT_EXCEPTION(s)                                                                \
  {                                                                                      \
    buflen = 0;                                                                          \
    buffer = nullptr;                                                                    \
    if (s == nullptr) {                                                                  \
      message = "No error message given!\n";                                             \
    } else {                                                                             \
      buflen = BoutException::BUFFER_LEN;                                                \
      buffer = new char[buflen];                                                         \
      bout_vsnprintf(buffer, buflen, s);                                                 \
      for (int i = 0; i < buflen; ++i) {                                                 \
        if (buffer[i] == 0) {                                                            \
          if (i > 0 && buffer[i - 1] == '\n') {                                          \
            buffer[i - 1] = 0;                                                           \
          }                                                                              \
          break;                                                                         \
        }                                                                                \
      }                                                                                  \
      message.assign(buffer);                                                            \
      delete[] buffer;                                                                   \
      buffer = nullptr;                                                                  \
    }                                                                                    \
    message = "====== Exception thrown ======\n" + message + "\n";                       \
                                                                                         \
    this->Backtrace();                                                                   \
  }

BoutException::BoutException(const char *s, ...) { INIT_EXCEPTION(s); }

BoutException::BoutException(const std::string &msg) {
  message = "====== Exception thrown ======\n" + msg + "\n";

  this->Backtrace();
}

const char *BoutException::what() const noexcept{
#ifdef BACKTRACE
  _tmp=message;
  _tmp+=BacktraceGenerate();
  return _tmp.c_str();
#else
  return message.c_str();
#endif
}

BoutRhsFail::BoutRhsFail(const char *s, ...) : BoutException::BoutException(nullptr) {
  INIT_EXCEPTION(s);
}

BoutIterationFail::BoutIterationFail(const char *s, ...)
    : BoutException::BoutException(nullptr) {
  INIT_EXCEPTION(s);
}
