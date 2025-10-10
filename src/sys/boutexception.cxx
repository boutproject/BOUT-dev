#include "bout/build_defines.hxx"

#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/msg_stack.hxx>
#include <bout/utils.hxx>
#include <mpi.h>

#if BOUT_USE_BACKTRACE
#include <dlfcn.h>
#include <execinfo.h>
#endif

#include <array>
#include <cstdlib>
#include <string>

#include <fmt/format.h>

namespace {
const std::string header{"====== Exception thrown ======\n"};
}

void BoutParallelThrowRhsFail(int status, const char* message) {
  int allstatus;
  MPI_Allreduce(&status, &allstatus, 1, MPI_INT, MPI_LOR, BoutComm::get());

  if (allstatus) {
    throw BoutRhsFail(message);
  }
}

BoutException::BoutException(std::string msg) : message(std::move(msg)) {
  makeBacktrace();
  if (std::getenv("BOUT_SHOW_BACKTRACE") != nullptr) {
    message = getBacktrace() + "\n" + message;
  }
}

BoutException::~BoutException() {
  // If an exception is thrown while a TRACE is active, we won't clear
  // up the msg_stack. We also won't know how many messages to pop, so
  // just clear everything
  msg_stack.clear();
#if BOUT_USE_BACKTRACE
  // Call required for memory allocated by `backtrace_symbols`
  free(messages); // NOLINT
#endif
}

std::string BoutException::getBacktrace() const {
  std::string backtrace_message;
#if BOUT_USE_BACKTRACE
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
    void* ptr = trace[i];
    if (dladdr(trace[i], &info)) {
      // Additionally, check whether this is the default offset for an executable
      if (info.dli_fbase != reinterpret_cast<void*>(0x400000)) {
        ptr = reinterpret_cast<void*>(reinterpret_cast<size_t>(trace[i])
                                      - reinterpret_cast<size_t>(info.dli_fbase));
      }
    }

    // Pipe stderr to /dev/null to avoid cluttering output
    // when addr2line fails or is not installed
    const auto syscom = fmt::format(
        FMT_STRING("addr2line {:p} -Cfpie {:.{}s} 2> /dev/null"), ptr, messages[i], p);
    // last parameter is the file name of the symbol
    FILE* file = popen(syscom.c_str(), "r");
    if (file != nullptr) {
      std::array<char, 1024> out{};
      char* retstr = nullptr;
      std::string buf;
      while ((retstr = fgets(out.data(), out.size() - 1, file)) != nullptr) {
        buf += retstr;
      }
      int const status = pclose(file);
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
#if BOUT_USE_BACKTRACE
  trace_size = backtrace(trace.data(), TRACE_MAX);
  messages = backtrace_symbols(trace.data(), trace_size);
#endif
}
