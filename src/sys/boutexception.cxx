#include <bout/boutcomm.hxx>
#include <bout/boutexception.hxx>
#include <bout/msg_stack.hxx>
#include <bout/utils.hxx>

#include <mpi.h>

#include <cpptrace/cpptrace.hpp> // IWYU pragma: keep
#include <cpptrace/formatting.hpp>
#include <cpptrace/utils.hpp>

#include <cstdlib>
#include <string>
#include <utility>

#include <fmt/format.h>

namespace {
const std::string header{"====== Exception thrown ======\n"};
}

bool BoutException::show_backtrace = true;

void BoutParallelThrowRhsFail(int status, const char* message) {
  int allstatus = 0;
  MPI_Allreduce(&status, &allstatus, 1, MPI_INT, MPI_LOR, BoutComm::get());

  if (allstatus != 0) {
    throw BoutRhsFail(message);
  }
}

BoutException::BoutException(std::string msg) : message(std::move(msg)) {
  const char* show_backtrace_env_var = std::getenv("BOUT_SHOW_BACKTRACE");
  const auto should_show_backtrace =
      show_backtrace
      and ((show_backtrace_env_var == nullptr)
           or (show_backtrace_env_var != nullptr
               and std::string{show_backtrace_env_var} != "0"));
  if (should_show_backtrace) {
    message = getBacktrace();
  }
}

BoutException::~BoutException() {
  // If an exception is thrown while a TRACE is active, we won't clear
  // up the msg_stack. We also won't know how many messages to pop, so
  // just clear everything
  msg_stack.clear();
}

std::string BoutException::getBacktrace() const {
  using namespace cpptrace;

  const auto colours = isatty(stdout_fileno) || isatty(stderr_fileno)
                           ? formatter::color_mode::always
                           : formatter::color_mode::none;

  auto formatter = cpptrace::formatter{}
                       .addresses(formatter::address_mode::none)
                       .break_before_filename(true)
                       .colors(colours)
                       .snippets(true)
                       .symbols(formatter::symbol_mode::pretty)
                       .filter([](const stacktrace_frame& frame) {
                         return (
                             // Don't include our exception machinery
                             (frame.symbol.find("BoutException::") == std::string::npos)
                             // Don't include pre-main functions
                             and (frame.symbol.find("__libc_start") == std::string::npos)
                             and (frame.symbol != "_start"));
                       })
                       .filtered_frame_placeholders(false);

  const std::string backtrace_message = formatter.format(generate_trace());

  return fmt::format("{}\n{}\n{}{}\n", backtrace_message, msg_stack.getDump(), header,
                     message);
}
