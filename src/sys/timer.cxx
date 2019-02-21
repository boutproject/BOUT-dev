#include <bout/sys/timer.hxx>
#include <mpi.h>

Timer::Timer() : timing(getInfo("")) {
  timing.started = MPI_Wtime();
  timing.running = true;
}

Timer::Timer(const std::string& label) : timing(getInfo(label)) {
  timing.started = MPI_Wtime();
  timing.running = true;
}

Timer::~Timer() {
  double finished = MPI_Wtime();
  timing.running = false;
  timing.time += finished - timing.started;
}

// Static method to clean up all memory
void Timer::cleanup() { info.clear(); }

std::map<std::string, Timer::timer_info> Timer::info;

Timer::timer_info& Timer::getInfo(const std::string& label) {
  auto it = info.find(label);
  if (it == info.end()) {
    // Not in map, so create it
    auto timer = info.emplace(label, timer_info{0.0, false, 0.0});
    // timer is a pair of an iterator and bool
    // The iterator is a pair of key, value
    return timer.first->second;
  }
  return it->second;
}

double Timer::getTime(const Timer::timer_info& info) {
  if (info.running) {
    return info.time + (MPI_Wtime() - info.started);
  }
  return info.time;
}

double Timer::resetTime(Timer::timer_info& info) {
  double val = info.time;
  info.time = 0.0;
  if (info.running) {
    double cur_time = MPI_Wtime();
    val += cur_time - info.started;
    info.started = cur_time;
  }
  return val;
}
