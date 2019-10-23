#include <mpi.h>
#include <bout/sys/timer.hxx>

using namespace std;

Timer::Timer() {
  timing = getInfo("");
  if (timing->counter == 0) {
    timing->started = MPI_Wtime();
    timing->running = true;
  }
  timing->counter++;
}

Timer::Timer(const std::string &label) {
  timing = getInfo(label);
  if (timing->counter == 0) {
    timing->started = MPI_Wtime();
    timing->running = true;
  }
  timing->counter++;
}

Timer::~Timer() {
  timing->counter--;
  if (timing->counter == 0) {
    double finished = MPI_Wtime();
    timing->running = false;
    timing->time += finished - timing->started;
  }
}

double Timer::getTime() {
  if(timing->running)
    return timing->time + (MPI_Wtime() - timing->started);
  return timing->time;
}

double Timer::resetTime() {
  double val = timing->time;
  timing->time = 0.0;
  if(timing->running) {
    double cur_time = MPI_Wtime();
    val += cur_time - timing->started;
    timing->started = cur_time;
  }
  return val;
}

double Timer::getTime(const std::string &label) {
  timer_info* t = getInfo(label);
  if(t->running)
    return t->time + (MPI_Wtime() - t->started);
  return t->time;
}

double Timer::resetTime(const std::string &label) {
  timer_info* t = getInfo(label);
  double val = t->time;
  t->time = 0.0;
  if(t->running) {
    double cur_time = MPI_Wtime();
    val += cur_time - t->started;
    t->started = cur_time;
  }
  return val;
}

// Static method to clean up all memory
void Timer::cleanup() {
  // Iterate over map
  for(const auto& it : info) {
    delete it.second;
  }
  info.clear();
}

map<std::string, Timer::timer_info*> Timer::info;

Timer::timer_info* Timer::getInfo(const std::string &label) {
  map<string, timer_info*>::iterator it(info.find(label));
  if(it == info.end()) {
    // Not in map, so create it
    timer_info *t = new timer_info;
    t->time = 0.0;
    t->running = false;
    t->counter = 0;
    info[label] = t;
    return t;
  }
  return it->second;
}
