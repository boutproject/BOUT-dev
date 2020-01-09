#include "bout/sys/timer.hxx"

Timer::Timer() : timing(getInfo("")) {
  if (timing.counter == 0) {
    timing.started = clock_type::now();
    timing.running = true;
    ++timing.hits;
  }
  timing.counter += 1;
}

Timer::Timer(const std::string& label) : timing(getInfo(label)) {
  if (timing.counter == 0) {
    timing.started = clock_type::now();
    timing.running = true;
    ++timing.hits;
  }
  timing.counter += 1;
}

Timer::~Timer() {
  timing.counter -= 1;
  if (timing.counter == 0) {
    const auto elapsed = clock_type::now() - timing.started;
    timing.running = false;
    timing.time += elapsed;
    timing.total_time += elapsed;
  }
}

void Timer::cleanup() { info.clear(); }

std::map<std::string, Timer::timer_info> Timer::info;

Timer::timer_info& Timer::getInfo(const std::string& label) {
  auto it = info.find(label);
  if (it == info.end()) {
    auto timer = info.emplace(
        label, timer_info{seconds{0}, seconds{0}, false, clock_type::now(), 0, 0});
    return timer.first->second;
  }
  return it->second;
}

double Timer::getTime(const Timer::timer_info& info) {
  if (info.running) {
    return seconds{info.time + (clock_type::now() - info.started)}.count();
  }
  return seconds{info.time}.count();
}

double Timer::getTotalTime(const Timer::timer_info& info) {
  if (info.running) {
    return seconds{info.total_time + (clock_type::now() - info.started)}.count();
  }
  return seconds{info.total_time}.count();
}

double Timer::resetTime(Timer::timer_info& info) {
  auto current_duration = info.time;
  info.time = clock_type::duration{0};
  if (info.running) {
    const auto current_time = clock_type::now();
    const auto elapsed = current_time - info.started;
    current_duration += elapsed;
    info.started = current_time;
    info.total_time += elapsed;
  }
  return seconds{current_duration}.count();
}

void Timer::listAllInfo() {
  const std::string headerOne = "Timer name";
  const std::string headerTime = "Current time (s)";
  const std::string headerHits = "Hits";
  const std::string headerTotal = "Total time (s)";
  const std::string separator = " | ";

  auto max_width = static_cast<unsigned int>(headerTime.length());

  output << "Timer report \n\n";
  output << std::setw(max_width) << headerOne << separator << headerTime
         << separator << headerHits << separator << headerTotal << "\n";
  output << std::setw(max_width) << std::string(max_width, '-') << separator
          << std::string(max_width, '-') << separator
          << std::string(max_width, '-') << separator
          << std::string(max_width, '-') << "\n";
  for (const auto& kv : info) {
    output << std::left << std::setw(max_width) << kv.first << separator
           << kv.second.time.count() << separator << << kv.second.hits 
           << separator << kv.second.total_time.count() << ")\n";
  }
  output << "\n";
}
