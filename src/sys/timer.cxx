#include "bout/sys/timer.hxx"

#include <fmt/format.h>

#include <algorithm>

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
  const std::string headerName = "Timer name";
  const std::string headerTime = "Current time (s)";
  const std::string headerTotal = "Total time (s)";
  const std::string headerHits = "Hits";

  const auto longest_width = [](const std::string& header, const auto& op) {
    return std::max(header.size(),
                    op(*std::max_element(std::begin(info), std::end(info),
                                         [&op](const auto& lhs, const auto& rhs) {
                                           return op(lhs) < op(rhs);
                                         })));
  };

  const auto max_width =
      longest_width(headerName, [](const auto& it) { return it.first.length(); });
  const auto time_width = longest_width(headerTime, [](const auto& it) {
    return fmt::format("{:.9g}", it.second.time.count()).length();
  });
  const auto total_time_width = longest_width(headerTotal, [](const auto& it) {
    return fmt::format("{:.9g}", it.second.total_time.count()).length();
  });
  const auto hit_width = longest_width(headerHits, [](const auto& it) {
    return fmt::format("{}", it.second.hits).length();
  });

  using namespace std::string_literals;

  output.write("{0:<{1}} | {2:<{3}} | {5:<{6}} | {4}\n", headerName, max_width,
               headerTime, time_width, headerHits, headerTotal, total_time_width);
  output.write("{0:-<{1}}-|-{0:-<{2}}-|-{0:-<{4}}-|-{0:-<{3}}\n", ""s, max_width,
               time_width, hit_width, total_time_width);
  for (const auto& kv : info) {
    output.write("{0:<{1}} | {2:<{3}.9g} | {6:<{7}.9g} | {4:<{5}}\n", kv.first, max_width,
                 kv.second.time.count(), time_width, kv.second.hits, hit_width,
                 kv.second.total_time.count(), total_time_width);
  }
}
