#include "bout/sys/timer.hxx"

#include "fmt/format.h"

#include <algorithm>

Timer::Timer() : timing(getInfo("")) {
  if (timing.counter == 0) {
    timing.started = clock_type::now();
    timing.running = true;
    timing.ntimes++;
  }
  timing.counter += 1;
}

Timer::Timer(const std::string& label) : timing(getInfo(label)) {
  if (timing.counter == 0) {
    timing.started = clock_type::now();
    timing.running = true;
    timing.ntimes++;
  }
  timing.counter += 1;
}

Timer::~Timer() {
  timing.counter -= 1;
  if (timing.counter == 0) {
    auto finished = clock_type::now();
    timing.running = false;
    timing.time += finished - timing.started;
  }
}

void Timer::cleanup() { info.clear(); }

std::map<std::string, Timer::timer_info> Timer::info;

Timer::timer_info& Timer::getInfo(const std::string& label) {
  auto it = info.find(label);
  if (it == info.end()) {
    auto timer = info.emplace(
                              label, timer_info{seconds{0}, false, clock_type::now(), 0, 0});
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

double Timer::resetTime(Timer::timer_info& info) {
  auto val = info.time;
  info.time = clock_type::duration{0};
  if (info.running) {
    auto cur_time = clock_type::now();
    val += cur_time - info.started;
    info.started = cur_time;
  }
  return seconds{val}.count();
}

void Timer::listAllInfo() {
  const std::string headerName = "Timer name";
  const std::string headerTime = "Time (s)";
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
  const auto hit_width = longest_width(headerHits, [](const auto& it) {
    return fmt::format("{}", it.second.ntimes).length();
  });

  using namespace std::string_literals;

  output.write("\nTimer report \n\n");
  output.write("{0:<{1}} | {2:<{3}} | {4}\n", headerName, max_width, headerTime,
               time_width, headerHits);
  output.write("{0:-<{1}}-|-{0:-<{2}}-|-{0:-<{3}}\n", ""s, max_width, time_width,
               hit_width);
  for (const auto& kv : info) {
    output.write("{0:<{1}} | {2:<{3}.9g} | {4:<{5}}\n", kv.first, max_width,
                 kv.second.time.count(), time_width, kv.second.ntimes, hit_width);
  }
  output.write("\n");
}
