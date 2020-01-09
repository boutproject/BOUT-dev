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
  using namespace std::string_literals;
  const auto header_name = "Timer name"s;
  const auto header_time = "Total time (s)"s;
  const auto header_hits = "Hits"s;
  const auto header_mean = "Mean time/hit (s)"s;

  // Helper lambda to get the longest item in info, after applying
  // another lambda (op), and also taking into account the header. op
  // transforms from a map key-value pair to the length of either the
  // key or a member of the value (Timer::timing_info). As
  // std::max_element returns an iterator, it needs to be dereferenced
  // and have op applied to it
  const auto longest_width = [](const std::string& header, const auto& op) {
    return std::max(header.length(),
                    op(*std::max_element(std::begin(info), std::end(info),
                                         [&op](const auto& lhs, const auto& rhs) {
                                           return op(lhs) < op(rhs);
                                         })));
  };

  // Get the widths of the four columns: timer name, total time,
  // number of hits, mean time per hit
  const auto max_width =
      longest_width(header_name, [](const auto& it) { return it.first.length(); });
  const auto total_time_width = longest_width(header_time, [](const auto& it) {
    return fmt::format("{:.9g}", it.second.total_time.count()).length();
  });
  const auto hit_width = longest_width(header_hits, [](const auto& it) {
    return fmt::format("{}", it.second.hits).length();
  });
  // Assume that the mean will be about the same length as the total
  const auto mean_time_width = longest_width(header_mean, [](const auto& it) {
    return fmt::format("{:.9g}", it.second.total_time.count()).length();
  });

  output.write("{0:<{1}} | {2:<{3}} | {4:<{5}} | {6:<{7}}\n", header_name, max_width,
               header_time, total_time_width, header_hits, hit_width, header_mean,
               mean_time_width);
  output.write("{0:-<{1}}-|-{0:-<{2}}-|-{0:-<{3}}-|-{0:-<{4}}\n", ""s, max_width,
               total_time_width, hit_width, mean_time_width);
  for (const auto& kv : info) {
    output.write("{0:<{1}} | {2:<{3}.9g} | {4:<{5}} | {6:<{7}.9g}\n", kv.first, max_width,
                 kv.second.total_time.count(), total_time_width, kv.second.hits,
                 hit_width, kv.second.total_time.count() / kv.second.hits,
                 mean_time_width);
  }
}
