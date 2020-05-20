#include "bout/sys/timer.hxx"

#include <fmt/core.h>

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

void Timer::printTimeReport() {
  using namespace std::string_literals;
  const auto header_name = "Timer name"s;
  const auto header_time = "Total time (s)"s;
  const auto header_frac = "% of top"s;
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
  // Only print percentages to two d.p.
  const auto frac_width = std::max(header_frac.length(), "100.00%"s.length());
  const auto hit_width = longest_width(header_hits, [](const auto& it) {
    return fmt::format("{}", it.second.hits).length();
  });
  // Assume that the mean will be about the same length as the total
  const auto mean_time_width = longest_width(header_mean, [](const auto& it) {
    return fmt::format("{:.9g}", it.second.total_time.count()).length();
  });

  output.write("{0:<{1}} | {2:<{3}} | {4:<{5}} | {6:<{7}} | {8:<{9}}\n", header_name,
               max_width, header_time, total_time_width, header_frac, frac_width,
               header_hits, hit_width, header_mean, mean_time_width);
  output.write("{0:-<{1}}-|-{0:-<{2}}-|-{0:-<{3}}-|-{0:-<{4}}-|-{0:-<{5}}\n", ""s,
               max_width, total_time_width, frac_width, hit_width, mean_time_width);

  // To sort the map by the total time, we need to copy the key-value
  // pairs into another container that we can sort ourselves
  std::vector<std::pair<std::string, timer_info>> sorted_info;
  sorted_info.reserve(info.size());
  std::transform(begin(info), end(info), std::back_inserter(sorted_info),
                 [](const auto& it) { return it; });
  // Sort so that the largest total time is first
  std::sort(begin(sorted_info), end(sorted_info), [](const auto& a, const auto& b) {
    return a.second.total_time > b.second.total_time;
  });

  for (const auto& kv : sorted_info) {
    output.write("{0:<{1}} | {2:<{3}.9g} | {8:<{9}.2%} | {4:<{5}} | {6:<{7}.9g}\n",
                 kv.first, max_width, kv.second.total_time.count(), total_time_width,
                 kv.second.hits, hit_width, kv.second.total_time.count() / kv.second.hits,
                 mean_time_width,
                 kv.second.total_time.count() / sorted_info[0].second.total_time.count(),
                 frac_width);
  }
}
