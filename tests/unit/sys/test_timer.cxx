#include "gtest/gtest.h"

#include "bout/sys/timer.hxx"

#include <chrono>
#include <thread>

namespace bout {
namespace testing {
using ms = std::chrono::duration<double, std::chrono::milliseconds::period>;
using seconds = std::chrono::duration<double, std::chrono::seconds::period>;
constexpr double TimerTolerance{1e-3};
constexpr auto sleep_length = ms(0.5);
} // namespace testing
} // namespace bout

TEST(TimerTest, GetTime) {
  using namespace bout::testing;

  auto start = std::chrono::high_resolution_clock::now();

  Timer timer{};

  std::this_thread::sleep_for(sleep_length);

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(timer.getTime(), elapsed.count(), TimerTolerance);
}

TEST(TimerTest, GetUnknownTimeLabel) {
  EXPECT_DOUBLE_EQ(Timer::getTime("Does not exist"), 0.0);
}

TEST(TimerTest, GetTimeLabelInScope) {
  using namespace bout::testing;

  auto start = std::chrono::high_resolution_clock::now();

  Timer timer{"GetTimeLabelInScope test"};

  std::this_thread::sleep_for(sleep_length);

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelInScope test"), elapsed.count(),
              TimerTolerance);
}

TEST(TimerTest, GetTimeLabelOutOfScope) {
  using namespace bout::testing;

  auto start = std::chrono::high_resolution_clock::now();

  {
    Timer timer{"GetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(sleep_length);
  }

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;
  
  EXPECT_NEAR(Timer::getTime("GetTimeLabelOutOfScope test"), elapsed.count(),
              TimerTolerance);
}

TEST(TimerTest, GetTimeLabelRepeat) {
  using namespace bout::testing;

  auto start = std::chrono::high_resolution_clock::now();

  {
    Timer timer{"GetTimeLabelRepeat test"};

    std::this_thread::sleep_for(sleep_length);
  }

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelRepeat test"), elapsed.count(),
              TimerTolerance);

  {
    Timer timer{"GetTimeLabelRepeat test"};

    std::this_thread::sleep_for(sleep_length);
  }

  auto end2 = std::chrono::high_resolution_clock::now();
  seconds elapsed2 = end2 - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelRepeat test"), elapsed2.count(),
              TimerTolerance);
}

TEST(TimerTest, ResetTime) {
  using namespace bout::testing;

  Timer timer{"ResetTime test"};

  std::this_thread::sleep_for(sleep_length);

  timer.resetTime();

  auto start = std::chrono::high_resolution_clock::now();

  std::this_thread::sleep_for(sleep_length);

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(timer.getTime(), elapsed.count(), TimerTolerance);
}

TEST(TimerTest, ResetTimeLabelInScope) {
  using namespace bout::testing;

  Timer timer{"ResetTimeLabelInScope test"};

  std::this_thread::sleep_for(sleep_length);

  Timer::resetTime("ResetTimeLabelInScope test");

  auto start = std::chrono::high_resolution_clock::now();

  std::this_thread::sleep_for(sleep_length);

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("ResetTimeLabelInScope test"), elapsed.count(),
              TimerTolerance);
}

TEST(TimerTest, ResetTimeLabelOutOfScope) {
  using namespace bout::testing;

  {
    Timer timer{"ResetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(sleep_length);
  }

  Timer::resetTime("ResetTimeLabelOutOfScope test");

  auto start = std::chrono::high_resolution_clock::now();

  {
    Timer timer{"ResetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(sleep_length);
  }

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("ResetTimeLabelOutOfScope test"),
              elapsed.count(), TimerTolerance);
}

TEST(TimerTest, Cleanup) {
  using namespace bout::testing;

  {
    Timer timer{"Cleanup test"};

    std::this_thread::sleep_for(sleep_length);
  }

  Timer::cleanup();

  auto start = std::chrono::high_resolution_clock::now();

  {
    Timer timer{"Cleanup test"};

    std::this_thread::sleep_for(sleep_length);
  }

  auto end = std::chrono::high_resolution_clock::now();
  seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("Cleanup test"), elapsed.count(),
              TimerTolerance);
}
