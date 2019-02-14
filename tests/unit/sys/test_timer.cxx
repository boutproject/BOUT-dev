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
constexpr auto sleep_length_seconds = seconds(sleep_length);
} // namespace testing
} // namespace bout

TEST(TimerTest, GetTime) {
  using namespace bout::testing;

  Timer timer{};

  std::this_thread::sleep_for(sleep_length);

  EXPECT_NEAR(timer.getTime(), sleep_length_seconds.count(), TimerTolerance);
}

TEST(TimerTest, GetUnknownTimeLabel) {
  EXPECT_DOUBLE_EQ(Timer::getTime("Does not exist"), 0.0);
}

TEST(TimerTest, GetTimeLabelInScope) {
  using namespace bout::testing;

  Timer timer{"GetTimeLabelInScope test"};

  std::this_thread::sleep_for(sleep_length);

  EXPECT_NEAR(Timer::getTime("GetTimeLabelInScope test"), sleep_length_seconds.count(),
              TimerTolerance);
}

TEST(TimerTest, GetTimeLabelOutOfScope) {
  using namespace bout::testing;

  {
    Timer timer{"GetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(sleep_length);
  }

  EXPECT_NEAR(Timer::getTime("GetTimeLabelOutOfScope test"), sleep_length_seconds.count(),
              TimerTolerance);
}

TEST(TimerTest, GetTimeLabelRepeat) {
  using namespace bout::testing;

  {
    Timer timer{"GetTimeLabelRepeat test"};

    std::this_thread::sleep_for(sleep_length);
  }

  EXPECT_NEAR(Timer::getTime("GetTimeLabelRepeat test"), sleep_length_seconds.count(),
              TimerTolerance);

  {
    Timer timer{"GetTimeLabelRepeat test"};

    std::this_thread::sleep_for(sleep_length);
  }

  EXPECT_NEAR(Timer::getTime("GetTimeLabelRepeat test"), 2 * sleep_length_seconds.count(),
              TimerTolerance);
}

TEST(TimerTest, ResetTime) {
  using namespace bout::testing;

  Timer timer{"ResetTime test"};

  std::this_thread::sleep_for(sleep_length);

  timer.resetTime();

  std::this_thread::sleep_for(sleep_length);

  EXPECT_NEAR(timer.getTime(), sleep_length_seconds.count(), TimerTolerance);
}

TEST(TimerTest, ResetTimeLabelInScope) {
  using namespace bout::testing;

  Timer timer{"ResetTimeLabelInScope test"};

  std::this_thread::sleep_for(sleep_length);

  Timer::resetTime("ResetTimeLabelInScope test");

  std::this_thread::sleep_for(sleep_length);

  EXPECT_NEAR(Timer::getTime("ResetTimeLabelInScope test"), sleep_length_seconds.count(),
              TimerTolerance);
}

TEST(TimerTest, ResetTimeLabelOutOfScope) {
  using namespace bout::testing;

  {
    Timer timer{"ResetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(sleep_length);
  }

  Timer::resetTime("ResetTimeLabelOutOfScope test");

  {
    Timer timer{"ResetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(sleep_length);
  }

  EXPECT_NEAR(Timer::getTime("ResetTimeLabelOutOfScope test"),
              sleep_length_seconds.count(), TimerTolerance);
}

TEST(TimerTest, Cleanup) {
  using namespace bout::testing;

  {
    Timer timer{"Cleanup test"};

    std::this_thread::sleep_for(sleep_length);
  }

  Timer::cleanup();

  {
    Timer timer{"Cleanup test"};

    std::this_thread::sleep_for(sleep_length);
  }

  EXPECT_NEAR(Timer::getTime("Cleanup test"), sleep_length_seconds.count(),
              TimerTolerance);
}
