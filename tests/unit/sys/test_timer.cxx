#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "bout/sys/timer.hxx"

#include <chrono>
#include <thread>

namespace bout {
namespace testing {
using ms = std::chrono::duration<double, std::chrono::milliseconds::period>;
constexpr double TimerTolerance{0.5e-3};
constexpr auto sleep_length = ms(1.);
} // namespace testing
} // namespace bout

TEST(TimerTest, GetTime) {
  auto start = Timer::clock_type::now();

  Timer timer{};

  std::this_thread::sleep_for(bout::testing::sleep_length);

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(timer.getTime(), elapsed.count(), bout::testing::TimerTolerance);
}

TEST(TimerTest, GetUnknownTimeLabel) {
  EXPECT_DOUBLE_EQ(Timer::getTime("Does not exist"), 0.0);
}

TEST(TimerTest, GetTimeLabelInScope) {
  auto start = Timer::clock_type::now();

  Timer timer{"GetTimeLabelInScope test"};

  std::this_thread::sleep_for(bout::testing::sleep_length);

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelInScope test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, GetTimeLabelOutOfScope) {
  auto start = Timer::clock_type::now();

  {
    Timer timer{"GetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelOutOfScope test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, GetTimeLabelSubScope) {
  auto start = Timer::clock_type::now();

  Timer timer{"GetTimeLabelSubScope test"};

  {
    Timer timer{"GetTimeLabelSubScope test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  std::this_thread::sleep_for(bout::testing::sleep_length);

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelSubScope test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, GetTimeLabelRepeat) {
  auto start = Timer::clock_type::now();

  {
    Timer timer{"GetTimeLabelRepeat test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelRepeat test"), elapsed.count(),
              bout::testing::TimerTolerance);

  {
    Timer timer{"GetTimeLabelRepeat test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  auto end2 = Timer::clock_type::now();
  Timer::seconds elapsed2 = end2 - start;

  EXPECT_NEAR(Timer::getTime("GetTimeLabelRepeat test"), elapsed2.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, ResetTime) {
  Timer timer{"ResetTime test"};

  std::this_thread::sleep_for(bout::testing::sleep_length);

  timer.resetTime();

  auto start = Timer::clock_type::now();

  std::this_thread::sleep_for(bout::testing::sleep_length);

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(timer.getTime(), elapsed.count(), bout::testing::TimerTolerance);
}

TEST(TimerTest, ResetTimeLabelInScope) {
  Timer timer{"ResetTimeLabelInScope test"};

  std::this_thread::sleep_for(bout::testing::sleep_length);

  Timer::resetTime("ResetTimeLabelInScope test");

  auto start = Timer::clock_type::now();

  std::this_thread::sleep_for(bout::testing::sleep_length);

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("ResetTimeLabelInScope test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, ResetTimeLabelOutOfScope) {
  {
    Timer timer{"ResetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  Timer::resetTime("ResetTimeLabelOutOfScope test");

  auto start = Timer::clock_type::now();

  {
    Timer timer{"ResetTimeLabelOutOfScope test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("ResetTimeLabelOutOfScope test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, GetTotalTime) {
  const auto start = Timer::clock_type::now();
  Timer timer{"GetTotalTime test"};

  std::this_thread::sleep_for(bout::testing::sleep_length);

  timer.resetTime();

  std::this_thread::sleep_for(bout::testing::sleep_length);

  const auto end = Timer::clock_type::now();
  const Timer::seconds elapsed = end - start;

  EXPECT_NEAR(timer.getTotalTime(), elapsed.count(), bout::testing::TimerTolerance);
}

TEST(TimerTest, GetTotalTimeFromLabel) {
  const auto start = Timer::clock_type::now();
  Timer timer{"GetTotalTimeFromLabel test"};

  std::this_thread::sleep_for(bout::testing::sleep_length);

  timer.resetTime();

  std::this_thread::sleep_for(bout::testing::sleep_length);

  const auto end = Timer::clock_type::now();
  const Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTotalTime("GetTotalTimeFromLabel test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, Cleanup) {
  {
    Timer timer{"Cleanup test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  Timer::cleanup();

  auto start = Timer::clock_type::now();

  {
    Timer timer{"Cleanup test"};

    std::this_thread::sleep_for(bout::testing::sleep_length);
  }

  auto end = Timer::clock_type::now();
  Timer::seconds elapsed = end - start;

  EXPECT_NEAR(Timer::getTime("Cleanup test"), elapsed.count(),
              bout::testing::TimerTolerance);
}

TEST(TimerTest, ListAllInfo) {
  Timer::cleanup();

  {
    Timer time1{"one"};
    std::this_thread::sleep_for(bout::testing::sleep_length);
    {
      Timer time2{"two"};
      std::this_thread::sleep_for(bout::testing::sleep_length);
    }
    {
      Timer time2{"two"};
      std::this_thread::sleep_for(bout::testing::sleep_length);
    }
    Timer time_long{"18 characters long"};
  }

  Timer::resetTime("two");

  std::streambuf* old_cout_rdbuf(std::cout.rdbuf());
  std::stringstream cout_capture;
  std::cout.rdbuf(cout_capture.rdbuf());

  Timer::listAllInfo();

  std::cout.rdbuf(old_cout_rdbuf);

  using namespace ::testing;
  EXPECT_THAT(cout_capture.str(), HasSubstr("Timer name         |"));
  EXPECT_THAT(cout_capture.str(), ContainsRegex("one *| 0\\.\\d+ | 1    | 0\\.\\d+"));
  EXPECT_THAT(cout_capture.str(), ContainsRegex("two *| 0 * | 2    | 0\\.\\d+"));
}
