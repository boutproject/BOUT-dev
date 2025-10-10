#include "gtest/gtest.h"

#include "test_fakesolver.hxx"
#include "bout/boutexception.hxx"
#include "bout/output.hxx"
#include "bout/solver.hxx"

#include <algorithm>
#include <iterator>

TEST(SolverFactoryTest, GetInstance) { EXPECT_NO_THROW(SolverFactory::getInstance()); }

TEST(SolverFactoryTest, GetDefaultSolverType) {
  EXPECT_NE(SolverFactory::getDefaultType(), "");
}

TEST(SolverFactoryTest, RegisterSolver) {
  auto available = SolverFactory::getInstance().listAvailable();

  auto found_fake = std::find(begin(available), end(available), "fake_solver");

  EXPECT_NE(found_fake, end(available));
}

TEST(SolverFactoryTest, Create) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "fake_solver";
  auto solver = SolverFactory::getInstance().create();

  solver->run();

  EXPECT_TRUE(static_cast<FakeSolver*>(solver.get())->run_called);

  Options::cleanup();
}

TEST(SolverFactoryTest, CreateDefault) {
  WithQuietOutput quiet{output_info};

  EXPECT_NO_THROW(SolverFactory::getInstance().create());

  Options::cleanup();
}

TEST(SolverFactoryTest, CreateFromOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  options["type"] = "fake_solver";
  auto solver = SolverFactory::getInstance().create(&options);

  solver->run();

  EXPECT_TRUE(static_cast<FakeSolver*>(solver.get())->run_called);
}

TEST(SolverFactoryTest, CreateFromName) {
  WithQuietOutput quiet{output_info};

  constexpr auto fail_run = 13;
  Options::root()["solver"]["fail_run"] = fail_run;
  auto solver = SolverFactory::getInstance().create("fake_solver");

  EXPECT_EQ(solver->run(), fail_run);

  Options::cleanup();
}

TEST(SolverFactoryTest, CreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  constexpr auto fail_run = 31;
  Options options;
  options["fail_run"] = fail_run;
  auto solver = SolverFactory::getInstance().create("fake_solver", &options);

  EXPECT_EQ(solver->run(), fail_run);
}

TEST(SolverFactoryTest, BadCreate) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "bad_solver";
  EXPECT_THROW(SolverFactory::getInstance().create(), BoutException);
  Options::cleanup();
}

TEST(SolverFactoryTest, BadCreateFromOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  options["type"] = "bad_solver";
  EXPECT_THROW(SolverFactory::getInstance().create(&options), BoutException);
}

TEST(SolverFactoryTest, BadCreateFromName) {
  WithQuietOutput quiet{output_info};

  EXPECT_THROW(SolverFactory::getInstance().create("bad_solver"), BoutException);
  Options::cleanup();
}

TEST(SolverFactoryTest, BadCreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  EXPECT_THROW(SolverFactory::getInstance().create("bad_solver", &options),
               BoutException);
}
