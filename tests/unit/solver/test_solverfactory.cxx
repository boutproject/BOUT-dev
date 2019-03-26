#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "test_extras.hxx"
#include "bout/solver.hxx"
#include "bout/solverfactory.hxx"

#include <algorithm>

namespace {
class FakeSolver : public Solver {
public:
  FakeSolver(Options* options) : Solver(options) {}
  ~FakeSolver() = default;
  int run() { return (*options)["number"].withDefault(42); }
};

RegisterSolver<FakeSolver> register_fake("fake_solver");
} // namespace

TEST(SolverFactoryTest, GetInstance) { EXPECT_NE(SolverFactory::getInstance(), nullptr); }

TEST(SolverFactoryTest, GetDefaultSolverType) {
  EXPECT_NE(SolverFactory::getDefaultSolverType(), "");
}

TEST(SolverFactoryTest, RegisterSolver) {
  auto available = SolverFactory::getInstance()->listAvailable();

  auto found_fake = std::find(begin(available), end(available), "fake_solver");

  EXPECT_NE(found_fake, end(available));
}

TEST(SolverFactoryTest, Create) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "fake_solver";
  auto solver = SolverFactory::getInstance()->createSolver();

  EXPECT_EQ(solver->run(), 42);

  Options::cleanup();
}

TEST(SolverFactoryTest, CreateDefault) {
  WithQuietOutput quiet{output_info};

  EXPECT_NO_THROW(SolverFactory::getInstance()->createSolver());

  Options::cleanup();
}

TEST(SolverFactoryTest, CreateFromOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  options["type"] = "fake_solver";
  auto solver = SolverFactory::getInstance()->createSolver(&options);

  EXPECT_EQ(solver->run(), 42);
}

TEST(SolverFactoryTest, CreateFromName) {
  WithQuietOutput quiet{output_info};

  constexpr auto number = 13;
  Options::root()["solver"]["number"] = number;
  auto solver = SolverFactory::getInstance()->createSolver("fake_solver");

  EXPECT_EQ(solver->run(), number);

  Options::cleanup();
}

TEST(SolverFactoryTest, CreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  constexpr auto number = 13;
  Options options;
  options["number"] = number;
  auto solver = SolverFactory::getInstance()->createSolver("fake_solver", &options);

  EXPECT_EQ(solver->run(), number);
}

TEST(SolverFactoryTest, BadCreate) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "bad_solver";
  EXPECT_THROW(SolverFactory::getInstance()->createSolver(), BoutException);
  Options::cleanup();
}

TEST(SolverFactoryTest, BadCreateFromOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  options["type"] = "bad_solver";
  EXPECT_THROW(SolverFactory::getInstance()->createSolver(&options), BoutException);
}

TEST(SolverFactoryTest, BadCreateFromName) {
  WithQuietOutput quiet{output_info};

  EXPECT_THROW(SolverFactory::getInstance()->createSolver("bad_solver"), BoutException);
  Options::cleanup();
}

TEST(SolverFactoryTest, BadCreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  EXPECT_THROW(SolverFactory::getInstance()->createSolver("bad_solver", &options), BoutException);
}
