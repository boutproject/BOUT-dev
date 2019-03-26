#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
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

class SolverTest : public FakeMeshFixture {
public:
  SolverTest() : FakeMeshFixture() {}
  virtual ~SolverTest() = default;

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_progress{output_progress};
};

TEST_F(SolverTest, Create) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "fake_solver";
  auto solver = Solver::create();

  EXPECT_EQ(solver->run(), 42);

  Options::cleanup();
}

TEST_F(SolverTest, CreateDefault) {
  WithQuietOutput quiet{output_info};

  EXPECT_NO_THROW(Solver::create());

  Options::cleanup();
}

TEST_F(SolverTest, CreateFromOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  options["type"] = "fake_solver";
  auto solver = Solver::create(&options);

  EXPECT_EQ(solver->run(), 42);
}

TEST_F(SolverTest, CreateFromName) {
  WithQuietOutput quiet{output_info};

  constexpr auto number = 13;
  Options::root()["solver"]["number"] = number;
  auto solver = Solver::create("fake_solver");

  EXPECT_EQ(solver->run(), number);

  Options::cleanup();
}

TEST_F(SolverTest, CreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  constexpr auto number = 13;
  Options options;
  options["number"] = number;
  auto solver = Solver::create("fake_solver", &options);

  EXPECT_EQ(solver->run(), number);
}

TEST_F(SolverTest, BadCreate) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "bad_solver";
  EXPECT_THROW(Solver::create(), BoutException);
  Options::cleanup();
}

TEST_F(SolverTest, BadCreateFromOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  options["type"] = "bad_solver";
  EXPECT_THROW(Solver::create(&options), BoutException);
}

TEST_F(SolverTest, BadCreateFromName) {
  WithQuietOutput quiet{output_info};

  EXPECT_THROW(Solver::create("bad_solver"), BoutException);
  Options::cleanup();
}

TEST_F(SolverTest, BadCreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  Options options;
  EXPECT_THROW(Solver::create("bad_solver", &options), BoutException);
}

TEST_F(SolverTest, AddField2D) {
  Options options;
  FakeSolver solver{&options};

  Field2D field{};
  EXPECT_NO_THROW(solver.add(field, "field"));
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);

#if CHECK > 0
  EXPECT_THROW(solver.add(field, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);
#endif

  EXPECT_NO_THROW(solver.add(field, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 2);
  EXPECT_EQ(solver.n3Dvars(), 0);
}

TEST_F(SolverTest, AddField3D) {
  Options options;
  FakeSolver solver{&options};

  Field3D field{};
  EXPECT_NO_THROW(solver.add(field, "field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);

#if CHECK > 0
  EXPECT_THROW(solver.add(field, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);
#endif

  EXPECT_NO_THROW(solver.add(field, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 2);
}

TEST_F(SolverTest, Init) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.init(0, 0));
  EXPECT_THROW(solver.init(0, 0), BoutException);

  Field2D field{};
  EXPECT_THROW(solver.add(field, "field"), BoutException);
}

TEST_F(SolverTest, SplitOperator) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_FALSE(solver.splitOperator());
}

TEST_F(SolverTest, ResetInternalFields) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_THROW(solver.resetInternalFields(), BoutException);
}

TEST_F(SolverTest, SetMaxTimestep) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.setMaxTimestep(4.5));
}

TEST_F(SolverTest, GetCurrentTimestep) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_EQ(solver.getCurrentTimestep(), 0.0);
}
