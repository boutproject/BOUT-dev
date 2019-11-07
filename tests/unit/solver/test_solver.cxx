#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "test_fakesolver.hxx"
#include "bout/solver.hxx"
#include "bout/solverfactory.hxx"

#include <algorithm>
#include <string>
#include <vector>

namespace {
/// A sentinel value for whether a `FakeMonitor` has been called or not
constexpr static int called_sentinel{-999};

/// A `Monitor` that returns a bad value when called past its \p
/// trigger_time_
class FakeMonitor : public Monitor {
public:
  FakeMonitor(BoutReal timestep = -1, BoutReal trigger_time_ = 0.)
      : Monitor(timestep), trigger_time(trigger_time_) {}
  ~FakeMonitor() = default;
  auto call(Solver*, BoutReal time, int iter, int) -> int {
    last_called = iter;
    return time > trigger_time ? -1 : 0;
  }
  auto getTimestepShim() -> BoutReal { return getTimestep(); }
  auto setTimestepShim(BoutReal timestep) -> void { setTimestep(timestep); }
  void cleanup() { cleaned = true; }

  int last_called{called_sentinel};
  bool cleaned{false};

private:
  BoutReal trigger_time{0.0};
};

} // namespace

class SolverTest : public FakeMeshFixture {
public:
  SolverTest() : FakeMeshFixture() {
    Options::root()["field"]["function"] = "1.0";
    Options::root()["field"]["solution"] = "2.0";
    Options::root()["another_field"]["function"] = "3.0";
    Options::root()["another_field"]["solution"] = "4.0";
    Options::root()["vector_x"]["function"] = "5.0";
    Options::root()["vector_y"]["function"] = "6.0";
    Options::root()["vector_z"]["function"] = "7.0";
    Options::root()["another_vectorx"]["function"] = "8.0";
    Options::root()["another_vectory"]["function"] = "9.0";
    Options::root()["another_vectorz"]["function"] = "10.0";
  }
  virtual ~SolverTest() { Options::cleanup(); }

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_progress{output_progress};
};

TEST_F(SolverTest, Create) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "fake_solver";
  auto solver = Solver::create();

  solver->run();

  EXPECT_TRUE(static_cast<FakeSolver*>(solver.get())->run_called);

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

  solver->run();

  EXPECT_TRUE(static_cast<FakeSolver*>(solver.get())->run_called);
}

TEST_F(SolverTest, CreateFromName) {
  WithQuietOutput quiet{output_info};

  constexpr auto fail_run = 13;
  Options::root()["solver"]["fail_run"] = fail_run;
  auto solver = Solver::create("fake_solver");

  EXPECT_EQ(solver->run(), fail_run);

  Options::cleanup();
}

TEST_F(SolverTest, CreateFromNameAndOptions) {
  WithQuietOutput quiet{output_info};

  constexpr auto fail_run = 13;
  Options options;
  options["fail_run"] = fail_run;
  auto solver = Solver::create("fake_solver", &options);

  EXPECT_EQ(solver->run(), fail_run);
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

  Field2D field1{}, field2{};
  EXPECT_NO_THROW(solver.add(field1, "field"));
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);
  EXPECT_TRUE(IsFieldEqual(field1, 1.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(field2, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);
#endif

  EXPECT_NO_THROW(solver.add(field2, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 2);
  EXPECT_EQ(solver.n3Dvars(), 0);
  EXPECT_TRUE(IsFieldEqual(field2, 3.0));

  const auto expected_names = std::vector<std::string>{"field", "another_field"};
  EXPECT_EQ(solver.listField2DNames(), expected_names);
}

TEST_F(SolverTest, AddField2DMMS) {
  Options options;
  options["mms"] = true;
  options["mms_initialise"] = true;
  FakeSolver solver{&options};

  Field2D field1{}, field2{};
  EXPECT_NO_THROW(solver.add(field1, "field"));
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);
  EXPECT_TRUE(IsFieldEqual(field1, 2.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(field2, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);
#endif

  EXPECT_NO_THROW(solver.add(field2, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 2);
  EXPECT_EQ(solver.n3Dvars(), 0);
  EXPECT_TRUE(IsFieldEqual(field2, 4.0));

  const auto expected_names = std::vector<std::string>{"field", "another_field"};
  EXPECT_EQ(solver.listField2DNames(), expected_names);
}

TEST_F(SolverTest, AddField3D) {
  Options options;
  FakeSolver solver{&options};

  Field3D field1{}, field2{};
  EXPECT_NO_THROW(solver.add(field1, "field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);
  EXPECT_TRUE(IsFieldEqual(field1, 1.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(field2, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);
#endif

  EXPECT_NO_THROW(solver.add(field2, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 2);
  EXPECT_TRUE(IsFieldEqual(field2, 3.0));

  const auto expected_names = std::vector<std::string>{"field", "another_field"};
  EXPECT_EQ(solver.listField3DNames(), expected_names);
}

TEST_F(SolverTest, AddField3DMMS) {
  Options options;
  options["mms"] = true;
  options["mms_initialise"] = true;
  FakeSolver solver{&options};

  Field3D field1{}, field2{};
  EXPECT_NO_THROW(solver.add(field1, "field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);
  EXPECT_TRUE(IsFieldEqual(field1, 2.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(field2, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);
#endif

  EXPECT_NO_THROW(solver.add(field2, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 2);
  EXPECT_TRUE(IsFieldEqual(field2, 4.0));

  const auto expected_names = std::vector<std::string>{"field", "another_field"};
  EXPECT_EQ(solver.listField3DNames(), expected_names);
}

TEST_F(SolverTest, AddVector2D) {
  Options options;
  FakeSolver solver{&options};

  Vector2D vector1{}, vector2{};
  EXPECT_NO_THROW(solver.add(vector1, "vector"));
  EXPECT_EQ(solver.n2Dvars(), 3);
  EXPECT_EQ(solver.n3Dvars(), 0);
  EXPECT_TRUE(IsFieldEqual(vector1.x, 5.0));
  EXPECT_TRUE(IsFieldEqual(vector1.y, 6.0));
  EXPECT_TRUE(IsFieldEqual(vector1.z, 7.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(vector2, "vector"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 3);
  EXPECT_EQ(solver.n3Dvars(), 0);
#endif

  vector2.covariant = false;
  EXPECT_NO_THROW(solver.add(vector2, "another_vector"));
  EXPECT_EQ(solver.n2Dvars(), 6);
  EXPECT_EQ(solver.n3Dvars(), 0);
  EXPECT_TRUE(IsFieldEqual(vector2.x, 8.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 9.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 10.0));

  const auto expected_names = std::vector<std::string>{"vector", "another_vector"};
  EXPECT_EQ(solver.listVector2DNames(), expected_names);
}

TEST_F(SolverTest, AddVector3D) {
  Options options;
  FakeSolver solver{&options};

  Vector3D vector1{}, vector2{};
  EXPECT_NO_THROW(solver.add(vector1, "vector"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 3);
  EXPECT_TRUE(IsFieldEqual(vector1.x, 5.0));
  EXPECT_TRUE(IsFieldEqual(vector1.y, 6.0));
  EXPECT_TRUE(IsFieldEqual(vector1.z, 7.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(vector2, "vector"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 3);
#endif

  vector2.covariant = false;
  EXPECT_NO_THROW(solver.add(vector2, "another_vector"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 6);
  EXPECT_TRUE(IsFieldEqual(vector2.x, 8.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 9.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 10.0));

  const auto expected_names = std::vector<std::string>{"vector", "another_vector"};
  EXPECT_EQ(solver.listVector3DNames(), expected_names);
}

TEST_F(SolverTest, ConstraintField2D) {
  Options options;
  FakeSolver solver{&options};

  Field2D field1{}, field2{};
  EXPECT_NO_THROW(solver.constraint(field1, field1, "field"));
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);

#if CHECK > 0
  EXPECT_THROW(solver.constraint(field2, field2, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);

  EXPECT_THROW(solver.constraint(field2, field2, ""), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);

  solver.changeHasConstraints(false);
  EXPECT_THROW(solver.constraint(field2, field2, "some_other_name"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 1);
  EXPECT_EQ(solver.n3Dvars(), 0);
  solver.changeHasConstraints(true);
#endif

  EXPECT_NO_THROW(solver.constraint(field2, field2, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 2);
  EXPECT_EQ(solver.n3Dvars(), 0);

  const auto expected_names = std::vector<std::string>{"field", "another_field"};
  EXPECT_EQ(solver.listField2DNames(), expected_names);
}

TEST_F(SolverTest, ConstraintField3D) {
  Options options;
  FakeSolver solver{&options};

  Field3D field1{}, field2{};
  EXPECT_NO_THROW(solver.constraint(field1, field1, "field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);

#if CHECK > 0
  EXPECT_THROW(solver.constraint(field2, field2, "field"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);

  EXPECT_THROW(solver.constraint(field2, field2, ""), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);

  solver.changeHasConstraints(false);
  EXPECT_THROW(solver.constraint(field2, field2, "some_other_name"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 1);
  solver.changeHasConstraints(true);
#endif

  EXPECT_NO_THROW(solver.constraint(field2, field2, "another_field"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 2);

  const auto expected_names = std::vector<std::string>{"field", "another_field"};
  EXPECT_EQ(solver.listField3DNames(), expected_names);
}

TEST_F(SolverTest, ConstraintVector2D) {
  Options options;
  FakeSolver solver{&options};

  Vector2D vector1{}, vector2{};
  EXPECT_NO_THROW(solver.constraint(vector1, vector1, "vector"));
  EXPECT_EQ(solver.n2Dvars(), 3);
  EXPECT_EQ(solver.n3Dvars(), 0);

#if CHECK > 0
  EXPECT_THROW(solver.constraint(vector2, vector2, "vector"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 3);
  EXPECT_EQ(solver.n3Dvars(), 0);

  EXPECT_THROW(solver.constraint(vector2, vector2, ""), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 3);
  EXPECT_EQ(solver.n3Dvars(), 0);

  solver.changeHasConstraints(false);
  EXPECT_THROW(solver.constraint(vector2, vector2, "some_other_name"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 3);
  EXPECT_EQ(solver.n3Dvars(), 0);
  solver.changeHasConstraints(true);
#endif

  vector2.covariant = false;
  EXPECT_NO_THROW(solver.constraint(vector2, vector2, "another_vector"));
  EXPECT_EQ(solver.n2Dvars(), 6);
  EXPECT_EQ(solver.n3Dvars(), 0);

  const auto expected_names = std::vector<std::string>{"vector", "another_vector"};
  EXPECT_EQ(solver.listVector2DNames(), expected_names);
}

TEST_F(SolverTest, ConstraintVector3D) {
  Options options;
  FakeSolver solver{&options};

  Vector3D vector1{}, vector2{};
  EXPECT_NO_THROW(solver.constraint(vector1, vector1, "vector"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 3);

#if CHECK > 0
  EXPECT_THROW(solver.constraint(vector2, vector2, "vector"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 3);

  EXPECT_THROW(solver.constraint(vector2, vector2, ""), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 3);

  solver.changeHasConstraints(false);
  EXPECT_THROW(solver.constraint(vector2, vector2, "some_other_name"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 3);
  solver.changeHasConstraints(true);
#endif

  vector2.covariant = false;
  EXPECT_NO_THROW(solver.constraint(vector2, vector2, "another_vector"));
  EXPECT_EQ(solver.n2Dvars(), 0);
  EXPECT_EQ(solver.n3Dvars(), 6);

  const auto expected_names = std::vector<std::string>{"vector", "another_vector"};
  EXPECT_EQ(solver.listVector3DNames(), expected_names);
}

TEST_F(SolverTest, NoInitTwice) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.init(0, 0));
  EXPECT_THROW(solver.init(0, 0), BoutException);
}

TEST_F(SolverTest, NoAddAfterInit) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.init(0, 0));

  Field2D field1{};
  EXPECT_THROW(solver.add(field1, "field"), BoutException);
  Field3D field2{};
  EXPECT_THROW(solver.add(field2, "field"), BoutException);
  Vector2D vector1{};
  EXPECT_THROW(solver.add(vector1, "vector"), BoutException);
  Vector3D vector2{};
  EXPECT_THROW(solver.add(vector2, "vector"), BoutException);
}

TEST_F(SolverTest, NoConstraintsAfterInit) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.init(0, 0));

  Field2D field1{};
  EXPECT_THROW(solver.constraint(field1, field1, "field"), BoutException);
  Field3D field2{};
  EXPECT_THROW(solver.constraint(field2, field2, "field"), BoutException);
  Vector2D vector1{};
  EXPECT_THROW(solver.constraint(vector1, vector1, "vector"), BoutException);
  Vector3D vector2{};
  EXPECT_THROW(solver.constraint(vector2, vector2, "vector"), BoutException);
}

TEST_F(SolverTest, SplitOperator) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_FALSE(solver.splitOperator());

  rhsfunc fake_rhs = [](BoutReal) -> int { return 0; };
  solver.setSplitOperator(fake_rhs, fake_rhs);

  EXPECT_TRUE(solver.splitOperator());
}

TEST_F(SolverTest, ResetInternalFields) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_THROW(solver.resetInternalFields(), BoutException);
}

TEST_F(SolverTest, SetMaxTimestep) {
  Options options;
  FakeSolver solver{&options};

  auto expected = 4.5;
  EXPECT_NO_THROW(solver.setMaxTimestep(expected));
  EXPECT_EQ(solver.getMaxTimestepShim(), expected);
}

TEST_F(SolverTest, GetCurrentTimestep) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_EQ(solver.getCurrentTimestep(), 0.0);
}

TEST_F(SolverTest, HasConstraints) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_TRUE(solver.constraints());

  solver.changeHasConstraints(false);

  EXPECT_FALSE(solver.constraints());
}

TEST_F(SolverTest, GetLocalN) {
  Options options;
  FakeSolver solver{&options};

  Options::root()["field2"]["evolve_bndry"] = true;
  Options::root()["field4"]["evolve_bndry"] = true;
  Options::root()["input"]["transform_from_field_aligned"] = false;

  constexpr auto localmesh_nx = 5;
  constexpr auto localmesh_ny = 7;
  constexpr auto localmesh_nz = 9;

  FakeMesh localmesh{localmesh_nx, localmesh_ny, localmesh_nz};
  localmesh.createDefaultRegions();
  localmesh.createBoundaryRegions();
  localmesh.setCoordinates(nullptr);

  Field2D field1{bout::globals::mesh};
  Field2D field2{&localmesh};
  Field3D field3{&localmesh};
  Field3D field4{bout::globals::mesh};

  solver.add(field1, "field1");
  solver.add(field2, "field2");
  solver.add(field3, "field3");
  solver.add(field4, "field4");

  solver.init(0, 0);

  static_cast<FakeMesh*>(field1.getMesh())->createBoundaryRegions();

  constexpr auto globalmesh_nx_no_boundry = SolverTest::nx - 2;
  constexpr auto globalmesh_ny_no_boundry = SolverTest::ny - 2;
  constexpr auto localmesh_nx_no_boundry = localmesh_nx - 2;
  constexpr auto localmesh_ny_no_boundry = localmesh_ny - 2;
  constexpr auto expected_total =
      (globalmesh_nx_no_boundry * globalmesh_ny_no_boundry)
      + (localmesh_nx * localmesh_ny)
      + (localmesh_nx_no_boundry * localmesh_ny_no_boundry * localmesh_nz)
      + (nx * ny * nz);

  EXPECT_EQ(solver.getLocalNShim(), expected_total);
}

TEST_F(SolverTest, HavePreconditioner) {
  PhysicsPrecon preconditioner = [](BoutReal time, BoutReal gamma,
                                    BoutReal delta) -> int {
    return static_cast<int>(time + gamma + delta);
  };

  Options options;
  FakeSolver solver{&options};

  EXPECT_FALSE(solver.haveUserPreconShim());

  solver.setPrecon(preconditioner);

  EXPECT_TRUE(solver.haveUserPreconShim());
}

TEST_F(SolverTest, RunPreconditioner) {
  PhysicsPrecon preconditioner = [](BoutReal time, BoutReal gamma,
                                    BoutReal delta) -> int {
    return static_cast<int>(time + gamma + delta);
  };

  Options options;
  FakeSolver solver{&options};

  solver.setPrecon(preconditioner);

  constexpr auto time = 1.0;
  constexpr auto gamma = 2.0;
  constexpr auto delta = 3.0;
  constexpr auto expected = time + gamma + delta;

  EXPECT_EQ(solver.runPreconShim(time, gamma, delta), expected);
}

TEST_F(SolverTest, AddMonitor) {
  Options options;
  FakeSolver solver{&options};

  FakeMonitor monitor;
  EXPECT_NO_THROW(monitor.setTimestepShim(10.0));
  EXPECT_EQ(monitor.getTimestepShim(), 10.0);

  EXPECT_NO_THROW(solver.addMonitor(&monitor));

  EXPECT_THROW(monitor.setTimestepShim(20.0), BoutException);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 0, 0));

  EXPECT_EQ(monitor.last_called, 0);
}

TEST_F(SolverTest, AddMonitorFront) {
  WithQuietOutput quiet{output_error};
  Options options;
  FakeSolver solver{&options};

  FakeMonitor monitor1;
  FakeMonitor monitor2;
  EXPECT_NO_THROW(solver.addMonitor(&monitor1, Solver::FRONT));
  EXPECT_NO_THROW(solver.addMonitor(&monitor2, Solver::FRONT));

  // Everything's fine
  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 0, 0));

  EXPECT_EQ(monitor1.last_called, 0);
  EXPECT_EQ(monitor2.last_called, 0);

  // One monitor signals to quit
  EXPECT_THROW(solver.callMonitorsShim(5.0, 1, 0), BoutException);

  EXPECT_EQ(monitor1.last_called, 0);
  EXPECT_EQ(monitor2.last_called, 1);

  // Last timestep
  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 9, 10));

  EXPECT_EQ(monitor1.last_called, 9);
  EXPECT_EQ(monitor2.last_called, 9);
  EXPECT_TRUE(monitor1.cleaned);
  EXPECT_TRUE(monitor2.cleaned);
}

TEST_F(SolverTest, AddMonitorBack) {
  WithQuietOutput quiet{output_error};
  Options options;
  FakeSolver solver{&options};

  FakeMonitor monitor1;
  FakeMonitor monitor2;
  EXPECT_NO_THROW(solver.addMonitor(&monitor1, Solver::BACK));
  EXPECT_NO_THROW(solver.addMonitor(&monitor2, Solver::BACK));

  // Everything's fine
  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 0, 0));

  EXPECT_EQ(monitor1.last_called, 0);
  EXPECT_EQ(monitor2.last_called, 0);

  // One monitor signals to quit
  EXPECT_THROW(solver.callMonitorsShim(5.0, 1, 0), BoutException);

  EXPECT_EQ(monitor1.last_called, 1);
  EXPECT_EQ(monitor2.last_called, 0);

  // Last timestep
  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 9, 10));

  EXPECT_EQ(monitor1.last_called, 9);
  EXPECT_EQ(monitor2.last_called, 9);
  EXPECT_TRUE(monitor1.cleaned);
  EXPECT_TRUE(monitor2.cleaned);
}

TEST_F(SolverTest, AddMonitorCheckFrequencies) {
  Options options;
  FakeSolver solver{&options};

  FakeMonitor default_timestep;
  FakeMonitor smaller_timestep{0.1};
  FakeMonitor even_smaller_timestep{0.01};
  FakeMonitor larger_timestep{2.};
  FakeMonitor incompatible_timestep{3.14259};

  EXPECT_NO_THROW(solver.addMonitor(&default_timestep));
  EXPECT_NO_THROW(solver.addMonitor(&smaller_timestep));
  EXPECT_NO_THROW(solver.addMonitor(&even_smaller_timestep));
  EXPECT_NO_THROW(solver.addMonitor(&larger_timestep));
  EXPECT_THROW(solver.addMonitor(&incompatible_timestep), BoutException);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, -1, 0));

  EXPECT_EQ(default_timestep.last_called, -1);
  EXPECT_EQ(smaller_timestep.last_called, -1);
  EXPECT_EQ(even_smaller_timestep.last_called, -1);
  EXPECT_EQ(larger_timestep.last_called, -1);
  EXPECT_EQ(incompatible_timestep.last_called, called_sentinel);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 9, 0));

  EXPECT_EQ(default_timestep.last_called, 0);
  EXPECT_EQ(smaller_timestep.last_called, 0);
  EXPECT_EQ(even_smaller_timestep.last_called, 9);
  EXPECT_EQ(larger_timestep.last_called, -1);
  EXPECT_EQ(incompatible_timestep.last_called, called_sentinel);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 10, 0));

  EXPECT_EQ(default_timestep.last_called, 0);
  EXPECT_EQ(smaller_timestep.last_called, 0);
  EXPECT_EQ(even_smaller_timestep.last_called, 10);
  EXPECT_EQ(larger_timestep.last_called, -1);
  EXPECT_EQ(incompatible_timestep.last_called, called_sentinel);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 199, 0));

  EXPECT_EQ(default_timestep.last_called, 19);
  EXPECT_EQ(smaller_timestep.last_called, 19);
  EXPECT_EQ(even_smaller_timestep.last_called, 199);
  EXPECT_EQ(larger_timestep.last_called, 0);
  EXPECT_EQ(incompatible_timestep.last_called, called_sentinel);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 399, 0));

  EXPECT_EQ(default_timestep.last_called, 39);
  EXPECT_EQ(smaller_timestep.last_called, 39);
  EXPECT_EQ(even_smaller_timestep.last_called, 399);
  EXPECT_EQ(larger_timestep.last_called, 1);
  EXPECT_EQ(incompatible_timestep.last_called, called_sentinel);

  solver.init(0, 0);

  FakeMonitor too_small_postinit_timestep{0.001};
  EXPECT_THROW(solver.addMonitor(&too_small_postinit_timestep), BoutException);
  FakeMonitor larger_postinit_timestep{4.};
  EXPECT_NO_THROW(solver.addMonitor(&larger_postinit_timestep, Solver::BACK));

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 399, 0));

  EXPECT_EQ(default_timestep.last_called, 39);
  EXPECT_EQ(smaller_timestep.last_called, 39);
  EXPECT_EQ(even_smaller_timestep.last_called, 399);
  EXPECT_EQ(larger_timestep.last_called, 1);
  EXPECT_EQ(larger_postinit_timestep.last_called, 0);
  EXPECT_EQ(incompatible_timestep.last_called, called_sentinel);
}

TEST_F(SolverTest, RemoveMonitor) {
  Options options;
  FakeSolver solver{&options};

  FakeMonitor monitor1;
  FakeMonitor monitor2;
  EXPECT_NO_THROW(solver.addMonitor(&monitor1, Solver::BACK));
  EXPECT_NO_THROW(solver.addMonitor(&monitor2, Solver::BACK));

  solver.removeMonitor(&monitor1);

  std::list<Monitor*> expected{&monitor2};
  EXPECT_EQ(solver.getMonitorsShim(), expected);

  // Removing same monitor again should be a no-op
  solver.removeMonitor(&monitor1);
  EXPECT_EQ(solver.getMonitorsShim(), expected);
}

namespace {
auto timestep_monitor1(Solver*, BoutReal simtime, BoutReal lastdt) -> int {
  return simtime * lastdt < 0. ? 1 : 0;
}
auto timestep_monitor2(Solver*, BoutReal simtime, BoutReal lastdt) -> int {
  return simtime + lastdt < 0. ? 2 : 0;
}
} // namespace

TEST_F(SolverTest, AddTimestepMonitor) {
  Options options;
  options["monitor_timestep"] = true;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.addTimestepMonitor(timestep_monitor1));
  EXPECT_NO_THROW(solver.addTimestepMonitor(timestep_monitor2));

  EXPECT_EQ(solver.callTimestepMonitorsShim(1., 1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(1., -1.), 1);
  EXPECT_EQ(solver.callTimestepMonitorsShim(-1., -1.), 2);
}

TEST_F(SolverTest, RemoveTimestepMonitor) {
  Options options;
  options["monitor_timestep"] = true;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.addTimestepMonitor(timestep_monitor1));
  EXPECT_NO_THROW(solver.addTimestepMonitor(timestep_monitor2));

  solver.removeTimestepMonitor(timestep_monitor1);

  EXPECT_EQ(solver.callTimestepMonitorsShim(1., 1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(1., -1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(-1., -1.), 2);

  solver.removeTimestepMonitor(timestep_monitor1);

  EXPECT_EQ(solver.callTimestepMonitorsShim(1., 1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(1., -1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(-1., -1.), 2);
}

TEST_F(SolverTest, DontCallTimestepMonitors) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.addTimestepMonitor(timestep_monitor1));
  EXPECT_NO_THROW(solver.addTimestepMonitor(timestep_monitor2));

  EXPECT_EQ(solver.callTimestepMonitorsShim(1., 1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(1., -1.), 0);
  EXPECT_EQ(solver.callTimestepMonitorsShim(-1., -1.), 0);
}

TEST_F(SolverTest, BasicSolve) {
  Options options;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  EXPECT_NO_THROW(solver.solve());

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);
}

TEST_F(SolverTest, SolveBadInit) {
  Options options;
  options["fail_init"] = -1;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  EXPECT_THROW(solver.solve(), BoutException);

  EXPECT_TRUE(solver.init_called);
  EXPECT_FALSE(solver.run_called);
}

TEST_F(SolverTest, SolveBadRun) {
  Options options;
  options["fail_run"] = -1;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  EXPECT_EQ(solver.solve(), -1);

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);
}

TEST_F(SolverTest, SolveThrowRun) {
  WithQuietOutput quiet_error{output_error};

  Options options;
  options["throw_run"] = true;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  EXPECT_THROW(solver.solve(), BoutException);

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);
}

TEST_F(SolverTest, SolveFixDefaultTimestep) {
  Options options;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  FakeMonitor default_timestep;
  FakeMonitor smaller_timestep{0.1};
  FakeMonitor even_smaller_timestep{0.01};
  FakeMonitor larger_timestep{2.};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);
  solver.addMonitor(&even_smaller_timestep);
  solver.addMonitor(&larger_timestep);

  EXPECT_NO_THROW(solver.solve(100, 1.));

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 99, 0));

  EXPECT_EQ(default_timestep.last_called, 0);
  EXPECT_EQ(smaller_timestep.last_called, 9);
  EXPECT_EQ(even_smaller_timestep.last_called, 99);
  EXPECT_EQ(larger_timestep.last_called, called_sentinel);
}

TEST_F(SolverTest, SolveFixDefaultTimestepBad) {
  Options options;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  FakeMonitor default_timestep;
  FakeMonitor smaller_timestep{0.1};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);

  EXPECT_THROW(solver.solve(100, 3.142), BoutException);

  EXPECT_FALSE(solver.init_called);
  EXPECT_FALSE(solver.run_called);
}

TEST_F(SolverTest, SolveFixDefaultTimestepSmaller) {
  Options options;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  FakeMonitor default_timestep;
  FakeMonitor smaller_timestep{0.1};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);

  EXPECT_NO_THROW(solver.solve(100, 0.01));

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 99, 0));

  EXPECT_EQ(default_timestep.last_called, 99);
  EXPECT_EQ(smaller_timestep.last_called, 9);
}

TEST_F(SolverTest, SolveFixDefaultTimestepLarger) {
  Options options;
  FakeSolver solver{&options};

  Options::root()["dump_on_restart"] = false;

  FakeMonitor default_timestep;
  FakeMonitor smaller_timestep{0.1};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);

  EXPECT_NO_THROW(solver.solve(100, 1.));

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);

  EXPECT_NO_THROW(solver.callMonitorsShim(0.0, 99, 0));

  EXPECT_EQ(default_timestep.last_called, 9);
  EXPECT_EQ(smaller_timestep.last_called, 99);
}
