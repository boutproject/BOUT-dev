#include "gmock/gmock.h"
#include "gtest/gtest.h"

#include "test_extras.hxx"
#include "test_fakesolver.hxx"
#include "bout/boutexception.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/physicsmodel.hxx"
#include "bout/solver.hxx"
#include "bout/sys/uuid.h"

#include <string>
#include <vector>

#include "fake_mesh_fixture.hxx"

namespace {
constexpr int nout = 100;
constexpr BoutReal timestep = 1.;

using testing::_;
using testing::NiceMock;
using testing::StrictMock;

class MockMonitor : public Monitor {
public:
  MockMonitor(BoutReal timestep = -1) : Monitor(timestep) {}
  MockMonitor(const MockMonitor&) = delete;
  MockMonitor(MockMonitor&&) = delete;
  MockMonitor& operator=(const MockMonitor&) = delete;
  MockMonitor& operator=(MockMonitor&&) = delete;
  ~MockMonitor() override = default;

  MOCK_METHOD(int, call, (Solver * solver, BoutReal time, int iter, int), (override));
  auto getTimestepShim() -> BoutReal { return getTimestep(); }
  auto setTimestepShim(BoutReal timestep) -> void { setTimestep(timestep); }
  MOCK_METHOD(void, cleanup, (), (override));
};

class MockPhysicsModel : public PhysicsModel {
public:
  // Don't enable the output/restart files
  MockPhysicsModel() : PhysicsModel(bout::globals::mesh, false, false) {}
  MOCK_METHOD(int, init, (bool restarting), (override));
  // Mock postInit even though it's not pure virtual because it does
  // stuff with files
  MOCK_METHOD(int, postInit, (bool restarting), (override));
  MOCK_METHOD(int, rhs, (BoutReal time), (override));
  MOCK_METHOD(int, convective, (BoutReal t), (override));
  MOCK_METHOD(int, diffusive, (BoutReal t), (override));
  MOCK_METHOD(int, diffusive, (BoutReal t, bool linear), (override));
  MOCK_METHOD(int, outputMonitor, (BoutReal simtime, int iter, int NOUT), (override));
  MOCK_METHOD(int, timestepMonitor, (BoutReal simtime, BoutReal dt), (override));

  int preconditioner(BoutReal time, BoutReal gamma, BoutReal delta) {
    return static_cast<int>(time + gamma + delta);
  }

  int jacobian(BoutReal time) { return static_cast<int>(time); }

  // Expose some protected methods to aid testing
  using PhysicsModel::setJacobian;
  using PhysicsModel::setPrecon;
  using PhysicsModel::setSplitOperator;
};

} // namespace

class SolverTest : public FakeMeshFixture {
public:
  SolverTest() {
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
    Options::root()["input"]["error_on_unused_options"] = false;
  }
  SolverTest(const SolverTest&) = delete;
  SolverTest(SolverTest&&) = delete;
  SolverTest& operator=(const SolverTest&) = delete;
  SolverTest& operator=(SolverTest&&) = delete;
  ~SolverTest() override { Options::cleanup(); }

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_progress{output_progress};
};

TEST_F(SolverTest, Create) {
  WithQuietOutput quiet{output_info};

  Options::root()["solver"]["type"] = "fake_solver";
  auto solver = Solver::create();

  solver->run();

  EXPECT_TRUE(dynamic_cast<FakeSolver*>(solver.get())->run_called);

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

  EXPECT_TRUE(dynamic_cast<FakeSolver*>(solver.get())->run_called);
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

TEST_F(SolverTest, SetModel) {
  Options options;
  FakeSolver solver{&options};

  MockPhysicsModel model{};
  EXPECT_CALL(model, init).Times(1);
  EXPECT_CALL(model, postInit).Times(1);

  solver.setModel(&model);

  // Can't set a second model
  EXPECT_THROW(solver.setModel(&model), BoutException);
}

TEST_F(SolverTest, SetModelAfterInit) {
  Options options;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};

  solver.init();

  EXPECT_THROW(solver.setModel(&model), BoutException);
}

TEST_F(SolverTest, AddField2D) {
  Options options;
  FakeSolver solver{&options};

  Field2D field1{};
  Field2D field2{};
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

  Field2D field1{};
  Field2D field2{};
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

  Field3D field1{};
  Field3D field2{};
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

  Field3D field1{};
  Field3D field2{};
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

  Vector2D vector1{};
  Vector2D vector2{};
  EXPECT_NO_THROW(solver.add(vector1, "vector"));
#if not(BOUT_USE_METRIC_3D)
  constexpr int n2d = 3;
  constexpr int n3d = 0;
#else
  constexpr int n2d = 0;
  constexpr int n3d = 3;
#endif
  EXPECT_EQ(solver.n2Dvars(), n2d);
  EXPECT_EQ(solver.n3Dvars(), n3d);
  EXPECT_TRUE(IsFieldEqual(vector1.x, 5.0));
  EXPECT_TRUE(IsFieldEqual(vector1.y, 6.0));
  EXPECT_TRUE(IsFieldEqual(vector1.z, 7.0));

#if CHECK > 0
  EXPECT_THROW(solver.add(vector2, "vector"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), n2d);
  EXPECT_EQ(solver.n3Dvars(), n3d);
#endif

  vector2.covariant = false;
  EXPECT_NO_THROW(solver.add(vector2, "another_vector"));
  EXPECT_EQ(solver.n2Dvars(), n2d * 2);
  EXPECT_EQ(solver.n3Dvars(), n3d * 2);
  EXPECT_TRUE(IsFieldEqual(vector2.x, 8.0));
  EXPECT_TRUE(IsFieldEqual(vector2.y, 9.0));
  EXPECT_TRUE(IsFieldEqual(vector2.z, 10.0));

  const auto expected_names = std::vector<std::string>{"vector", "another_vector"};
  EXPECT_EQ(solver.listVector2DNames(), expected_names);
}

TEST_F(SolverTest, AddVector3D) {
  Options options;
  FakeSolver solver{&options};

  Vector3D vector1{};
  Vector3D vector2{};
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

  Field2D field1{};
  Field2D field2{};
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

  Field3D field1{};
  Field3D field2{};
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

  Vector2D vector1{};
  Vector2D vector2{};
  EXPECT_NO_THROW(solver.constraint(vector1, vector1, "vector"));
#if not(BOUT_USE_METRIC_3D)
  constexpr int n2d = 3;
  constexpr int n3d = 0;
#else
  constexpr int n2d = 0;
  constexpr int n3d = 3;
#endif
  EXPECT_EQ(solver.n2Dvars(), n2d);
  EXPECT_EQ(solver.n3Dvars(), n3d);

#if CHECK > 0
  EXPECT_THROW(solver.constraint(vector2, vector2, "vector"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), n2d);
  EXPECT_EQ(solver.n3Dvars(), n3d);

  EXPECT_THROW(solver.constraint(vector2, vector2, ""), BoutException);
  EXPECT_EQ(solver.n2Dvars(), n2d);
  EXPECT_EQ(solver.n3Dvars(), n3d);

  solver.changeHasConstraints(false);
  EXPECT_THROW(solver.constraint(vector2, vector2, "some_other_name"), BoutException);
  EXPECT_EQ(solver.n2Dvars(), n2d);
  EXPECT_EQ(solver.n3Dvars(), n3d);
  solver.changeHasConstraints(true);
#endif

  vector2.covariant = false;
  EXPECT_NO_THROW(solver.constraint(vector2, vector2, "another_vector"));
  EXPECT_EQ(solver.n2Dvars(), n2d * 2);
  EXPECT_EQ(solver.n3Dvars(), n3d * 2);

  const auto expected_names = std::vector<std::string>{"vector", "another_vector"};
  EXPECT_EQ(solver.listVector2DNames(), expected_names);
}

TEST_F(SolverTest, ConstraintVector3D) {
  Options options;
  FakeSolver solver{&options};

  Vector3D vector1{};
  Vector3D vector2{};
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

  EXPECT_NO_THROW(solver.init());
  EXPECT_THROW(solver.init(), BoutException);
}

TEST_F(SolverTest, NoAddAfterInit) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_NO_THROW(solver.init());

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

  EXPECT_NO_THROW(solver.init());

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

  NiceMock<MockPhysicsModel> model{};

  solver.setModel(&model);

  EXPECT_FALSE(solver.splitOperator());

  model.setSplitOperator();

  EXPECT_TRUE(solver.splitOperator());
}

TEST_F(SolverTest, ResetInternalFields) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_THROW(solver.resetInternalFields(), BoutException);
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

  solver.init();

  constexpr auto globalmesh_nx_no_boundry = SolverTest::nx - 2;
  constexpr auto globalmesh_ny_no_boundry = SolverTest::ny - 2;
  constexpr auto localmesh_nx_no_boundry = localmesh_nx - 2;
  constexpr auto localmesh_ny_no_boundry = localmesh_ny - 2;
  constexpr auto expected_total =
      (globalmesh_nx_no_boundry * globalmesh_ny_no_boundry)
      + (localmesh_nx * localmesh_ny)
      + (localmesh_nx_no_boundry * localmesh_ny_no_boundry * localmesh_nz)
      + (nx * ny * nz);

  EXPECT_EQ(solver.getLocalN(), expected_total);
}

TEST_F(SolverTest, HavePreconditioner) {
  Options options;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};

  solver.setModel(&model);

  EXPECT_FALSE(solver.hasPreconditioner());

  model.setPrecon(&MockPhysicsModel::preconditioner);

  EXPECT_TRUE(solver.hasPreconditioner());
}

TEST_F(SolverTest, RunPreconditioner) {
  Options options;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};

  solver.setModel(&model);
  model.setPrecon(&MockPhysicsModel::preconditioner);

  constexpr auto time = 1.0;
  constexpr auto gamma = 2.0;
  constexpr auto delta = 3.0;
  constexpr auto expected = time + gamma + delta;

  EXPECT_EQ(solver.runPreconditioner(time, gamma, delta), expected);
}

TEST_F(SolverTest, HasJacobian) {
  Options options;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  EXPECT_FALSE(solver.hasJacobian());

  model.setJacobian(&MockPhysicsModel::jacobian);

  EXPECT_TRUE(solver.hasJacobian());
}

TEST_F(SolverTest, RunJacobian) {
  Options options;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};

  solver.setModel(&model);
  model.setJacobian(&MockPhysicsModel::jacobian);

  constexpr auto time = 4.0;
  constexpr auto expected = 4;

  EXPECT_EQ(solver.runJacobian(time), expected);
}

TEST_F(SolverTest, AddMonitor) {
  Options options;
  FakeSolver solver{&options};
  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  StrictMock<MockMonitor> monitor;
  EXPECT_NO_THROW(monitor.setTimestepShim(10.0));
  EXPECT_EQ(monitor.getTimestepShim(), 10.0);

  EXPECT_NO_THROW(solver.addMonitor(&monitor));

  EXPECT_THROW(monitor.setTimestepShim(20.0), BoutException);

  EXPECT_CALL(monitor, call(_, 0.0, 10, 10));
  EXPECT_CALL(monitor, cleanup());
  EXPECT_NO_THROW(solver.call_monitors(0.0, 10, 10));
}

TEST_F(SolverTest, AddMonitorFront) {
  WithQuietOutput quiet{output_error};
  Options options;
  FakeSolver solver{&options};
  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  StrictMock<MockMonitor> monitor1;
  StrictMock<MockMonitor> monitor2;
  EXPECT_NO_THROW(solver.addMonitor(&monitor1, Solver::FRONT));
  EXPECT_NO_THROW(solver.addMonitor(&monitor2, Solver::FRONT));

  // Everything's fine
  EXPECT_CALL(monitor1, call(_, 0.0, 0, nout));
  EXPECT_CALL(monitor2, call(_, 0.0, 0, nout));
  solver.call_monitors(0.0, 0, nout);

  // One monitor signals to quit
  ON_CALL(monitor2, call(_, 5.0, _, _)).WillByDefault(testing::Return(-1));
  EXPECT_CALL(monitor2, call(_, 5.0, 1, nout));
  EXPECT_CALL(monitor1, cleanup());
  EXPECT_CALL(monitor2, cleanup());

  EXPECT_THROW(solver.call_monitors(5.0, 1, nout), BoutException);

  // Last timestep
  EXPECT_CALL(monitor1, call(_, 0.0, nout, nout));
  EXPECT_CALL(monitor2, call(_, 0.0, nout, nout));
  EXPECT_CALL(monitor1, cleanup());
  EXPECT_CALL(monitor2, cleanup());
  solver.call_monitors(0.0, nout, nout);
}

TEST_F(SolverTest, AddMonitorBack) {
  WithQuietOutput quiet{output_error};
  Options options;
  FakeSolver solver{&options};
  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  StrictMock<MockMonitor> monitor1;
  StrictMock<MockMonitor> monitor2;
  EXPECT_NO_THROW(solver.addMonitor(&monitor1, Solver::BACK));
  EXPECT_NO_THROW(solver.addMonitor(&monitor2, Solver::BACK));

  // Everything's fine
  EXPECT_CALL(monitor1, call(_, 0.0, 0, nout));
  EXPECT_CALL(monitor2, call(_, 0.0, 0, nout));
  solver.call_monitors(0.0, 0, nout);

  // One monitor signals to quit
  ON_CALL(monitor1, call(_, 5.0, _, _)).WillByDefault(testing::Return(-1));
  EXPECT_CALL(monitor1, call(_, 5.0, 1, nout));
  EXPECT_CALL(monitor1, cleanup());
  EXPECT_CALL(monitor2, cleanup());
  EXPECT_THROW(solver.call_monitors(5.0, 1, nout), BoutException);

  // Last timestep
  EXPECT_CALL(monitor1, call(_, 0.0, nout, nout));
  EXPECT_CALL(monitor2, call(_, 0.0, nout, nout));
  EXPECT_CALL(monitor1, cleanup());
  EXPECT_CALL(monitor2, cleanup());
  solver.call_monitors(0.0, nout, nout);
}

TEST_F(SolverTest, AddMonitorCheckFrequencies) {
  constexpr int nout_total = nout * 100;

  Options options;
  FakeSolver solver{&options};
  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  StrictMock<MockMonitor> default_timestep;
  StrictMock<MockMonitor> smaller_timestep{timestep / 10};
  StrictMock<MockMonitor> even_smaller_timestep{timestep / 100};
  StrictMock<MockMonitor> larger_timestep{timestep * 2};
  StrictMock<MockMonitor> incompatible_timestep{3.14259};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);
  solver.addMonitor(&even_smaller_timestep);
  solver.addMonitor(&larger_timestep);
  EXPECT_THROW(solver.addMonitor(&incompatible_timestep), BoutException);

  // Everything should trigger on initial timestep
  EXPECT_CALL(default_timestep, call(_, _, 0, nout));
  EXPECT_CALL(smaller_timestep, call(_, _, 0, nout * 10));
  EXPECT_CALL(even_smaller_timestep, call(_, _, 0, nout * 100));
  EXPECT_CALL(larger_timestep, call(_, _, 0, nout / 2));

  // This call is *required* to resolve all the monitor timesteps
  solver.solve(nout, timestep);

  // Everything should trigger when the periods line up (here, at 10%
  // of total internal timesteps)
  EXPECT_CALL(default_timestep, call(_, _, (nout / 10), nout));
  EXPECT_CALL(smaller_timestep, call(_, _, ((nout * 10) / 10), nout * 10));
  EXPECT_CALL(even_smaller_timestep, call(_, _, ((nout * 100) / 10), nout * 100));
  EXPECT_CALL(larger_timestep, call(_, _, ((nout / 2) / 10), nout / 2));
  solver.call_monitors(0.0, (nout_total / 10), nout_total);

  // But on the next internal timestep, only the fastest monitor should trigger.
  EXPECT_CALL(even_smaller_timestep, call(_, _, ((nout * 100) / 10) + 1, nout * 100));
  solver.call_monitors(0.0, (nout_total / 10) + 1, nout_total);

  // It's not possible to add a new monitor with a timestep shorter
  // than the current internal timestep, once we've initialised the solver
  MockMonitor too_small_postinit_timestep{timestep / 1000};
  EXPECT_THROW(solver.addMonitor(&too_small_postinit_timestep), BoutException);

  // But we can add monitors with *longer* timesteps, provided it's
  // multiple of the internal timestep
  StrictMock<MockMonitor> larger_postinit_timestep{timestep * 4};
  solver.addMonitor(&larger_postinit_timestep, Solver::BACK);

  // The *penultimate* timestep should only trigger the fastest monitor
  EXPECT_CALL(even_smaller_timestep, call(_, _, (nout * 100) - 1, nout * 100));
  solver.call_monitors(0.0, nout_total - 1, nout_total);

  // Everything should trigger on the final timestep
  EXPECT_CALL(default_timestep, call(_, _, nout, nout));
  EXPECT_CALL(smaller_timestep, call(_, _, (nout * 10), nout * 10));
  EXPECT_CALL(even_smaller_timestep, call(_, _, (nout * 100), nout * 100));
  EXPECT_CALL(larger_timestep, call(_, _, (nout / 2), nout / 2));
  EXPECT_CALL(larger_postinit_timestep, call(_, _, (nout / 4), nout / 4));

  // The final timestep should also trigger all monitors to cleanup
  EXPECT_CALL(default_timestep, cleanup());
  EXPECT_CALL(smaller_timestep, cleanup());
  EXPECT_CALL(even_smaller_timestep, cleanup());
  EXPECT_CALL(larger_timestep, cleanup());
  EXPECT_CALL(larger_postinit_timestep, cleanup());

  solver.call_monitors(0.0, nout_total, nout_total);
}

TEST_F(SolverTest, RemoveMonitor) {
  Options options;
  FakeSolver solver{&options};

  StrictMock<MockMonitor> monitor1;
  StrictMock<MockMonitor> monitor2;
  solver.addMonitor(&monitor1, Solver::BACK);
  solver.addMonitor(&monitor2, Solver::BACK);

  solver.removeMonitor(&monitor1);

  std::list<FakeSolver::MonitorInfo> expected{{&monitor2, ""}};
  EXPECT_EQ(solver.getMonitors(), expected);

  // Removing same monitor again should be a no-op
  solver.removeMonitor(&monitor1);
  EXPECT_EQ(solver.getMonitors(), expected);
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

  solver.addTimestepMonitor(timestep_monitor1);
  solver.addTimestepMonitor(timestep_monitor2);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);
  EXPECT_CALL(model, timestepMonitor).Times(1);

  EXPECT_EQ(solver.call_timestep_monitors(1., 1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(1., -1.), 1);
  EXPECT_EQ(solver.call_timestep_monitors(-1., -1.), 2);
}

TEST_F(SolverTest, RemoveTimestepMonitor) {
  Options options;
  options["monitor_timestep"] = true;
  FakeSolver solver{&options};

  solver.addTimestepMonitor(timestep_monitor1);
  solver.addTimestepMonitor(timestep_monitor2);

  solver.removeTimestepMonitor(timestep_monitor1);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);
  EXPECT_CALL(model, timestepMonitor).Times(2);

  EXPECT_EQ(solver.call_timestep_monitors(1., 1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(1., -1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(-1., -1.), 2);

  solver.removeTimestepMonitor(timestep_monitor1);

  EXPECT_CALL(model, timestepMonitor).Times(2);

  EXPECT_EQ(solver.call_timestep_monitors(1., 1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(1., -1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(-1., -1.), 2);
}

TEST_F(SolverTest, DontCallTimestepMonitors) {
  Options options;
  FakeSolver solver{&options};

  solver.addTimestepMonitor(timestep_monitor1);
  solver.addTimestepMonitor(timestep_monitor2);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);
  EXPECT_CALL(model, timestepMonitor).Times(0);

  EXPECT_EQ(solver.call_timestep_monitors(1., 1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(1., -1.), 0);
  EXPECT_EQ(solver.call_timestep_monitors(-1., -1.), 0);
}

TEST_F(SolverTest, BasicSolve) {
  Options options;
  FakeSolver solver{&options};

  MockPhysicsModel model{};
  EXPECT_CALL(model, init(false)).Times(1);
  EXPECT_CALL(model, postInit(false)).Times(1);
  solver.setModel(&model);

  EXPECT_CALL(model, rhs(0)).Times(1);

  Options::cleanup();
  solver.solve();

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);
}

TEST_F(SolverTest, GetRunID) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_THROW(solver.getRunID(), BoutException);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);
  solver.solve();

  EXPECT_NO_THROW(solver.getRunID());
  EXPECT_TRUE(uuids::uuid::is_valid_uuid(solver.getRunID()));
}

TEST_F(SolverTest, GetRunRestartFrom) {
  Options options;
  FakeSolver solver{&options};

  EXPECT_THROW(solver.getRunRestartFrom(), BoutException);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);
  solver.solve();

  EXPECT_NO_THROW(solver.getRunRestartFrom());
  // It would be valid if this was a restart
  // But hard to check that case without mocking DataFile
  EXPECT_FALSE(uuids::uuid::is_valid_uuid(solver.getRunRestartFrom()));
}

TEST_F(SolverTest, SolveBadInit) {
  Options options;
  options["fail_init"] = -1;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  EXPECT_THROW(solver.solve(), BoutException);

  EXPECT_TRUE(solver.init_called);
  EXPECT_FALSE(solver.run_called);
}

TEST_F(SolverTest, SolveBadRun) {
  Options options;
  options["fail_run"] = -1;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  Options::cleanup();
  EXPECT_EQ(solver.solve(), -1);

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);
}

TEST_F(SolverTest, SolveThrowRun) {
  WithQuietOutput quiet_error{output_error};

  Options options;
  options["throw_run"] = true;
  FakeSolver solver{&options};

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);
  Options::cleanup();

  EXPECT_THROW(solver.solve(), BoutException);

  EXPECT_TRUE(solver.init_called);
  EXPECT_TRUE(solver.run_called);
}

TEST_F(SolverTest, SolveFixDefaultTimestep) {
  constexpr int nout_total = nout * 100;

  Options options;
  FakeSolver solver{&options};

  StrictMock<MockMonitor> default_timestep;
  StrictMock<MockMonitor> smaller_timestep{timestep / 10};
  StrictMock<MockMonitor> even_smaller_timestep{timestep / 100};
  StrictMock<MockMonitor> larger_timestep{timestep * 2};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);
  solver.addMonitor(&even_smaller_timestep);
  solver.addMonitor(&larger_timestep);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  Options::cleanup();

  EXPECT_CALL(default_timestep, call(_, _, 0, nout));
  EXPECT_CALL(smaller_timestep, call(_, _, 0, nout * 10));
  EXPECT_CALL(even_smaller_timestep, call(_, _, 0, nout * 100));
  EXPECT_CALL(larger_timestep, call(_, _, 0, nout / 2));

  solver.solve(nout, timestep);

  EXPECT_CALL(default_timestep, call(_, _, nout, nout));
  EXPECT_CALL(smaller_timestep, call(_, _, (nout * 10), nout * 10));
  EXPECT_CALL(even_smaller_timestep, call(_, _, (nout * 100), nout * 100));
  EXPECT_CALL(larger_timestep, call(_, _, (nout / 2), nout / 2));

  EXPECT_CALL(default_timestep, cleanup());
  EXPECT_CALL(smaller_timestep, cleanup());
  EXPECT_CALL(even_smaller_timestep, cleanup());
  EXPECT_CALL(larger_timestep, cleanup());
  solver.call_monitors(0.0, nout_total, nout_total);
}

TEST_F(SolverTest, SolveFixDefaultTimestepBad) {
  Options options;
  FakeSolver solver{&options};

  StrictMock<MockMonitor> default_timestep;
  StrictMock<MockMonitor> smaller_timestep{0.1};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  EXPECT_THROW(solver.solve(100, 3.142), BoutException);

  EXPECT_FALSE(solver.init_called);
  EXPECT_FALSE(solver.run_called);
}

TEST_F(SolverTest, SolveFixDefaultTimestepSmaller) {
  Options options;
  FakeSolver solver{&options};

  StrictMock<MockMonitor> default_timestep;
  StrictMock<MockMonitor> larger_timestep{timestep * 10};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&larger_timestep);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  Options::cleanup();

  EXPECT_CALL(default_timestep, call(_, _, 0, nout));
  EXPECT_CALL(larger_timestep, call(_, _, 0, nout / 10));
  solver.solve(nout, timestep);

  EXPECT_CALL(default_timestep, call(_, _, nout, nout));
  EXPECT_CALL(larger_timestep, call(_, _, (nout / 10), nout / 10));
  EXPECT_CALL(default_timestep, cleanup());
  EXPECT_CALL(larger_timestep, cleanup());
  solver.call_monitors(0.0, nout, nout);
}

TEST_F(SolverTest, SolveFixDefaultTimestepLarger) {
  Options options;
  FakeSolver solver{&options};

  StrictMock<MockMonitor> default_timestep;
  StrictMock<MockMonitor> smaller_timestep{timestep / 10};

  solver.addMonitor(&default_timestep);
  solver.addMonitor(&smaller_timestep);

  NiceMock<MockPhysicsModel> model{};
  solver.setModel(&model);

  Options::cleanup();

  EXPECT_CALL(default_timestep, call(_, _, 0, nout));
  EXPECT_CALL(smaller_timestep, call(_, _, 0, nout * 10));
  solver.solve(nout, timestep);

  EXPECT_CALL(default_timestep, call(_, _, nout, nout));
  EXPECT_CALL(smaller_timestep, call(_, _, (nout * 10), (nout * 10)));
  EXPECT_CALL(default_timestep, cleanup());
  EXPECT_CALL(smaller_timestep, cleanup());
  solver.call_monitors(0.0, (nout * 10), (nout * 10));
}
