#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "test_extras.hxx"
#include "bout/solver.hxx"
#include "bout/solverfactory.hxx"

#include <algorithm>
#include <vector>
#include <string>

namespace {
class FakeSolver : public Solver {
public:
  FakeSolver(Options* options) : Solver(options) {
    has_constraints = true;
  }
  ~FakeSolver() = default;
  int run() { return (*options)["number"].withDefault(42); }

  void changeHasConstraints(bool new_value) { has_constraints = new_value; }

  auto listField2DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(f2d), end(f2d), std::back_inserter(result),
                   [](const VarStr<Field2D>& f) { return f.name; });
    return result;
  }

  auto listField3DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(f3d), end(f3d), std::back_inserter(result),
                   [](const VarStr<Field3D>& f) { return f.name; });
    return result;
  }

  auto listVector2DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(v2d), end(v2d), std::back_inserter(result),
                   [](const VarStr<Vector2D>& f) { return f.name; });
    return result;
  }

  auto listVector3DNames() -> std::vector<std::string> {
    std::vector<std::string> result{};
    std::transform(begin(v3d), end(v3d), std::back_inserter(result),
                   [](const VarStr<Vector3D>& f) { return f.name; });
    return result;
  }

  int getLocalNHelper() { return getLocalN(); }
};

RegisterSolver<FakeSolver> register_fake("fake_solver");
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

  rhsfunc fake_rhs = [](BoutReal) -> int { return 0;};
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

  EXPECT_NO_THROW(solver.setMaxTimestep(4.5));
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

  Field2D field1{}, field2{};
  Field3D field3{}, field4{};

  solver.add(field1, "field1");
  solver.add(field2, "field2");
  solver.add(field3, "field3");
  solver.add(field4, "field4");

  solver.init(0, 0);

  static_cast<FakeMesh*>(field1.getMesh())->createBoundaryRegions();

  constexpr auto nx_no_boundry = nx - 2;
  constexpr auto ny_no_boundry = ny - 2;
  constexpr auto expected_total = (nx_no_boundry * ny_no_boundry) + (nx * ny)
                                  + (nx_no_boundry * ny_no_boundry * nz) + (nx * ny * nz);

  EXPECT_EQ(solver.getLocalNHelper(), expected_total);
}
