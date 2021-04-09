#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "field_factory.hxx"
#include "output.hxx"
#include "test_extras.hxx"
#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/traits.hxx"
#include "bout/paralleltransform.hxx"

/// Global mesh
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;
using bout::generator::Context;

// Reuse the "standard" fixture for FakeMesh
template <typename T>
class FieldFactoryCreationTest : public FakeMeshFixture {
public:
  FieldFactoryCreationTest() : FakeMeshFixture() {
    Options options;
    options["input"]["transform_from_field_aligned"] = false;
    factory = FieldFactory{mesh, &options};
  }

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_warn{output_warn};

  FieldFactory factory;

  // We can't just decide which FieldFactory::create?D function to
  // call with
  //
  //     if (bout::utils::is_Field3D<T>::value) {
  //       return factory.create3D(...);
  //     } else {
  //       return factory.create2D(...);
  //     }
  //
  // as this is *runtime* so the compiler still needs to evaluate both
  // branches -- until C++17 when we can use `if constexpr ...`
  template <class... Args>
  T createDispatch(std::true_type, Args&&... args) {
    return factory.create3D(std::forward<Args>(args)...);
  }

  template <class... Args>
  T createDispatch(std::false_type, Args&&... args) {
    return factory.create2D(std::forward<Args>(args)...);
  }

  // Generic way of calling either FieldFactory::create2D or
  // FieldFactory::create3D
  template <class... Args>
  T create(Args&&... args) {
    return createDispatch(bout::utils::is_Field3D<T>{}, std::forward<Args>(args)...);
  }
};

using Fields = ::testing::Types<Field2D, Field3D>;

TYPED_TEST_SUITE(FieldFactoryCreationTest, Fields);

TYPED_TEST(FieldFactoryCreationTest, CreateFromValueGenerator) {
  auto value = BoutReal{4.};
  auto output = this->create(generator(value));

  EXPECT_TRUE(IsFieldEqual(output, value));
}

TYPED_TEST(FieldFactoryCreationTest, CreateFromPointerGenerator) {
  auto value = BoutReal{5.};
  auto output = this->create(generator(&value));

  EXPECT_TRUE(IsFieldEqual(output, value));
}

TYPED_TEST(FieldFactoryCreationTest, CreateFromFunction) {
  FuncPtr function = [](BoutReal, BoutReal x, BoutReal, BoutReal) -> BoutReal {
    return x + 1.;
  };

  auto generator = std::make_shared<FieldFunction>(FieldFunction{function});

  auto output = this->create(generator);

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return index.x() + 1.; },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateNull) {
  FieldNull null{};
  auto output = this->create(null.clone({}));

  EXPECT_TRUE(IsFieldEqual(output, 0.0));
}

TYPED_TEST(FieldFactoryCreationTest, CreatePi) {
  auto output = this->create("pi");

  EXPECT_TRUE(IsFieldEqual(output, PI));
}

TYPED_TEST(FieldFactoryCreationTest, CreatePiSymbol) {
  auto output = this->create("π");

  EXPECT_TRUE(IsFieldEqual(output, PI));
}

TYPED_TEST(FieldFactoryCreationTest, CreateAddition) {
  auto output = this->create("1 + 1");

  EXPECT_TRUE(IsFieldEqual(output, 2.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateSubtraction) {
  auto output = this->create("11 - 1");

  EXPECT_TRUE(IsFieldEqual(output, 10.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateMultiplication) {
  auto output = this->create("3 * 1");

  EXPECT_TRUE(IsFieldEqual(output, 3.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateDivision) {
  auto output = this->create("12 / 2");

  EXPECT_TRUE(IsFieldEqual(output, 6.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateNested) {
  auto output = this->create("1 + (12 / (2 * 3))");

  EXPECT_TRUE(IsFieldEqual(output, 3.));
}

TYPED_TEST(FieldFactoryCreationTest, ParseThenCreateNested) {
  auto generator = this->factory.parse("1 + (12 / (2 * 3))");
  auto output = this->create(generator);

  EXPECT_TRUE(IsFieldEqual(output, 3.));
}

TYPED_TEST(FieldFactoryCreationTest, ParseNull) {
  EXPECT_THROW(this->create(FieldGeneratorPtr{nullptr}), BoutException);
}

TYPED_TEST(FieldFactoryCreationTest, CreateX) {
  auto output = this->create("x");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateY) {
  auto output = this->create("y");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return TWOPI * index.y(); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateZ) {
  auto output = this->create("z");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return TWOPI * index.z() / FieldFactoryCreationTest<TypeParam>::nz;
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateXStaggered) {
  // Need this->mesh_staggered to access member of base FakeMeshFixture because
  // derived FieldFactoryCreationTest is a template clas
  auto output = this->create("x", nullptr, this->mesh_staggered, CELL_XLOW);

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return index.x() - 0.5; },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
  EXPECT_EQ(output.getLocation(), CELL_XLOW);
}

TYPED_TEST(FieldFactoryCreationTest, CreateYStaggered) {
  // Need this->mesh_staggered to access member of base FakeMeshFixture because
  // derived FieldFactoryCreationTest is a template clas
  auto output = this->create("y", nullptr, this->mesh_staggered, CELL_YLOW);

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return TWOPI * (index.y() - 0.5);
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
  EXPECT_EQ(output.getLocation(), CELL_YLOW);
}

TYPED_TEST(FieldFactoryCreationTest, CreateZStaggered) {
  // Need this->mesh_staggered to access member of base FakeMeshFixture because
  // derived FieldFactoryCreationTest is a template clas
  auto output = this->create("z", nullptr, this->mesh_staggered, CELL_ZLOW);

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {

        auto offset = BoutReal{0.0};
        if (bout::utils::is_Field3D<TypeParam>::value) {
          offset = 0.5;
        }

        return TWOPI * (index.z() - offset) / FieldFactoryCreationTest<TypeParam>::nz;
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
  EXPECT_EQ(output.getLocation(), CELL_ZLOW);
}

TYPED_TEST(FieldFactoryCreationTest, CreateT) {
  auto time = BoutReal{77.0};
  auto output = this->create("t", nullptr, nullptr, CELL_CENTRE, time);

  EXPECT_TRUE(IsFieldEqual(output, time));
}

TYPED_TEST(FieldFactoryCreationTest, CreateSinX) {
  auto output = this->create("sin(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::sin(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateCosX) {
  auto output = this->create("cos(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::cos(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateTanX) {
  auto output = this->create("tan(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::tan(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateAsinX) {
  auto output = this->create("asin(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::asin(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateAcosX) {
  auto output = this->create("acos(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::acos(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateAtanX) {
  auto output = this->create("atan(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::atan(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateAtanX2) {
  auto output = this->create("atan(x, 2)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::atan2(index.x(), 2.);
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateSinhX) {
  auto output = this->create("sinh(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::sinh(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateCoshX) {
  auto output = this->create("cosh(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::cosh(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateTanhX) {
  auto output = this->create("tanh(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::tanh(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateExpX) {
  auto output = this->create("exp(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::exp(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateLogX) {
  auto output = this->create("log(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::log(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateGaussX) {
  auto output = this->create("gauss(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::exp(-std::pow(index.x(), 2) / 2.) / std::sqrt(TWOPI);
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateGaussWithWidthX) {
  auto output = this->create("gauss(x, 4.)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        constexpr auto width = BoutReal{4.};
        return std::exp(-std::pow(index.x() / width, 2) / 2.)
               / (std::sqrt(TWOPI) * width);
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateAbsX) {
  auto output = this->create("abs(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::abs(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateSqrtX) {
  auto output = this->create("sqrt(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::sqrt(index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateHeavisideXPi) {
  auto output = this->create("h(x - pi)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        if (index.x() > PI) {
          return 1.0;
        } else {
          return 0.0;
        }
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateErfX) {
  auto output = this->create("erf(x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return std::erf(index.x()); },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateFmodX) {
  auto output = this->create("fmod(pi, x)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::fmod(PI, index.x());
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateMinX) {
  auto output = this->create("min(2, 3, 1, 4)");

  EXPECT_TRUE(IsFieldEqual(output, 1.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateMaxX) {
  auto output = this->create("max(2, 3, 1, 4)");

  EXPECT_TRUE(IsFieldEqual(output, 4.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateClampX) {
  // Check that the first argument is within low and high limits
  // Also check that each input can be an expression
  
  auto output = this->create("clamp(1 + 1, 1, 3)");

  EXPECT_TRUE(IsFieldEqual(output, 2.));

  output = this->create("clamp(-1, 2 - 1, 3)");

  EXPECT_TRUE(IsFieldEqual(output, 1.));
  
  output = this->create("clamp(5, 1, 6 / 2)");

  EXPECT_TRUE(IsFieldEqual(output, 3.));
}

TYPED_TEST(FieldFactoryCreationTest, CreatePowX) {
  auto output = this->create("power(x, 2)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        return std::pow(index.x(), 2);
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateTanhHatX) {
  auto output = this->create("tanhhat(x, 1, 2, 3)");

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal {
        constexpr auto width = BoutReal{1.};
        constexpr auto centre = BoutReal{2.};
        constexpr auto steepness = BoutReal{3.};

        return 0.5
               * (std::tanh(steepness * (index.x() - (centre - 0.5 * width)))
                  - std::tanh(steepness * (index.x() - (centre + 0.5 * width))));
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateRoundX) {
  auto output = this->create("round(pi)");

  EXPECT_TRUE(IsFieldEqual(output, 3.));
}

TYPED_TEST(FieldFactoryCreationTest, CreateWithLookup) {
  auto a_value = int{6};

  auto options = Options{};
  options["a"] = a_value;
  auto output = this->create("x + a", &options);

  auto expected = makeField<TypeParam>(
      [&a_value](typename TypeParam::ind_type& index) -> BoutReal {
        return index.x() + a_value;
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, ParseFromCache) {
  auto a_value = int{6};

  auto options = Options{};
  options["a"] = a_value;
  this->factory.parse("x + a", &options);

  auto output = this->create("x + a");

  auto expected = makeField<TypeParam>(
      [&a_value](typename TypeParam::ind_type& index) -> BoutReal {
        return index.x() + a_value;
      },
      mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TYPED_TEST(FieldFactoryCreationTest, CreateOnMesh) {
  constexpr auto nx = int{1};
  constexpr auto ny = int{1};
  constexpr auto nz = int{1};

  FakeMesh localmesh{nx, ny, nz};
  localmesh.createDefaultRegions();
  localmesh.setCoordinates(std::make_shared<Coordinates>(
      &localmesh, Field2D{1.0}, Field2D{1.0}, BoutReal{1.0}, Field2D{1.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{1.0}, Field2D{1.0}, Field2D{1.0}, Field2D{0.0}, Field2D{0.0}, Field2D{0.0},
      Field2D{0.0}, Field2D{0.0}));
  // No call to Coordinates::geometry() needed here

  localmesh.getCoordinates()->setParallelTransform(
      bout::utils::make_unique<ParallelTransformIdentity>(localmesh));

  auto output = this->create("x", nullptr, &localmesh);

  auto expected = makeField<TypeParam>(
      [](typename TypeParam::ind_type& index) -> BoutReal { return index.x(); },
      &localmesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
  EXPECT_EQ(output.getNx(), nx);
  EXPECT_EQ(output.getNy(), ny);
  EXPECT_EQ(output.getNz(), nz);
}

// The following tests still use the FieldFactory, but don't need to
// be typed and make take longer as they check that exceptions get
// thrown. Doing these twice will slow down the test unnecessarily

class FieldFactoryTest : public FakeMeshFixture {
public:
  FieldFactoryTest() : FakeMeshFixture{}, factory{mesh} {}
  virtual ~FieldFactoryTest() = default;

  WithQuietOutput quiet_info{output_info}, quiet{output}, quiet_error{output_error};
  WithQuietOutput quiet_warn{output_warn};

  FieldFactory factory;
};

TEST_F(FieldFactoryTest, RequireMesh) {
  delete bout::globals::mesh;
  bout::globals::mesh = nullptr;

  FieldFactory local_factory{nullptr, nullptr};

  EXPECT_THROW(local_factory.create2D("x", nullptr, nullptr), BoutException);
  EXPECT_THROW(local_factory.create3D("x", nullptr, nullptr), BoutException);
}

TEST_F(FieldFactoryTest, CreateOnMeshWithoutCoordinates) {
  static_cast<FakeMesh*>(mesh)->setCoordinates(nullptr);
  EXPECT_NO_THROW(factory.create3D("x"));
}

TEST_F(FieldFactoryTest, CleanCache) {
  auto a_value = int{6};

  auto options = Options{};
  options["a"] = a_value;
  factory.parse("x + a", &options);

  factory.cleanCache();

  EXPECT_THROW(factory.parse("x + a"), ParseException);
}

TEST_F(FieldFactoryTest, ParseSelfReference) {
  auto options = Options{};
  options["a"] = "a";

  EXPECT_THROW(factory.parse("a", &options), BoutException);
}

TEST_F(FieldFactoryTest, SinArgs) {
  EXPECT_THROW(factory.parse("sin()"), ParseException);
  EXPECT_THROW(factory.parse("sin(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, CosArgs) {
  EXPECT_THROW(factory.parse("cos()"), ParseException);
  EXPECT_THROW(factory.parse("cos(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, TanArgs) {
  EXPECT_THROW(factory.parse("tan()"), ParseException);
  EXPECT_THROW(factory.parse("tan(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, SinhArgs) {
  EXPECT_THROW(factory.parse("sinh()"), ParseException);
  EXPECT_THROW(factory.parse("sinh(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, CoshArgs) {
  EXPECT_THROW(factory.parse("cosh()"), ParseException);
  EXPECT_THROW(factory.parse("cosh(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, TanhArgs) {
  EXPECT_THROW(factory.parse("tanh()"), ParseException);
  EXPECT_THROW(factory.parse("tanh(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, ASinArgs) {
  EXPECT_THROW(factory.parse("asin()"), ParseException);
  EXPECT_THROW(factory.parse("asin(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, ACosArgs) {
  EXPECT_THROW(factory.parse("acos()"), ParseException);
  EXPECT_THROW(factory.parse("acos(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, ATanArgs) {
  EXPECT_THROW(factory.parse("atan()"), ParseException);
  EXPECT_THROW(factory.parse("atan(x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, ExpArgs) {
  EXPECT_THROW(factory.parse("exp()"), ParseException);
  EXPECT_THROW(factory.parse("exp(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, LogArgs) {
  EXPECT_THROW(factory.parse("log()"), ParseException);
  EXPECT_THROW(factory.parse("log(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, GaussArgs) {
  EXPECT_THROW(factory.parse("gauss()"), ParseException);
  EXPECT_THROW(factory.parse("gauss(x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, AbsArgs) {
  EXPECT_THROW(factory.parse("abs()"), ParseException);
  EXPECT_THROW(factory.parse("abs(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, SqrtArgs) {
  EXPECT_THROW(factory.parse("sqrt()"), ParseException);
  EXPECT_THROW(factory.parse("sqrt(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, HeavisideArgs) {
  EXPECT_THROW(factory.parse("h()"), ParseException);
  EXPECT_THROW(factory.parse("h(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, ErfArgs) {
  EXPECT_THROW(factory.parse("erf()"), ParseException);
  EXPECT_THROW(factory.parse("erf(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, FmodArgs) {
  EXPECT_THROW(factory.parse("fmod()"), ParseException);
  EXPECT_THROW(factory.parse("fmod(x)"), ParseException);
  EXPECT_THROW(factory.parse("fmod(x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, MinArgs) {
  EXPECT_THROW(factory.parse("min()"), ParseException);
}

TEST_F(FieldFactoryTest, MaxArgs) {
  EXPECT_THROW(factory.parse("max()"), ParseException);
}

TEST_F(FieldFactoryTest, ClampArgs) {
  EXPECT_THROW(factory.parse("clamp()"), ParseException);
  EXPECT_THROW(factory.parse("clamp(1)"), ParseException);
  EXPECT_THROW(factory.parse("clamp(1,2)"), ParseException);
}

TEST_F(FieldFactoryTest, PowerArgs) {
  EXPECT_THROW(factory.parse("power()"), ParseException);
  EXPECT_THROW(factory.parse("power(x)"), ParseException);
  EXPECT_THROW(factory.parse("power(x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, RoundArgs) {
  EXPECT_THROW(factory.parse("round()"), ParseException);
  EXPECT_THROW(factory.parse("round(x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, BallooningArgs) {
  EXPECT_THROW(factory.parse("ballooning()"), ParseException);
  EXPECT_THROW(factory.parse("ballooning(x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, MixmodeArgs) {
  EXPECT_THROW(factory.parse("mixmode()"), ParseException);
  EXPECT_THROW(factory.parse("mixmode(x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, TanhhatArgs) {
  EXPECT_THROW(factory.parse("tanhhat()"), ParseException);
  EXPECT_THROW(factory.parse("tanhhat(x, x, x, x, x)"), ParseException);
}

TEST_F(FieldFactoryTest, Where) {
  auto fieldgen = factory.parse("where({val}, 3, 5)");

  EXPECT_DOUBLE_EQ(fieldgen->generate(Context().set("val", 1.0)), 3);
  EXPECT_DOUBLE_EQ(fieldgen->generate(Context().set("val", -1.0)), 5);
}

TEST_F(FieldFactoryTest, Recursion) {
  // Need to enable recursion
  Options opt;
  opt["input"]["max_recursion_depth"] = 4; // Should be sufficient for n=6

  // Create a factory with a max_recursion_depth != 0
  FieldFactory factory_rec(nullptr, &opt);
  
  // Fibonacci sequence: 1 1 2 3 5 8
  opt["fib"] = "where({n} - 2.5, [n={n}-1](fib) + [n={n}-2](fib), 1)";
  
  auto gen = factory_rec.parse("fib", &opt);
  EXPECT_DOUBLE_EQ(gen->generate(Context().set("n", 3)), 2);
  EXPECT_DOUBLE_EQ(gen->generate(Context().set("n", 4)), 3);
  EXPECT_DOUBLE_EQ(gen->generate(Context().set("n", 5)), 5);
  EXPECT_DOUBLE_EQ(gen->generate(Context().set("n", 6)), 8);
  EXPECT_THROW(gen->generate(Context().set("n", 7)), BoutException); // Max recursion exceeded
}
// A mock ParallelTransform to test transform_from_field_aligned
// property of FieldFactory. For now, the transform just returns the
// negative of the input. Ideally, this will get moved to GoogleMock
// when we start using it.
//
// Can turn off the ability to do the transform. Should still be valid
class MockParallelTransform : public ParallelTransform {
public:
  MockParallelTransform(Mesh& mesh, bool allow_transform_)
      : ParallelTransform(mesh), allow_transform(allow_transform_) {}
  ~MockParallelTransform() = default;

  void calcParallelSlices(Field3D&) override {}

  bool canToFromFieldAligned() override { return allow_transform; }

  bool requiresTwistShift(bool, YDirectionType) override { return false; }

  void checkInputGrid() override {}

  Field3D fromFieldAligned(const Field3D& f, const std::string&) override {
    if (f.getDirectionY() != YDirectionType::Aligned) {
      throw BoutException("Unaligned field passed to fromFieldAligned");
    }
    return -f;
  }

  FieldPerp fromFieldAligned(const FieldPerp& f, const std::string&) override {
    if (f.getDirectionY() != YDirectionType::Aligned) {
      throw BoutException("Unaligned field passed to fromFieldAligned");
    }
    return -f;
  }

  Field3D toFieldAligned(const Field3D& f, const std::string&) override {
    if (f.getDirectionY() != YDirectionType::Standard) {
      throw BoutException("Aligned field passed to toFieldAligned");
    }
    return -f;
  }
  FieldPerp toFieldAligned(const FieldPerp& f, const std::string&) override {
    if (f.getDirectionY() != YDirectionType::Standard) {
      throw BoutException("Aligned field passed to toFieldAligned");
    }
    return -f;
  }

private:
  const bool allow_transform;
};

class FieldFactoryCreateAndTransformTest : public FakeMeshFixture {
public:
  FieldFactoryCreateAndTransformTest() : FakeMeshFixture{} {
    // We need Coordinates so a parallel transform is available as
    // FieldFactory::create3D wants to un-field-align the result
    static_cast<FakeMesh*>(mesh)->setCoordinates(test_coords);
  }

  WithQuietOutput quiet_info{output_info};
  WithQuietOutput quiet_warn{output_warn};
};

TEST_F(FieldFactoryCreateAndTransformTest, Create2D) {
  mesh->getCoordinates()->setParallelTransform(
    bout::utils::make_unique<MockParallelTransform>(*mesh, true));

  FieldFactory factory;

  auto output = factory.create2D("x");

  // Field2Ds can't be transformed, so expect no change
  auto expected = makeField<Field2D>(
      [](typename Field2D::ind_type& index) -> BoutReal { return index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TEST_F(FieldFactoryCreateAndTransformTest, Create3D) {
  mesh->getCoordinates()->setParallelTransform(
    bout::utils::make_unique<MockParallelTransform>(*mesh, true));

  FieldFactory factory;

  auto output = factory.create3D("x");

  auto expected = makeField<Field3D>(
      [](typename Field3D::ind_type& index) -> BoutReal { return -index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TEST_F(FieldFactoryCreateAndTransformTest, Create2DNoTransform) {
  mesh->getCoordinates()->setParallelTransform(
    bout::utils::make_unique<MockParallelTransform>(*mesh, true));

  Options options;
  options["input"]["transform_from_field_aligned"] = false;
  FieldFactory factory{mesh, &options};

  auto output = factory.create2D("x");

  // Field2Ds can't be transformed, so expect no change
  auto expected = makeField<Field2D>(
      [](typename Field2D::ind_type& index) -> BoutReal { return index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TEST_F(FieldFactoryCreateAndTransformTest, Create3DNoTransform) {
  mesh->getCoordinates()->setParallelTransform(
    bout::utils::make_unique<MockParallelTransform>(*mesh, true));

  Options options;
  options["input"]["transform_from_field_aligned"] = false;
  FieldFactory factory{mesh, &options};

  auto output = factory.create3D("x");

  auto expected = makeField<Field3D>(
      [](typename Field3D::ind_type& index) -> BoutReal { return index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TEST_F(FieldFactoryCreateAndTransformTest, Create2DCantTransform) {
  mesh->getCoordinates()->setParallelTransform(
    bout::utils::make_unique<MockParallelTransform>(*mesh, false));

  FieldFactory factory{mesh};

  auto output = factory.create2D("x");

  // Field2Ds can't be transformed, so expect no change
  auto expected = makeField<Field2D>(
      [](typename Field2D::ind_type& index) -> BoutReal { return index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}

TEST_F(FieldFactoryCreateAndTransformTest, Create3DCantTransform) {
  mesh->getCoordinates()->setParallelTransform(
    bout::utils::make_unique<MockParallelTransform>(*mesh, false));

  FieldFactory factory{mesh};

  auto output = factory.create3D("x");

  auto expected = makeField<Field3D>(
      [](typename Field3D::ind_type& index) -> BoutReal { return index.x(); }, mesh);

  EXPECT_TRUE(IsFieldEqual(output, expected));
}
