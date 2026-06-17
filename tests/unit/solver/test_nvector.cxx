#include "bout/build_defines.hxx"

#if BOUT_HAS_SUNDIALS_MANYVECTOR

#include <cmath>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutcomm.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/sundials_backports.hxx"
#include "gtest/gtest.h"

#include "../../../src/solver/nvector.hxx"

namespace {

struct NVectorDestroy {
  void operator()(N_Vector vec) const {
    if (vec != nullptr) {
      N_VDestroy(vec);
    }
  }
};

using UniqueNVector = std::unique_ptr<std::remove_pointer_t<N_Vector>, NVectorDestroy>;

auto makeNVector(N_Vector vec) -> UniqueNVector { return UniqueNVector(vec); }

template <typename FieldType>
auto makePositiveField(BoutReal offset = 0.0) -> FieldType {
  return makeField<FieldType>(
      [offset](typename FieldType::ind_type& i) {
        return offset + 1.0 + static_cast<BoutReal>((i.ind % 7) + 1);
      },
      bout::globals::mesh);
}

template <typename FieldType>
auto makeSignedField() -> FieldType {
  return makeField<FieldType>(
      [](typename FieldType::ind_type& i) {
        const auto magnitude = 1.0 + static_cast<BoutReal>(i.ind % 5);
        return (i.ind % 2 == 0) ? magnitude : -magnitude;
      },
      bout::globals::mesh);
}

template <typename FieldType>
auto regionSize(const FieldType& field, const std::string& region) -> sunindextype {
  return static_cast<sunindextype>(field.getRegion(region).size());
}

template <typename FieldType, typename UnaryOp>
auto sumField(const FieldType& field, const std::string& region, UnaryOp&& op)
    -> BoutReal {
  BoutReal result = 0.0;
  for (const auto& i : field.getRegion(region)) {
    result += op(field[i]);
  }
  return result;
}

template <typename FieldType, typename BinaryOp>
auto sumFields(const FieldType& lhs, const FieldType& rhs, const std::string& region,
               BinaryOp&& op) -> BoutReal {
  BoutReal result = 0.0;
  for (const auto& i : lhs.getRegion(region)) {
    result += op(lhs[i], rhs[i]);
  }
  return result;
}

template <typename T>
class BoutNVectorTest : public FakeMeshFixture {
public:
  BoutNVectorTest() : sunctx(createSUNContext(BoutComm::get())) {}

protected:
  sundials::Context sunctx;
};

using FieldTypes = ::testing::Types<Field2D, Field3D>;
TYPED_TEST_SUITE(BoutNVectorTest, FieldTypes);

TYPED_TEST(BoutNVectorTest, ConstructorAliasesFieldAndReportsLength) {
  auto field = makePositiveField<TypeParam>();
  auto vec = makeNVector(BoutNVector::create(this->sunctx, field, true));

  ASSERT_NE(vec.get(), nullptr);
  EXPECT_EQ(N_VGetVectorID(vec.get()), SUNDIALS_NVEC_CUSTOM);
  EXPECT_EQ(N_VGetLength(vec.get()), regionSize(field, "RGN_ALL"));
  EXPECT_EQ(&BoutNVector::get<TypeParam>(vec.get()), &field);

  field = 3.25;
  EXPECT_TRUE(IsFieldEqual(BoutNVector::get<TypeParam>(vec.get()), 3.25));
}

TYPED_TEST(BoutNVectorTest, ConstructorWithoutBoundaryUsesNoBoundaryLength) {
  TypeParam field{bout::globals::mesh};
  field.allocate();
  field = 1.0;
  for (const auto& i : field.getRegion("RGN_BNDRY")) {
    field[i] = 100.0;
  }

  auto vec = makeNVector(BoutNVector::create(this->sunctx, field, false));

  ASSERT_NE(vec.get(), nullptr);
  EXPECT_EQ(N_VGetLength(vec.get()), regionSize(field, "RGN_NOBNDRY"));
  EXPECT_DOUBLE_EQ(N_VL1Norm(vec.get()),
                   static_cast<BoutReal>(regionSize(field, "RGN_NOBNDRY")));
}

TYPED_TEST(BoutNVectorTest, LinearAndPointwiseOperations) {
  auto x = makePositiveField<TypeParam>(0.0);
  auto y = makePositiveField<TypeParam>(20.0);
  TypeParam z{bout::globals::mesh};
  z.allocate();
  z = 0.0;

  auto vx = makeNVector(BoutNVector::create(this->sunctx, x, true));
  auto vy = makeNVector(BoutNVector::create(this->sunctx, y, true));
  auto vz = makeNVector(BoutNVector::create(this->sunctx, z, true));

  N_VLinearSum(2.0, vx.get(), -1.0, vy.get(), vz.get());
  auto expected = 2.0 * x - y;
  EXPECT_TRUE(IsFieldEqual(z, expected));

  N_VProd(vx.get(), vy.get(), vz.get());
  expected = x * y;
  EXPECT_TRUE(IsFieldEqual(z, expected));

  N_VDiv(vy.get(), vx.get(), vz.get());
  expected = y / x;
  EXPECT_TRUE(IsFieldEqual(z, expected));

  N_VScale(-0.5, vx.get(), vz.get());
  expected = -0.5 * x;
  EXPECT_TRUE(IsFieldEqual(z, expected));

  auto signed_field = makeSignedField<TypeParam>();
  auto vsigned = makeNVector(BoutNVector::create(this->sunctx, signed_field, true));

  N_VAbs(vsigned.get(), vz.get());
  expected = abs(signed_field);
  EXPECT_TRUE(IsFieldEqual(z, expected));

  N_VInv(vx.get(), vz.get());
  expected = 1.0 / x;
  EXPECT_TRUE(IsFieldEqual(z, expected));

  N_VAddConst(vx.get(), 3.5, vz.get());
  expected = 3.5 + x;
  EXPECT_TRUE(IsFieldEqual(z, expected));

  N_VConst(-2.0, vz.get());
  EXPECT_TRUE(IsFieldEqual(z, -2.0));
}

TYPED_TEST(BoutNVectorTest, ReductionOperations) {
  auto x = makeSignedField<TypeParam>();
  auto weights = makePositiveField<TypeParam>(5.0);
  auto vx = makeNVector(BoutNVector::create(this->sunctx, x, true));
  auto vweights = makeNVector(BoutNVector::create(this->sunctx, weights, true));

  const auto dot_expected = sumFields(
      x, weights, "RGN_ALL", [](BoutReal lhs, BoutReal rhs) { return lhs * rhs; });
  EXPECT_DOUBLE_EQ(N_VDotProd(vx.get(), vweights.get()), dot_expected);

  const auto l1_expected =
      sumField(x, "RGN_ALL", [](BoutReal value) { return std::abs(value); });
  EXPECT_DOUBLE_EQ(N_VL1Norm(vx.get()), l1_expected);

  const auto wl2_sq_expected =
      sumFields(x, weights, "RGN_ALL",
                [](BoutReal lhs, BoutReal rhs) { return std::pow(lhs * rhs, 2); });
  const auto wl2_expected = std::sqrt(wl2_sq_expected);
  EXPECT_NEAR(N_VWL2Norm(vx.get(), vweights.get()), wl2_expected, BoutRealTolerance);

  const auto wrms_expected =
      wl2_expected / std::sqrt(static_cast<BoutReal>(regionSize(x, "RGN_ALL")));
  EXPECT_NEAR(N_VWrmsNorm(vx.get(), vweights.get()), wrms_expected, BoutRealTolerance);

  EXPECT_DOUBLE_EQ(N_VMaxNorm(vx.get()), max(abs(x), true));
  EXPECT_DOUBLE_EQ(N_VMin(vx.get()), min(x, true));
}

TYPED_TEST(BoutNVectorTest, CloneCreatesIndependentStorage) {
  auto field = makePositiveField<TypeParam>();
  auto vec = makeNVector(BoutNVector::create(this->sunctx, field, true));
  auto clone = makeNVector(N_VClone(vec.get()));

  ASSERT_NE(clone.get(), nullptr);
  EXPECT_EQ(N_VGetVectorID(clone.get()), SUNDIALS_NVEC_CUSTOM);
  EXPECT_EQ(N_VGetLength(clone.get()), N_VGetLength(vec.get()));
  EXPECT_NE(&BoutNVector::get<TypeParam>(clone.get()),
            &BoutNVector::get<TypeParam>(vec.get()));

  N_VConst(42.0, clone.get());
  EXPECT_TRUE(IsFieldEqual(BoutNVector::get<TypeParam>(clone.get()), 42.0));
  EXPECT_FALSE(IsFieldEqual(field, 42.0));
}

TYPED_TEST(BoutNVectorTest, SwapExchangesFieldStorage) {
  auto field = makePositiveField<TypeParam>(0.0);
  auto other = makePositiveField<TypeParam>(100.0);
  const auto original_field = field;
  const auto original_other = other;

  auto vec = makeNVector(BoutNVector::create(this->sunctx, field, true));

  BoutNVector::swap(vec.get(), other);

  EXPECT_TRUE(IsFieldEqual(field, original_other));
  EXPECT_TRUE(IsFieldEqual(other, original_field));
  EXPECT_TRUE(IsFieldEqual(BoutNVector::get<TypeParam>(vec.get()), original_other));
}

#if BOUT_HAS_SUNDIALS_MANYVECTOR
TEST_F(FakeMeshFixture, ManyVectorGetAndSwapSubvectors) {
  const auto sunctx = createSUNContext(BoutComm::get());
  auto field2d = makePositiveField<Field2D>(0.0);
  auto field3d = makePositiveField<Field3D>(50.0);
  auto replacement2d = makePositiveField<Field2D>(100.0);
  const auto original2d = field2d;
  const auto original3d = field3d;
  const auto original_replacement2d = replacement2d;

  std::vector<N_Vector> subvectors;
  subvectors.emplace_back(BoutNVector::create(sunctx, field2d, true));
  subvectors.emplace_back(BoutNVector::create(sunctx, field3d, true));

  auto manyvector = makeNVector(BoutNVector::create(sunctx, subvectors));

  ASSERT_NE(manyvector.get(), nullptr);
  EXPECT_EQ(&BoutNVector::get<Field2D>(manyvector.get(), 0), &field2d);
  EXPECT_EQ(&BoutNVector::get<Field3D>(manyvector.get(), 1), &field3d);
  EXPECT_EQ(N_VGetLength(manyvector.get()),
            regionSize(field2d, "RGN_ALL") + regionSize(field3d, "RGN_ALL"));

  BoutNVector::swap(manyvector.get(), replacement2d, 0);

  EXPECT_TRUE(IsFieldEqual(field2d, original_replacement2d));
  EXPECT_TRUE(IsFieldEqual(replacement2d, original2d));
  EXPECT_TRUE(IsFieldEqual(BoutNVector::get<Field2D>(manyvector.get(), 0),
                           original_replacement2d));
  EXPECT_TRUE(IsFieldEqual(BoutNVector::get<Field3D>(manyvector.get(), 1), original3d));
}
#endif

} // namespace

#endif
