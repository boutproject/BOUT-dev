#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include <vector>

#include "gtest/gtest.h"

#include <bout/petsc_interface.hxx>

#include "fake_mesh_fixture.hxx"

#include "../../../src/solver/impls/snes/snes.hxx"

namespace {

auto makeVec(const std::vector<PetscScalar>& values) -> bout::petsc::UniqueVec {
  bout::petsc::UniqueVec vec{new Vec{nullptr}};
  EXPECT_EQ(
      VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(values.size()), vec.get()),
      PETSC_SUCCESS);
  PetscScalar* data = nullptr;
  EXPECT_EQ(VecGetArray(*vec, &data), PETSC_SUCCESS);
  for (PetscInt i = 0; i < static_cast<PetscInt>(values.size()); ++i) {
    data[i] = values[static_cast<std::size_t>(i)];
  }
  EXPECT_EQ(VecRestoreArray(*vec, &data), PETSC_SUCCESS);
  return vec;
}

auto getVecValues(Vec vec) -> std::vector<PetscScalar> {
  const PetscScalar* data = nullptr;
  PetscInt size = 0;
  EXPECT_EQ(VecGetLocalSize(vec, &size), PETSC_SUCCESS);
  EXPECT_EQ(VecGetArrayRead(vec, &data), PETSC_SUCCESS);
  std::vector<PetscScalar> values(static_cast<std::size_t>(size));
  for (PetscInt i = 0; i < size; ++i) {
    values[static_cast<std::size_t>(i)] = data[static_cast<std::size_t>(i)];
  }
  EXPECT_EQ(VecRestoreArrayRead(vec, &data), PETSC_SUCCESS);
  return values;
}

void expectVecValues(Vec vec, const std::vector<PetscScalar>& expected) {
  EXPECT_EQ(getVecValues(vec), expected);
}

void expectVecValuesNear(Vec vec, const std::vector<PetscScalar>& expected,
                         PetscReal tolerance = 1e-12) {
  const auto actual = getVecValues(vec);
  ASSERT_EQ(actual.size(), expected.size());
  for (std::size_t i = 0; i < actual.size(); ++i) {
    EXPECT_NEAR(actual[i], expected[i], tolerance);
  }
}

class PredictorTest : public FakeMeshFixture {
public:
  PredictorTest() : predictor(options) {}

protected:
  Options options;
  bout::Predictor predictor;
};

TEST_F(PredictorTest, ConstantPredictReturnsMostRecentState) {
  auto state = makeVec({1.0, -2.0, 3.5});
  auto out = makeVec({0.0, 0.0, 0.0});
  Vec out_vec = *out;

  predictor.push_state(1.5, *state);
  predictor.predict(bout::SnesPredictor::constant, 2.0, out_vec);

  expectVecValues(*out, {1.0, -2.0, 3.5});
}

TEST_F(PredictorTest, LinearPredictFallsBackToConstantWithOneState) {
  auto state = makeVec({2.0, 4.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = *out;

  predictor.push_state(3.0, *state);
  predictor.predict(bout::SnesPredictor::linear, 5.0, out_vec);

  expectVecValues(*out, {2.0, 4.0});
}

TEST_F(PredictorTest, LinearPredictExtrapolatesFromTwoStates) {
  auto older = makeVec({1.0, 5.0});
  auto newer = makeVec({3.0, 9.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = *out;

  predictor.push_state(1.0, *older);
  predictor.push_state(2.0, *newer);
  predictor.predict(bout::SnesPredictor::linear, 2.5, out_vec);

  expectVecValues(*out, {4.0, 11.0});
}

TEST_F(PredictorTest, LinearPredictInterpolatesInsideInterval) {
  auto older = makeVec({2.0, 10.0});
  auto newer = makeVec({6.0, 14.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *older);
  predictor.push_state(2.0, *newer);
  predictor.predict(bout::SnesPredictor::linear, 1.0, out_vec);

  expectVecValues(*out, {4.0, 12.0});
}

TEST_F(PredictorTest, PushStateUsesTwoMostRecentStatesForLinearPrediction) {
  auto oldest = makeVec({0.0});
  auto middle = makeVec({1.0});
  auto newest = makeVec({4.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *oldest);
  predictor.push_state(1.0, *middle);
  predictor.push_state(2.0, *newest);
  predictor.predict(bout::SnesPredictor::linear, 3.0, out_vec);

  expectVecValues(*out, {7.0});
}

TEST_F(PredictorTest, SetDefaultChangesPredictBehavior) {
  auto older = makeVec({2.0});
  auto newer = makeVec({5.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *older);
  predictor.push_state(1.0, *newer);

  predictor.setDefault(bout::SnesPredictor::constant);
  predictor.predict(2.0, out_vec);
  expectVecValues(*out, {5.0});

  predictor.setDefault(bout::SnesPredictor::linear);
  predictor.predict(2.0, out_vec);
  expectVecValues(*out, {8.0});
}

TEST_F(PredictorTest, RescaleSkipsUnallocatedHistorySlots) {
  auto norms = makeVec({2.0, 4.0});

  EXPECT_NO_THROW(predictor.rescale(*norms));
}

TEST_F(PredictorTest, RescaleAppliesToAllAllocatedHistoryStates) {
  auto older = makeVec({4.0, 12.0});
  auto newer = makeVec({8.0, 20.0});
  auto norms = makeVec({2.0, 4.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *older);
  predictor.push_state(1.0, *newer);
  predictor.rescale(*norms);

  predictor.predict(bout::SnesPredictor::constant, 2.0, out_vec);
  expectVecValues(*out, {4.0, 5.0});

  predictor.predict(bout::SnesPredictor::linear, 0.0, out_vec);
  expectVecValues(*out, {2.0, 3.0});
}

TEST_F(PredictorTest, RescalePreservesSubsequentLinearPrediction) {
  auto older = makeVec({2.0, 8.0});
  auto newer = makeVec({6.0, 12.0});
  auto norms = makeVec({2.0, 4.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = *out;

  predictor.push_state(1.0, *older);
  predictor.push_state(3.0, *newer);
  predictor.rescale(*norms);
  predictor.predict(bout::SnesPredictor::linear, 4.0, out_vec);

  expectVecValues(*out, {4.0, 3.5});
}

TEST_F(PredictorTest, BeRationalFallsBackToLinearWithTwoStates) {
  auto older = makeVec({8.0});
  auto newer = makeVec({4.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *older);
  predictor.push_state(1.0, *newer);
  predictor.predict(bout::SnesPredictor::be_rational, 2.0, out_vec);

  expectVecValues(*out, {0.0});
}

TEST_F(PredictorTest, BeRationalPredictsMonotoneBackwardEulerDecay) {
  auto oldest = makeVec({8.0});
  auto middle = makeVec({4.0});
  auto newest = makeVec({2.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *oldest);
  predictor.push_state(1.0, *middle);
  predictor.push_state(2.0, *newest);
  predictor.predict(bout::SnesPredictor::be_rational, 3.0, out_vec);

  expectVecValuesNear(*out, {1.0});
}

TEST_F(PredictorTest, BeRationalFallsBackToConstantOnSignChange) {
  auto oldest = makeVec({2.0});
  auto middle = makeVec({1.0});
  auto newest = makeVec({2.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *oldest);
  predictor.push_state(1.0, *middle);
  predictor.push_state(2.0, *newest);
  predictor.predict(bout::SnesPredictor::be_rational, 3.0, out_vec);

  expectVecValues(*out, {2.0});
}

TEST_F(PredictorTest, BeRationalFallsBackToConstantBelowNoiseThreshold) {
  auto oldest = makeVec({0.9999990});
  auto middle = makeVec({0.9999995});
  auto newest = makeVec({1.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *oldest);
  predictor.push_state(1.0, *middle);
  predictor.push_state(2.0, *newest);
  predictor.predict(bout::SnesPredictor::be_rational, 3.0, out_vec);

  expectVecValuesNear(*out, {1.0});
}

TEST(PredictorOptionsTest, BlendUpdatesSubsequentPrediction) {
  Options options;
  options["blend"] = true;
  options["blend_weight"] = 1.0;
  options["blend_initial"] = 1.0;
  bout::Predictor predictor(options);

  auto older = makeVec({4.0});
  auto newer = makeVec({2.0});
  auto converged = makeVec({1.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(1.0, *older);
  predictor.push_state(2.0, *newer);
  predictor.predict(bout::SnesPredictor::linear, 3.0, out_vec);
  predictor.push_state(3.0, *converged);
  predictor.predict(bout::SnesPredictor::linear, 4.0, out_vec);

  expectVecValuesNear(*out, {0.5});
}

TEST(PredictorOptionsTest, ConstantPredictionDoesNotUpdateBlend) {
  Options options;
  options["blend"] = true;
  options["blend_weight"] = 1.0;
  options["blend_initial"] = 0.25;
  bout::Predictor predictor(options);

  auto older = makeVec({4.0});
  auto newer = makeVec({2.0});
  auto converged = makeVec({1.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(1.0, *older);
  predictor.push_state(2.0, *newer);
  predictor.predict(bout::SnesPredictor::constant, 3.0, out_vec);
  predictor.push_state(3.0, *converged);
  predictor.predict(bout::SnesPredictor::linear, 4.0, out_vec);

  expectVecValuesNear(*out, {0.75});
}

TEST(PredictorOptionsTest, BlendUpdateIsClampedByBlendMax) {
  Options options;
  options["blend"] = true;
  options["blend_weight"] = 1.0;
  options["blend_initial"] = 1.0;
  options["blend_max"] = 0.25;
  bout::Predictor predictor(options);

  auto older = makeVec({4.0});
  auto newer = makeVec({2.0});
  auto converged = makeVec({1.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(1.0, *older);
  predictor.push_state(2.0, *newer);
  predictor.predict(bout::SnesPredictor::linear, 3.0, out_vec);
  predictor.push_state(3.0, *converged);
  predictor.predict(bout::SnesPredictor::linear, 4.0, out_vec);

  expectVecValuesNear(*out, {0.75});
}

namespace {

auto nextPredictionWithShadow(bool enable_shadow) -> PetscScalar {
  Options options;
  options["blend"] = true;
  options["blend_weight"] = 1.0;
  options["blend_initial"] = 1.0;
  options["shadow"] = enable_shadow;
  bout::Predictor predictor(options);

  auto oldest = makeVec({8.0});
  auto middle = makeVec({4.0});
  auto newest = makeVec({2.0});
  auto converged = makeVec({1.0});
  auto out = makeVec({0.0});
  Vec out_vec = *out;

  predictor.push_state(0.0, *oldest);
  predictor.push_state(1.0, *middle);
  predictor.push_state(2.0, *newest);
  predictor.predict(bout::SnesPredictor::linear, 3.0, out_vec);
  predictor.push_state(3.0, *converged);
  predictor.predict(bout::SnesPredictor::linear, 4.0, out_vec);

  return getVecValues(*out).front();
}

} // namespace

TEST(PredictorOptionsTest, ShadowPredictionsDoNotInterfereWithBlendUpdate) {
  EXPECT_NEAR(nextPredictionWithShadow(false), 0.5, 1e-12);
  EXPECT_NEAR(nextPredictionWithShadow(true), 0.5, 1e-12);
}

} // namespace

#endif
