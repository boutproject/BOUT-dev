#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"

#include "../../../src/solver/impls/snes/snes.hxx"

namespace {

struct PetscVecDestroy {
  void operator()(std::remove_pointer_t<Vec>* vec) const {
    if (vec != nullptr) {
      ::VecDestroy(&vec);
    }
  }
};

using UniqueVec = std::unique_ptr<std::remove_pointer_t<Vec>, PetscVecDestroy>;

auto makeVec(const std::vector<PetscScalar>& values) -> UniqueVec {
  Vec vec{nullptr};
  EXPECT_EQ(VecCreateSeq(PETSC_COMM_SELF, static_cast<PetscInt>(values.size()), &vec),
            PETSC_SUCCESS);
  PetscScalar* data = nullptr;
  EXPECT_EQ(VecGetArray(vec, &data), PETSC_SUCCESS);
  for (PetscInt i = 0; i < static_cast<PetscInt>(values.size()); ++i) {
    data[i] = values[static_cast<std::size_t>(i)];
  }
  EXPECT_EQ(VecRestoreArray(vec, &data), PETSC_SUCCESS);
  return UniqueVec(vec);
}

auto getVecValues(Vec vec) -> std::vector<PetscScalar> {
  const PetscScalar* data = nullptr;
  PetscInt size = 0;
  EXPECT_EQ(VecGetLocalSize(vec, &size), PETSC_SUCCESS);
  EXPECT_EQ(VecGetArrayRead(vec, &data), PETSC_SUCCESS);
  std::vector<PetscScalar> values(static_cast<std::size_t>(size));
  for (PetscInt i = 0; i < size; ++i) {
    values[static_cast<std::size_t>(i)] = data[i];
  }
  EXPECT_EQ(VecRestoreArrayRead(vec, &data), PETSC_SUCCESS);
  return values;
}

void expectVecValues(Vec vec, const std::vector<PetscScalar>& expected) {
  EXPECT_EQ(getVecValues(vec), expected);
}

class PredictorTest : public FakeMeshFixture {
public:
  PredictorTest() : predictor(options) {}

protected:
  Options options;
  Predictor predictor;
};

TEST_F(PredictorTest, ConstantPredictReturnsMostRecentState) {
  auto state = makeVec({1.0, -2.0, 3.5});
  auto out = makeVec({0.0, 0.0, 0.0});
  Vec out_vec = out.get();

  predictor.push_state(1.5, state.get());
  predictor.predict(BoutSnesPredictor::constant, 2.0, out_vec);

  expectVecValues(out.get(), {1.0, -2.0, 3.5});
}

TEST_F(PredictorTest, LinearPredictFallsBackToConstantWithOneState) {
  auto state = makeVec({2.0, 4.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = out.get();

  predictor.push_state(3.0, state.get());
  predictor.predict(BoutSnesPredictor::linear, 5.0, out_vec);

  expectVecValues(out.get(), {2.0, 4.0});
}

TEST_F(PredictorTest, LinearPredictExtrapolatesFromTwoStates) {
  auto older = makeVec({1.0, 5.0});
  auto newer = makeVec({3.0, 9.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = out.get();

  predictor.push_state(1.0, older.get());
  predictor.push_state(2.0, newer.get());
  predictor.predict(BoutSnesPredictor::linear, 2.5, out_vec);

  expectVecValues(out.get(), {4.0, 11.0});
}

TEST_F(PredictorTest, LinearPredictInterpolatesInsideInterval) {
  auto older = makeVec({2.0, 10.0});
  auto newer = makeVec({6.0, 14.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = out.get();

  predictor.push_state(0.0, older.get());
  predictor.push_state(2.0, newer.get());
  predictor.predict(BoutSnesPredictor::linear, 1.0, out_vec);

  expectVecValues(out.get(), {4.0, 12.0});
}

TEST_F(PredictorTest, PushStateKeepsOnlyTwoMostRecentStates) {
  auto oldest = makeVec({0.0});
  auto middle = makeVec({1.0});
  auto newest = makeVec({4.0});
  auto out = makeVec({0.0});
  Vec out_vec = out.get();

  predictor.push_state(0.0, oldest.get());
  predictor.push_state(1.0, middle.get());
  predictor.push_state(2.0, newest.get());
  predictor.predict(BoutSnesPredictor::linear, 3.0, out_vec);

  expectVecValues(out.get(), {7.0});
}

TEST_F(PredictorTest, SetDefaultChangesPredictBehavior) {
  auto older = makeVec({2.0});
  auto newer = makeVec({5.0});
  auto out = makeVec({0.0});
  Vec out_vec = out.get();

  predictor.push_state(0.0, older.get());
  predictor.push_state(1.0, newer.get());

  predictor.setDefault(BoutSnesPredictor::constant);
  predictor.predict(2.0, out_vec);
  expectVecValues(out.get(), {5.0});

  predictor.setDefault(BoutSnesPredictor::linear);
  predictor.predict(2.0, out_vec);
  expectVecValues(out.get(), {8.0});
}

TEST_F(PredictorTest, RescaleSkipsUnallocatedHistorySlots) {
  auto norms = makeVec({2.0, 4.0});

  EXPECT_NO_THROW(predictor.rescale(norms.get()));
}

TEST_F(PredictorTest, RescaleAppliesToAllAllocatedHistoryStates) {
  auto older = makeVec({4.0, 12.0});
  auto newer = makeVec({8.0, 20.0});
  auto norms = makeVec({2.0, 4.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = out.get();

  predictor.push_state(0.0, older.get());
  predictor.push_state(1.0, newer.get());
  predictor.rescale(norms.get());

  predictor.predict(BoutSnesPredictor::constant, 2.0, out_vec);
  expectVecValues(out.get(), {4.0, 5.0});

  predictor.predict(BoutSnesPredictor::linear, 0.0, out_vec);
  expectVecValues(out.get(), {2.0, 3.0});
}

TEST_F(PredictorTest, RescalePreservesSubsequentLinearPrediction) {
  auto older = makeVec({2.0, 8.0});
  auto newer = makeVec({6.0, 12.0});
  auto norms = makeVec({2.0, 4.0});
  auto out = makeVec({0.0, 0.0});
  Vec out_vec = out.get();

  predictor.push_state(1.0, older.get());
  predictor.push_state(3.0, newer.get());
  predictor.rescale(norms.get());
  predictor.predict(BoutSnesPredictor::linear, 4.0, out_vec);

  expectVecValues(out.get(), {4.0, 3.5});
}

} // namespace

#endif
