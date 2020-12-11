#include "bout/build_config.hxx"

#include <utility>

#include "test_extras.hxx"
#include "gtest/gtest.h"

#include "field2d.hxx"
#include "field3d.hxx"
#include "fieldperp.hxx"
#include "bout/operatorstencil.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#if BOUT_HAS_PETSC

#include <petscconf.h>

namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals
} // namespace bout

// The unit tests use the global mesh
using namespace bout::globals;

// Reuse the "standard" fixture for FakeMesh
template <typename F>
class PetscVectorTest : public FakeMeshFixture {
public:
  using ind_type = typename F::ind_type;
  WithQuietOutput all{output};
  F field;
  OperatorStencil<ind_type> stencil;
  IndexerPtr<F> indexer;

  PetscVectorTest()
      : FakeMeshFixture(), field(bout::globals::mesh),
        stencil(squareStencil<ind_type>(bout::globals::mesh)),
        indexer(std::make_shared<GlobalIndexer<F>>(bout::globals::mesh, stencil)) {
    field.allocate();
    field = 1.5;
    PetscErrorPrintf = PetscErrorPrintfNone;
  }

  virtual ~PetscVectorTest() {
    PetscErrorPrintf = PetscErrorPrintfDefault;
  }
};

using FieldTypes = ::testing::Types<Field3D, Field2D, FieldPerp>;
TYPED_TEST_SUITE(PetscVectorTest, FieldTypes);

void testArraysEqual(PetscScalar* s1, PetscScalar* s2, PetscInt n) {
  for (int i = 0; i < n; i++) {
    EXPECT_DOUBLE_EQ(s1[i], s2[i]);
  }
}

void testVectorsEqual(Vec* v1, Vec* v2) {
  PetscScalar *v1Contents, *v2Contents;
  PetscInt n1, n2;
  VecGetArray(*v1, &v1Contents);
  VecGetArray(*v2, &v2Contents);
  VecGetLocalSize(*v1, &n1);
  VecGetLocalSize(*v2, &n2);
  ASSERT_EQ(n1, n2);
  testArraysEqual(v1Contents, v2Contents, n1);
}

// Test constructor from field
TYPED_TEST(PetscVectorTest, FieldConstructor) {
  BOUT_FOR(i, this->field.getRegion("RGN_ALL")) {
    this->field[i] = static_cast<BoutReal>(i.ind);
  }
  PetscVector<TypeParam> vector(this->field, this->indexer);
  Vec* vectorPtr = vector.get();
  PetscScalar* vecContents;
  PetscInt n;
  VecGetArray(*vectorPtr, &vecContents);
  VecGetLocalSize(*vectorPtr, &n);
  ASSERT_EQ(n, this->field.getNx() * this->field.getNy() * this->field.getNz());
  TypeParam result = vector.toField();
  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) { EXPECT_EQ(result[i], this->field[i]); }
}

// Test copy constructor
TYPED_TEST(PetscVectorTest, CopyConstructor) {
  SCOPED_TRACE("CopyConstructor");
  PetscVector<TypeParam> vector(this->field, this->indexer);
  PetscVector<TypeParam> copy(vector);
  Vec *vectorPtr = vector.get(), *copyPtr = copy.get();
  EXPECT_NE(vectorPtr, copyPtr);
  testVectorsEqual(vectorPtr, copyPtr);
}

// Test move constructor
TYPED_TEST(PetscVectorTest, MoveConstructor) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  Vec vectorPtr = *vector.get();
  EXPECT_NE(vectorPtr, nullptr);
  PetscVector<TypeParam> moved(std::move(vector));
  Vec movedPtr = *moved.get();
  EXPECT_EQ(vectorPtr, movedPtr);
}

// Test assignment from field
TYPED_TEST(PetscVectorTest, FieldAssignment) {
  SCOPED_TRACE("FieldAssignment");
  PetscVector<TypeParam> vector(this->field, this->indexer);
  const TypeParam val(-10.);
  vector = val;
  Vec* vectorPtr = vector.get();
  PetscScalar* vecContents;
  PetscInt n;
  VecGetArray(*vectorPtr, &vecContents);
  VecGetLocalSize(*vectorPtr, &n);
  ASSERT_EQ(n, this->field.getNx() * this->field.getNy() * this->field.getNz());
  TypeParam result = vector.toField();
  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) { EXPECT_EQ(result[i], val[i]); }
}

// Test copy assignment
TYPED_TEST(PetscVectorTest, CopyAssignment) {
  SCOPED_TRACE("CopyAssignment");
  PetscVector<TypeParam> vector(this->field, this->indexer);
  PetscVector<TypeParam> copy = vector;
  Vec *vectorPtr = vector.get(), *copyPtr = copy.get();
  EXPECT_NE(vectorPtr, copyPtr);
  testVectorsEqual(vectorPtr, copyPtr);
}

// Test move assignment
TYPED_TEST(PetscVectorTest, MoveAssignment) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  Vec vectorPtr = *vector.get();
  EXPECT_NE(vectorPtr, nullptr);
  PetscVector<TypeParam> moved = std::move(vector);
  Vec movedPtr = *moved.get();
  EXPECT_EQ(vectorPtr, movedPtr);
}

// Test getting elements
TYPED_TEST(PetscVectorTest, TestGetElements) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    vector(i) = (2.5 * this->field[i] - 1.0);
  }
  Vec* rawvec = vector.get();
  PetscScalar* vecContents;
  VecAssemblyBegin(*rawvec);
  VecAssemblyEnd(*rawvec);
  VecGetArray(*rawvec, &vecContents);
  TypeParam result = vector.toField();
  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    EXPECT_EQ(result[i], 2.5 * this->field[i] - 1.0);
  }
}

// Test getting constant elements
TYPED_TEST(PetscVectorTest, TestGetElementsConst) {
  const PetscVector<TypeParam> vector(this->field, this->indexer);
  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    const BoutReal element = vector(i);
    EXPECT_EQ(element, this->field[i]);
  }
}

// Test assemble
TYPED_TEST(PetscVectorTest, TestAssemble) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  Vec* rawvec = vector.get();
  const PetscInt i = 4;
  const PetscScalar r = 3.141592;
  VecSetValues(*rawvec, 1, &i, &r, INSERT_VALUES);
  vector.assemble();
  PetscScalar* vecContents;
  VecGetArray(*rawvec, &vecContents);
  ASSERT_EQ(vecContents[i], r);
}

#ifdef PETSC_USE_DEBUG

// Test trying to get an element from an uninitialised vector
TYPED_TEST(PetscVectorTest, TestGetUninitialised) {
  PetscVector<TypeParam> vector;
  typename TypeParam::ind_type index(0);
  EXPECT_THROW(vector(index), BoutException);
}

#if CHECKLEVEL >= 3
// Test trying to get an element that is out of bounds
TYPED_TEST(PetscVectorTest, TestGetOutOfBounds) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  typename TypeParam::ind_type index1(this->field.getNx() * this->field.getNy()
                                      * this->field.getNz());
  EXPECT_THROW(vector(index1), BoutException);
  typename TypeParam::ind_type index2(-1);
  EXPECT_THROW(vector(index2), BoutException);
  typename TypeParam::ind_type index3(10000000);
  EXPECT_THROW(vector(index3), BoutException);
}
#endif // CHECKLEVEL >= 3

// Test trying to use both INSERT_VALUES and ADD_VALUES
TYPED_TEST(PetscVectorTest, TestMixedSetting) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  typename TypeParam::ind_type i = *(this->field.getRegion("RGN_NOBNDRY").begin());
  typename TypeParam::ind_type j(i.ind + 1);
  const PetscScalar r = 3.141592;
  vector(i) = r;
  EXPECT_THROW(vector(j) += r, BoutException);
}

// Test destroy
TYPED_TEST(PetscVectorTest, TestDestroy) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  Vec oldVec = *vector.get();
  Vec newVec;
  PetscErrorCode err;
  vector.destroy();
  err = VecDuplicate(oldVec, &newVec);
  ASSERT_NE(err, 0); // If original vector was destroyed, should not
                     // be able to duplicate it.
}

#endif // PETSC_USE_DEBUG

// Test swap
TYPED_TEST(PetscVectorTest, TestSwap) {
  PetscVector<TypeParam> lhs(this->field, this->indexer), rhs(this->field, this->indexer);
  Vec l0 = *lhs.get(), r0 = *rhs.get();
  EXPECT_NE(l0, nullptr);
  EXPECT_NE(r0, nullptr);
  swap(lhs, rhs);
  Vec l1 = *lhs.get(), r1 = *rhs.get();
  EXPECT_NE(l0, l1);
  EXPECT_NE(r0, r1);
  EXPECT_EQ(l0, r1);
  EXPECT_EQ(r0, l1);
}

#endif // BOUT_HAS_PETSC
