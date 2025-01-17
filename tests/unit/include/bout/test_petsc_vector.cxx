#include "bout/build_defines.hxx"

#include <utility>

#include "gtest/gtest.h"

#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/fieldperp.hxx"
#include "bout/operatorstencil.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/region.hxx"

#if BOUT_HAS_PETSC

#include <petscconf.h>

#include "fake_mesh_fixture.hxx"

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
      : FakeMeshFixture(), field(1.5, bout::globals::mesh),
        stencil(squareStencil<ind_type>(bout::globals::mesh)),
        indexer(std::make_shared<GlobalIndexer<F>>(bout::globals::mesh, stencil)) {
    PetscErrorPrintf = PetscErrorPrintfNone;
  }

  virtual ~PetscVectorTest() { PetscErrorPrintf = PetscErrorPrintfDefault; }
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
  Vec* vectorPtr = vector.get();
  Vec* copyPtr = copy.get();
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

TYPED_TEST(PetscVectorTest, SetElement) {
  SCOPED_TRACE("FieldAssignment");
  PetscVector<TypeParam> vector(this->field, this->indexer);
  const TypeParam val(-10.);

  BOUT_FOR(index, val.getRegion("RGN_ALL")) { vector(index) = val[index]; }
  vector.assemble();

  Vec* vectorPtr = vector.get();
  PetscScalar* vecContents = nullptr;
  PetscInt size = 0;
  VecGetArray(*vectorPtr, &vecContents);
  VecGetLocalSize(*vectorPtr, &size);
  ASSERT_EQ(size, this->field.size());
  TypeParam result = vector.toField();
  BOUT_FOR(i, this->field.getRegion("RGN_NOY")) { EXPECT_EQ(result[i], val[i]); }
}

// Test getting elements
TYPED_TEST(PetscVectorTest, TestGetElements) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  BOUT_FOR(i, this->field.getRegion("RGN_NOBNDRY")) {
    vector(i) = (2.5 * this->field[i] - 1.0);
  }
  vector.assemble();
  TypeParam result = vector.toField();
  Vec* rawvec = vector.get();
  PetscScalar* vecContents = nullptr;
  VecAssemblyBegin(*rawvec);
  VecAssemblyEnd(*rawvec);
  VecGetArray(*rawvec, &vecContents);
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

TYPED_TEST(PetscVectorTest, TestMixedSetting) {
  PetscVector<TypeParam> vector(this->field, this->indexer);
  typename TypeParam::ind_type index1 = *(this->field.getRegion("RGN_NOBNDRY").begin());
  typename TypeParam::ind_type index2(index1.ind + 1);
  const PetscScalar r = 3.141592;
  vector(index1) = r;
  vector(index2) += r;
  vector.assemble();
  PetscScalar* vecContents = nullptr;
  VecGetArray(*(vector.get()), &vecContents);
  ASSERT_EQ(vecContents[index1.ind], r);
  ASSERT_EQ(vecContents[index2.ind], this->field[index2] + r);
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
