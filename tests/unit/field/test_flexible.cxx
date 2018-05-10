#include "gtest/gtest.h"

#include "field2d.hxx"
#include "field3d.hxx"
#include "flexible.hxx"

#include "test_extras.hxx"
#include "unused.hxx"

#include <algorithm>
#include <list>
#include <vector>

/// Global mesh
extern Mesh *mesh;

class InterpMesh : public FakeMesh {
public:
  InterpMesh(int nx, int ny, int nz) : FakeMesh(nx, ny, nz) {
    // Mesh only on one process, so global and local indices are the
    // same
    GlobalNx = nx;
    GlobalNy = ny;
    GlobalNz = nz;
    LocalNx = nx;
    LocalNy = ny;
    LocalNz = nz;
    // Small "inner" region
    xstart = 2;
    xend = nx - 3;
    ystart = 2;
    yend = ny - 3;

    transform = std::unique_ptr<ParallelTransform>(new ParallelTransformIdentity());

    // Unused variables
    periodicX = false;
    NXPE = 1;
    PE_XIND = 0;
    StaggerGrids = true;
    IncIntShear = false;
    maxregionblocksize = MAXREGIONBLOCKSIZE;
  }
};

/// Test fixture to make sure the global mesh is our fake one
template <typename T> class FlexibleFieldTest : public ::testing::Test {
public:
  typedef std::list<T> List;
  static T shared_;
  T value_;
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new InterpMesh(nx, ny, nz);
    mesh->createDefaultRegions();
  }

  static void TearDownTestCase() {
    delete mesh;
    mesh = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

class FlexibleFieldTest2 : public ::testing::Test {
public:
  static void SetUpTestCase() {
    // Delete any existing mesh
    if (mesh != nullptr) {
      delete mesh;
      mesh = nullptr;
    }
    mesh = new InterpMesh(nx, ny, nz);
    mesh->createDefaultRegions();
  }

  static void TearDownTestCase() {
    delete mesh;
    mesh = nullptr;
  }

public:
  static const int nx;
  static const int ny;
  static const int nz;
};

template <typename T> const int FlexibleFieldTest<T>::nx = 7;
template <typename T> const int FlexibleFieldTest<T>::ny = 6;
template <typename T> const int FlexibleFieldTest<T>::nz = 7;

const int FlexibleFieldTest2::nx = 7;
const int FlexibleFieldTest2::ny = 6;
const int FlexibleFieldTest2::nz = 2;

typedef ::testing::Types<Field2D, Field3D> FlexibleFieldTypes;
TYPED_TEST_CASE(FlexibleFieldTest, FlexibleFieldTypes);

TYPED_TEST(FlexibleFieldTest, ConstructorField) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  Flexible<TypeParam> flex(field);
  TypeParam copy = flex;
  for (auto i : field) {
    EXPECT_DOUBLE_EQ(field[i], copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, ConstructorPass) {
  Flexible<TypeParam> flex(3., mesh);
  TypeParam copy = flex;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(3., copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, ConstructorMove) {
  Flexible<TypeParam> flex(TypeParam(3., mesh));
  TypeParam copy = flex;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(3., copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, ConstructorReal) {
  Flexible<TypeParam> flex(3.);
  TypeParam copy = flex;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(3., copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, AssignMove) {
  Flexible<TypeParam> flex;
  flex = TypeParam(3., mesh);
  TypeParam copy = flex;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(3., copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, AssignCopy) {
  Flexible<TypeParam> flex;
  TypeParam fu(3., mesh);
  TypeParam copy = flex = fu;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(3., copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, AssignReal) {
  Flexible<TypeParam> flex;
  flex = 3.;
  TypeParam copy = flex;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(3., copy[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, set) {
  TypeParam src;
  Flexible<TypeParam> flex;
  src = 1;
  // Set makes a copy
  flex.set(src);
  // Changing the source should not change the flexible field
  src = 2;
  TypeParam copy = flex;
  for (auto i : copy) {
    EXPECT_DOUBLE_EQ(1., copy[i]);
  }
}

#define ALL_CELL_LOCS                                                                    \
  { CELL_CENTRE, CELL_XLOW, CELL_YLOW, CELL_ZLOW}

TYPED_TEST(FlexibleFieldTest, MultiplactionRhsSame) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  Flexible<TypeParam> flex(field);
  TypeParam one = 1.;
  for (CELL_LOC loc : ALL_CELL_LOCS) {
    one.setLocation(loc);
    TypeParam test = one * flex;
    for (auto i : field.region(RGN_NOBNDRY)) {
      EXPECT_DOUBLE_EQ(field[i], test[i]);
    }
  }
}

TYPED_TEST(FlexibleFieldTest, MultiplactionLhsSame) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  Flexible<TypeParam> flex(field);
  TypeParam one = 1.;
  for (CELL_LOC loc : ALL_CELL_LOCS) {
    one.setLocation(loc);
    TypeParam test = flex * one;
    for (auto i : field.region(RGN_NOBNDRY)) {
      EXPECT_DOUBLE_EQ(field[i], test[i]);
    }
  }
}

TYPED_TEST(FlexibleFieldTest, InterpolationX) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  for (auto i : field) {
    field[i] = i.x;
  }
  Flexible<TypeParam> flex(field);
  TypeParam test = 1.;
  test.setLocation(CELL_XLOW);
  test = test * flex;
  for (auto i : test.region(RGN_NOBNDRY)) {
    EXPECT_DOUBLE_EQ(test[i], i.x - .5);
  }
}

TYPED_TEST(FlexibleFieldTest, InterpolationY) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  for (auto i : field) {
    field[i] = i.y;
  }
  Flexible<TypeParam> flex(field);
  TypeParam test = 1.;
  test.setLocation(CELL_YLOW);
  test = test * flex;
  for (auto i : test.region(RGN_NOBNDRY)) {
    EXPECT_DOUBLE_EQ(test[i], i.y - .5);
  }
}

TYPED_TEST(FlexibleFieldTest, applyTDerivBoundary) {
  Flexible<TypeParam> flex(1);
  // Is not implemented yet - thus it should throw
  // If someone implemnts it, a matching test should be written.
  EXPECT_THROW(flex.applyTDerivBoundary(), BoutException);
}

// Disabled: XLOW->YLOW interpolation is not implemented
// TYPED_TEST(FlexibleFieldTest, InterpolationXY) {
//   TypeParam field(3.);
//   field.setLocation(CELL_XLOW);
//   for (auto i : field) {
//     field[i] = i.x + 3 * i.y;
//   }
//   Flexible<TypeParam> flex(field);
//   TypeParam test = 1.;
//   test.setLocation(CELL_YLOW);
//   test = test * flex;
//   for (auto i : test.region(RGN_NOBNDRY)) {
//     EXPECT_DOUBLE_EQ(test[i], i.x + .5 + 3 * i.y - 1.5);
//   }
// }

TYPED_TEST(FlexibleFieldTest, InterpolationZ) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  int nz = field.getNz();
  int breakval = 3;
  for (auto i:field){
    field[i]=((i.z + breakval)%nz - breakval)*2-1; // Continuous through zero, discontinuity at nz-breakval
  }
  Flexible<TypeParam> flex(field);
  TypeParam test=1.;
  test.setLocation(CELL_ZLOW);
  test = test * flex;
  if (field.is3D()) {
    //for (auto i:test.region(RGN_NOZ)){
    //  EXPECT_DOUBLE_EQ(test[i] , 0);
    //}
    for (int i=0; i<nz-breakval-2; i++) {
      EXPECT_DOUBLE_EQ(test(2,2,i), i*2-2);
    }
    for (int i=nz-breakval+2; i<nz; i++) {
      EXPECT_DOUBLE_EQ(test(2,2,i), (i-nz)*2-2);
    }
  } else {
    for (auto i:test.region(RGN_NOBNDRY)){
      EXPECT_DOUBLE_EQ(test[i] , field[i]);
    }
  }
}

TYPED_TEST(FlexibleFieldTest, MultiplactionRhsField3D) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  Flexible<TypeParam> flex(field);
  for (CELL_LOC loc : ALL_CELL_LOCS) {
    for (CELL_LOC loc : ALL_CELL_LOCS) {
      field = 1.;
      field.setLocation(loc);
      flex.set(field);
    }
    field = 0.;
    field.setLocation(loc);
    flex.set(field);
    Field3D test = 1.;
    test.setLocation(loc);
    test = test * flex;
    for (auto i : test.region(RGN_NOBNDRY)) {
      EXPECT_DOUBLE_EQ(test[i], 0);
    }
  }
}

TYPED_TEST(FlexibleFieldTest, MultiplactionLhsField3D) {
  TypeParam field(3.);
  field.setLocation(CELL_CENTRE);
  Flexible<TypeParam> flex(field);
  for (CELL_LOC loc : ALL_CELL_LOCS) {
    for (CELL_LOC loc : ALL_CELL_LOCS) {
      field = 1.;
      field.setLocation(loc);
      flex.set(field);
    }
    field = 0.;
    field.setLocation(loc);
    flex.set(field);
    Field3D test = 1.;
    test.setLocation(loc);
    test = flex * test;
    for (auto i : test.region(RGN_NOBNDRY)) {
      EXPECT_DOUBLE_EQ(test[i], 0);
    }
  }
}

#define FlexTestFieldOp(name, op, field, in1, in2, res1, res2)                           \
  TYPED_TEST(FlexibleFieldTest, name) {                                                  \
    TypeParam input(3.);                                                                 \
    input.setLocation(CELL_CENTRE);                                                      \
    Flexible<TypeParam> flex(input);                                                     \
    for (CELL_LOC loc : ALL_CELL_LOCS) {                                                 \
      for (CELL_LOC loc : ALL_CELL_LOCS) {                                               \
        input = 1.;                                                                      \
        input.setLocation(loc);                                                          \
        flex.set(input);                                                                 \
      }                                                                                  \
      input = in1;                                                                       \
      input.setLocation(loc);                                                            \
      flex.set(input);                                                                   \
      field tmp((BoutReal)in2);                                                          \
      tmp.setLocation(loc);                                                              \
      auto test = flex op tmp;                                                           \
      for (auto i : test.region(RGN_NOBNDRY)) {                                                              \
        EXPECT_DOUBLE_EQ(test[i], res1);                                                 \
      }                                                                                  \
      test = tmp op flex;                                                                \
      for (auto i : test) {                                                              \
        EXPECT_DOUBLE_EQ(test[i], res2);                                                 \
      }                                                                                  \
    }                                                                                    \
  }

#define FlexTestFieldOpReal(name, op, in1, in2, res1, res2)                              \
  TYPED_TEST(FlexibleFieldTest, name) {                                                  \
    TypeParam input(3.);                                                                 \
    input.setLocation(CELL_CENTRE);                                                      \
    Flexible<TypeParam> flex(input);                                                     \
    input = in1;                                                                         \
    flex.set(input);                                                                     \
    BoutReal tmp(in2);                                                                   \
    auto test = flex op tmp;                                                             \
    for (auto i : test.region(RGN_NOBNDRY)) {                                                                \
      EXPECT_DOUBLE_EQ(test[i], res1);                                                   \
    }                                                                                    \
    test = tmp op flex;                                                                  \
    for (auto i : test.region(RGN_NOBNDRY)) {                                                                \
      EXPECT_DOUBLE_EQ(test[i], res2);                                                   \
    }                                                                                    \
  }

#define FlexTestFieldUpdate(name, op, field, in1, in2, res1, res2)                       \
  TYPED_TEST(FlexibleFieldTest, name) {                                                  \
    for (CELL_LOC loc1 : ALL_CELL_LOCS) {                                                \
      TypeParam input((BoutReal)in1);                                                    \
      input.setLocation(loc1);                                                           \
      Flexible<TypeParam> flex(input);                                                   \
      for (CELL_LOC loc2 : ALL_CELL_LOCS) {                                              \
        input = in2;                                                                     \
        input.setLocation(loc2);                                                         \
        if (loc1 == loc2) {                                                              \
          flex op input;                                                                 \
          TypeParam tmp = flex;                                                          \
          EXPECT_EQ(tmp.getLocation(), loc1);                                            \
          for (auto i : tmp.region(RGN_NOBNDRY)) {                                                           \
            EXPECT_DOUBLE_EQ(tmp[i], res1);                                              \
          }                                                                              \
        } else {                                                                         \
          EXPECT_THROW(flex op input, BoutException);                                    \
        }                                                                                \
      }                                                                                  \
    }                                                                                    \
  }

#define FlexTestFieldUpdateReal(name, op, in1, in2, res1, res2)                          \
  TYPED_TEST(FlexibleFieldTest, name) {                                                  \
    for (CELL_LOC loc1 : ALL_CELL_LOCS) {                                                \
      TypeParam input((BoutReal)in1);                                                    \
      input.setLocation(loc1);                                                           \
      Flexible<TypeParam> flex(input);                                                   \
      BoutReal tmp(in2);                                                                 \
      flex op tmp;                                                                       \
      for (auto i : input.region(RGN_NOBNDRY)) {                                                             \
        EXPECT_DOUBLE_EQ(flex[i], res1);                                                 \
      }                                                                                  \
    }                                                                                    \
  }

FlexTestFieldOp(MultiplicationField2D, *, Field2D, 1, 0, 0, 0);
FlexTestFieldOp(MultiplicationField3D, *, Field3D, 1, 0, 0, 0);
FlexTestFieldOpReal(MultiplicationReal, *, 1, 0, 0, 0);
FlexTestFieldOp(DivisionField2D, /, Field2D, 2, 1, 2, .5);
FlexTestFieldOp(DivisionField3D, /, Field3D, 2, 1, 2, .5);
FlexTestFieldOpReal(DivisionReal, /, 2, 1, 2, .5);
FlexTestFieldOp(AdditionField2D, +, Field2D, 2, 1, 3, 3);
FlexTestFieldOp(AdditionField3D, +, Field3D, 2, 1, 3, 3);
FlexTestFieldOpReal(AdditionReal, +, 2, 1, 3, 3);
FlexTestFieldOp(SubstructionField2D, -, Field2D, 2, 1, 1, -1);
FlexTestFieldOp(SubstructionField3D, -, Field3D, 2, 1, 1, -1);
FlexTestFieldOpReal(SubstructionReal, -, 2, 1, 1, -1);

FlexTestFieldUpdate(UpdateMultiplicationField2D, *=, Field2D, 1, 2, 2, 2);
FlexTestFieldUpdate(UpdateMultiplicationField3D, *=, Field3D, 1, 2, 2, 2);
FlexTestFieldUpdateReal(UpdateMultiplicationReal, *=, 1, 2, 2, 2);
FlexTestFieldUpdate(UpdateDivisionField2D, /=, Field2D, 2, 1, 2, .5);
FlexTestFieldUpdate(UpdateDivisionField3D, /=, Field3D, 2, 1, 2, .5);
FlexTestFieldUpdateReal(UpdateDivisionReal, /=, 2, 1, 2, .5);
FlexTestFieldUpdate(UpdateAdditionField2D, +=, Field2D, 2, 1, 3, 3);
FlexTestFieldUpdate(UpdateAdditionField3D, +=, Field3D, 2, 1, 3, 3);
FlexTestFieldUpdateReal(UpdateAdditionReal, +=, 2, 1, 3, 3);
FlexTestFieldUpdate(UpdateSubstructionField2D, -=, Field2D, 2, 1, 1, -1);
FlexTestFieldUpdate(UpdateSubstructionField3D, -=, Field3D, 2, 1, 1, -1);
FlexTestFieldUpdateReal(UpdateSubstructionReal, -=, 2, 1, 1, -1);

TYPED_TEST(FlexibleFieldTest, AccessDataIterator) {
  Flexible<TypeParam> flex(3.);
  TypeParam test = flex;
  for (auto i : test) {
    EXPECT_DOUBLE_EQ(test[i], flex[i]);
  }
}

TYPED_TEST(FlexibleFieldTest, AccessIndices) {
  Flexible<TypeParam> flex(3.);
  TypeParam test = flex;
  for (auto i : test) {
    Indices ii{i.x, i.y, i.z};
    EXPECT_DOUBLE_EQ(test[ii], flex[ii]);
  }
}

TYPED_TEST(FlexibleFieldTest, AccessIndicesConst) {
  const TypeParam test(3.);
  const Flexible<TypeParam> flex(test);
  for (auto i : test) {
    Indices ii{i.x, i.y, i.z};
    EXPECT_DOUBLE_EQ(test[ii], flex[ii]);
  }
}

// TYPED_TEST(FlexibleFieldTest, AccessDataIterator2) {
//   Flexible<TypeParam> flex(3.);
//   TypeParam test=flex;
//   for (auto i:flex.region()){
//     EXPECT_DOUBLE_EQ(test[i] , flex[i]);
//   }
// }

TYPED_TEST(FlexibleFieldTest, isReal) {
  Flexible<TypeParam> flex(3.);
  TypeParam test = flex;
  EXPECT_EQ(test.isReal(), flex.isReal());
}

TYPED_TEST(FlexibleFieldTest, is3D) {
  Flexible<TypeParam> flex(3.);
  TypeParam test = flex;
  EXPECT_EQ(test.is3D(), flex.is3D());
}

TYPED_TEST(FlexibleFieldTest, byteSize) {
  Flexible<TypeParam> flex(3.);
  TypeParam test = flex;
  EXPECT_EQ(test.byteSize(), flex.byteSize());
}

TYPED_TEST(FlexibleFieldTest, BoutRealSize) {
  Flexible<TypeParam> flex(3.);
  TypeParam test = flex;
  EXPECT_EQ(test.BoutRealSize(), flex.BoutRealSize());
}

class FakeField : public Field {
public:
  FakeField(int i=1) { location = (CELL_LOC)i; }
  virtual CELL_LOC getLocation() const override { return location; }
  bool is3D() const { return dim; }
  bool isReal() const { return real; }
  bool isAllocated() const { return alloc; }
  int byteSize() const { return bs; }
  int BoutRealSize() const { return brs; }
  virtual const IndexRange region(REGION UNUSED(rgn)) const {
    return IndexRange{0, 0, 0, 0, 0, 0};
  }
  virtual const BoutReal &operator[](const Indices &) const override { return val; }
  virtual BoutReal &operator[](const Indices &) { return val; }
  virtual BoutReal& operator[](const DataIterator &) { return val; }
  void allocate() { alloc = true; };
  void accept(FieldVisitor &UNUSED(v)){};
  void doneComms() { coms = true; };
  void applyBoundary(bool) { bndry = true; };
  bool accepted, real, dim;
  int brs, bs;
  bool coms, bndry, alloc;
  CELL_LOC location;
  BoutReal val;
  const DataIterator iterator() const { return DataIterator(0, 0, 0, 0, 0, 0); }
  const DataIterator begin() const { return iterator(); }
  const DataIterator end() const { return iterator(); }

private:
};

FakeField interp_to(const FakeField &f, CELL_LOC loc, REGION UNUSED(region)) {
  FakeField tmp(f);
  tmp.location = loc;
  return tmp;
}

TEST_F(FlexibleFieldTest2, FakeField) {
  FakeField(2);
  EXPECT_NO_THROW(Flexible<FakeField>(2));
  EXPECT_THROW(Flexible<FakeField>(10), BoutException);
  EXPECT_THROW(Flexible<FakeField>(0), BoutException);
  EXPECT_THROW(Flexible<FakeField>(-10), BoutException);
}

#define PASS_ON_TEST(isReal, real, true, false)                                          \
  TEST_F(FlexibleFieldTest2, isReal) {                                                   \
    FakeField test(2);                                                                   \
    Flexible<FakeField> flex(test);                                                      \
    test.real = true;                                                                    \
    flex = test;                                                                         \
    EXPECT_EQ(test.isReal(), flex.isReal());                                             \
    test.real = false;                                                                   \
    flex = test;                                                                         \
    EXPECT_EQ(test.isReal(), flex.isReal());                                             \
  }

PASS_ON_TEST(isReal, real, true, false);
PASS_ON_TEST(is3D, dim, true, false);
PASS_ON_TEST(byteSize, bs, 128, 77);
PASS_ON_TEST(BoutRealSize, brs, 33, 74);

#undef PASS_ON_TEST

TEST_F(FlexibleFieldTest2, AccessIntField2D) {
  Field2D test(3.);
  for (auto i : test) {
    test[i] = i.x + i.y;
  }
  Flexible<Field2D> flex(test);
  for (auto i : test) {
    EXPECT_DOUBLE_EQ(test(i.x, i.y), flex(i.x, i.y));
  }
}

#define PASS_ON_TEST(doneComms, coms)                                                    \
  TEST_F(FlexibleFieldTest2, doneComms) {                                                \
    FakeField f(2);                                                                      \
    f.coms = false;                                                                      \
    Flexible<FakeField> flex(f);                                                         \
    EXPECT_FALSE(flex.get(CELL_CENTRE).coms);                                            \
    flex.doneComms();                                                                    \
    EXPECT_TRUE(flex.get(CELL_CENTRE).coms);                                             \
  }

PASS_ON_TEST(doneComms, coms);
PASS_ON_TEST(applyBoundary, bndry);
// PASS_ON_TEST(allocate,alloc);

#undef PASS_ON_TEST
