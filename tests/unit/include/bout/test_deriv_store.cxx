#include "gtest/gtest.h"

#include <bout/deriv_store.hxx>
#include <bout_types.hxx>
#include <boutexception.hxx>
#include <output.hxx>
#include <unused.hxx>

#include <typeindex>

using FieldType = std::vector<BoutReal>;

using standardType = DerivativeStore<FieldType>::standardFunc;
using flowType = DerivativeStore<FieldType>::upwindFunc;

void standardReturnTenSetToOne(const FieldType& UNUSED(inp), FieldType& out,
                               const REGION = RGN_ALL) {
  out.resize(10, 1.0);
}

void flowReturnSixSetToTwo(const FieldType& UNUSED(vel), const FieldType& UNUSED(inp),
                           FieldType& out, const REGION = RGN_ALL) {
  out.resize(6, 2.0);
}

TEST(DerivativeStoreTest, CanGetInstance) {
  EXPECT_NO_THROW(DerivativeStore<FieldType>::getInstance());
}

TEST(DerivativeStoreTest, IsAllEmpty) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  EXPECT_TRUE(store.isEmpty());
}

TEST(DerivativeStoreTest, IsEmptySpecificType) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  EXPECT_TRUE(store.isEmpty(DERIV::Standard, DIRECTION::X, STAGGER::None));
}

TEST(DerivativeStoreTest, GetAvailableMethodsEmpty) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  EXPECT_TRUE(store.isEmpty());
  const auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 0);
}

TEST(DerivativeStoreTest, RegisterStandardMethod) {
  auto& store = DerivativeStore<FieldType>::getInstance();

  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "FirstStandard");
  const auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);

  store.reset();
}

TEST(DerivativeStoreTest, RegisterStandardMethodAndClear) {
  auto& store = DerivativeStore<FieldType>::getInstance();

  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "FirstStandard");

  auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());

  store.reset();
  methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 0);
}

TEST(DerivativeStoreTest, RegisterMatchingMethodStandardMethodTwice) {
  auto& store = DerivativeStore<FieldType>::getInstance();

  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "FirstStandard");
  auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());

  // Try to register another method with the same key
  EXPECT_THROW(store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X,
                                        STAGGER::None, "FirstStandard"),
               BoutException);

  methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());

  store.reset();
}

TEST(DerivativeStoreTest, RegisterMatchingMethodStandardMethodTwo) {
  auto& store = DerivativeStore<FieldType>::getInstance();

  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "FirstStandard");
  auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);

  // Register another method with a different key
  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "SecondStandard");

  methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 2);

  EXPECT_NE(methods.find("FirstStandard"), methods.end());
  EXPECT_NE(methods.find("SecondStandard"), methods.end());
  store.reset();
}

TEST(DerivativeStoreTest, RegisterStandardSecond) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  const DERIV type = DERIV::StandardSecond;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(standardType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
  store.reset();
}

TEST(DerivativeStoreTest, RegisterMatchingMethodStandardSecondMethodTwice) {
  auto& store = DerivativeStore<FieldType>::getInstance();

  store.registerDerivative(standardType{}, DERIV::StandardSecond, DIRECTION::X,
                           STAGGER::None, "SecondStandard");

  // Try to register another method with the same key
  EXPECT_THROW(store.registerDerivative(standardType{}, DERIV::StandardSecond,
                                        DIRECTION::X, STAGGER::None, "SecondStandard"),
               BoutException);

  auto methods = store.getAvailableMethods(DERIV::StandardSecond, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("SecondStandard"), methods.end());

  store.reset();
}

TEST(DerivativeStoreTest, RegisterStandardFourth) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  const DERIV type = DERIV::StandardFourth;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(standardType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
  store.reset();
}

TEST(DerivativeStoreTest, RegisterMatchingMethodStandardFourthMethodTwice) {
  auto& store = DerivativeStore<FieldType>::getInstance();

  store.registerDerivative(standardType{}, DERIV::StandardFourth, DIRECTION::X,
                           STAGGER::None, "FourthStandard");

  // Try to register another method with the same key
  EXPECT_THROW(store.registerDerivative(standardType{}, DERIV::StandardFourth,
                                        DIRECTION::X, STAGGER::None, "FourthStandard"),
               BoutException);

  auto methods = store.getAvailableMethods(DERIV::StandardFourth, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FourthStandard"), methods.end());

  store.reset();
}

TEST(DerivativeStoreTest, RegisterUpwind) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  const DERIV type = DERIV::Upwind;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
  store.reset();
}

TEST(DerivativeStoreTest, RegisterFlux) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  const DERIV type = DERIV::Flux;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
  store.reset();
}

TEST(DerivativeStoreTest, RegisterStandardAndGetBack) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  const DERIV type = DERIV::Standard;
  const DIRECTION dir = DIRECTION::X;
  const std::string firstName = "FIRSTSTANDARD";

  store.forceDefaultMethod("UNKNOWN", type, dir, STAGGER::None);

  store.registerDerivative(standardReturnTenSetToOne, type, dir, STAGGER::None,
                           firstName);
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find(firstName), methods.end());

  // Note pass functions around by value etc. so we can't just compare the address
  // of the returned function against the original method.
  standardType returned = store.getStandardDerivative(firstName, dir);

  FieldType inOrig, outOrig;
  standardReturnTenSetToOne(inOrig, outOrig);

  FieldType outRet;
  returned(inOrig, outRet, RGN_ALL);

  EXPECT_EQ(outOrig.size(), 10);
  ASSERT_EQ(outRet.size(), outOrig.size());

  for (unsigned int i = 0; i < outOrig.size(); i++) {
    EXPECT_EQ(outOrig[i], outRet[i]);
  };

  store.reset();
}

TEST(DerivativeStoreTest, RegisterFlowAndGetBack) {
  auto& store = DerivativeStore<FieldType>::getInstance();
  const DERIV type = DERIV::Upwind;
  const DIRECTION dir = DIRECTION::X;
  const std::string firstName = "FIRSTSTANDARD";

  store.forceDefaultMethod("UNKNOWN", type, dir, STAGGER::None);

  store.registerDerivative(flowReturnSixSetToTwo, type, dir, STAGGER::None, firstName);
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find(firstName), methods.end());

  // Note pass functions around by value etc. so we can't just compare the address
  // of the returned function against the original method.
  flowType returned = store.getUpwindDerivative(firstName, dir);

  FieldType inOrig, outOrig;
  flowReturnSixSetToTwo(inOrig, inOrig, outOrig);

  FieldType outRet;
  returned(inOrig, inOrig, outRet, RGN_ALL);

  EXPECT_EQ(outOrig.size(), 6);
  ASSERT_EQ(outRet.size(), outOrig.size());

  for (unsigned int i = 0; i < outOrig.size(); i++) {
    EXPECT_EQ(outOrig[i], outRet[i]);
  };

  store.reset();
}
