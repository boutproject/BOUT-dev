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
                               const std::string& = "RGN_ALL") {
  out.resize(10, 1.0);
}

void flowReturnSixSetToTwo(const FieldType& UNUSED(vel), const FieldType& UNUSED(inp),
                           FieldType& out, const std::string& = "RGN_ALL") {
  out.resize(6, 2.0);
}

class DerivativeStoreTest : public ::testing::Test {
public:
  DerivativeStoreTest() : store(DerivativeStore<FieldType>::getInstance()) {}
  ~DerivativeStoreTest() override { store.reset(); }

  DerivativeStore<FieldType>& store;
};

TEST_F(DerivativeStoreTest, CanGetInstance) {
  EXPECT_NO_THROW(DerivativeStore<FieldType>::getInstance());
}

TEST_F(DerivativeStoreTest, IsAllEmpty) { EXPECT_TRUE(store.isEmpty()); }

TEST_F(DerivativeStoreTest, IsEmptySpecificType) {
  EXPECT_TRUE(store.isEmpty(DERIV::Standard, DIRECTION::X, STAGGER::None));
}

TEST_F(DerivativeStoreTest, GetAvailableMethodsEmpty) {
  EXPECT_TRUE(store.isEmpty());
  const auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 0);
}

TEST_F(DerivativeStoreTest, RegisterStandardMethod) {
  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "FirstStandard");
  const auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
}

TEST_F(DerivativeStoreTest, RegisterStandardMethodAndClear) {

  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "FirstStandard");

  auto methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());

  store.reset();

  methods = store.getAvailableMethods(DERIV::Standard, DIRECTION::X);
  EXPECT_EQ(methods.size(), 0);
}

TEST_F(DerivativeStoreTest, RegisterMatchingMethodStandardMethodTwice) {

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
}

TEST_F(DerivativeStoreTest, RegisterMatchingMethodStandardMethodTwo) {

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
}

TEST_F(DerivativeStoreTest, RegisterStandardSecond) {
  const DERIV type = DERIV::StandardSecond;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(standardType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterMatchingMethodStandardSecondMethodTwice) {

  store.registerDerivative(standardType{}, DERIV::StandardSecond, DIRECTION::X,
                           STAGGER::None, "SecondStandard");

  // Try to register another method with the same key
  EXPECT_THROW(store.registerDerivative(standardType{}, DERIV::StandardSecond,
                                        DIRECTION::X, STAGGER::None, "SecondStandard"),
               BoutException);

  auto methods = store.getAvailableMethods(DERIV::StandardSecond, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("SecondStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterStandardFourth) {
  const DERIV type = DERIV::StandardFourth;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(standardType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterMatchingMethodStandardFourthMethodTwice) {

  store.registerDerivative(standardType{}, DERIV::StandardFourth, DIRECTION::X,
                           STAGGER::None, "FourthStandard");

  // Try to register another method with the same key
  EXPECT_THROW(store.registerDerivative(standardType{}, DERIV::StandardFourth,
                                        DIRECTION::X, STAGGER::None, "FourthStandard"),
               BoutException);

  auto methods = store.getAvailableMethods(DERIV::StandardFourth, DIRECTION::X);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FourthStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterUpwindAsStandard) {
  const DERIV type = DERIV::Standard;
  const DIRECTION dir = DIRECTION::X;

  EXPECT_THROW(
      store.registerDerivative(flowType{}, type, dir, STAGGER::None, "BadSignature"),
      BoutException);
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 0);
  EXPECT_EQ(methods.find("BadSignature"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterStandardAsUpwind) {
  const DERIV type = DERIV::Upwind;
  const DIRECTION dir = DIRECTION::X;

  EXPECT_THROW(
      store.registerDerivative(standardType{}, type, dir, STAGGER::None, "BadSignature"),
      BoutException);
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 0);
  EXPECT_EQ(methods.find("BadSignature"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterUpwind) {
  const DERIV type = DERIV::Upwind;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterUpwindTwice) {
  const DERIV type = DERIV::Upwind;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard");

  EXPECT_THROW(
      store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard"),
      BoutException);

  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterFlux) {
  const DERIV type = DERIV::Flux;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard");
  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterFluxTwice) {
  const DERIV type = DERIV::Flux;
  const DIRECTION dir = DIRECTION::X;

  store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard");

  EXPECT_THROW(
      store.registerDerivative(flowType{}, type, dir, STAGGER::None, "FirstStandard"),
      BoutException);

  auto methods = store.getAvailableMethods(type, dir);
  EXPECT_EQ(methods.size(), 1);
  EXPECT_NE(methods.find("FirstStandard"), methods.end());
}

TEST_F(DerivativeStoreTest, RegisterStandardAndGetBack) {
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
  standardReturnTenSetToOne({}, outOrig);

  FieldType outRet;
  returned(inOrig, outRet, "RGN_ALL");

  EXPECT_EQ(outOrig.size(), 10);
  ASSERT_EQ(outRet.size(), outOrig.size());

  EXPECT_EQ(outOrig, outRet);
}

TEST_F(DerivativeStoreTest, RegisterFlowAndGetBack) {
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
  returned(inOrig, inOrig, outRet, "RGN_ALL");

  EXPECT_EQ(outOrig.size(), 6);
  ASSERT_EQ(outRet.size(), outOrig.size());

  EXPECT_EQ(outOrig, outRet);
}

TEST_F(DerivativeStoreTest, GetUnknownDerivativeFromStandard) {
  // Register something we're not just throwing because the store is empty
  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "something");
  EXPECT_THROW(store.getStandardDerivative("unknown", DIRECTION::X), BoutException);
}

TEST_F(DerivativeStoreTest, GetUnknownDerivativeFromFlow) {
  // Register something we're not just throwing because the store is empty
  store.registerDerivative(flowType{}, DERIV::Flux, DIRECTION::X, STAGGER::None,
                           "something");
  EXPECT_THROW(store.getFlowDerivative("unknown", DIRECTION::X), BoutException);
}

TEST_F(DerivativeStoreTest, GetUpwindFromStandardDerivative) {
  // Register something we're not just throwing because the store is empty
  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "bad type");
  store.registerDerivative(flowType{}, DERIV::Flux, DIRECTION::X, STAGGER::None,
                           "bad type");
  EXPECT_THROW(
      store.getStandardDerivative("bad type", DIRECTION::X, STAGGER::None, DERIV::Upwind),
      BoutException);
}

TEST_F(DerivativeStoreTest, GetFluxFromStandardDerivative) {
  // Register something we're not just throwing because the store is empty
  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "bad type");
  store.registerDerivative(flowType{}, DERIV::Flux, DIRECTION::X, STAGGER::None,
                           "bad type");
  EXPECT_THROW(
      store.getStandardDerivative("bad type", DIRECTION::X, STAGGER::None, DERIV::Flux),
      BoutException);
}

TEST_F(DerivativeStoreTest, GetStandardFromFlowDerivative) {
  // Register something we're not just throwing because the store is empty
  store.registerDerivative(standardType{}, DERIV::Standard, DIRECTION::X, STAGGER::None,
                           "bad type");
  store.registerDerivative(flowType{}, DERIV::Flux, DIRECTION::X, STAGGER::None,
                           "bad type");
  EXPECT_THROW(
      store.getFlowDerivative("bad type", DIRECTION::X, STAGGER::None, DERIV::Standard),
      BoutException);
}
