#include "gtest/gtest.h"
#include "boutexception.hxx"

#include <iostream>
#include <string>

TEST(BoutExceptionTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutException("test"), BoutException);
}

TEST(BoutExceptionTest, WhatTest) {
  try {
    throw BoutException(std::string("Test message"));
  } catch (BoutException &e) {
    std::string message(e.what());
    EXPECT_NE(message.find("Test message"), std::string::npos);
  }
  try {
    throw BoutException("%s", "second");
  } catch (BoutException &e) {
    std::string message(e.what());
    EXPECT_NE(message.find("second"), std::string::npos);
  }
}

TEST(BoutRhsFailTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutRhsFail("RHS Fail test"), BoutRhsFail);
}

TEST(BoutParallelRhsFailTest, ThrowCorrect) {
  EXPECT_THROW(BoutParallelThrowRhsFail(1, "RHS Fail test"), BoutRhsFail);
}

TEST(BoutParallelRhsFailTest, ThrowNoError) {
  EXPECT_NO_THROW(BoutParallelThrowRhsFail(0, "RHS Fail test"));
}

TEST(BoutIterationFailTest, ThrowCorrect) {
  EXPECT_THROW(throw BoutIterationFail("Iteration fail test"), BoutIterationFail);
}
