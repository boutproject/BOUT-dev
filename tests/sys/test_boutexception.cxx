#include "gtest/gtest.h"
#include "boutexception.hxx"

#include <iostream>
#include <string>

TEST(BoutExceptionTest, WhatTest) {
  try {
    throw BoutException("Test message");
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
