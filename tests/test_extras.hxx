#ifndef TEST_EXTRAS_H__
#define TEST_EXTRAS_H__

#include "gtest/gtest.h"

::testing::AssertionResult IsSubString(const std::string &str,
                                       const std::string &substring);

#endif //  TEST_EXTRAS_H__
