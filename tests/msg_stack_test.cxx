#include "msg_stack.hxx"
#include "gtest/gtest.h"

#include <iostream>
#include <string>

TEST(MsgStackTest, BasicTest) {

  // msg_stack is global variable defined in msg_stack.hxx
  msg_stack.push("First");
  auto first = msg_stack.getDump();
  auto first_dump = "====== Back trace ======\n -> First\n";

  EXPECT_EQ(first_dump, first);
  {
    // The TRACE macro MUST follow this line to get the line number correct
    std::string line = std::to_string(__LINE__ + 1);
    TRACE("Second");
    auto second = msg_stack.getDump();
    auto second_dump = "====== Back trace ======\n -> Second on line " + line +
                       " of 'msg_stack_test.cxx'\n -> First\n";

    EXPECT_EQ(second_dump, second);
  }

  // Should now contain only the first message
  auto third = msg_stack.getDump();

  EXPECT_EQ(first_dump, third);
}
