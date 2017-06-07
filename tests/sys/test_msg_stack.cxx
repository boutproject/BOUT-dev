#include "gtest/gtest.h"
#include "msg_stack.hxx"

#include <iostream>
#include <string>

TEST(MsgStackTest, BasicTest) {

  MsgStack msg_stack;
  msg_stack.push("First");
  auto first = msg_stack.getDump();
  auto first_dump = "====== Back trace ======\n -> First\n";

  EXPECT_EQ(first_dump, first);
}

TEST(MsgStackTest, PopTest) {

  MsgStack msg_stack;
  msg_stack.push("First");
  msg_stack.pop();

  auto first = msg_stack.getDump();
  auto first_dump = "====== Back trace ======\n -> First\n";

  EXPECT_NE(first_dump, first);
}

TEST(MsgStackTest, CanUseGlobal) {

  // msg_stack is global variable defined in msg_stack.hxx
  msg_stack.push("First");
  auto first = msg_stack.getDump();
  auto first_dump = "====== Back trace ======\n -> First\n";

  msg_stack.pop();
  EXPECT_EQ(first_dump, first);
}

TEST(MsgStackTest, TraceMacroTest) {

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
                       " of 'sys/test_msg_stack.cxx'\n -> First\n";

    EXPECT_EQ(second_dump, second);
  }

  // Should now contain only the first message
  auto third = msg_stack.getDump();

  EXPECT_EQ(first_dump, third);
}
