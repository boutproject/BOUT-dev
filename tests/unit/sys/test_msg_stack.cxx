// These tests rely on MsgStack::getDump, and so won't work without it
#if CHECK > 1

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

  auto dump = msg_stack.getDump();
  auto expected_dump = "====== Back trace ======\n";

  EXPECT_EQ(dump, expected_dump);
}

TEST(MsgStackTest, ClearTest) {
  MsgStack msg_stack;

  msg_stack.push("First");
  msg_stack.push("Second");
  msg_stack.clear();

  auto dump = msg_stack.getDump();
  auto expected_dump = "====== Back trace ======\n";

  EXPECT_EQ(dump, expected_dump);
}

TEST(MsgStackTest, CanUseGlobal) {
  msg_stack.clear();

  // msg_stack is global variable defined in msg_stack.hxx
  msg_stack.push("First");
  auto first = msg_stack.getDump();
  auto first_dump = "====== Back trace ======\n -> First\n";

  msg_stack.pop();
  EXPECT_EQ(first_dump, first);
}

TEST(MsgStackTest, TraceMacroTest) {
  msg_stack.clear();

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
#endif
