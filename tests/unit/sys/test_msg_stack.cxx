// These tests rely on MsgStack::getDump, and so won't work without it
#if BOUT_USE_MSGSTACK

#include "gtest/gtest.h"
#include "msg_stack.hxx"
#include "test_extras.hxx"

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

TEST(MsgStackTest, PopValueTest) {
  MsgStack msg_stack;

  msg_stack.push("First");
  msg_stack.push("Second");
  msg_stack.push("Third");
  msg_stack.pop(2);

  auto dump = msg_stack.getDump();
  auto first_dump = "====== Back trace ======\n -> Second\n -> First\n";
  EXPECT_EQ(dump, first_dump);

  msg_stack.pop(-5); // Points to first message only
  dump = msg_stack.getDump();
  auto second_dump = "====== Back trace ======\n";
  EXPECT_EQ(dump, second_dump);
}

TEST(MsgStackTest, ReallocStorageTest) {
  MsgStack msg_stack;

  std::string expected_dump = "";

  for (int i = 0; i < 20; i++) {
    msg_stack.push("Message {:d}", i);
    expected_dump = " -> Message " + std::to_string(i) + "\n" + expected_dump;
  }
  expected_dump = "====== Back trace ======\n" + expected_dump;
  auto dump = msg_stack.getDump();
  EXPECT_EQ(dump, expected_dump);
}

TEST(MsgStackTest, PushReturnTest) {
  MsgStack msg_stack;

  for (int i = 0; i < 6; i++) {
    EXPECT_EQ(msg_stack.push("Message {:d}", i), i);
  }
}

TEST(MsgStackTest, NoMessageTest) {
  MsgStack msg_stack;
  msg_stack.push();
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
    auto second_dump = "====== Back trace ======\n -> Second on line " + line;

    EXPECT_TRUE(IsSubString(second, second_dump));
  }

  // Should now contain only the first message
  auto third = msg_stack.getDump();

  EXPECT_EQ(first_dump, third);
}

TEST(MsgStackTest, DumpTest) {
  // Code to capture output -- see test_output.cxx
  // Write cout to buffer instead of stdout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf(std::cout.rdbuf());
  std::cout.rdbuf(buffer.rdbuf());

  MsgStack msg_stack;

  msg_stack.push("First");
  msg_stack.dump();
  auto first_dump = "====== Back trace ======\n -> First\n";

  EXPECT_EQ(first_dump, buffer.str());

  // Clear buffer
  buffer.str("");
  // When done redirect cout to its old self
  std::cout.rdbuf(sbuf);
}

#endif
