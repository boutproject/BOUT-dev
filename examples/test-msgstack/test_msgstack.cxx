#include "bout/assert.hxx"
#include "msg_stack.hxx"

#include <iostream>
#include <string>

int main() {

  // msg_stack is global variable defined in msg_stack.hxx
  msg_stack.push("First");
  auto first = msg_stack.getDump();

  std::cout << first;
  ASSERT0(first == "====== Back trace ======\n -> First\n");

  {
    // The TRACE macro MUST follow this line to get the line number correct
    std::string line = std::to_string(__LINE__ + 1);
    TRACE("Second");
    auto second = msg_stack.getDump();
    auto second_dump = "====== Back trace ======\n -> Second on line " + line +
                       " of 'test_msgstack.cxx'\n -> First\n";
    std::cout << second;
    ASSERT0(second == second_dump);
  }

  // Should now contain only the first message
  auto third = msg_stack.getDump();

  std::cout << third;
  ASSERT0(third == "====== Back trace ======\n -> First\n");
}
