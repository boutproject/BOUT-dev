#include "gtest/gtest.h"
#include "output.hxx"

#include <string>

// stdout redirection code from https://stackoverflow.com/a/4043813/2043465
class OutputTest : public ::testing::Test {
public:
  OutputTest() : sbuf(std::cout.rdbuf()) {
    // Redirect cout to our stringstream buffer or any other ostream
    std::cout.rdbuf(buffer.rdbuf());
  }

  ~OutputTest() {
    // Clear buffer
    buffer.str("");
    // When done redirect cout to its old self
    std::cout.rdbuf(sbuf);
  }

  // Write cout to buffer instead of stdout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf;
};

TEST_F(OutputTest, JustStdOutCpp) {
  Output local_output;
  local_output << "Hello, world!" << 1 << "\n";

  EXPECT_EQ(buffer.str(), "Hello, world!1\n");
}

TEST_F(OutputTest, JustStdOutPrintf) {
  Output local_output;
  local_output.write("%s%d\n", "Hello, world!", 2);

  EXPECT_EQ(buffer.str(), "Hello, world!2\n");
}

TEST_F(OutputTest, JustStdOutGlobalInstance) {
  output << "Hello, world!\n";

  EXPECT_EQ(buffer.str(), "Hello, world!\n");
}

TEST_F(OutputTest, DisableEnableStdout) {
  Output local_output;

  local_output.disable();
  local_output << "disabled\n";

  EXPECT_NE(buffer.str(), "disabled\n");

  local_output.enable();
  local_output << "enabled\n";

  EXPECT_EQ(buffer.str(), "enabled\n");
}
