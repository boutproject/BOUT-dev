#include "gtest/gtest.h"
#include "output.hxx"

#include <cstdio>
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

TEST_F(OutputTest, OpenFile) {
  Output local_output;

  // Get a filename for a temporary file
  char *filename = std::tmpnam(nullptr);

  std::string test_output = "To stdout and file\n";

  local_output.open(filename);
  local_output << test_output;

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  EXPECT_EQ(test_output, test_buffer.str());
  EXPECT_EQ(test_output, buffer.str());

  std::remove(filename);
}

TEST_F(OutputTest, JustPrint) {
  Output local_output;

  // Get a filename for a temporary file
  char *filename = std::tmpnam(nullptr);

  std::string test_output = "To stdout only\n";

  local_output.open(filename);
  local_output.print(test_output.c_str());

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  EXPECT_EQ("", test_buffer.str());
  EXPECT_EQ(test_output, buffer.str());

  std::remove(filename);
}

TEST_F(OutputTest, DisableEnableStdout) {
  Output local_output;

  // Get a filename for a temporary file
  char *filename = std::tmpnam(nullptr);

  std::string file_only = "To file only\n";
  std::string file_and_stdout = "To stdout and file\n";

  // Open temporary file and close stdout
  local_output.open(filename);
  local_output.disable();

  local_output << file_only;

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();

  EXPECT_EQ(file_only, test_buffer.str());
  EXPECT_EQ("", buffer.str());

  // Enable stdout again
  local_output.enable();
  local_output << file_and_stdout;

  test_buffer << test_file.rdbuf();

  // File should contain both outputs, stdout only latter
  EXPECT_EQ(file_only + file_and_stdout, test_buffer.str());
  EXPECT_EQ(file_and_stdout, buffer.str());

  test_file.close();
  std::remove(filename);
}

TEST_F(OutputTest, CleanupAndGetInstance) {
  Output *local_output = Output::getInstance();
  EXPECT_NE(local_output, nullptr);

  local_output->cleanup();

  // Get a new instance
  Output *new_output = Output::getInstance();

  *new_output << "Hello, world!\n";

  EXPECT_EQ(buffer.str(), "Hello, world!\n");
}
