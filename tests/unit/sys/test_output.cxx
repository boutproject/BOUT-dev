#include "gtest/gtest.h"
#include "boutexception.hxx"
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

  virtual ~OutputTest() {
    // Clear buffer
    buffer.str("");
    // When done redirect cout to its old self
    std::cout.rdbuf(sbuf);

    std::remove(filename.c_str());
  }

  // Write cout to buffer instead of stdout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf;
  // A temporary filename
  std::string filename{std::tmpnam(nullptr)};
};

TEST_F(OutputTest, JustStdOutCpp) {
  Output local_output;
  local_output << "Hello, world!" << 1 << std::endl;

  EXPECT_EQ(buffer.str(), "Hello, world!1\n");
}

TEST_F(OutputTest, JustStdOutPrintf) {
  Output local_output;
  local_output.write("{:s}{:d}\n", "Hello, world!", 2);

  EXPECT_EQ(buffer.str(), "Hello, world!2\n");
}

TEST_F(OutputTest, JustStdOutGlobalInstance) {
  output << "Hello, world!" << 3 << std::endl;

  EXPECT_EQ(buffer.str(), "Hello, world!3\n");
}

TEST_F(OutputTest, OpenFile) {
  Output local_output;

  std::string test_output = "To stdout and file\n";

  local_output.open(filename);
  local_output << test_output;

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  EXPECT_EQ(test_output, buffer.str());
  EXPECT_EQ(test_output, test_buffer.str());
}

TEST_F(OutputTest, JustPrint) {
  Output local_output;

  std::string test_output = "To stdout only\n";

  local_output.open(filename);
  local_output.print(test_output);

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  EXPECT_EQ("", test_buffer.str());
  EXPECT_EQ(test_output, buffer.str());
}

TEST_F(OutputTest, DisableEnableStdout) {
  Output local_output;

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
}

TEST_F(OutputTest, GetInstance) {
  Output *local_output = Output::getInstance();

  EXPECT_NE(local_output, nullptr);

  Output *new_output = Output::getInstance();

  EXPECT_EQ(local_output, new_output);
}

TEST_F(OutputTest, ConditionalGetBase) {
  Output local_output_base;
  ConditionalOutput local_output(&local_output_base);

  EXPECT_EQ(local_output.getBase(), &local_output_base);
}

TEST_F(OutputTest, ConditionalCheckIsEnabled) {
  Output local_output_base;
  ConditionalOutput local_output(&local_output_base);

  EXPECT_TRUE(local_output.isEnabled());
  local_output.disable();
  EXPECT_FALSE(local_output.isEnabled());
  local_output.enable();
  EXPECT_TRUE(local_output.isEnabled());
  local_output.enable(false);
  EXPECT_FALSE(local_output.isEnabled());
  local_output.enable(true);
  EXPECT_TRUE(local_output.isEnabled());
}

TEST_F(OutputTest, ConditionalJustStdOutCpp) {
  Output local_output_base;
  ConditionalOutput local_output(&local_output_base);

  local_output << "Hello, world!" << 4 << std::endl;

  EXPECT_EQ(buffer.str(), "Hello, world!4\n");
}

TEST_F(OutputTest, ConditionalJustStdOutPrintf) {
  Output local_output_base;
  ConditionalOutput local_output(&local_output_base);

  local_output.write("{:s}{:d}\n", "Hello, world!", 5);

  EXPECT_EQ(buffer.str(), "Hello, world!5\n");
}

TEST_F(OutputTest, ConditionalDisable) {
  Output local_output_base;
  ConditionalOutput local_output(&local_output_base);

  local_output.disable();
  local_output << "Hello, world!" << 6;
  local_output << std::endl;

  EXPECT_EQ(buffer.str(), "");
}

TEST_F(OutputTest, ConditionalJustStdOutGlobalInstances) {

  output_warn.enable();
  output_warn << "warn output\n";
  EXPECT_EQ(buffer.str(), "warn output\n");

  buffer.str("");
  output_info.enable();
  output_info << "info output\n";
  EXPECT_EQ(buffer.str(), "info output\n");

  buffer.str("");
  output_progress.enable();
  output_progress << "progress output\n";
  EXPECT_EQ(buffer.str(), "progress output\n");

  buffer.str("");
  output_error.enable();
  output_error << "error output\n";
  EXPECT_EQ(buffer.str(), "error output\n");

  buffer.str("");
  output_debug.enable();
  output_debug << "debug output\n";
#ifdef DEBUG_ENABLED
  EXPECT_EQ(buffer.str(), "debug output\n");
#else
  EXPECT_EQ(buffer.str(), "");
#endif
}

TEST_F(OutputTest, ConditionalJustPrint) {
  Output local_output_base;
  ConditionalOutput local_output(&local_output_base);

  std::string test_output = "To stdout only\n";

  local_output.open(filename);
  local_output.print(test_output);

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  EXPECT_EQ("", test_buffer.str());
  EXPECT_EQ(test_output, buffer.str());
}

TEST_F(OutputTest, ConditionalMultipleLayersGetBase) {
  Output local_output_base;
  ConditionalOutput local_output_first(&local_output_base);
  ConditionalOutput local_output_second(&local_output_first);

  EXPECT_EQ(local_output_first.getBase(), &local_output_base);
  EXPECT_NE(local_output_second.getBase(), &local_output_first);
  EXPECT_EQ(local_output_second.getBase(), &local_output_base);
}

TEST_F(OutputTest, ConditionalMultipleLayersJustStdOut) {
  Output local_output_base;
  ConditionalOutput local_output_first(&local_output_base);
  ConditionalOutput local_output_second(&local_output_first);

  local_output_second << "Hello, world!" << 7 << std::endl;

  EXPECT_EQ(buffer.str(), "Hello, world!7\n");
}

TEST_F(OutputTest, ConditionalMultipleLayersJustStdOutPrintf) {
  Output local_output_base;
  ConditionalOutput local_output_first(&local_output_base);
  ConditionalOutput local_output_second(&local_output_first);

  local_output_second.write("{:s}{:d}\n", "Hello, world!", 8);

  EXPECT_EQ(buffer.str(), "Hello, world!8\n");
}

TEST_F(OutputTest, DummyCheckEnableDoesntWork) {
  DummyOutput dummy;

  EXPECT_FALSE(dummy.isEnabled());
  dummy.enable();
  EXPECT_FALSE(dummy.isEnabled());
  dummy.enable(true);
  EXPECT_FALSE(dummy.isEnabled());
  dummy.enable(false);
  EXPECT_FALSE(dummy.isEnabled());
  dummy.disable();
  EXPECT_FALSE(dummy.isEnabled());
}

TEST_F(OutputTest, DummyOutputStdOut) {
  DummyOutput dummy;
  dummy << "Vanish to the void" << std::endl;
  dummy.write("Vanish to the void\n");

  EXPECT_EQ(buffer.str(), "");
}

TEST_F(OutputTest, DummyJustPrint) {
  DummyOutput dummy;

  std::string test_output = "To stdout only\n";

  dummy.open(filename);
  dummy.print(test_output);

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  EXPECT_EQ("", test_buffer.str());
  EXPECT_EQ("", buffer.str());
}
