#include "gtest/gtest.h"
#include "test_extras.hxx"
#include "optionsreader.hxx"

#include "boutexception.hxx"
#include "output.hxx"
#include "utils.hxx"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// stdout redirection code from https://stackoverflow.com/a/4043813/2043465
class OptionsReaderTest : public ::testing::Test {
public:
  OptionsReaderTest() : sbuf(std::cout.rdbuf()) {
    // Redirect cout to our stringstream buffer or any other ostream
    std::cout.rdbuf(buffer.rdbuf());
    output_info.disable();
  }

  ~OptionsReaderTest() {
    // Clear buffer
    buffer.str("");
    // When done redirect cout to its old self
    std::cout.rdbuf(sbuf);

    // Make sure options singleton is clean
    Options::cleanup();

    output_info.enable();
  }

  // Write cout to buffer instead of stdout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf;
};

TEST_F(OptionsReaderTest, BadFilename) {
  OptionsReader reader;
  EXPECT_THROW(reader.read(nullptr, NULL), BoutException);
}

TEST_F(OptionsReaderTest, BadCommandLineMultipleEquals) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=42=31");

  EXPECT_THROW(reader.parseCommandLine(options, 2, argv), BoutException);

  free(argv[0]);
  free(argv[1]);
  free(argv);
}

TEST_F(OptionsReaderTest, BadCommandLineEmptyKey) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("=42");

  EXPECT_THROW(reader.parseCommandLine(options, 2, argv), BoutException);

  free(argv[0]);
  free(argv[1]);
  free(argv);
}

TEST_F(OptionsReaderTest, BadCommandLineEmptyValue) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=");

  EXPECT_THROW(reader.parseCommandLine(options, 2, argv), BoutException);

  free(argv[0]);
  free(argv[1]);
  free(argv);
}

TEST_F(OptionsReaderTest, ParseCommandLine) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=42");

  reader.parseCommandLine(options, 2, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);

  free(argv[0]);
  free(argv[1]);
  free(argv);
}

TEST_F(OptionsReaderTest, ParseCommandLineGlobalInstance) {
  OptionsReader *reader = OptionsReader::getInstance();
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=42");

  reader->parseCommandLine(options, 2, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);

  reader->cleanup();

  free(argv[0]);
  free(argv[1]);
  free(argv);
}

TEST_F(OptionsReaderTest, ParseCommandLineWithSpaces) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(4*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key");
  argv[2] = copy_string("=");
  argv[3] = copy_string("42");

  reader.parseCommandLine(options, 4, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);

  free(argv[0]);
  free(argv[1]);
  free(argv[2]);
  free(argv[3]);
  free(argv);
}

TEST_F(OptionsReaderTest, ParseCommandLineWithTrailingSpace) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(3*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=");
  argv[2] = copy_string("42");

  reader.parseCommandLine(options, 3, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);

  free(argv[0]);
  free(argv[1]);
  free(argv[2]);
  free(argv);
}

TEST_F(OptionsReaderTest, ParseCommandLineFlag) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(3*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("-flag");
  argv[2] = copy_string("command");

  reader.parseCommandLine(options, 3, argv);

  ASSERT_TRUE(options->isSet("flag"));

  bool flag;
  options->get("flag", flag, false, false);

  EXPECT_EQ(flag, true);

  ASSERT_TRUE(options->isSet("command"));

  bool command;
  options->get("command", command, false, false);

  EXPECT_EQ(command, true);

  free(argv[0]);
  free(argv[1]);
  free(argv[2]);
  free(argv);
}

TEST_F(OptionsReaderTest, ParseCommandLineWithSection) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2*(sizeof(char*))));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("subsection1:int_key=42");

  reader.parseCommandLine(options, 2, argv);

  EXPECT_FALSE(options->isSet("int_key"));

  Options *section1 = options->getSection("subsection1");

  ASSERT_TRUE(section1->isSet("int_key"));

  int sub_value;
  section1->get("int_key", sub_value, -1, false);

  EXPECT_EQ(sub_value, 42);

  free(argv[0]);
  free(argv[1]);
  free(argv);
}

TEST_F(OptionsReaderTest, ReadFile) {
  const std::string text = R"(
flag
[section1]
int_key = 34
real_key = 42.34e-67
[section1:subsection2]
bool_key = false
)";

  char *filename = std::tmpnam(nullptr);
  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  reader.read(options, filename);

  ASSERT_TRUE(options->isSet("flag"));

  bool flag;
  options->get("flag", flag, false);

  EXPECT_TRUE(flag);

  EXPECT_FALSE(options->isSet("int_key"));

  Options *section1 = options->getSection("section1");

  ASSERT_TRUE(section1->isSet("int_key"));

  int int_value;
  section1->get("int_key", int_value, -1, false);

  EXPECT_EQ(int_value, 34);

  ASSERT_TRUE(section1->isSet("real_key"));

  BoutReal real_value;
  section1->get("real_key", real_value, -1, false);

  EXPECT_DOUBLE_EQ(real_value, 42.34e-67);

  Options *subsection2 = section1->getSection("subsection2");

  ASSERT_TRUE(subsection2->isSet("bool_key"));

  bool bool_value;
  subsection2->get("bool_key", bool_value, true);

  EXPECT_FALSE(bool_value);

  std::remove(filename);
}

TEST_F(OptionsReaderTest, ReadBadFile) {
  char *filename = std::tmpnam(nullptr);
  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, filename), BoutException);
}

TEST_F(OptionsReaderTest, ReadBadFileSectionIncomplete) {
  const std::string text = R"(
[section1
int_key = 34
)";

  char *filename = std::tmpnam(nullptr);
  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, filename), BoutException);
};

TEST_F(OptionsReaderTest, ReadBadFileSectionEmptyName) {
  const std::string text = R"(
[]
int_key = 34
)";

  char *filename = std::tmpnam(nullptr);
  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, filename), BoutException);
};

TEST_F(OptionsReaderTest, WriteFile) {
  char *filename = std::tmpnam(nullptr);
  OptionsReader reader;
  Options *options = Options::getRoot();

  options->set("bool_key", true, "test");
  Options *section1 = options->getSection("section1");
  section1->set("int_key", 17, "test");
  section1->set("real_key", 6.17e23, "test");
  Options *subsection2 = section1->getSection("subsection2");
  subsection2->set("string_key", "BOUT++", "test");

  reader.write(options, filename);

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  std::vector<std::string> expected = {"bool_key = true",        "[section1]",
                                       "int_key = 17",           "real_key = 6.17000000000000006e+23",
                                       "[section1:subsection2]", "string_key = BOUT++"};

  for (auto &result : expected) {
    EXPECT_TRUE(IsSubString(test_buffer.str(), result));
  }

  std::remove(filename);
}

TEST_F(OptionsReaderTest, WriteBadFile) {
  std::string filename1 = std::tmpnam(nullptr);
  std::string filename = filename1 + std::tmpnam(nullptr);
  OptionsReader reader;
  Options *options = Options::getRoot();

  options->set("bool_key", true, "test");
  Options *section1 = options->getSection("section1");
  section1->set("int_key", 17, "test");

  EXPECT_THROW(reader.write(options, filename.c_str()), BoutException);

  std::remove(filename.c_str());
}
