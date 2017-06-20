#include "gtest/gtest.h"
#include "test_extras.hxx"
#include "optionsreader.hxx"

#include "boutexception.hxx"
#include "utils.hxx"

#include <fstream>
#include <cstdio>

// stdout redirection code from https://stackoverflow.com/a/4043813/2043465
class OptionsReaderTest : public ::testing::Test {
public:
  OptionsReaderTest() : sbuf(std::cout.rdbuf()) {
    // Redirect cout to our stringstream buffer or any other ostream
    std::cout.rdbuf(buffer.rdbuf());
  }

  ~OptionsReaderTest() {
    // Clear buffer
    buffer.str("");
    // When done redirect cout to its old self
    std::cout.rdbuf(sbuf);

    // Make sure options singleton is clean
    Options::cleanup();
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

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=42=31");

  EXPECT_THROW(reader.parseCommandLine(options, 2, argv), BoutException);
}

TEST_F(OptionsReaderTest, BadCommandLineEmptyKey) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("=42");

  EXPECT_THROW(reader.parseCommandLine(options, 2, argv), BoutException);
}

TEST_F(OptionsReaderTest, BadCommandLineEmptyValue) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=");

  EXPECT_THROW(reader.parseCommandLine(options, 2, argv), BoutException);
}

TEST_F(OptionsReaderTest, ParseCommandLine) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=42");

  reader.parseCommandLine(options, 2, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);
}

TEST_F(OptionsReaderTest, ParseCommandLineGlobalInstance) {
  OptionsReader *reader = OptionsReader::getInstance();
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=42");

  reader->parseCommandLine(options, 2, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);
}

TEST_F(OptionsReaderTest, ParseCommandLineWithSpaces) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(4));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key");
  argv[2] = copy_string("=");
  argv[3] = copy_string("42");

  reader.parseCommandLine(options, 4, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);
}

TEST_F(OptionsReaderTest, ParseCommandLineWithTrailingSpace) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(3));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("int_key=");
  argv[2] = copy_string("42");

  reader.parseCommandLine(options, 3, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);
}

TEST_F(OptionsReaderTest, ParseCommandLineFlag) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(3));
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
}

TEST_F(OptionsReaderTest, ParseCommandLineWithSection) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("subsection1:int_key=42");

  reader.parseCommandLine(options, 2, argv);

  EXPECT_FALSE(options->isSet("int_key"));

  Options *section1 = options->getSection("subsection1");

  ASSERT_TRUE(section1->isSet("int_key"));

  int sub_value;
  section1->get("int_key", sub_value, -1, false);

  EXPECT_EQ(sub_value, 42);
}

TEST_F(OptionsReaderTest, ReadFile) {
  const std::string text = R"(
[section1]
int_key = 34
real_key = 42.34e-67
[section1:subsection2]
bool_key = false
flag
)";

  char *filename = std::tmpnam(nullptr);
  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  reader.read(options, filename);

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

  ASSERT_TRUE(subsection2->isSet("flag"));

  bool flag;
  subsection2->get("flag", flag, false);

  EXPECT_TRUE(flag);
}

TEST_F(OptionsReaderTest, ReadBadFile) {
  char *filename = std::tmpnam(nullptr);
  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, filename), BoutException);
}
