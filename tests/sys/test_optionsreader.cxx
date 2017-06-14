#include "gtest/gtest.h"
#include "optionsreader.hxx"

#include "boutexception.hxx"
#include "utils.hxx"

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

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("-flag");

  reader.parseCommandLine(options, 2, argv);

  ASSERT_TRUE(options->isSet("flag"));

  bool value;
  options->get("flag", value, false, false);

  EXPECT_EQ(value, true);
}

TEST_F(OptionsReaderTest, ParseCommandLineWithSection) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  char **argv = static_cast<char **>(malloc(2));
  argv[0] = copy_string("prog");
  argv[1] = copy_string("subsection1:int_key=42");

  reader.parseCommandLine(options, 2, argv);

  ASSERT_TRUE(options->isSet("int_key"));

  int value;
  options->get("int_key", value, -1, false);

  EXPECT_EQ(value, 42);

  Options *section1 = options->getSection("subsection1");

  int sub_value;
  section1->get("int_key", sub_value, -1, false);

  EXPECT_EQ(sub_value, 42);
}
