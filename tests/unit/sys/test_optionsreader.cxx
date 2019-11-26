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
  }

  ~OptionsReaderTest() override {
    // Clear buffer
    buffer.str("");
    // When done redirect cout to its old self
    std::cout.rdbuf(sbuf);

    // Make sure options singleton is clean
    Options::cleanup();

    std::remove(filename.c_str());
  }

  // Write cout to buffer instead of stdout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf *sbuf;

  WithQuietOutput quiet{output_info};
  // A temporary filename
  std::string filename{std::tmpnam(nullptr)};
};

TEST_F(OptionsReaderTest, BadFilename) {
  OptionsReader reader;
  EXPECT_THROW(reader.read(nullptr, nullptr), BoutException);
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

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  reader.read(options, "%s", filename.c_str());

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
}

TEST_F(OptionsReaderTest, ReadBadFile) {
  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, "%s", filename.c_str()), BoutException);
}

TEST_F(OptionsReaderTest, ReadBadFileSectionIncomplete) {
  const std::string text = R"(
[section1
int_key = 34
)";

  std::ofstream test_file(filename.c_str(), std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, "%s", filename.c_str()), BoutException);
};

TEST_F(OptionsReaderTest, ReadBadFileSectionEmptyName) {
  const std::string text = R"(
[]
int_key = 34
)";

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  Options *options = Options::getRoot();
  EXPECT_THROW(reader.read(options, "%s", filename.c_str()), BoutException);
};

TEST_F(OptionsReaderTest, WriteFile) {
  OptionsReader reader;
  Options *options = Options::getRoot();

  options->set("bool_key", true, "test");
  Options *section1 = options->getSection("section1");
  section1->set("int_key", 17, "test");
  section1->set("real_key", 6.17e23, "test");
  Options *subsection2 = section1->getSection("subsection2");
  subsection2->set("string_key", "BOUT++", "test");

  reader.write(options, "%s", filename.c_str());

  std::ifstream test_file(filename);
  std::stringstream test_buffer;
  test_buffer << test_file.rdbuf();
  test_file.close();

  std::vector<std::string> expected = {"bool_key = true",        "[section1]",
                                       "int_key = 17",           "real_key = 6.17e+23",
                                       "[section1:subsection2]", "string_key = BOUT++"};

  for (auto &result : expected) {
    EXPECT_TRUE(IsSubString(test_buffer.str(), result));
  }
}

TEST_F(OptionsReaderTest, WriteBadFile) {
  std::string filename1 = filename + std::tmpnam(nullptr);
  OptionsReader reader;
  Options *options = Options::getRoot();

  options->set("bool_key", true, "test");
  Options *section1 = options->getSection("section1");
  section1->set("int_key", 17, "test");

  EXPECT_THROW(reader.write(options, "%s", filename1.c_str()), BoutException);

  std::remove(filename1.c_str());
}

TEST_F(OptionsReaderTest, ReadEmptyString) {
const std::string text = R"(
value =
)";

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();
  
  Options opt;
  OptionsReader reader;

  reader.read(&opt, "%s", filename.c_str());

  std::string val = opt["value"];
  EXPECT_TRUE(val.empty());
}

TEST_F(OptionsReaderTest, ReadFileEscapedChars) {
  const std::string text = R"(
some-value = 3  # Names can contain symbols, but need to be escaped

[h2+]           # Sections can contain symbols
another = 5
one-more = 2

[tests]
test1 = some\-value                  # Escaping of single characters
test2 = h2\+:another * some\-value   # Including in section names
test3 = `some-value` + 1             # Escaping character sequence
test4 = `h2+:one-more` * 2           # Escape sequence including :
test5 = `h2+`:another                # Escape only start of sequence
test6 = h2`+`:on`e-`more             # Escape sequences in the middle
)";

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  reader.read(Options::getRoot(), "%s", filename.c_str());

  auto options = Options::root()["tests"];
  
  EXPECT_EQ(options["test1"].as<int>(), 3);
  EXPECT_EQ(options["test2"].as<int>(), 15);
  EXPECT_EQ(options["test3"].as<int>(), 4);
  EXPECT_EQ(options["test4"].as<int>(), 4);
  EXPECT_EQ(options["test5"].as<int>(), 5);
  EXPECT_EQ(options["test6"].as<int>(), 2);
}

// Variable names must not contain colons, escaped or otherwise
// Sections can contain colons, where they indicate subsections.
// That is tested elsewhere.
TEST_F(OptionsReaderTest, ReadFileVariablesNoColons) {
  const std::string text = R"(
some:value = 3 
)";

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  
  EXPECT_THROW(reader.read(Options::getRoot(), "%s", filename.c_str()), BoutException);
}

TEST_F(OptionsReaderTest, ReadUnicodeNames) {
  const std::string text = R"(

α = 1.3
重要的數字 = 3

[tests]
結果 = 重要的數字 + 5
value = α*(1 + 重要的數字)
twopi = 2 * π   # Unicode symbol defined for pi
)";

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  reader.read(Options::getRoot(), "%s", filename.c_str());

  auto options = Options::root()["tests"];
  
  EXPECT_EQ(options["結果"].as<int>(), 8);
  EXPECT_DOUBLE_EQ(options["value"].as<BoutReal>(), 1.3*(1+3));
  EXPECT_DOUBLE_EQ(options["twopi"].as<BoutReal>(), 2 * 3.141592653589793);
}

TEST_F(OptionsReaderTest, ReadMultiLine) {
  const std::string text = R"(

result = (1 +
          2 + 
          3)

value = [a = 1,
         b = 2 * 2](
           {a} + {b})

)";

  std::ofstream test_file(filename, std::ios::out);
  test_file << text;
  test_file.close();

  OptionsReader reader;
  reader.read(Options::getRoot(), filename.c_str());

  auto options = Options::root();
  
  EXPECT_EQ(options["result"].as<int>(), 6);
  EXPECT_EQ(options["value"].as<int>(), 5);
}
