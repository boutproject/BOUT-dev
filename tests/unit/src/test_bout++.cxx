#include "bout/build_config.hxx"

#include "gtest/gtest.h"

#include "bout.hxx"
#include "boutexception.hxx"
#include "test_extras.hxx"
#include "utils.hxx"
#include "bout/version.hxx"

#include <algorithm>
#include <csignal>
#include <iostream>
#include <string>
#include <vector>

std::vector<char*> get_c_string_vector(std::vector<std::string>& vec_args) {
  std::vector<char*> c_args{};
  std::transform(begin(vec_args), end(vec_args), std::back_inserter(c_args),
                 [](std::string& arg) { return &arg.front(); });
  return c_args;
}

TEST(ParseCommandLineArgsDeathTest, HelpShortOption) {
  std::vector<std::string> v_args{"test", "-h"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  auto cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(std::cerr.rdbuf());

  EXPECT_EXIT(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
              ::testing::ExitedWithCode(0), _("Usage:"));

  std::cout.rdbuf(cout_buf);
}

TEST(ParseCommandLineArgsDeathTest, HelpLongOption) {
  std::vector<std::string> v_args{"test", "--help"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  auto cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(std::cerr.rdbuf());

  EXPECT_EXIT(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
              ::testing::ExitedWithCode(0), _("Usage:"));

  std::cout.rdbuf(cout_buf);
}

TEST(ParseCommandLineArgs, DataDir) {
  std::vector<std::string> v_args{"test", "-d", "test_data_directory"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.data_dir, "test_data_directory");
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, DataDirBad) {
  std::vector<std::string> v_args{"test", "-d"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  EXPECT_THROW(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
               BoutException);
}

TEST(ParseCommandLineArgs, OptionsFile) {
  std::vector<std::string> v_args{"test", "-f", "test_options_file"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.opt_file, "test_options_file");
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, OptionsFileBad) {
  std::vector<std::string> v_args{"test", "-f"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  EXPECT_THROW(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
               BoutException);
}

TEST(ParseCommandLineArgs, SettingsFile) {
  std::vector<std::string> v_args{"test", "-o", "test_settings_file"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.set_file, "test_settings_file");
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, SettingsFileBad) {
  std::vector<std::string> v_args{"test", "-o"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  EXPECT_THROW(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
               BoutException);
}

TEST(ParseCommandLineArgs, OptionsAndSettingsFilesSame) {
  std::vector<std::string> v_args{"test", "-o", "same", "-f", "same"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  EXPECT_THROW(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
               BoutException);
}

TEST(ParseCommandLineArgs, OptionsAndSettingsFilesDifferent) {
  std::vector<std::string> v_args{"test", "-o", "same", "-f", "different"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  EXPECT_NO_THROW(bout::experimental::parseCommandLineArgs(c_args.size(), argv));
}

TEST(ParseCommandLineArgs, LogFile) {
  std::vector<std::string> v_args{"test", "-l", "test_log_file"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.log_file, "test_log_file");
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, LogFileBad) {
  std::vector<std::string> v_args{"test", "-l"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  EXPECT_THROW(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
               BoutException);
}

TEST(ParseCommandLineArgs, VerbosityShort) {
  std::vector<std::string> v_args{"test", "-v"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 5);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, VerbosityShortMultiple) {
  std::vector<std::string> v_args{"test", "-v", "-v"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 6);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, VerbosityLong) {
  std::vector<std::string> v_args{"test", "--verbose"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 5);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, VerbosityLongMultiple) {
  std::vector<std::string> v_args{"test", "--verbose", "--verbose"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 6);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, QuietShort) {
  std::vector<std::string> v_args{"test", "-q"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 3);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, QuietShortMultiple) {
  std::vector<std::string> v_args{"test", "-q", "-q"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 2);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, QuietLong) {
  std::vector<std::string> v_args{"test", "--quiet"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 3);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, QuietLongMultiple) {
  std::vector<std::string> v_args{"test", "--quiet", "--quiet"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_EQ(args.verbosity, 2);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, ColorShort) {
  std::vector<std::string> v_args{"test", "-c"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_TRUE(args.color_output);
  EXPECT_EQ(args.original_argv, v_args);
}

TEST(ParseCommandLineArgs, ColorLong) {
  std::vector<std::string> v_args{"test", "--color"};
  auto v_args_copy = v_args;
  auto c_args = get_c_string_vector(v_args_copy);
  char** argv = c_args.data();

  auto args = bout::experimental::parseCommandLineArgs(c_args.size(), argv);

  EXPECT_TRUE(args.color_output);
  EXPECT_EQ(args.original_argv, v_args);
}

class PrintStartupTest : public ::testing::Test {
public:
  PrintStartupTest() : sbuf(std::cout.rdbuf()) {
    // Redirect cout to our stringstream buffer or any other ostream
    std::cout.rdbuf(buffer.rdbuf());
  }

  virtual ~PrintStartupTest() {
    // Clear buffer
    buffer.str("");
    // When done redirect cout to its old self
    std::cout.rdbuf(sbuf);
  }

  // Write cout to buffer instead of stdout
  std::stringstream buffer;
  // Save cout's buffer here
  std::streambuf* sbuf;
};

TEST_F(PrintStartupTest, Header) {
  bout::experimental::printStartupHeader(4, 8);

  EXPECT_TRUE(IsSubString(buffer.str(), bout::version::full));
  EXPECT_TRUE(IsSubString(buffer.str(), _("4 of 8")));
}

TEST_F(PrintStartupTest, CompileTimeOptions) {
  bout::experimental::printCompileTimeOptions();

  EXPECT_TRUE(IsSubString(buffer.str(), _("Compile-time options:\n")));
  EXPECT_TRUE(IsSubString(buffer.str(), _("Signal")));
  EXPECT_TRUE(IsSubString(buffer.str(), "netCDF"));
  EXPECT_TRUE(IsSubString(buffer.str(), "OpenMP"));
  EXPECT_TRUE(IsSubString(buffer.str(), _("Compiled with flags")));
}

TEST_F(PrintStartupTest, CommandLineArguments) {
  std::vector<std::string> args{"-d", "test1", "test2", "test3"};
  bout::experimental::printCommandLineArguments(args);

  for (auto& arg : args) {
    EXPECT_TRUE(IsSubString(buffer.str(), arg));
  }
}

#if BOUT_USE_SIGNAL

#if BOUT_USE_SIGFPE
#include <fenv.h>
#endif

class SignalHandlerTest : public ::testing::Test {
public:
  SignalHandlerTest() = default;
  virtual ~SignalHandlerTest() {
    std::signal(SIGUSR1, SIG_DFL);
    std::signal(SIGFPE, SIG_DFL);
    std::signal(SIGSEGV, SIG_DFL);
#if BOUT_USE_SIGFPE
    std::signal(SIGFPE, SIG_DFL);
    fedisableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#endif
  }
};

using SignalHandlerTestDeathTest = SignalHandlerTest;

TEST_F(SignalHandlerTestDeathTest, SegFault) {
  bout::experimental::setupSignalHandler(bout::experimental::defaultSignalHandler);
  // This test is *incredibly* expensive, maybe as much as 1s, so only test the one signal
  EXPECT_DEATH(std::raise(SIGSEGV), "SEGMENTATION FAULT");
}
#endif

TEST(BoutInitialiseFunctions, SetRunStartInfo) {
  WithQuietOutput quiet{output_info};

  Options options;

  bout::experimental::setRunStartInfo(options);

  auto run_section = options["run"];

  ASSERT_TRUE(run_section.isSection());
  EXPECT_TRUE(run_section.isSet("version"));
#ifdef REVISION
  EXPECT_TRUE(run_section.isSet("revision"));
#endif
  EXPECT_TRUE(run_section.isSet("started"));
}

TEST(BoutInitialiseFunctions, SetRunFinishInfo) {
  WithQuietOutput quiet{output_info};

  Options options;

  bout::experimental::setRunFinishInfo(options);

  ASSERT_TRUE(options["run"].isSection());
  EXPECT_TRUE(options["run"].isSet("finished"));
}

TEST(BoutInitialiseFunctions, CheckDataDirectoryIsAccessible) {
  using namespace bout::experimental;
  EXPECT_THROW(checkDataDirectoryIsAccessible("./bad/non/existent/directory"),
               BoutException);
  EXPECT_THROW(checkDataDirectoryIsAccessible(__FILE__), BoutException);
  EXPECT_NO_THROW(checkDataDirectoryIsAccessible("."));
}

TEST(BoutInitialiseFunctions, SavePIDtoFile) {
  WithQuietOutput quiet{output_info};

  EXPECT_NO_THROW(bout::experimental::savePIDtoFile("/tmp", 1));

  std::string filename{"/tmp/.BOUT.pid.1"};
  std::ifstream pid_file;
  pid_file.open(filename);

  EXPECT_TRUE(pid_file.good());

  std::stringstream contents;
  contents << pid_file.rdbuf();

  EXPECT_GT(contents.str().length(), 0);

  std::remove(filename.c_str());

  EXPECT_THROW(bout::experimental::savePIDtoFile("/does/likely/not/exists", 2), BoutException);
}
