#include "gtest/gtest.h"

#include "bout.hxx"
#include "boutexception.hxx"
#include "utils.hxx"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

std::vector<char*> get_c_string_vector(std::vector<std::string>& vec_args) {
  std::vector<char*> c_args{};
  std::transform(begin(vec_args), end(vec_args), std::back_inserter(c_args),
                 [](std::string& arg) { return &arg.front(); });
  return c_args;
}

TEST(ParseCommandLineArgs, HelpShortOption) {
  std::vector<std::string> v_args{"test", "-h"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  auto cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(std::cerr.rdbuf());

  EXPECT_EXIT(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
              ::testing::ExitedWithCode(0), "Usage:");

  std::cout.rdbuf(cout_buf);
}

TEST(ParseCommandLineArgs, HelpLongOption) {
  std::vector<std::string> v_args{"test", "--help"};
  auto c_args = get_c_string_vector(v_args);
  char** argv = c_args.data();

  auto cout_buf = std::cout.rdbuf();
  std::cout.rdbuf(std::cerr.rdbuf());

  EXPECT_EXIT(bout::experimental::parseCommandLineArgs(c_args.size(), argv),
              ::testing::ExitedWithCode(0), "Usage:");

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
