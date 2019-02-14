#include "gtest/gtest.h"

#include <bout/template_combinations.hxx>
#include <boutexception.hxx>
#include <output.hxx>

#include <typeindex>

using innerType = std::vector<std::string>;
using registrationType = std::vector<innerType>;

struct test1 {};
struct test2 {};

const std::map<std::type_index, std::string> humanName{
    {std::type_index(typeid(int)), "int"},
    {std::type_index(typeid(long)), "long"},
    {std::type_index(typeid(float)), "float"},
    {std::type_index(typeid(double)), "double"},
    {std::type_index(typeid(std::string)), "std::string"},
    {std::type_index(typeid(bool)), "boolean"},
    {std::type_index(typeid(test1)), "test1 struct"},
    {std::type_index(typeid(test2)), "test2 struct"},

};

template <typename T>
std::string typeName() {
  // Note because this routine is templated this is evaluated at compile time
  // and as such
  // we need to make sure that there is an entry for type T in humanName at
  // compile time
  // -- in other words we can't expand the map in the below routines.
  return humanName.at(std::type_index(typeid(T)));
}

// Helper struct that acts as the bottom level object that deals
// with the provided types
struct registerItems {
  template <typename T>
  void operator()(T) {
    innerType thisCase;
    thisCase.push_back(typeName<T>());
    registrations.push_back(thisCase);
  };
  template <typename T, typename T2>
  void operator()(T, T2) {
    innerType thisCase;
    thisCase.push_back(typeName<T>());
    thisCase.push_back(typeName<T2>());
    registrations.push_back(thisCase);
  };
  template <typename T, typename T2, typename T3>
  void operator()(T, T2, T3) {
    innerType thisCase;
    thisCase.push_back(typeName<T>());
    thisCase.push_back(typeName<T2>());
    thisCase.push_back(typeName<T3>());
    registrations.push_back(thisCase);
  };
  template <typename T, typename T2, typename T3, typename T4>
  void operator()(T, T2, T3, T4) {
    innerType thisCase;
    thisCase.push_back(typeName<T>());
    thisCase.push_back(typeName<T2>());
    thisCase.push_back(typeName<T3>());
    thisCase.push_back(typeName<T4>());
    registrations.push_back(thisCase);
  };

  registrationType& registrations;
};

TEST(TemplateCombinationsTest, SingleSetSingleTypeItem) {
  registrationType registrations;
  registerItems reg{registrations};

  produceCombinations<Set<int>> junk(reg);

  // Check number of registrations
  ASSERT_EQ(reg.registrations.size(), 1);
  // Check size of each registration
  for (const auto& entry : registrations) {
    ASSERT_EQ(entry.size(), 1);
  }

  EXPECT_EQ(reg.registrations[0][0], typeName<int>());
}

TEST(TemplateCombinationsTest, SingleSetTwoTypesItem) {
  registrationType registrations;
  registerItems reg{registrations};

  produceCombinations<Set<int, float>> junk(reg);

  // Check number of registrations
  ASSERT_EQ(reg.registrations.size(), 2);
  // Check size of each registration
  for (const auto& entry : registrations) {
    ASSERT_EQ(entry.size(), 1);
  }

  EXPECT_EQ(reg.registrations[0][0], typeName<int>());
  EXPECT_EQ(reg.registrations[1][0], typeName<float>());
}

TEST(TemplateCombinationsTest, TwoSetsOneTypeItem) {
  registrationType registrations;
  registerItems reg{registrations};

  produceCombinations<Set<int>, Set<float>> junk(reg);

  // Check number of registrations
  ASSERT_EQ(reg.registrations.size(), 1);
  // Check size of each registration
  for (const auto& entry : registrations) {
    ASSERT_EQ(entry.size(), 2);
  }

  EXPECT_EQ(reg.registrations[0][0], typeName<int>());
  EXPECT_EQ(reg.registrations[0][1], typeName<float>());
}

TEST(TemplateCombinationsTest, TwoSetsMultiTypeItem) {
  registrationType registrations;
  registerItems reg{registrations};

  produceCombinations<Set<int, long>, Set<float, double, std::string>> junk(reg);

  innerType firstSet{typeName<int>(), typeName<long>()};
  innerType secondSet{typeName<float>(), typeName<double>(), typeName<std::string>()};

  // Check number of registrations
  ASSERT_EQ(reg.registrations.size(), firstSet.size() * secondSet.size());
  // Check size of each registration
  for (const auto& entry : registrations) {
    ASSERT_EQ(entry.size(), 2);
  }

  int count = 0;
  for (const auto& first : firstSet) {
    EXPECT_EQ(reg.registrations[count][0], first);
    for (const auto& second : secondSet) {
      EXPECT_EQ(reg.registrations[count][1], second);
      count++;
    }
  };
}

TEST(TemplateCombinationsTest, FourSetsMultiTypeItem) {
  registrationType registrations;
  registerItems reg{registrations};

  produceCombinations<Set<int, long>, Set<float, double, std::string>, Set<bool>,
                      Set<test1, test2>>
      junk(reg);

  innerType firstSet{typeName<int>(), typeName<long>()};
  innerType secondSet{typeName<float>(), typeName<double>(), typeName<std::string>()};
  innerType thirdSet{typeName<bool>()};
  innerType fourthSet{typeName<test1>(), typeName<test2>()};

  // Check number of registrations
  ASSERT_EQ(reg.registrations.size(),
            firstSet.size() * secondSet.size() * thirdSet.size() * fourthSet.size());
  // Check size of each registration
  for (const auto& entry : registrations) {
    ASSERT_EQ(entry.size(), 4);
  }

  int count = 0;
  for (const auto& first : firstSet) {
    EXPECT_EQ(reg.registrations[count][0], first);
    for (const auto& second : secondSet) {
      EXPECT_EQ(reg.registrations[count][1], second);
      for (const auto& third : thirdSet) {
        EXPECT_EQ(reg.registrations[count][2], third);
        for (const auto& fourth : fourthSet) {
          EXPECT_EQ(reg.registrations[count][3], fourth);
          count++;
        }
      }
    }
  };
}
