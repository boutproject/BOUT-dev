#include "gtest/gtest.h"

#include "bout/generic_factory.hxx"

#include <string>
#include <vector>

class Base {
public:
  Base() {}
  virtual std::string foo() { return "Base"; }
};

class Derived1 : public Base {
public:
  Derived1() {}
  std::string foo() override { return "Derived1"; }
};

class Derived2 : public Base {
public:
  Derived2() {}
  std::string foo() override { return "Derived2"; }
};

namespace {
RegisterInFactory<Base, Base> registerme("base");
RegisterInFactory<Base, Derived1> registerme1("derived1");
RegisterInFactory<Base, Derived2> registerme2("derived2");
} // namespace

TEST(GenericFactory, RegisterAndCreate) {

  auto base_ = Factory<Base>::getInstance().create("base");
  EXPECT_EQ(base_->foo(), "Base");

  auto derived1_ = Factory<Base>::getInstance().create("derived1");
  EXPECT_EQ(derived1_->foo(), "Derived1");

  auto derived2_ = Factory<Base>::getInstance().create("derived2");
  EXPECT_EQ(derived2_->foo(), "Derived2");
}

TEST(GenericFactory, ListAvailable) {
  auto available = Factory<Base>::getInstance().listAvailable();
  std::vector<std::string> expected{"base", "derived1", "derived2"};

  EXPECT_EQ(available, expected);
}
