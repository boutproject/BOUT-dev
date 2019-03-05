#include "gtest/gtest.h"

#include "boutexception.hxx"
#include "bout/generic_factory.hxx"

#include <exception>
#include <memory>
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

class BaseComplicated {
public:
  std::string name;
  BaseComplicated(std::string name) : name(name) {}
  virtual std::string foo() { return name; }
};

class DerivedComplicated1 : public BaseComplicated {
public:
  DerivedComplicated1(std::string name) : BaseComplicated(name) {}
};

class DerivedComplicated2 : public BaseComplicated {
public:
  DerivedComplicated2(std::string name) : BaseComplicated(name) {}
};

// Save some typing later
typedef Factory<BaseComplicated, std::function<BaseComplicated *(const std::string &)>>
    ComplicatedFactory;

// We need to specialise the helper class to pass arguments to the constructor
template<typename DerivedType>
class RegisterInFactory<BaseComplicated, DerivedType> {
public:
  RegisterInFactory(const std::string &name) {
    ComplicatedFactory::getInstance().add(
        name, [](std::string name) -> BaseComplicated * { return new DerivedType(name); });
  }
};

namespace {
RegisterInFactory<BaseComplicated, BaseComplicated>
    registerme3("basecomplicated");
RegisterInFactory<BaseComplicated, DerivedComplicated1>
    registerme4("derivedcomplicated1");
RegisterInFactory<BaseComplicated, DerivedComplicated2>
    registerme5("derivedcomplicated2");
} // namespace

TEST(GenericFactory, RegisterAndCreate) {

  std::unique_ptr<Base> base_{Factory<Base>::getInstance().create("base")};
  EXPECT_EQ(base_->foo(), "Base");

  std::unique_ptr<Base> derived1_{Factory<Base>::getInstance().create("derived1")};
  EXPECT_EQ(derived1_->foo(), "Derived1");

  std::unique_ptr<Base> derived2_{Factory<Base>::getInstance().create("derived2")};
  EXPECT_EQ(derived2_->foo(), "Derived2");
}

TEST(GenericFactory, ListAvailable) {
  auto available = Factory<Base>::getInstance().listAvailable();
  std::vector<std::string> expected{"base", "derived1", "derived2"};

  EXPECT_EQ(available, expected);
}

TEST(GenericFactory, Remove) {
  auto available = Factory<Base>::getInstance().listAvailable();
  std::vector<std::string> expected{"base", "derived1", "derived2"};

  EXPECT_EQ(available, expected);

  EXPECT_TRUE(Factory<Base>::getInstance().remove("derived2"));

  std::vector<std::string> expected2{"base", "derived1"};
  available = Factory<Base>::getInstance().listAvailable();

  EXPECT_EQ(available, expected2);

  // Better add this back in as this object is shared between tests!
  RegisterInFactory<Base, Derived2> registerme("derived2");
}

TEST(GenericFactory, RemoveNonexistant) {
  // Try a remove for something that shouldn't be registered
  EXPECT_FALSE(Factory<Base>::getInstance().remove("derived83"));
}

TEST(GenericFactory, GetUnknownType) {
  EXPECT_THROW(Factory<Base>::getInstance().create("unknown"), BoutException);
}

TEST(GenericFactory, Complicated) {

  std::unique_ptr<BaseComplicated> base_{
      ComplicatedFactory::getInstance().create("basecomplicated", "BaseComplicated")};
  EXPECT_EQ(base_->foo(), "BaseComplicated");

  std::unique_ptr<BaseComplicated> derived1_{ComplicatedFactory::getInstance().create(
      "derivedcomplicated1", "DerivedComplicated1")};
  EXPECT_EQ(derived1_->foo(), "DerivedComplicated1");

  std::unique_ptr<BaseComplicated> derived2_{ComplicatedFactory::getInstance().create(
      "derivedcomplicated2", "DerivedComplicated2")};
  EXPECT_EQ(derived2_->foo(), "DerivedComplicated2");
}
