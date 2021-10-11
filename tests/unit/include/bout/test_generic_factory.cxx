#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "boutexception.hxx"
#include "bout/generic_factory.hxx"

#include <exception>
#include <memory>
#include <string>
#include <vector>

class Options;

class Base {
public:
  Base(Options*) {}
  virtual ~Base() = default;
  virtual std::string foo() { return "Base"; }
};

class Derived1 : public Base {
public:
  Derived1(Options*) : Base({}) {}
  std::string foo() override { return "Derived1"; }
};

class Derived2 : public Base {
public:
  Derived2(Options*) : Base({}) {}
  std::string foo() override { return "Derived2"; }
};

namespace bout {
template <>
struct ArgumentHelper<Base> : ArgumentHelperBase {
  ArgumentHelper(Options& options)
      : something(options["something"].doc("Test option").withDefault(2)) {}
  int something;
  static PreconditionResult checkPreconditions(MAYBE_UNUSED(Options* options)) {
    return {true, ""};
  }
};

template <>
struct ArgumentHelper<Derived1> : ArgumentHelper<Base> {
  ArgumentHelper(Options& options)
      : ArgumentHelper<Base>(options),
        another_thing(options["another_thing"].doc("Another test option").withDefault(4)) {}
  int another_thing;
  static PreconditionResult checkPreconditions(Options* options) {
    if (options->isSet("nothing")) {
      return {false, "'Nothing' was set"};
    }
    return {true, ""};
  }
};
} // namespace bout

namespace {
class BaseFactory : public Factory<Base, BaseFactory, Options*> {
public:
  static constexpr auto type_name = "Base";
  static constexpr auto section_name = "base";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "base";
};
constexpr decltype(BaseFactory::type_name) BaseFactory::type_name;
constexpr decltype(BaseFactory::section_name) BaseFactory::section_name;
constexpr decltype(BaseFactory::option_name) BaseFactory::option_name;
constexpr decltype(BaseFactory::default_type) BaseFactory::default_type;

BaseFactory::RegisterInFactory<Base> registerme("base");
BaseFactory::RegisterInFactory<Derived1> registerme1("derived1");
BaseFactory::RegisterInFactory<Derived2> registerme2("derived2");
BaseFactory::RegisterUnavailableInFactory dontregisterme("not here",
                                                         "this is only a test");
} // namespace

class BaseComplicated {
public:
  std::string name;
  BaseComplicated(std::string name) : name(std::move(name)) {}
  virtual ~BaseComplicated() = default;
  virtual std::string foo() { return name; }
};

class DerivedComplicated1 : public BaseComplicated {
public:
  DerivedComplicated1(std::string name) : BaseComplicated(std::move(name)) {}
};

class DerivedComplicated2 : public BaseComplicated {
public:
  DerivedComplicated2(std::string name) : BaseComplicated(std::move(name)) {}
};

class ComplicatedFactory
    : public Factory<BaseComplicated, ComplicatedFactory, const std::string&> {
public:
  static constexpr auto type_name = "Complicated";
  static constexpr auto section_name = "complicated";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "basecomplicated";
};

constexpr decltype(ComplicatedFactory::type_name) ComplicatedFactory::type_name;
constexpr decltype(ComplicatedFactory::section_name) ComplicatedFactory::section_name;
constexpr decltype(ComplicatedFactory::option_name) ComplicatedFactory::option_name;
constexpr decltype(ComplicatedFactory::default_type) ComplicatedFactory::default_type;

namespace {
ComplicatedFactory::RegisterInFactory<BaseComplicated> registerme3("basecomplicated");
ComplicatedFactory::RegisterInFactory<DerivedComplicated1>
    registerme4("derivedcomplicated1");
ComplicatedFactory::RegisterInFactory<DerivedComplicated2>
    registerme5("derivedcomplicated2");
} // namespace

TEST(GenericFactory, RegisterAndCreate) {

  auto base_{BaseFactory::getInstance().create("base")};
  EXPECT_EQ(base_->foo(), "Base");

  auto derived1_{BaseFactory::getInstance().create("derived1")};
  EXPECT_EQ(derived1_->foo(), "Derived1");

  auto derived2_{BaseFactory::getInstance().create("derived2")};
  EXPECT_EQ(derived2_->foo(), "Derived2");
}

TEST(GenericFactory, ListAvailable) {
  auto available = BaseFactory::getInstance().listAvailable();
  std::vector<std::string> expected{"base", "derived1", "derived2"};

  EXPECT_EQ(available, expected);
}

TEST(GenericFactory, ListUnavailable) {
  auto available = BaseFactory::getInstance().listUnavailableReasons();
  std::vector<std::string> expected{"not here (this is only a test)"};

  EXPECT_EQ(available, expected);
}

TEST(GenericFactory, Remove) {
  auto available = BaseFactory::getInstance().listAvailable();
  std::vector<std::string> expected{"base", "derived1", "derived2"};

  EXPECT_EQ(available, expected);

  EXPECT_TRUE(BaseFactory::getInstance().remove("derived2"));

  std::vector<std::string> expected2{"base", "derived1"};
  available = BaseFactory::getInstance().listAvailable();

  EXPECT_EQ(available, expected2);

  // Better add this back in as this object is shared between tests!
  BaseFactory::RegisterInFactory<Derived2> registerme("derived2");
}

TEST(GenericFactory, RemoveNonexistant) {
  // Try a remove for something that shouldn't be registered
  EXPECT_FALSE(BaseFactory::getInstance().remove("derived83"));
}

TEST(GenericFactory, GetUnknownType) {
  EXPECT_THROW(BaseFactory::getInstance().create("unknown"), BoutException);
}

TEST(GenericFactory, GetUnavailableType) {
  EXPECT_THROW(BaseFactory::getInstance().create("not here"), BoutException);
}

TEST(GenericFactory, Help) {
  const auto help = BaseFactory::getInstance().help("derived1");
  EXPECT_THAT(help, testing::ContainsRegex("base:something"));
  EXPECT_THAT(help, testing::ContainsRegex("base:another_thing"));
}

TEST(GenericFactory, HelpForNonExistent) {
  EXPECT_THROW(BaseFactory::getInstance().help("not here"), BoutException);
}

TEST(GenericFactory, CheckPreconditionTrue) {
  Options options;
  const auto preconditions = BaseFactory::getInstance().checkPreconditions("derived1")(&options);
  EXPECT_TRUE(preconditions);
}

TEST(GenericFactory, CheckPreconditionFalse) {
  Options options;
  options["nothing"] = true;
  const auto preconditions = BaseFactory::getInstance().checkPreconditions("derived1")(&options);
  EXPECT_FALSE(preconditions);
  EXPECT_EQ(preconditions.reason, "'Nothing' was set");
}

TEST(GenericFactory, AssertPrecondition) {
  Options options;
  EXPECT_NO_THROW(bout::assertPreconditions<Derived1>(&options));
  options["nothing"] = true;
  EXPECT_THROW(bout::assertPreconditions<Derived1>(&options), BoutException);
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
