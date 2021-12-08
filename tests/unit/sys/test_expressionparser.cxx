#include "gtest/gtest.h"

#include "bout/sys/expressionparser.hxx"
#include "bout_types.hxx"
#include "unused.hxx"
#include "test_extras.hxx"

#include <vector>

using bout::generator::Context;

// Need to inherit from ExpressionParser in order to expose the
// protected parseString as a public method
class ExpressionParserSubClass : public ExpressionParser {
public:
  using ExpressionParser::parseString;
  using ExpressionParser::fuzzyFind;
};

class ExpressionParserTest : public ::testing::Test {
public:
  ~ExpressionParserTest() override = default;
  ExpressionParserSubClass parser;
  std::vector<double> x_array = {-1., 0., 1., 5., 10., 3.14e8};
  std::vector<double> y_array = {-1., 0., 1., 5., 10., 3.14e8};
  std::vector<double> z_array = {-1., 0., 1., 5., 10., 3.14e8};
  std::vector<double> t_array = {-1., 0., 1., 5., 10., 3.14e8};
  WithQuietOutput quiet_warn{output_warn};
};

Context LegacyContext(BoutReal x, BoutReal y, BoutReal z, BoutReal t) {
  return Context().set("x", x, "y", y, "z", z, "t", t);
}

/// For testing, a generator function of two inputs
class BinaryGenerator : public FieldGenerator {
public:
  BinaryGenerator(std::shared_ptr<FieldGenerator> a = nullptr,
                  std::shared_ptr<FieldGenerator> b = nullptr)
      : a(std::move(a)), b(std::move(b)) {}

  std::shared_ptr<FieldGenerator>
  clone(const std::list<std::shared_ptr<FieldGenerator>> args) override {
    if (args.size() != 2) {
      throw ParseException(
          "Incorrect number of arguments to increment function. Expecting 2, got {:d}",
          args.size());
    }

    return std::make_shared<BinaryGenerator>(args.front(), args.back());
  }
  BoutReal generate(const Context& ctx) override {
    return a->generate(ctx) + b->generate(ctx);
  }
  std::string str() const override {
    return std::string{"add(" + a->str() + ", " + b->str() + ")"};
  }

private:
  std::shared_ptr<FieldGenerator> a, b;
};

class IncrementGenerator : public FieldGenerator {
public:
  IncrementGenerator(std::shared_ptr<FieldGenerator> gen = nullptr)
      : gen(std::move(gen)) {}

  std::shared_ptr<FieldGenerator>
  clone(const std::list<std::shared_ptr<FieldGenerator>> args) override {
    if (args.size() != 1) {
      throw ParseException(
          "Incorrect number of arguments to increment function. Expecting 1, got {:d}",
          args.size());
    }

    return std::make_shared<IncrementGenerator>(args.front());
  }
  BoutReal generate(const Context& ctx) override {
    return gen->generate(ctx) + 1;
  }
  std::string str() const override {
    return std::string{"increment(" + gen->str() + ")"};
  }

private:
  std::shared_ptr<FieldGenerator> gen;
};

// Function that takes no arguments and returns 4.0
class NullaryGenerator : public FieldGenerator {
public:
  NullaryGenerator() = default;

  std::shared_ptr<FieldGenerator>
  clone(const std::list<std::shared_ptr<FieldGenerator>> args) override {
    if (args.size() != 0) {
      throw ParseException(
          "Incorrect number of arguments to nullary function. Expecting 0, got {:d}",
          args.size());
    }

    return std::make_shared<NullaryGenerator>();
  }

  BoutReal generate(const Context &) override {
    return 4.0;
  }
};

TEST_F(ExpressionParserTest, Parse2) {
  auto fieldgen = parser.parseString("2");
  EXPECT_EQ(fieldgen->str(), "2");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), 2);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseX) {
  auto fieldgen = parser.parseString("x");
  EXPECT_EQ(fieldgen->str(), "x");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseY) {
  auto fieldgen = parser.parseString("y");
  EXPECT_EQ(fieldgen->str(), "y");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
	  EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), y);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseZ) {
  auto fieldgen = parser.parseString("z");
  EXPECT_EQ(fieldgen->str(), "z");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), z);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseT) {
  auto fieldgen = parser.parseString("t");
  EXPECT_EQ(fieldgen->str(), "t");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), t);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseXPlus2) {
  auto fieldgen = parser.parseString("x + 2");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x + 2);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseXTimesMinus4) {
  auto fieldgen = parser.parseString("x*(-4)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x * (-4));
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseXDividedBy3e8) {
  auto fieldgen = parser.parseString("x / 3.e8");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x / 3.e8);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ParseXSquared) {
  auto fieldgen = parser.parseString("x^2");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x * x);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, UnknownExpression) {
  EXPECT_THROW(parser.parseString("foo"), ParseException);
}

TEST_F(ExpressionParserTest, UnknownGenerator) {
  EXPECT_THROW(parser.parseString("aaa(1)"), ParseException);
}

TEST_F(ExpressionParserTest, BadNumbers) {
  EXPECT_THROW(parser.parseString("1.1.4"), ParseException);
  EXPECT_THROW(parser.parseString("2.1e4.5."), ParseException);
  EXPECT_THROW(parser.parseString("3.e8e4"), ParseException);
  EXPECT_THROW(parser.parseString("4()"), ParseException);
  EXPECT_THROW(parser.parseString("5_000"), ParseException);
}

TEST_F(ExpressionParserTest, BadFunctions) {
  parser.addGenerator("increment", std::make_shared<IncrementGenerator>());

  EXPECT_THROW(parser.parseString("increment(--)"), ParseException);
  EXPECT_THROW(parser.parseString("increment(x, )"), ParseException);
  EXPECT_THROW(parser.parseString("increment(x = 1)"), ParseException);
}

TEST_F(ExpressionParserTest, BadExpressions) {
  parser.addGenerator("increment", std::make_shared<IncrementGenerator>());

  EXPECT_THROW(parser.parseString("x = 1"), ParseException);
  EXPECT_THROW(parser.parseString("increment(x"), ParseException);
  EXPECT_THROW(parser.parseString("increment"), ParseException);
  EXPECT_THROW(parser.parseString("2]"), ParseException);
  EXPECT_THROW(parser.parseString("2)"), ParseException);
  EXPECT_THROW(parser.parseString("4+"), ParseException);
  EXPECT_THROW(parser.parseString("+4"), ParseException);
  EXPECT_THROW(parser.parseString("\n"), ParseException);
  EXPECT_THROW(parser.parseString("(3"), ParseException);
  EXPECT_THROW(parser.parseString("2-3[4"), ParseException);
  EXPECT_THROW(parser.parseString("[val = 42]{val}"), ParseException);
}

TEST_F(ExpressionParserTest, AddGenerator) {
  parser.addGenerator("increment", std::make_shared<IncrementGenerator>());

  auto fieldgen = parser.parseString("increment(x)");
  EXPECT_EQ(fieldgen->str(), "increment(x)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x + 1);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, AddFieldValue) {
  parser.addGenerator("fourty_two", std::make_shared<FieldValue>(42.0));

  auto fieldgen = parser.parseString("fourty_two(x)");

  // This next EXPECT depends on how stringstream renders doubles. We
  // may need to get rid of it if it proves fragile
  EXPECT_EQ(fieldgen->str(), "42");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), 42.0);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, AddNullaryFunction) {
  parser.addGenerator("nullary", std::make_shared<NullaryGenerator>());

  auto fieldgen = parser.parseString("nullary()");
  EXPECT_EQ(fieldgen->str(), "?");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), 4.0);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, CloneBinaryOp) {
  FieldBinary goodFieldBinary(nullptr, nullptr, '*');
  std::list<FieldGeneratorPtr> args;
  args.push_front(nullptr);
  args.push_front(nullptr);
  auto clonedFieldBinary = goodFieldBinary.clone(args);
  parser.addBinaryOp('&', clonedFieldBinary, 10);
  auto clonedFieldgen = parser.parseString("2 & (x + 3)");
  auto actualFieldgen = parser.parseString("2 * (x + 3)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
  	  auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(actualFieldgen->generate(ctx),
                           clonedFieldgen->generate(ctx));
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, BadCloneBinaryOp) {
  FieldBinary goodFieldBinary(nullptr, nullptr, '*');
  std::list<FieldGeneratorPtr> args;
  EXPECT_THROW(auto badFieldBinary = goodFieldBinary.clone(args), ParseException);
}

TEST_F(ExpressionParserTest, BadBinaryOp) {
  // Refers to an unrecognised binary operator "?"
  parser.addBinaryOp('&', std::make_shared<FieldBinary>(nullptr, nullptr, '?'), 5);
  auto fieldgen = parser.parseString("2 & x + 3");
  auto ctx = LegacyContext(0,0,0,0);

  EXPECT_THROW(fieldgen->generate(ctx), ParseException);
}

TEST_F(ExpressionParserTest, AddBinaryOp) {
  // Add a synonym for multiply with a lower precedence than addition
  parser.addBinaryOp('&', std::make_shared<FieldBinary>(nullptr, nullptr, '*'), 5);

  auto fieldgen = parser.parseString("2 & x + 3");
  EXPECT_EQ(fieldgen->str(), "(2*(x+3))");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), 2 * (x + 3));
        }
      }
    }
  }
}

TEST(ParseExceptionTest, WhatTest) {
  try {
    throw ParseException("test message");
  } catch (ParseException &e) {
    std::string message{e.what()};
    EXPECT_NE(message.find("test message"), std::string::npos);
  }
}

TEST_F(ExpressionParserTest, EscapeSymbol) {
  auto fieldgen = parser.parseString("`x`");
  EXPECT_EQ(fieldgen->str(), "x");
  
  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x);
        }
      }
    }
  }
}

// Backslash can be used to escape a single character
TEST_F(ExpressionParserTest, GeneratorNameEscape) {
  parser.addGenerator("one+", std::make_shared<IncrementGenerator>());

  auto fieldgen = parser.parseString("one\\+(x)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x + 1);
        }
      }
    }
  }
}

// Back-ticks can be used to escape sequences of characters
TEST_F(ExpressionParserTest, GeneratorNameLongEscape) {
  parser.addGenerator("++", std::make_shared<IncrementGenerator>());

  auto fieldgen = parser.parseString("`++`(x)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x + 1);
        }
      }
    }
  }
}

// Back-ticks can be used for part of a symbol
TEST_F(ExpressionParserTest, GeneratorNamePartEscape) {
  parser.addGenerator("one+this", std::make_shared<IncrementGenerator>());

  auto fieldgen = parser.parseString("one`+`this(x)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x + 1);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, AddBinaryGenerator) {
  parser.addGenerator("add", std::make_shared<BinaryGenerator>());

  auto fieldgen = parser.parseString("add(x,y)");
  EXPECT_EQ(fieldgen->str(), "add(x, y)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), x + y);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ImplicitMultiply) {

  auto fieldgen = parser.parseString("2x + 3y");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), 2*x + 3*y);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ImplicitMultiplyBracket) {

  auto fieldgen = parser.parseString("2(x + 3y)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx), 2*(x + 3*y));
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, BadImplicitMultiply) {
  EXPECT_THROW(parser.parseString("x2"), ParseException);
  EXPECT_THROW(parser.parseString("(1+x)2"), ParseException);
  EXPECT_THROW(parser.parseString("2 2"), ParseException);
}

TEST_F(ExpressionParserTest, PassParameter) {

  auto fieldgen = parser.parseString("{value}");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx.set("value", y + z)), y + z);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, MissingBrace) {
  EXPECT_THROW(parser.parseString("{"), ParseException);
  EXPECT_THROW(parser.parseString("2 + 3 * {something + 2"), ParseException);
  EXPECT_THROW(parser.parseString("2 + 3 * {something + {another}"), ParseException);
  EXPECT_THROW(parser.parseString("2 + 3 * something} + 2"), ParseException);
}

TEST_F(ExpressionParserTest, PassParameterImplicitMultiply) {

  auto fieldgen = parser.parseString("x - 3{value}");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(fieldgen->generate(ctx.set("value", 1 + y + z)), x - 3 * (1 + y + z));
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, PassMultipleParameters) {

  auto fieldgen = parser.parseString("x + {value} - 2*{other}");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          auto ctx = LegacyContext(x, y, z, t);
          EXPECT_DOUBLE_EQ(
              fieldgen->generate(ctx.set("value", 1 + y + z).set("other", x - y)),
              x + (1 + y + z) - 2 * (x - y));
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, ContextValue) {
  auto fieldgen = parser.parseString("[val = 42]({val})");
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 42);
}

TEST_F(ExpressionParserTest, ContextValueReplace) {
  auto fieldgen = parser.parseString("[val = 42]({val})");
  EXPECT_DOUBLE_EQ(fieldgen->generate(Context().set("val", 21)), 42);
}

TEST_F(ExpressionParserTest, ContextValueExpr) {
  auto fieldgen = parser.parseString("[val = 21]({val} + {val})");
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 42);
}

TEST_F(ExpressionParserTest, ContextValueExprTwoArgs) {
  auto fieldgen = parser.parseString("[val = 21, val2 = 13]({val} + {val2})");
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 34);
}

class GeneratorCloneCopy : public FieldGenerator {
public:
  explicit GeneratorCloneCopy(FieldGeneratorPtr expr) : expr(expr) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<GeneratorCloneCopy>(expr);
  }
  double generate(const Context& ctx) override {
    return expr->generate(ctx);
  }
private:
  FieldGeneratorPtr expr;
};

TEST_F(ExpressionParserTest, ContextFunction) {
  parser.addGenerator("func",
                      std::make_shared<GeneratorCloneCopy>(parser.parseString("2 * {x}")));

  auto fieldgen = parser.parseString("[x=3](func)");
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 6);
}

TEST_F(ExpressionParserTest, ContextLocal) {
  auto gen = parser.parseString("[x={x}-1]({x}) + {x}");
  EXPECT_DOUBLE_EQ(gen->generate(Context().set("x", 5)), 9); // 4 + 5
}

TEST_F(ExpressionParserTest, ContextTwice) {
  auto gen = parser.parseString("[x={x}-1]({x}) + [x={x}-2]({x})");
  EXPECT_DOUBLE_EQ(gen->generate(Context().set("x", 5)), 7); // 4 + 3
}

TEST_F(ExpressionParserTest, SumNothing) {
  auto fieldgen = parser.parseString("sum(i, 0, 42)");
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 0.0);
}

TEST_F(ExpressionParserTest, SumOne) {
  auto fieldgen = parser.parseString("sum(i, 1, 42)");
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 42);
}

TEST_F(ExpressionParserTest, SumExpr) {
  auto fieldgen = parser.parseString("sum(i, 2 + 1, {i}^2)"); // => 0^2 + 1^2 + 2^2
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 5);
}

TEST_F(ExpressionParserTest, SumNestedScope) {
  auto fieldgen = parser.parseString("sum(i, 3, sum(i, 2*{i}, {i}+1))"); // => (0) + (1 + 2) + (1 + 2 + 3 + 4) = 13
  EXPECT_DOUBLE_EQ(fieldgen->generate({}), 13);
}

TEST_F(ExpressionParserTest, FuzzyFind) {
  // We need some generators to lookup, but we don't care what they are
  parser.addGenerator("increment", {});
  parser.addGenerator("decrement", {});
  parser.addGenerator("multiply", {});
  parser.addGenerator("divide", {});

  auto matches = parser.fuzzyFind("recrement");
  EXPECT_EQ(matches.size(), 2);
  auto first_match = matches.begin();
  EXPECT_EQ(first_match->name, "decrement");
  EXPECT_EQ(first_match->distance, 1);

  auto CAPS_matches = parser.fuzzyFind("MULTIPLY");
  EXPECT_EQ(CAPS_matches.size(), 1);
  auto first_CAPS_match = CAPS_matches.begin();
  EXPECT_EQ(first_CAPS_match->name, "multiply");
  EXPECT_EQ(first_CAPS_match->distance, 1);
}
