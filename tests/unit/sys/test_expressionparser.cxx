#include "gtest/gtest.h"

#include "bout/sys/expressionparser.hxx"
#include "bout_types.hxx"
#include "unused.hxx"

#include <vector>

// Need to inherit from ExpressionParser in order to expose the
// protected parseString as a public method
class ExpressionParserSubClass : public ExpressionParser {
public:
  using ExpressionParser::parseString;
};

class ExpressionParserTest : public ::testing::Test {
public:
  ExpressionParserSubClass parser;
  std::vector<double> x_array = {-1., 0., 1., 5., 10., 3.14e8};
  std::vector<double> y_array = {-1., 0., 1., 5., 10., 3.14e8};
  std::vector<double> z_array = {-1., 0., 1., 5., 10., 3.14e8};
  std::vector<double> t_array = {-1., 0., 1., 5., 10., 3.14e8};
};

class IncrementGenerator : public FieldGenerator {
public:
  IncrementGenerator(std::shared_ptr<FieldGenerator> gen = nullptr) : gen(gen) {}

  std::shared_ptr<FieldGenerator>
  clone(const std::list<std::shared_ptr<FieldGenerator>> args) {
    if (args.size() != 1) {
      throw ParseException(
          "Incorrect number of arguments to increment function. Expecting 1, got %d",
          args.size());
    }

    return std::make_shared<IncrementGenerator>(args.front());
  }

  BoutReal generate(BoutReal x, BoutReal y, BoutReal z, BoutReal t) {
    return gen->generate(x, y, z, t) + 1;
  }
  const std::string str() { return std::string{"increment(" + gen->str() + ")"}; }

private:
  std::shared_ptr<FieldGenerator> gen;
};

// Function that takes no arguments and returns 4.0
class NullaryGenerator : public FieldGenerator {
public:
  NullaryGenerator() {}

  std::shared_ptr<FieldGenerator>
  clone(const std::list<std::shared_ptr<FieldGenerator>> args) {
    if (args.size() != 0) {
      throw ParseException(
          "Incorrect number of arguments to nullary function. Expecting 0, got %d",
          args.size());
    }

    return std::make_shared<NullaryGenerator>();
  }

  BoutReal generate(BoutReal UNUSED(x), BoutReal UNUSED(y), BoutReal UNUSED(z),
                    BoutReal UNUSED(t)) {
    return 4.0;
  }
  const std::string str() { return std::string{"nullary()"}; }
};

TEST_F(ExpressionParserTest, Parse2) {
  auto fieldgen = parser.parseString("2");
  EXPECT_EQ(fieldgen->str(), "2");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), 2);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), x);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), y);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), z);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), t);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), x + 2);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), x * (-4));
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), x / 3.e8);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), x * x);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, UnknownExpression) {
  EXPECT_THROW(parser.parseString("foo"), ParseException);
}

TEST_F(ExpressionParserTest, BadNumbers) {
  EXPECT_THROW(parser.parseString("1.1.4"), ParseException);
  EXPECT_THROW(parser.parseString("2.1e4.5."), ParseException);
  EXPECT_THROW(parser.parseString("3.e8e4"), ParseException);
  EXPECT_THROW(parser.parseString("4ee"), ParseException);
  EXPECT_THROW(parser.parseString("5G"), ParseException);
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
  EXPECT_THROW(parser.parseString("4+"), ParseException);
  EXPECT_THROW(parser.parseString("+4"), ParseException);
  EXPECT_THROW(parser.parseString("\n"), ParseException);
  EXPECT_THROW(parser.parseString("(3"), ParseException);
  EXPECT_THROW(parser.parseString("2-3[4"), ParseException);
}

TEST_F(ExpressionParserTest, AddGenerator) {
  parser.addGenerator("increment", std::make_shared<IncrementGenerator>());

  auto fieldgen = parser.parseString("increment(x)");
  EXPECT_EQ(fieldgen->str(), "increment(x)");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), x + 1);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), 42.0);
        }
      }
    }
  }
}

TEST_F(ExpressionParserTest, AddNullaryFunction) {
  parser.addGenerator("nullary", std::make_shared<NullaryGenerator>());

  auto fieldgen = parser.parseString("nullary()");
  EXPECT_EQ(fieldgen->str(), "nullary()");

  for (auto x : x_array) {
    for (auto y : y_array) {
      for (auto z : z_array) {
        for (auto t : t_array) {
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), 4.0);
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
          EXPECT_DOUBLE_EQ(actualFieldgen->generate(x, y, z, t),
                           clonedFieldgen->generate(x, y, z, t));
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
  EXPECT_THROW(fieldgen->generate(0., 0., 0., 0.), ParseException);
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
          EXPECT_DOUBLE_EQ(fieldgen->generate(x, y, z, t), 2 * (x + 3));
        }
      }
    }
  }
}

TEST(ParseExceptionTest, WhatTest) {
  try {
    throw ParseException("%s", "test message");
  } catch (ParseException &e) {
    std::string message{e.what()};
    EXPECT_NE(message.find("test message"), std::string::npos);
  }
}
