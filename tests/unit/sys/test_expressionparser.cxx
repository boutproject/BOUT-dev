#include "gtest/gtest.h"
#include "bout/sys/expressionparser.hxx"

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

TEST_F(ExpressionParserTest, MissingBracket) {
  EXPECT_THROW(parser.parseString("sin(x"), ParseException);
}

TEST_F(ExpressionParserTest, BadNumbers) {
  EXPECT_THROW(parser.parseString("1.1.4"), ParseException);
  EXPECT_THROW(parser.parseString("2.1e4.5."), ParseException);
  EXPECT_THROW(parser.parseString("3.e8e4"), ParseException);
  EXPECT_THROW(parser.parseString("4ee"), ParseException);
  EXPECT_THROW(parser.parseString("5G"), ParseException);
}
