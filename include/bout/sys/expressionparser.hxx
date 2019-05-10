/*!************************************************************************
 * \file expressionparser.hxx
 *
 * Parses strings containing expressions, returning a tree of generators
 *
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __EXPRESSION_PARSER_H__
#define __EXPRESSION_PARSER_H__

#include "bout/format.hxx"
#include "unused.hxx"

#include <exception>
#include <list>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

#include "position.hxx"

class FieldGenerator;
using FieldGeneratorPtr = std::shared_ptr<FieldGenerator>;

//////////////////////////////////////////////////////////

/*!
 * Represents an operation which generates a value at a given (x,y,z) location,
 * perhaps using other generators passed to clone()
 */
class FieldGenerator {
public:
  virtual ~FieldGenerator() = default;

  /// Virtual constructor. Makes a copy of this FieldGenerator,
  /// initialised with the given list of arguments. It is up to the implementations
  /// to test whether the correct number of arguments is passed.
  ///
  /// @param[in] args   A (possibly empty) list of arguments to the generator function
  virtual FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> UNUSED(args)) {
    return nullptr;
  }

  /// Generate a value at the given coordinates (x,y,z,t)
  /// This should be deterministic, always returning the same value given the same inputs
  virtual double generate(Position pos) = 0;

  /// Create a string representation of the generator, for debugging output
  virtual std::string str() const { return std::string("?"); }
};

/*!
 * Parses expressions, turning strings into FieldGenerator
 * objects which can be used to perform calculations.
 *
 * This class is not intended to be used directly, but should be inherited
 * to add additional functionality
 */
class ExpressionParser {
public:
  ExpressionParser();
  virtual ~ExpressionParser() = default;

  /// Add a generator to the parser, which can then be recognised and used
  /// in expressions.
  ///
  /// @param[in] name  The name to be recognised in expressions. This should
  ///                  start with a letter and contain no whitespace, only
  ///                  alphanumeric letters and underscores.
  /// @param[in] g     The class inheriting from FieldGenerator. When recognised
  ///                  in an expression, the clone() function will be called
  ///                  to build a tree of generators
  void addGenerator(const std::string& name, FieldGeneratorPtr g);

  /// Add a binary operator such as +,-,*,/,^
  ///
  /// @param[in] sym  The operator symbol. This must be a single character
  /// @param[in] b    The FieldGenerator to use. When the symbol is recognised,
  ///                 b->clone() will be called with two input arguments
  /// @param[in] precedence  The precedence of the operator, which decides which order
  ///                        operators are performed in. Higher precedence operations
  ///                        are done before low precedence operations.
  /// Binary operators already defined are
  ///  +, -  precedence = 10
  ///  *, /  precedence = 20
  ///  ^     precedence = 30
  ///
  void addBinaryOp(char sym, FieldGeneratorPtr b, int precedence);

protected:
  /// This will be called to resolve any unknown symbols
  virtual FieldGeneratorPtr resolve(std::string& UNUSED(name)) const { return nullptr; }

  /// Parses a given string into a tree of FieldGenerator objects
  FieldGeneratorPtr parseString(const std::string& input) const;

  /// Characters which cannot be used in symbols; all other allowed
  /// In addition, whitespace cannot be used
  /// Adding a binary operator adds its symbol to this string
  std::string reserved_chars = "+-*/^[](){},";

private:
  std::map<std::string, FieldGeneratorPtr> gen; ///< Generators, addressed by name
  std::map<char, std::pair<FieldGeneratorPtr, int>> bin_op; ///< Binary operations

  /// Lexing info, used when splitting input into tokens
  struct LexInfo {

    LexInfo(const std::string& input, const std::string& reserved_chars = "");

    /// Current token. -1 for number, -2 for symbol, -3 for {string}, 0 for "end of input"
    signed char curtok = 0;
    double curval; ///< Value if a number
    std::string curident;       ///< Identifier, variable or function name
    signed char LastChar;       ///< The last character read from the string
    std::stringstream ss;       ///< Used to read values from the input string
    std::string reserved_chars; ///< Reserved characters, not in symbols
    char nextToken();           ///< Get the next token in the string
  };

  FieldGeneratorPtr parseIdentifierExpr(LexInfo& lex) const;
  FieldGeneratorPtr parseParenExpr(LexInfo& lex) const;

  /// Parse a primary expression, one of:
  ///   - number
  ///   - identifier
  ///   - ( ... )
  ///   - [ ... ]
  ///   - a unary '-', which is converted to '0 -'
  ///   A ParseException is thrown if none of these is found
  FieldGeneratorPtr parsePrimary(LexInfo& lex) const;
  FieldGeneratorPtr parseBinOpRHS(LexInfo& lex, int prec, FieldGeneratorPtr lhs) const;
  FieldGeneratorPtr parseExpression(LexInfo& lex) const;
};

//////////////////////////////////////////////////////

/// Binary operators
class FieldBinary : public FieldGenerator {
public:
  FieldBinary(FieldGeneratorPtr l, FieldGeneratorPtr r, char o)
      : lhs(std::move(l)), rhs(std::move(r)), op(o) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override;
  double generate(Position pos) override;

  std::string str() const override {
    return std::string("(") + lhs->str() + std::string(1, op) + rhs->str()
           + std::string(")");
  }

private:
  FieldGeneratorPtr lhs, rhs;
  char op;
};

/// Represent fixed values
class FieldValue : public FieldGenerator {
public:
  FieldValue(double val) : value(val) {}

  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldValue>(value);
  }

  double generate(Position UNUSED(pos)) override {
    return value;
  }
  std::string str() const override {
    std::stringstream ss;
    ss << value;
    return ss.str();
  }

private:
  double value;
};

//////////////////////////////////////////////////////

class ParseException : public std::exception {
public:
  ParseException(const char*, ...) BOUT_FORMAT_ARGS(2, 3);
  ~ParseException() override = default;

  const char* what() const noexcept override;

protected:
  std::string message;
};

#endif // __EXPRESSION_PARSER_H__
