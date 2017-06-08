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

class FieldGenerator;
class ExpressionParser;
class ParseException;

#ifndef __EXPRESSION_PARSER_H__
#define __EXPRESSION_PARSER_H__

#include "unused.hxx"
#include "utils.hxx"

#include <string>
#include <map>
#include <list>
#include <utility>
#include <sstream>

#include <exception>

//////////////////////////////////////////////////////////

/*!
 * Represents an operation which generates a value at a given (x,y,z) location,
 * perhaps using other generators passed to clone()
 */
class FieldGenerator {
public:
  virtual ~FieldGenerator() { }

  /// Virtual constructor. Makes a copy of this FieldGenerator,
  /// initialised with the given list of arguments. It is up to the implementations
  /// to test whether the correct number of arguments is passed.
  ///
  /// @param[in] args   A (possibly empty) list of arguments to the generator function
  virtual FieldGenerator* clone(const std::list<FieldGenerator*> UNUSED(args)) {return NULL;}

  /// Generate a value at the given coordinates (x,y,z,t)
  /// This should be deterministic, always returning the same value given the same inputs
  virtual double generate(double x, double y, double z, double t) = 0;

  /// Create a string representation of the generator, for debugging output
  virtual const std::string str() {return std::string("?");}
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
  virtual ~ExpressionParser();

  /// Add a generator to the parser, which can then be recognised and used
  /// in expressions.
  ///
  /// @param[in] name  The name to be recognised in expressions. This should
  ///                  start with a letter and contain no whitespace, only
  ///                  alphanumeric letters and underscores.
  /// @param[in] g     The class inheriting from FieldGenerator. When recognised
  ///                  in an expression, the clone() function will be called
  ///                  to build a tree of generators
  void addGenerator(std::string name, FieldGenerator* g);

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
  void addBinaryOp(char sym, FieldGenerator* b, int precedence);
  
protected:
  /// This will be called to resolve any unknown symbols
  virtual FieldGenerator* resolve(std::string &UNUSED(name)) {return NULL;}

  /// Parses a given string into a tree of FieldGenerator objects
  FieldGenerator* parseString(const std::string &input);
  
private:
  
  std::map<std::string, FieldGenerator*> gen;  ///< Generators, addressed by name
  std::map<char, std::pair<FieldGenerator*, int> > bin_op; ///< Binary operations
  
  /// List of allocated generators
  std::list<FieldGenerator*> genheap;

  /// Lexing info, used when splitting input into tokens
  struct LexInfo {
    
    LexInfo(std::string input);
    
    signed char curtok;  ///< Current token. -1 for number, -2 for string, 0 for "end of input"
    double curval; ///< Value if a number
    std::string curident; ///< Identifier, variable or function name
    signed char LastChar;   ///< The last character read from the string
    std::stringstream ss; ///< Used to read values from the input string
    char nextToken(); ///< Get the next token in the string
    
    int getPos(); ///< Return position in the input
  };
  
  FieldGenerator* parseIdentifierExpr(LexInfo &lex);
  FieldGenerator* parseParenExpr(LexInfo &lex);
  FieldGenerator* parsePrimary(LexInfo &lex);
  FieldGenerator* parseBinOpRHS(LexInfo &lex, int prec, FieldGenerator* lhs);
  FieldGenerator* parseExpression(LexInfo &lex);
  
  /// Record generator in list, and return it
  FieldGenerator* record( FieldGenerator* g) {
    genheap.push_back(g);
    return g;
  }
};

//////////////////////////////////////////////////////

/// Binary operators
class FieldBinary : public FieldGenerator {
public:
  FieldBinary(FieldGenerator* l, FieldGenerator* r, char o) : lhs(l), rhs(r), op(o) {}
  FieldGenerator* clone(const std::list<FieldGenerator*> args);
  double generate(double x, double y, double z, double t);

  const std::string str() {return std::string("(")+lhs->str()+std::string(1,op)+rhs->str()+std::string(")");}
private:
  FieldGenerator *lhs, *rhs;
  char op;
};

/// Represent fixed values
class FieldValue : public FieldGenerator {
public:
  FieldValue(double val) : value(val) {}
  FieldGenerator* clone(const std::list<FieldGenerator*> UNUSED(args)) { return new FieldValue(value); }
  double generate(double UNUSED(x), double UNUSED(y), double UNUSED(z), double UNUSED(t)) { return value; }
  const std::string str() {
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
  template <typename... Args>
  ParseException(const std::string &format, Args... args)
      : message(string_format(format, args...)) {}
  virtual ~ParseException() throw() {}

  const char *what() const throw() { return message.c_str(); }

protected:
  std::string message;
};

#endif // __EXPRESSION_PARSER_H__
