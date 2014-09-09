/**************************************************************************
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

#include <string>
#include <map>
#include <list>
#include <utility>
#include <sstream>

#include <exception>

//////////////////////////////////////////////////////////
// Generates a value at a given (x,y,z) location,
// perhaps using other generators passed to clone()

class FieldGenerator {
public:
  virtual ~FieldGenerator() { }
  virtual FieldGenerator* clone(const std::list<FieldGenerator*> args) {return NULL;}
  virtual double generate(double x, double y, double z, double t) = 0;
  virtual const std::string str() {return std::string("?");}
};

class ExpressionParser {
public:
  ExpressionParser();
  virtual ~ExpressionParser();
  
  void addGenerator(std::string name, FieldGenerator* g);
  void addBinaryOp(char sym, FieldGenerator* b, int precedence);
  
protected:
  /// This will be called to resolve any unknown symbols
  virtual FieldGenerator* resolve(std::string &name) {return NULL;}

  /// Parses a given string into a tree of FieldGenerator objects
  FieldGenerator* parseString(const std::string &input);
  
private:
  
  std::map<std::string, FieldGenerator*> gen;  // Generators
  std::map<char, std::pair<FieldGenerator*, int> > bin_op; // Binary operations
  
  // List of allocated generators
  std::list<FieldGenerator*> genheap;
  
  struct LexInfo {
    // Lexing info
    
    LexInfo(std::string input);
    
    char curtok;  // Current token. -1 for number, -2 for string, 0 for "end of input"
    double curval; // Value if a number
    std::string curident; // Identifier
    char LastChar;
    std::stringstream ss;
    char nextToken();
    
    int getPos(); // Return position in the input
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
  FieldGenerator* clone(const std::list<FieldGenerator*> args) { return new FieldValue(value); }
  double generate(double x, double y, double z, double t) { return value; }
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
  ParseException(const char *, ...);
  virtual ~ParseException() throw() {}
  
  const char* what() const throw();
  
protected:
  std::string message;
};


#endif // __EXPRESSION_PARSER_H__
