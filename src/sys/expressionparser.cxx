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

#include <bout/sys/expressionparser.hxx>

#include <utils.hxx> // for lowercase

#include <stdarg.h>
#include <stdio.h>

using std::string;
using std::list;
using std::stringstream;

#include <stdlib.h>

/////////////////////////////////////////////
namespace { // These classes only visible in this file
  
  // Basic generators: Numerical value, 'x', 'y' and 'z'
  
  class FieldX : public FieldGenerator {
  public:
    FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
      return std::make_shared<FieldX>();
    }
    double generate(double x, double UNUSED(y), double UNUSED(z),
                    double UNUSED(t)) override {
      return x;
    }
    const std::string str() override { return std::string("x"); }
  };
  
  class FieldY : public FieldGenerator {
  public:
    FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
      return std::make_shared<FieldY>();
    }
    double generate(double UNUSED(x), double y, double UNUSED(z),
                    double UNUSED(t)) override {
      return y;
    }
    const std::string str() override { return std::string("y"); }
  };

  class FieldZ : public FieldGenerator {
  public:
    FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
      return std::make_shared<FieldZ>();
    }
    double generate(double UNUSED(x), double UNUSED(y), double z,
                    double UNUSED(t)) override {
      return z;
    }
    const std::string str() override { return std::string("z"); }
  };
  
  class FieldT : public FieldGenerator {
  public:
    FieldGeneratorPtr clone(const list<FieldGeneratorPtr> UNUSED(args)) override {
      return std::make_shared<FieldT>();
    }
    double generate(double UNUSED(x), double UNUSED(y), double UNUSED(z),
                    double t) override {
      return t;
    }
    const std::string str() override { return std::string("t"); }
  };
}

FieldGeneratorPtr FieldBinary::clone(const list<FieldGeneratorPtr> args) {
  if(args.size() != 2)
    throw ParseException("Binary operator expecting 2 arguments. Got '%d'", args.size());
  
  return std::make_shared<FieldBinary>(args.front(), args.back(), op);
}

BoutReal FieldBinary::generate(double x, double y, double z, double t) {
  BoutReal lval = lhs->generate(x,y,z,t);
  BoutReal rval = rhs->generate(x,y,z,t);
  switch(op) {
  case '+': return lval + rval;
  case '-': return lval - rval;
  case '*': return lval * rval;
  case '/': return lval / rval;
  case '^': return pow(lval, rval);
  }
  // Unknown operator.
  throw ParseException("Unknown binary operator '%c'", op);
}

/////////////////////////////////////////////

ExpressionParser::ExpressionParser() {
  // Add standard binary operations
  addBinaryOp('+', std::make_shared<FieldBinary>(nullptr, nullptr, '+'), 10);
  addBinaryOp('-', std::make_shared<FieldBinary>(nullptr, nullptr, '-'), 10);
  addBinaryOp('*', std::make_shared<FieldBinary>(nullptr, nullptr, '*'), 20);
  addBinaryOp('/', std::make_shared<FieldBinary>(nullptr, nullptr, '/'), 20);
  addBinaryOp('^', std::make_shared<FieldBinary>(nullptr, nullptr, '^'), 30);
  
  // Add standard generators
  addGenerator("x", std::make_shared<FieldX>());
  addGenerator("y", std::make_shared<FieldY>());
  addGenerator("z", std::make_shared<FieldZ>());
  addGenerator("t", std::make_shared<FieldT>());
}

void ExpressionParser::addGenerator(const string &name, FieldGeneratorPtr g) {
  gen[name] = g;
}

void ExpressionParser::addBinaryOp(char sym, FieldGeneratorPtr b, int precedence) {
  bin_op[sym] = std::make_pair(b, precedence);
}

FieldGeneratorPtr ExpressionParser::parseString(const string &input) {
  // Allocate a new lexer
  LexInfo lex(input);

  // Parse
  return parseExpression(lex);
}

//////////////////////////////////////////////////////////
// Private functions

FieldGeneratorPtr ExpressionParser::parseIdentifierExpr(LexInfo &lex) {
  string name = lowercase(lex.curident);
  lex.nextToken();
  
  if(lex.curtok == '(') {
    // Argument list. Find if a generator or function
    
    auto it = gen.find(name);
    if(it == gen.end())
      throw ParseException("Couldn't find generator '%s'", name.c_str());
    
    // Parse arguments (if any)
    list<FieldGeneratorPtr> args;
    
    lex.nextToken();
    if(lex.curtok == ')') {
      // Empty list
      lex.nextToken();
      return it->second->clone(args);
    }
    do{
      // Should be an expression
      args.push_back(parseExpression(lex));
      
      // Now either a comma or ')'
      
      if(lex.curtok == ')') {
        // Finished list
        lex.nextToken();
        return it->second->clone(args);
      }
      if(lex.curtok != ',') {
        throw ParseException("Expecting ',' or ')' in function argument list (%s)\n",
                             name.c_str());
      }
      lex.nextToken();
    }while(true);
    
  }else {
    // No arguments. Search in generator list
    auto it = gen.find(name);
    if(it == gen.end()) {
      // Not in internal map. Try to resolve
      FieldGeneratorPtr g = resolve(name);
      if(g == nullptr)
        throw ParseException("Couldn't find generator '%s'", name.c_str());
      return g;
    }
    list<FieldGeneratorPtr> args;
    return it->second->clone(args);
  }
}

FieldGeneratorPtr ExpressionParser::parseParenExpr(LexInfo &lex) {
  lex.nextToken(); // eat '('
  
  FieldGeneratorPtr g = parseExpression(lex);

  if ((lex.curtok != ')') && (lex.curtok != ']'))
    throw ParseException("Expecting ')' or ']' but got curtok=%d (%c)",
                         static_cast<int>(lex.curtok), lex.curtok);

  lex.nextToken(); // eat ')'
  return g;
}

FieldGeneratorPtr ExpressionParser::parsePrimary(LexInfo &lex) {
  switch(lex.curtok) {
  case -1: { // a number
    lex.nextToken(); // Eat number
    return std::make_shared<FieldValue>(lex.curval);
  }
  case -2: {
    return parseIdentifierExpr(lex);
  }
  case '-': {
    // Unary minus
    // Don't eat the minus, and return an implicit zero
    return std::make_shared<FieldValue>(0.0);
  }
  case '(':
  case '[':
    return parseParenExpr(lex);
  }
  throw ParseException("Unexpected token %d (%c)",
                       static_cast<int>(lex.curtok), lex.curtok);
}

FieldGeneratorPtr ExpressionParser::parseBinOpRHS(LexInfo &lex, int ExprPrec, FieldGeneratorPtr lhs) {
  
  while(true) {
    // Check for end of input
    if((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ','))
      return lhs;
    
    // Next token should be a binary operator
    auto it = bin_op.find(lex.curtok);
    
    if(it == bin_op.end())
      throw ParseException("Unexpected binary operator '%c'", lex.curtok);
    
    FieldGeneratorPtr op = it->second.first;
    int TokPrec = it->second.second;
  
    if (TokPrec < ExprPrec)
      return lhs;
    
    lex.nextToken(); // Eat binop
    
    FieldGeneratorPtr rhs = parsePrimary(lex);
    
    if((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ',')) {
      // Done
    
      list<FieldGeneratorPtr> args;
      args.push_front(lhs);
      args.push_back(rhs);
      return op->clone(args);
    }
    
    // Find next binop
    it = bin_op.find(lex.curtok);
    
    if(it == bin_op.end())
      throw ParseException("Unexpected character '%c'", lex.curtok);
    
    int NextPrec = it->second.second;
    if (TokPrec < NextPrec) {
      rhs = parseBinOpRHS(lex, TokPrec+1, rhs);
    }
    
    // Merge lhs and rhs into new lhs
    list<FieldGeneratorPtr> args;
    args.push_front(lhs);
    args.push_back(rhs);
    lhs = op->clone(args);
  }
}

FieldGeneratorPtr ExpressionParser::parseExpression(LexInfo &lex) {
  FieldGeneratorPtr lhs = parsePrimary(lex);
  return parseBinOpRHS(lex, 0, lhs);
}

//////////////////////////////////////////////////////////
// LexInfo

ExpressionParser::LexInfo::LexInfo(const std::string &input) {
  ss.clear();
  ss.str(input); // Set the input stream
  ss.seekg(0, std::ios_base::beg);
  
  LastChar = static_cast<signed char>(ss.get()); // First char from stream
  nextToken(); // Get first token
}

char ExpressionParser::LexInfo::nextToken() {
  while(isspace(LastChar))
    LastChar = static_cast<signed char>(ss.get());
  
  if(!ss.good()) {
    curtok = 0;
    return 0;
  }
  
  if (isalpha(LastChar)) { // identifier: [a-zA-Z][a-zA-Z0-9_:]*
    curident.clear();
    do {
      curident += LastChar;
      LastChar = static_cast<signed char>(ss.get());
    }while(isalnum(LastChar) || (LastChar == '_') || (LastChar == ':'));
    curtok = -2;
    return curtok;
  }
  
  // Handle numbers

  if (isdigit(LastChar) || (LastChar == '.')) {   // Number: [0-9.]+
    bool gotdecimal = false, gotexponent = false;
    std::string NumStr;
    
    while(true) {
      if(LastChar == '.') {
        if(gotdecimal || gotexponent) {
          throw ParseException("Unexpected '.' in number expression");
        }
        gotdecimal = true;
      }else if((LastChar == 'E') || (LastChar == 'e')) {
        if(gotexponent) {
          throw ParseException("ExpressionParser error: Unexpected extra 'e' in number expression");
        }
        gotexponent = true;
        // Next character should be a '+' or '-' or digit
        NumStr += 'e';
        LastChar = static_cast<signed char>(ss.get());
        if((LastChar != '+') && (LastChar != '-') && !isdigit(LastChar)) {
          throw ParseException("ExpressionParser error: Expecting '+', '-' or number after 'e'");
        }
      }else if(!isdigit(LastChar))
        break;
      
      NumStr += LastChar;
      LastChar = static_cast<signed char>(ss.get());
    }
    
    curval = std::stod(NumStr);
    curtok = -1;
    return curtok;
  }

  // LastChar is unsigned, explicitly cast
  curtok = LastChar;
  LastChar = static_cast<signed char>(ss.get());
  return curtok;
}

//////////////////////////////////////////////////////////
// ParseException


ParseException::ParseException(const char *s, ...) {
  if(s == nullptr)
    return;

  int buf_len=1024;
  char * buffer= new char[buf_len];
  bout_vsnprintf(buffer,buf_len, s);
  
  message.assign(buffer);
  delete[] buffer;
}

const char* ParseException::what() const noexcept {
  return message.c_str();
}

