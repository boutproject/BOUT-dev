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

using std::map;
using std::string;
using std::pair;
using std::list;
using std::stringstream;

#include <stdlib.h>

/////////////////////////////////////////////
namespace { // These classes only visible in this file
  
  // Basic generators: Numerical value, 'x', 'y' and 'z'
  
  class FieldX : public FieldGenerator {
  public:
    FieldGenerator* clone(const list<FieldGenerator*> UNUSED(args)) { return new FieldX(); }
    double generate(double x, double UNUSED(y), double UNUSED(z), double UNUSED(t)) {
      return x;
    }
    const std::string str() {return std::string("x");}
  };
  
  class FieldY : public FieldGenerator {
  public:
    FieldGenerator* clone(const list<FieldGenerator*> UNUSED(args)) { return new FieldY(); }
    double generate(double UNUSED(x), double y, double UNUSED(z), double UNUSED(t)) {
      return y;
    }
    const std::string str() {return std::string("y");}
  };

  class FieldZ : public FieldGenerator {
  public:
    FieldGenerator* clone(const list<FieldGenerator*> UNUSED(args)) { return new FieldZ(); }
    double generate(double UNUSED(x), double UNUSED(y), double z, double UNUSED(t)) {
      return z;
    }
    const std::string str() {return std::string("z");}
  };
  
  class FieldT : public FieldGenerator {
  public:
    FieldGenerator* clone(const list<FieldGenerator*> UNUSED(args)) { return new FieldT(); }
    double generate(double UNUSED(x), double UNUSED(y), double UNUSED(z), double t) {
      return t;
    }
    const std::string str() {return std::string("t");}
  };

  /// Unary minus
  class FieldUnary : public FieldGenerator {
  public:
    FieldUnary(FieldGenerator* g) : gen(g) {}
  
    FieldGenerator* clone(const list<FieldGenerator*> args) {
      if(args.size() != 1) {
        throw ParseException("Incorrect number of arguments to unary minus. Expecting 1, got %d", args.size());
      }
      return new FieldUnary(args.front());
    }
    double generate(double x, double y, double z, double t) {
      return -gen->generate(x,y,z,t);
    }
    const std::string str() {return std::string("(-")+gen->str()+std::string(")");}
  private:
    FieldGenerator *gen;
  };
}

FieldGenerator* FieldBinary::clone(const list<FieldGenerator*> args) {
  if(args.size() != 2)
    throw ParseException("Binary operator expecting 2 arguments. Got '%d'", args.size());
  
  return new FieldBinary(args.front(), args.back(), op);
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
  // Unknown operator. Throw an error?
  return 0.;
}

/////////////////////////////////////////////

ExpressionParser::ExpressionParser() {
  // Add standard binary operations
  addBinaryOp('+', new FieldBinary(NULL, NULL, '+'), 10);
  addBinaryOp('-', new FieldBinary(NULL, NULL, '-'), 10);
  addBinaryOp('*', new FieldBinary(NULL, NULL, '*'), 20);
  addBinaryOp('/', new FieldBinary(NULL, NULL, '/'), 20);
  addBinaryOp('^', new FieldBinary(NULL, NULL, '^'), 30);
  
  // Add standard generators
  addGenerator("x", new FieldX());
  addGenerator("y", new FieldY());
  addGenerator("z", new FieldZ());
  addGenerator("t", new FieldT());
}

ExpressionParser::~ExpressionParser() {
  // Free memory
  for(const auto& it : gen)
    delete it.second;
  
  for(const auto& it : bin_op)
    delete it.second.first;
  
  // Delete allocated generators
  for(const auto& it : genheap)
    delete it;
}

void ExpressionParser::addGenerator(string name, FieldGenerator* g) {
  gen[name] = g;
}

void ExpressionParser::addBinaryOp(char sym, FieldGenerator* b, int precedence) {
  bin_op[sym] = pair<FieldGenerator*, int>(b, precedence);
}

FieldGenerator* ExpressionParser::parseString(const string &input) {
  // Allocate a new lexer
  LexInfo lex(input);
  
  // Parse
  FieldGenerator *expr = parseExpression(lex);
  
  return expr;
}

//////////////////////////////////////////////////////////
// Private functions

FieldGenerator* ExpressionParser::parseIdentifierExpr(LexInfo &lex) {
  string name = lowercase(lex.curident);
  lex.nextToken();
  
  if(lex.curtok == '(') {
    // Argument list. Find if a generator or function
    
    map<string, FieldGenerator*>::iterator it = gen.find(name);
    if(it == gen.end())
      throw ParseException("Couldn't find generator '%s'", name.c_str());
    
    // Parse arguments (if any)
    list<FieldGenerator*> args;
    
    lex.nextToken();
    if(lex.curtok == ')') {
      // Empty list
      lex.nextToken();
      return record( it->second->clone(args) );
    }
    do{
      // Should be an expression
      FieldGenerator *a = parseExpression(lex);
      if(!a) {
	throw ParseException("Couldn't parse argument %d to function '%s'", 
                             args.size()+1, name.c_str());
      }
      args.push_back(a);
      
      // Now either a comma or ')'
      
      if(lex.curtok == ')') {
        // Finished list
        lex.nextToken();
        return record( it->second->clone(args) );
      }
      if(lex.curtok != ',') {
        throw ParseException("Expecting ',' or ')' in function argument list (%s)\n",
                             name.c_str());
      }
      lex.nextToken();
    }while(true);
    
  }else {
    // No arguments. Search in generator list
    map<string, FieldGenerator*>::iterator it = gen.find(name);
    if(it == gen.end()) {
      // Not in internal map. Try to resolve
      FieldGenerator* g = resolve(name);
      if(g == NULL)
        throw ParseException("Couldn't find generator '%s'", name.c_str());
      return g;
    }
    list<FieldGenerator*> args;
    return record( it->second->clone(args) );
  }
}

FieldGenerator* ExpressionParser::parseParenExpr(LexInfo &lex) {
  lex.nextToken(); // eat '('
  
  FieldGenerator* g = parseExpression(lex);
  if(!g)
    return NULL;
  
  if((lex.curtok != ')') && (lex.curtok != ']'))
    return NULL;
  lex.nextToken(); // eat ')'
  return g;
}

FieldGenerator* ExpressionParser::parsePrimary(LexInfo &lex) {
  switch(lex.curtok) {
  case -1: { // a number
    lex.nextToken(); // Eat number
    return record( new FieldValue(lex.curval) );
  }
  case -2: {
    return parseIdentifierExpr(lex);
  }
  case '-': {
    // Unary minus
    //lex.nextToken(); // Eat '-'
    //return record( new FieldUnary(parsePrimary(lex)) );
    // Don't eat the minus, and return an implicit zero
    return record( new FieldValue(0.0) );
  }
  case '(':
  case '[':
    return parseParenExpr(lex);
  }
  return NULL;
}

FieldGenerator* ExpressionParser::parseBinOpRHS(LexInfo &lex, int ExprPrec, FieldGenerator* lhs) {
  
  while(true) {
    // Check for end of input
    if((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ','))
      return lhs;
    
    // Next token should be a binary operator
    map<char, pair<FieldGenerator*, int> >::iterator it = bin_op.find(lex.curtok);
    
    if(it == bin_op.end())
      throw ParseException("Unexpected binary operator '%c'", lex.curtok);
    
    FieldGenerator* op = it->second.first;
    int TokPrec = it->second.second;
  
    if (TokPrec < ExprPrec)
      return lhs;
    
    lex.nextToken(); // Eat binop
    
    FieldGenerator* rhs = parsePrimary(lex);
    if(!rhs)
      return NULL;
    
    if((lex.curtok == 0) || (lex.curtok == ')') || (lex.curtok == ',')) {
      // Done
    
      list<FieldGenerator*> args;
      args.push_front(lhs);
      args.push_back(rhs);
      return record( op->clone(args) );
    }
    
    // Find next binop
    it = bin_op.find(lex.curtok);
    
    if(it == bin_op.end())
      throw ParseException("Unexpected character '%c'", lex.curtok);
    
    int NextPrec = it->second.second;
    if (TokPrec < NextPrec) {
      rhs = parseBinOpRHS(lex, TokPrec+1, rhs);
      if(!rhs)
	return 0;
    }
    
    // Merge lhs and rhs into new lhs
    list<FieldGenerator*> args;
    args.push_front(lhs);
    args.push_back(rhs);
    lhs = record( op->clone(args) );
  }
}

FieldGenerator* ExpressionParser::parseExpression(LexInfo &lex) {
  FieldGenerator* lhs = parsePrimary(lex);
  if(!lhs)
    return NULL;
  return parseBinOpRHS(lex, 0, lhs);
}

//////////////////////////////////////////////////////////
// LexInfo

ExpressionParser::LexInfo::LexInfo(string input) {
  ss.clear();
  ss.str(input); // Set the input stream
  ss.seekg(0, ios_base::beg);
  
  LastChar = ss.get(); // First char from stream
  nextToken(); // Get first token
}

char ExpressionParser::LexInfo::nextToken() {
  while(isspace(LastChar))
    LastChar = ss.get();
  
  if(!ss.good()) {
    curtok = 0;
    return 0;
  }
  
  if (isalpha(LastChar)) { // identifier: [a-zA-Z][a-zA-Z0-9_:]*
    curident.clear();
    do {
      curident += LastChar;
      LastChar = ss.get();
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
        LastChar = ss.get();
        if((LastChar != '+') && (LastChar != '-') && !isdigit(LastChar)) {
          throw ParseException("ExpressionParser error: Expecting '+', '-' or number after 'e'");
        }
      }else if(!isdigit(LastChar))
        break;
      
      NumStr += LastChar;
      LastChar = ss.get();
    }
    
    curval = strtod(NumStr.c_str(), 0);
    curtok = -1;
    return curtok;
  }
  
  curtok = LastChar;
  LastChar = ss.get();
  return curtok;
}

int ExpressionParser::LexInfo::getPos() {
  return (int) ss.tellg();
}

//////////////////////////////////////////////////////////
// ParseException


ParseException::ParseException(const char *s, ...) {
  if(s == (const char*) NULL)
    return;

  int buf_len=1024;
  char * buffer= new char[buf_len];
  bout_vsnprintf(buffer,buf_len, s);
  
  message.assign(buffer);
  delete[] buffer;
}

const char* ParseException::what() const throw() {
  return message.c_str();
}

