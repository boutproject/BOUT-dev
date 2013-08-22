/**************************************************************************
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

#include <field_factory.hxx>
#include <utils.hxx>

#include <stdlib.h>
#include <cmath>

#include <output.hxx>
#include <bout/constants.hxx>

#include "fieldgenerators.hxx"

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory(Mesh *m) : fieldmesh(m) {
  
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
  
  // Useful values
  addGenerator("pi", new FieldValue(PI));

  // Some standard functions
  addGenerator("sin", new FieldSin(NULL));
  addGenerator("cos", new FieldCos(NULL));
  addGenerator("tan", new FieldGenOneArg<tan>(NULL));
  
  addGenerator("acos", new FieldGenOneArg<acos>(NULL));
  addGenerator("asin", new FieldGenOneArg<asin>(NULL));
  addGenerator("atan", new FieldATan(NULL));
  
  addGenerator("sinh", new FieldSinh(NULL));
  addGenerator("cosh", new FieldCosh(NULL));
  addGenerator("tanh", new FieldTanh());

  addGenerator("exp", new FieldGenOneArg<exp>(NULL));
  addGenerator("log", new FieldGenOneArg<log>(NULL));
  addGenerator("gauss", new FieldGaussian(NULL, NULL));
  addGenerator("abs", new FieldAbs(NULL));
  addGenerator("sqrt", new FieldSqrt(NULL));
  addGenerator("h", new FieldHeaviside(NULL));

  addGenerator("min", new FieldMin());
  addGenerator("max", new FieldMax());
}

FieldFactory::~FieldFactory() {
  // Free memory
  for(map<string, FieldGenerator*>::iterator it = gen.begin(); it != gen.end(); it++)
    delete it->second;
  
  for(map<char, pair<FieldGenerator*, int> >::iterator it = bin_op.begin(); it != bin_op.end(); it++)
    delete it->second.first;
}

const Field2D FieldFactory::create2D(const string &value) {
  Field2D result = 0.;

  FieldGenerator* gen = parse(value);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 2D field from '"
           << value
           << "'" << endl;
    return result;
  }
  
  for(int x=0;x<fieldmesh->ngx;x++)
    for(int y=0;y<fieldmesh->ngy;y++)
      result[x][y] = gen->generate(fieldmesh, x,y,0);
  
  delete gen;

  return result;
}

const Field3D FieldFactory::create3D(const string &value) {
  Field3D result = 0.;

  FieldGenerator* gen = parse(value);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 3D field from '"
           << value
           << "'" << endl;
    return result;
  }
  
  for(int x=0;x<fieldmesh->ngx;x++)
    for(int y=0;y<fieldmesh->ngy;y++)
      for(int z=0;z<fieldmesh->ngz;z++)
        result[x][y][z] = gen->generate(fieldmesh, x,y,z);
  
  delete gen;

  return result;
}

void FieldFactory::addGenerator(string name, FieldGenerator* g) {
  gen[name] = g;
}

void FieldFactory::addBinaryOp(char sym, FieldGenerator* b, int precedence) {
  bin_op[sym] = pair<FieldGenerator*, int>(b, precedence);
}

//////////////////////////////////////////////////////////
// FieldFactory private functions

char FieldFactory::nextToken() {
  while(isspace(LastChar))
    LastChar = ss.get();
  
  if(!ss.good()) {
    curtok = 0;
    return 0;
  }
  
  if (isalpha(LastChar)) { // identifier: [a-zA-Z][a-zA-Z0-9_]*
    curident.clear();
    do {
      curident += LastChar;
      LastChar = ss.get();
    }while(isalnum(LastChar) || (LastChar == '_'));
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
          output << "FieldFactory error: Unexpected '.' in number expression" << endl;
	  curtok = 0;
          return 0;
        }
        gotdecimal = true;
      }else if((LastChar == 'E') || (LastChar == 'e')) {
        if(gotexponent) {
          output << "FieldFactory error: Unexpected extra 'e' in number expression" << endl;
	  curtok = 0;
          return 0;
        }
        gotexponent = true;
        // Next character should be a '+' or '-' or digit
        NumStr += 'e';
        LastChar = ss.get();
        if((LastChar != '+') && (LastChar != '-') && !isdigit(LastChar)) {
          output << "FieldFactory error: Expecting '+', '-' or number after 'e'"  << endl;
	  curtok = 0;
          return 0;
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

//////////////////////////////////////////////////////////

FieldGenerator* FieldFactory::parseIdentifierExpr() {
  string name = lowercase(curident);
  nextToken();
  
  if(curtok == '(') {
    // Argument list. Find if a generator or function
    
    map<string, FieldGenerator*>::iterator it = gen.find(name);
    if(it == gen.end()) {
      output << "FieldFactory error: Couldn't find generator '"
             << name << "'" << endl;
      return NULL;
    }
    
    // Parse arguments (if any)
    list<FieldGenerator*> args;
    
    nextToken();
    if(curtok == ')') {
      // Empty list
      nextToken();
      return it->second->clone(args);
    }
    do{
      // Should be an expression
      FieldGenerator *a = parseExpression();
      if(!a) {
	output << "FieldFactory error: Couldn't parse argument " << args.size()+1 
	       << " to " << name << " function" << endl;
        return NULL;
      }
      args.push_back(a);
      
      // Now either a comma or ')'
      
      if(curtok == ')') {
        // Finished list
        nextToken();
        return it->second->clone(args);
      }
      if(curtok != ',') {
        output << "FieldFactory error: Expecting ',' or ')' in function argument list (" << name << ")" << endl;
        return NULL;
      }
      nextToken();
    }while(true);
    
  }else {
    // No arguments. Search in generator list
    map<string, FieldGenerator*>::iterator it = gen.find(name);
    if(it == gen.end()) {
      output << "FieldFactory error: Can't find generator '" << name << "'" << endl;
      return NULL;
    }
    list<FieldGenerator*> args;
    return it->second->clone(args);
  }
}

FieldGenerator* FieldFactory::parseParenExpr() {
  nextToken(); // eat '('
  
  FieldGenerator* g = parseExpression();
  if(!g)
    return NULL;
  
  if((curtok != ')') && (curtok != ']'))
    return NULL;
  nextToken(); // eat ')'
  return g;
}

FieldGenerator* FieldFactory::parsePrimary() {
  switch(curtok) {
  case -1: { // a number
    nextToken(); // Eat number
    return new FieldValue(curval);
  }
  case -2: {
    return parseIdentifierExpr();
  }
  case '-': {
    // Unary minus
    nextToken(); // Eat '-'
    return new FieldUnary(parsePrimary());
  }
  case '(':
  case '[':
    return parseParenExpr();
  }
  return NULL;
}

FieldGenerator* FieldFactory::parseBinOpRHS(int ExprPrec, FieldGenerator* lhs) {
  // Check for end of input
  if((curtok == 0) || (curtok == ')') || (curtok == ','))
    return lhs;

  // Next token should be a binary operator
  map<char, pair<FieldGenerator*, int> >::iterator it = bin_op.find(curtok);
  
  if(it == bin_op.end()) {
    output << "FieldFactory error: Unexpected binary operator '" << curtok << "'" << endl;
    return NULL;
  }
  
  FieldGenerator* op = it->second.first;
  int TokPrec = it->second.second;
  
  if (TokPrec < ExprPrec)
    return lhs;
  
  nextToken(); // Eat binop
  
  FieldGenerator* rhs = parsePrimary();
  if(!rhs)
    return NULL;

  if((curtok == 0) || (curtok == ')') || (curtok == ',')) {
    // Done
    
    list<FieldGenerator*> args;
    args.push_front(lhs);
    args.push_back(rhs);
    return op->clone(args);
  }
    
  // Find next binop
  it = bin_op.find(curtok);
  
  if(it == bin_op.end()) {
    output << "FieldFactory error: Unexpected character '" << curtok << "'" << endl;
    return NULL;
  }
  
  int NextPrec = it->second.second;
  if (TokPrec < NextPrec) {
    rhs = parseBinOpRHS(TokPrec+1, rhs);
    if(!rhs)
      return 0;
  }
  
  // Merge lhs and rhs into new lhs
  list<FieldGenerator*> args;
  args.push_front(lhs);
  args.push_back(rhs);
  lhs = op->clone(args);
  
  return parseBinOpRHS(0, lhs);
}

FieldGenerator* FieldFactory::parseExpression() {
  FieldGenerator* lhs = parsePrimary();
  if(!lhs)
    return NULL;
  return parseBinOpRHS(0, lhs);
}

FieldGenerator* FieldFactory::parse(const string &input) {
  
  ss.clear();
  ss.str(input); // Set the input stream
  ss.seekg(0, ios_base::beg);
  
  LastChar = ss.get(); // First char from stream
  nextToken(); // Get first token
  
  return parseExpression();
}
