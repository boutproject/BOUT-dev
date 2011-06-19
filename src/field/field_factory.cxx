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

#include <globals.hxx>
#include <field_factory.hxx>
#include <utils.hxx>

#include <stdlib.h>
#include <cmath>

BoutReal FieldX::generate(int x, int y, int z) {
  return mesh->GlobalX(x);
}

BoutReal FieldY::generate(int x, int y, int z) {
  return TWOPI*mesh->GlobalY(y);
}

BoutReal FieldZ::generate(int x, int y, int z) {
  return TWOPI*((BoutReal) z) / ((BoutReal) (mesh->ngz-1));
}

FieldBinary::~FieldBinary() {
  if(lhs)
    delete lhs;
  if(rhs)
    delete rhs;
}

FieldGenerator* FieldBinary::clone(const list<FieldGenerator*> args) {
  if(args.size() != 2)
    return NULL;
  
  return new FieldBinary(args.front(), args.back(), op);
}

BoutReal FieldBinary::generate(int x, int y, int z) {
  BoutReal lval = lhs->generate(x,y,z);
  BoutReal rval = rhs->generate(x,y,z);
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

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory() {
  
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
}

FieldFactory::~FieldFactory() {
  
}

const Field3D FieldFactory::create2D(const string &value) {
  Field2D result = 0.;

  FieldGenerator* gen = parse(value);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 2D field from '"
           << value
           << "'" << endl;
    return result;
  }
  
  for(int x=0;x<mesh->ngx;x++)
    for(int y=0;y<mesh->ngy;y++)
      result[x][y] = gen->generate(x,y,0);
  
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
  
  for(int x=0;x<mesh->ngx;x++)
    for(int y=0;y<mesh->ngy;y++)
      for(int z=0;z<mesh->ngz;z++)
        result[x][y][z] = gen->generate(x,y,z);
  
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
          output << "Unexpected '.' in number expression" << endl;
          return 0;
        }
        gotdecimal = true;
      }else if((LastChar == 'E') || (LastChar == 'e')) {
        if(gotexponent) {
          output << "Unexpected extra 'e' in number expression" << endl;
          return 0;
        }
        gotexponent = true;
        // Next character should be a '+' or '-' or digit
        NumStr += 'e';
        LastChar = ss.get();
        if((LastChar != '+') && (LastChar != '-') && !isdigit(LastChar)) {
          output << "Expecting '+', '-' or number after 'e'"  << endl;
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
      if(!a) 
        return NULL;
      args.push_back(a);
      
      // Now either a comma or ')'
      
      if(curtok == ')') {
        // Finished list
        nextToken();
        return it->second->clone(args);
      }
      if(curtok != ',') {
        // error
        return NULL;
      }
    }while(true);
    
  }else {
    // No arguments. Search in generator list
    map<string, FieldGenerator*>::iterator it = gen.find(name);
    if(it == gen.end())
      return NULL;
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
  case -1: // a number
    return new FieldValue(curval);
  case -2: {
    return parseIdentifierExpr();
  }
  case '(':
  case '[':
    return parseParenExpr();
  }
  return NULL;
}

FieldGenerator* FieldFactory::parseBinOpRHS(int ExprPrec, FieldGenerator* lhs) {
  // Check for end of input
  if(curtok == 0)
    return lhs;

  // Next token should be a binary operator
  map<char, pair<FieldGenerator*, int> >::iterator it = bin_op.find(curtok);
  
  if(it == bin_op.end())
    return NULL;
  
  FieldGenerator* op = it->second.first;
  int TokPrec = it->second.second;
  
  if (TokPrec < ExprPrec)
    return lhs;
  
  nextToken(); // Eat binop
  
  FieldGenerator* rhs = parsePrimary();
  if(!rhs)
    return NULL;
  
  if(curtok == 0) {
    // Done
    
    list<FieldGenerator*> args;
    args.push_front(lhs);
    args.push_back(rhs);
    return op->clone(args);
  }
    
  // Find next binop
  it = bin_op.find(curtok);
  
  if(it == bin_op.end())
    return NULL;
  
  int NextPrec = it->second.second;
  if (TokPrec < NextPrec) {
    rhs = parseBinOpRHS(TokPrec+1, rhs);
    if(!rhs)
      return 0;
  }
  
  // Merge lhs and rhs
  list<FieldGenerator*> args;
  args.push_front(lhs);
  args.push_back(rhs);
  return op->clone(args);
}

FieldGenerator* FieldFactory::parseExpression() {
  FieldGenerator* lhs = parsePrimary();
  if(!lhs)
    return NULL;
  return parseBinOpRHS(0, lhs);
}

FieldGenerator* FieldFactory::parse(const string &input) {
  
  ss.str(input); // Set the input stream
  LastChar = ss.get(); // First char from stream
  nextToken(); // Get first token
  
  return parseExpression();
}
