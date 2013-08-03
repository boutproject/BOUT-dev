/**************************************************************************
 * Generate a field with specified values, mainly for creating
 * initial perturbations
 * 
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
class FieldFactory;

#ifndef __FIELD_FACTORY_H__
#define __FIELD_FACTORY_H__

#include "field2d.hxx"
#include "field3d.hxx"

#include "bout/mesh.hxx"

#include "bout/constants.hxx"

#include <string>
#include <map>
#include <sstream>
#include <list>
#include <utility>

#include "output.hxx"

using std::string;
using std::map;
using std::stringstream;
using std::list;
using std::pair;

//////////////////////////////////////////////////////////
// Generates a value at a given (x,y,z) location,
// perhaps using other generators passed to clone()

class FieldGenerator {
public:
  virtual ~FieldGenerator() { }
  virtual FieldGenerator* clone(const list<FieldGenerator*> args) {return NULL;}
  virtual BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) = 0;
};

//////////////////////////////////////////////////////////
// Create a tree of generators from an input string

class FieldFactory {
public:
  FieldFactory(Mesh *m);
  ~FieldFactory();
  
  const Field2D create2D(const string &value);
  const Field3D create3D(const string &value);
  
  void addGenerator(string name, FieldGenerator* g);
  void addBinaryOp(char sym, FieldGenerator* b, int precedence);
private:
  Mesh *fieldmesh;

  map<string, FieldGenerator*> gen;
  map<char, pair<FieldGenerator*, int> > bin_op; // Binary operations

  // Lexing info
  char curtok;  // Current token. -1 for number, -2 for string, 0 for "end of input"
  BoutReal curval; // Value if a number
  string curident; // Identifier
  char LastChar;
  stringstream ss;
  char nextToken();
  
  FieldGenerator* parseIdentifierExpr();
  FieldGenerator* parseParenExpr();
  FieldGenerator* parsePrimary();
  FieldGenerator* parseBinOpRHS(int prec, FieldGenerator* lhs);
  FieldGenerator* parseExpression();
  
  FieldGenerator* parse(const string &input);
};

#endif // __FIELD_FACTORY_H__
