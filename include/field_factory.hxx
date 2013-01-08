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
// Basic generators: Numerical value, 'x', 'y' and 'z'

class FieldValue : public FieldGenerator {
public:
  FieldValue(BoutReal val) : value(val) {}
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) { return value; }
private:
  BoutReal value;
};

class FieldX : public FieldGenerator {
public:
  FieldGenerator* clone(const list<FieldGenerator*> args) { return new FieldX(); }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
};

class FieldY : public FieldGenerator {
public:
  FieldGenerator* clone(const list<FieldGenerator*> args) { return new FieldY(); }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
};

class FieldZ : public FieldGenerator {
public:
  FieldGenerator* clone(const list<FieldGenerator*> args) { return new FieldZ(); }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
};

//////////////////////////////////////////////////////////
// Values

class FieldPI : public FieldGenerator {
public:
  FieldGenerator* clone(const list<FieldGenerator*> args) { return new FieldPI(); }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) { return PI;}
};

//////////////////////////////////////////////////////////
// Functions

class FieldSin : public FieldGenerator {
public:
  FieldSin(FieldGenerator* g) : gen(g) {}
  ~FieldSin() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

class FieldCos : public FieldGenerator {
public:
  FieldCos(FieldGenerator* g) : gen(g) {}
  ~FieldCos() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

class FieldSinh : public FieldGenerator {
public:
  FieldSinh(FieldGenerator* g) : gen(g) {}
  ~FieldSinh() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

class FieldCosh : public FieldGenerator {
public:
  FieldCosh(FieldGenerator* g) : gen(g) {}
  ~FieldCosh() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

class FieldGaussian : public FieldGenerator {
public:
  FieldGaussian(FieldGenerator *xin, FieldGenerator *sin) : X(xin), s(sin) {}
  ~FieldGaussian() {if(X) delete X; if(s) delete s;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *X, *s;
};

class FieldAbs : public FieldGenerator {
public:
  FieldAbs(FieldGenerator* g) : gen(g) {}
  ~FieldAbs() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

class FieldSqrt : public FieldGenerator {
public:
  FieldSqrt(FieldGenerator* g) : gen(g) {}
  ~FieldSqrt() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

class FieldHeaviside : public FieldGenerator {
public:
  FieldHeaviside(FieldGenerator* g) : gen(g) {}
  ~FieldHeaviside() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *gen;
};

/// Unary minus
class FieldUnary : public FieldGenerator {
public:
  FieldUnary(FieldGenerator* g) : gen(g) {}
  ~FieldUnary() {if(gen) delete gen;}
  
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() != 1) {
      output << "FieldFactory error: Incorrect number of arguments to unary minus. Expecting 1, got " << args.size() << endl;
      return NULL;
    }
    return new FieldUnary(args.front());
  }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) {
    return -gen->generate(fieldmesh, x,y,z);
  }
private:
  FieldGenerator *gen;
};

/// Binary operators
class FieldBinary : public FieldGenerator {
public:
  FieldBinary(FieldGenerator* l, FieldGenerator* r, char o) : lhs(l), rhs(r), op(o) {}
  ~FieldBinary();
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z);
private:
  FieldGenerator *lhs, *rhs;
  char op;
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
