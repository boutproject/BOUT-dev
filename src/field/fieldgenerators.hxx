/*
 * These classes are used by FieldFactory
 */

#ifndef __FIELDGENERATORS_H__
#define __FIELDGENERATORS_H__

#include <field_factory.hxx>

//////////////////////////////////////////////////////////
// Basic generators: Numerical value, 'x', 'y' and 'z'

class FieldValue : public FieldGenerator {
public:
  FieldValue(BoutReal val) : value(val) {}
  FieldGenerator* clone(const list<FieldGenerator*> args) { return new FieldValue(value); }
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

class FieldTanh : public FieldGenerator {
public:
  FieldTanh(FieldGenerator* g=NULL) : gen(g) {}
  ~FieldTanh() {if(gen) delete gen;}
  
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

#endif // __FIELDGENERATORS_H__
