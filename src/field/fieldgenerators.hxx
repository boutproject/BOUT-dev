/*
 * These classes are used by FieldFactory
 */

#ifndef __FIELDGENERATORS_H__
#define __FIELDGENERATORS_H__

#include <field_factory.hxx>

#include <cmath>

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

// Template class to define generators around a C function
typedef BoutReal(*single_arg_op)(BoutReal);
template<single_arg_op Op>
class FieldGenOneArg : public FieldGenerator { ///< Template for single-argument function
public:
  FieldGenOneArg(FieldGenerator* g) : gen(g) {}
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() != 1) {
      output << "FieldFactory error: Incorrect number of arguments to function. Expecting 1, got " << args.size() << endl;
      return NULL;
    }
    return new FieldGenOneArg<Op>(args.front());
  }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) {
    return Op(gen->generate(fieldmesh, x,y,z));
  }
private:
  FieldGenerator *gen;
};

typedef BoutReal(*double_arg_op)(BoutReal, BoutReal);
template<double_arg_op Op>
class FieldGenTwoArg : public FieldGenerator { ///< Template for two-argument function
public:
  FieldGenTwoArg(FieldGenerator* a, FieldGenerator* b) : A(a), B(b) {}
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() != 2) {
      output << "FieldFactory error: Incorrect number of arguments to function. Expecting 2, got " << args.size() << endl;
      return NULL;
    }
    return new FieldGenTwoArg<Op>(args.front(), args.back());
  }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) {
    return Op(A->generate(fieldmesh, x,y,z), B->generate(fieldmesh, x,y,z));
  }
private:
  FieldGenerator *A, *B;
};

class FieldATan : public FieldGenerator { // Arc Tangent
public:
  FieldATan(FieldGenerator* a, FieldGenerator* b=NULL) : A(a), B(b) {}
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() == 1) {
      return new FieldATan(args.front());
    }else if(args.size() == 2) {
      return new FieldATan(args.front(), args.back());
    }
    output << "FieldFactory error: Incorrect number of arguments to atan function. Expecting 1 or 2, got " << args.size() << endl;
    return NULL;
  }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) {
    if(B == NULL)
      return atan(A->generate(fieldmesh, x,y,z));
    return atan2(A->generate(fieldmesh, x,y,z), B->generate(fieldmesh, x,y,z));
  }
private:
  FieldGenerator *A, *B;
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

/// Minimum
class FieldMin : public FieldGenerator {
public:
  FieldMin() {}
  FieldMin(const list<FieldGenerator*> args) : input(args) {}
  ~FieldMin() {
    for(list<FieldGenerator*>::iterator it=input.begin(); it != input.end(); it++)
      delete (*it);
  }
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() == 0) {
      output << "FieldFactory error: min function must have some inputs\n";
      return NULL;
    }
    return new FieldMin(args);
  }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) {
    list<FieldGenerator*>::iterator it=input.begin();
    BoutReal result = (*it)->generate(fieldmesh, x,y,z);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(fieldmesh, x,y,z);
      if(val < result)
        result = val;
    }
    return result;
  }
private:
  list<FieldGenerator*> input;
};

/// Maximum
class FieldMax : public FieldGenerator {
public:
  FieldMax() {}
  FieldMax(const list<FieldGenerator*> args) : input(args) {}
  ~FieldMax() {
    for(list<FieldGenerator*>::iterator it=input.begin(); it != input.end(); it++)
      delete (*it);
  }
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() == 0) {
      output << "FieldFactory error: max function must have some inputs\n";
      return NULL;
    }
    return new FieldMax(args);
  }
  BoutReal generate(const Mesh *fieldmesh, int x, int y, int z) {
    list<FieldGenerator*>::iterator it=input.begin();
    BoutReal result = (*it)->generate(fieldmesh, x,y,z);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(fieldmesh, x,y,z);
      if(val > result)
        result = val;
    }
    return result;
  }
private:
  list<FieldGenerator*> input;
};

#endif // __FIELDGENERATORS_H__
