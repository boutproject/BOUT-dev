/*
 * These classes are used by FieldFactory
 */

#ifndef __FIELDGENERATORS_H__
#define __FIELDGENERATORS_H__

#include <field_factory.hxx>

#include <cmath>

using std::list;

//////////////////////////////////////////////////////////
// Generators from values

// Creates a Field Generator using a pointer to value
// WARNING: The value pointed to must remain in scope until this generator is finished
class FieldValuePtr : public FieldGenerator {
public:
  FieldValuePtr(BoutReal *val) : ptr(val) {}
  FieldGenerator* clone(const list<FieldGenerator*> args) { return new FieldValuePtr(ptr); }
  BoutReal generate(double x, double y, double z, double t) { return *ptr; }
private:
  BoutReal *ptr;
};

//////////////////////////////////////////////////////////
// Functions

class FieldSin : public FieldGenerator {
public:
  FieldSin(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
  const std::string str() {return std::string("sin(")+gen->str()+std::string(")");}
private:
  FieldGenerator *gen;
};

class FieldCos : public FieldGenerator {
public:
  FieldCos(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);

  const std::string str() {return std::string("cos(")+gen->str()+std::string(")");}
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
      throw ParseException("Incorrect number of arguments to function. Expecting 1, got %d", args.size());
    }
    return new FieldGenOneArg<Op>(args.front());
  }
  BoutReal generate(double x, double y, double z, double t) {
    return Op(gen->generate(x,y,z,t));
  }
  const std::string str() {return std::string("func(")+gen->str()+std::string(")");}
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
      throw ParseException("Incorrect number of arguments to function. Expecting 2, got %d", args.size());
    }
    return new FieldGenTwoArg<Op>(args.front(), args.back());
  }
  BoutReal generate(double x, double y, double z, double t) {
    return Op(A->generate(x,y,z,t), B->generate(x,y,z,t));
  }
  const std::string str() {return std::string("cos(")+A->str()+","+B->str()+std::string(")");}
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
    throw ParseException("Incorrect number of arguments to atan function. Expecting 1 or 2, got %d", args.size());
  }
  BoutReal generate(double x, double y, double z, double t) {
    if(B == NULL)
      return atan(A->generate(x,y,z,t));
    return atan2(A->generate(x,y,z,t), B->generate(x,y,z,t));
  }
private:
  FieldGenerator *A, *B;
};

class FieldSinh : public FieldGenerator {
public:
  FieldSinh(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGenerator *gen;
};

class FieldCosh : public FieldGenerator {
public:
  FieldCosh(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGenerator *gen;
};

class FieldTanh : public FieldGenerator {
public:
  FieldTanh(FieldGenerator* g=NULL) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGenerator *gen;
};

class FieldGaussian : public FieldGenerator {
public:
  FieldGaussian(FieldGenerator *xin, FieldGenerator *sin) : X(xin), s(sin) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGenerator *X, *s;
};

class FieldAbs : public FieldGenerator {
public:
  FieldAbs(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGenerator *gen;
};

class FieldSqrt : public FieldGenerator {
public:
  FieldSqrt(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGenerator *gen;
};

class FieldHeaviside : public FieldGenerator {
public:
  FieldHeaviside(FieldGenerator* g) : gen(g) {}

  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
  const std::string str() {return std::string("H(")+gen->str()+std::string(")");}
private:
  FieldGenerator *gen;
};

/// Minimum
class FieldMin : public FieldGenerator {
public:
  FieldMin() {}
  FieldMin(const list<FieldGenerator*> args) : input(args) {}
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() == 0) {
      throw ParseException("min function must have some inputs");
    }
    return new FieldMin(args);
  }
  BoutReal generate(double x, double y, double z, double t) {
    list<FieldGenerator*>::iterator it=input.begin();
    BoutReal result = (*it)->generate(x,y,z,t);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(x,y,z,t);
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
  FieldGenerator* clone(const list<FieldGenerator*> args) {
    if(args.size() == 0) {
      throw ParseException("max function must have some inputs");
    }
    return new FieldMax(args);
  }
  BoutReal generate(double x, double y, double z, double t) {
    list<FieldGenerator*>::iterator it=input.begin();
    BoutReal result = (*it)->generate(x,y,z,t);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(x,y,z,t);
      if(val > result)
        result = val;
    }
    return result;
  }
private:
  list<FieldGenerator*> input;
};


//////////////////////////////////////////////////////////
// Ballooning transform
// Use a truncated Ballooning transform to enforce periodicity
// in doubly periodic domains

#include <bout/mesh.hxx>

class FieldBallooning : public FieldGenerator {
public:
  FieldBallooning(Mesh *m, FieldGenerator* a = NULL, int n = 3) : mesh(m), arg(a), ball_n(n) {}
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  Mesh *mesh;
  FieldGenerator *arg;
  int ball_n;   // How many times around in each direction
};

//////////////////////////////////////////////////////////
// Mix of mode numbers
// This is similar to BOUT initialisation option 3

class FieldMixmode : public FieldGenerator {
public:
  FieldMixmode(FieldGenerator* a = NULL, BoutReal seed = 0.5);
  FieldGenerator* clone(const list<FieldGenerator*> args);
  BoutReal generate(double x, double y, double z, double t);
private:
  /// Generate a random number between 0 and 1
  /// given an arbitrary seed value
  BoutReal genRand(BoutReal seed);

  FieldGenerator *arg;
  BoutReal phase[14];
};

#endif // __FIELDGENERATORS_H__
