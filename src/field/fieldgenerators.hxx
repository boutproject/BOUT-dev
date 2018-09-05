/*!
 * \file fieldgenerators.hxx
 * 
 * These classes are used by FieldFactory
 */

#ifndef __FIELDGENERATORS_H__
#define __FIELDGENERATORS_H__

#include <field_factory.hxx>
#include <boutexception.hxx>
#include <unused.hxx>

#include <cmath>

using std::list;

//////////////////////////////////////////////////////////
// Generators from values

/// Creates a Field Generator using a pointer to value
/// WARNING: The value pointed to must remain in scope until this generator is finished
class FieldValuePtr : public FieldGenerator {
public:
  FieldValuePtr(BoutReal *val) : ptr(val) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > UNUSED(args)) { return std::make_shared<FieldValuePtr>(ptr); }
  BoutReal generate(double UNUSED(x), double UNUSED(y), double UNUSED(z), double UNUSED(t)) { return *ptr; }
private:
  BoutReal *ptr; 
};

//////////////////////////////////////////////////////////
// Functions

/// Sine function field generator
class FieldSin : public FieldGenerator {
public:
  FieldSin(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
  const std::string str() {return std::string("sin(")+gen->str()+std::string(")");}
private:
  FieldGeneratorPtr gen;
};

/// Cosine function field generator
class FieldCos : public FieldGenerator {
public:
  FieldCos(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);

  const std::string str() {return std::string("cos(")+gen->str()+std::string(")");}
private:
  FieldGeneratorPtr gen;
};

/// Template class to define generators around a C function
typedef BoutReal(*single_arg_op)(BoutReal);
template<single_arg_op Op>
class FieldGenOneArg : public FieldGenerator { ///< Template for single-argument function
public:
  FieldGenOneArg(FieldGeneratorPtr g) : gen(g) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args) {
    if(args.size() != 1) {
      throw ParseException("Incorrect number of arguments to function. Expecting 1, got %d", args.size());
    }
    return std::make_shared<FieldGenOneArg<Op>>(args.front());
  }
  BoutReal generate(double x, double y, double z, double t) {
    return Op(gen->generate(x,y,z,t));
  }
  const std::string str() {return std::string("func(")+gen->str()+std::string(")");}
private:
  FieldGeneratorPtr gen;
};

/// Template for a FieldGenerator with two input arguments
typedef BoutReal(*double_arg_op)(BoutReal, BoutReal);
template<double_arg_op Op>
class FieldGenTwoArg : public FieldGenerator { ///< Template for two-argument function
public:
  FieldGenTwoArg(FieldGeneratorPtr a, FieldGeneratorPtr b) : A(a), B(b) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args) {
    if(args.size() != 2) {
      throw ParseException("Incorrect number of arguments to function. Expecting 2, got %d", args.size());
    }
    return std::make_shared<FieldGenTwoArg<Op>>(args.front(), args.back());
  }
  BoutReal generate(double x, double y, double z, double t) {
    return Op(A->generate(x,y,z,t), B->generate(x,y,z,t));
  }
  const std::string str() {return std::string("cos(")+A->str()+","+B->str()+std::string(")");}
private:
  FieldGeneratorPtr A, B;
};

/// Arc (Inverse) tangent. Either one or two argument versions
class FieldATan : public FieldGenerator { 
public:
  FieldATan(FieldGeneratorPtr a, FieldGeneratorPtr b=nullptr) : A(a), B(b) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args) {
    if(args.size() == 1) {
      return std::make_shared<FieldATan>(args.front());
    }else if(args.size() == 2) {
      return std::make_shared<FieldATan>(args.front(), args.back());
    }
    throw ParseException("Incorrect number of arguments to atan function. Expecting 1 or 2, got %d", args.size());
  }
  BoutReal generate(double x, double y, double z, double t) {
    if(B == nullptr)
      return atan(A->generate(x,y,z,t));
    return atan2(A->generate(x,y,z,t), B->generate(x,y,z,t));
  }
private:
  FieldGeneratorPtr A, B;
};

/// Hyperbolic sine function
class FieldSinh : public FieldGenerator {
public:
  FieldSinh(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr gen;
};

/// Hyperbolic cosine
class FieldCosh : public FieldGenerator {
public:
  FieldCosh(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr gen;
};

/// Hyperbolic tangent
class FieldTanh : public FieldGenerator {
public:
  FieldTanh(FieldGeneratorPtr g=nullptr) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr gen;
};

/// Gaussian distribution, taking mean and width arguments
class FieldGaussian : public FieldGenerator {
public:
  FieldGaussian(FieldGeneratorPtr xin, FieldGeneratorPtr sin) : X(xin), s(sin) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr X, s;
};

/// Absolute value
class FieldAbs : public FieldGenerator {
public:
  FieldAbs(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr gen;
};

/// Square root function
class FieldSqrt : public FieldGenerator {
public:
  FieldSqrt(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr gen;
};

/// Heaviside function, switches between 0 and 1
class FieldHeaviside : public FieldGenerator {
public:
  FieldHeaviside(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
  const std::string str() {return std::string("H(")+gen->str()+std::string(")");}
private:
  FieldGeneratorPtr gen;
};

/// Generator for the error function erf
class FieldErf : public FieldGenerator {
public:
  FieldErf(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  FieldGeneratorPtr gen;
};

/// Minimum
class FieldMin : public FieldGenerator {
public:
  FieldMin() {}
  FieldMin(const list<FieldGeneratorPtr > args) : input(args) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args) {
    if(args.size() == 0) {
      throw ParseException("min function must have some inputs");
    }
    return std::make_shared<FieldMin>(args);
  }
  BoutReal generate(double x, double y, double z, double t) {
    list<FieldGeneratorPtr >::iterator it=input.begin();
    BoutReal result = (*it)->generate(x,y,z,t);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(x,y,z,t);
      if(val < result)
        result = val;
    }
    return result;
  }
private:
  list<FieldGeneratorPtr > input;
};

/// Maximum
class FieldMax : public FieldGenerator {
public:
  FieldMax() {}
  FieldMax(const list<FieldGeneratorPtr > args) : input(args) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args) {
    if(args.size() == 0) {
      throw ParseException("max function must have some inputs");
    }
    return std::make_shared<FieldMax>(args);
  }
  BoutReal generate(double x, double y, double z, double t) {
    list<FieldGeneratorPtr >::iterator it=input.begin();
    BoutReal result = (*it)->generate(x,y,z,t);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(x,y,z,t);
      if(val > result)
        result = val;
    }
    return result;
  }
private:
  list<FieldGeneratorPtr > input;
};

/// Generator to round to the nearest integer
class FieldRound : public FieldGenerator {
public:
  FieldRound(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args) {
    if(args.size() != 1)
      throw BoutException("round function must have one input");
    return std::make_shared<FieldRound>(args.front());
  }
  BoutReal generate(double x, double y, double z, double t) {
    BoutReal val = gen->generate(x,y,z,t);
    if(val > 0.0) {
      return static_cast<int>(val + 0.5);
    }
    return static_cast<int>(val - 0.5);
  }
private:
  FieldGeneratorPtr gen;
};

//////////////////////////////////////////////////////////
// Ballooning transform
// Use a truncated Ballooning transform to enforce periodicity
// in doubly periodic domains

#include <bout/mesh.hxx>

class FieldBallooning : public FieldGenerator {
public:
  FieldBallooning(Mesh *m, FieldGeneratorPtr a = nullptr, int n = 3) : mesh(m), arg(a), ball_n(n) {}
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  Mesh *mesh;
  FieldGeneratorPtr arg;
  int ball_n;   // How many times around in each direction
};

//////////////////////////////////////////////////////////
// Mix of mode numbers
// This is similar to BOUT initialisation option 3

class FieldMixmode : public FieldGenerator {
public:
  FieldMixmode(FieldGeneratorPtr a = nullptr, BoutReal seed = 0.5);
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  /// Generate a random number between 0 and 1 (exclusive)
  /// given an arbitrary seed value
  ///
  /// This PRNG has no memory, i.e. you need to call it 
  /// with a different seed each time.
  BoutReal genRand(BoutReal seed);

  FieldGeneratorPtr arg;
  BoutReal phase[14];
};

//////////////////////////////////////////////////////////
// TanhHat
class FieldTanhHat : public FieldGenerator {
public:
  // Constructor
  FieldTanhHat(FieldGeneratorPtr xin,
               FieldGeneratorPtr widthin,
               FieldGeneratorPtr centerin,
               FieldGeneratorPtr steepnessin)
        : X(xin), width(widthin), center(centerin), steepness(steepnessin) {};
  // Clone containing the list of arguments
  FieldGeneratorPtr clone(const list<FieldGeneratorPtr > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  // The (x,y,z,t) field
  FieldGeneratorPtr X;
  // Other input arguments to the FieldTanhHat
  FieldGeneratorPtr width, center, steepness;
};

#endif // __FIELDGENERATORS_H__
