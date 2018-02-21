/*!
 * \file fieldgenerators.hxx
 * 
 * These classes are used by FieldFactory
 */

#ifndef __FIELDGENERATORS_H__
#define __FIELDGENERATORS_H__

#include <bout/field_factory.hxx>
#include <bout/boutexception.hxx>
#include <bout/unused.hxx>

#include <cmath>

using std::list;

//////////////////////////////////////////////////////////
// Generators from values

/// Creates a Field Generator using a pointer to value
/// WARNING: The value pointed to must remain in scope until this generator is finished
class FieldValuePtr : public FieldGenerator {
public:
  FieldValuePtr(BoutReal *val) : ptr(val) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > UNUSED(args)) { return std::shared_ptr<FieldGenerator>( new  FieldValuePtr(ptr)); }
  BoutReal generate(double UNUSED(x), double UNUSED(y), double UNUSED(z), double UNUSED(t)) { return *ptr; }
private:
  BoutReal *ptr; 
};

//////////////////////////////////////////////////////////
// Functions

/// Sine function field generator
class FieldSin : public FieldGenerator {
public:
  FieldSin(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
  const std::string str() {return std::string("sin(")+gen->str()+std::string(")");}
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Cosine function field generator
class FieldCos : public FieldGenerator {
public:
  FieldCos(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);

  const std::string str() {return std::string("cos(")+gen->str()+std::string(")");}
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Template class to define generators around a C function
typedef BoutReal(*single_arg_op)(BoutReal);
template<single_arg_op Op>
class FieldGenOneArg : public FieldGenerator { ///< Template for single-argument function
public:
  FieldGenOneArg(std::shared_ptr<FieldGenerator> g) : gen(g) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args) {
    if(args.size() != 1) {
      throw ParseException("Incorrect number of arguments to function. Expecting 1, got %d", args.size());
    }
    return std::shared_ptr<FieldGenerator>( new  FieldGenOneArg<Op>(args.front()));
  }
  BoutReal generate(double x, double y, double z, double t) {
    return Op(gen->generate(x,y,z,t));
  }
  const std::string str() {return std::string("func(")+gen->str()+std::string(")");}
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Template for a FieldGenerator with two input arguments
typedef BoutReal(*double_arg_op)(BoutReal, BoutReal);
template<double_arg_op Op>
class FieldGenTwoArg : public FieldGenerator { ///< Template for two-argument function
public:
  FieldGenTwoArg(std::shared_ptr<FieldGenerator> a, std::shared_ptr<FieldGenerator> b) : A(a), B(b) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args) {
    if(args.size() != 2) {
      throw ParseException("Incorrect number of arguments to function. Expecting 2, got %d", args.size());
    }
    return std::shared_ptr<FieldGenerator>( new  FieldGenTwoArg<Op>(args.front(), args.back()));
  }
  BoutReal generate(double x, double y, double z, double t) {
    return Op(A->generate(x,y,z,t), B->generate(x,y,z,t));
  }
  const std::string str() {return std::string("cos(")+A->str()+","+B->str()+std::string(")");}
private:
  std::shared_ptr<FieldGenerator> A, B;
};

/// Arc (Inverse) tangent. Either one or two argument versions
class FieldATan : public FieldGenerator { 
public:
  FieldATan(std::shared_ptr<FieldGenerator> a, std::shared_ptr<FieldGenerator> b=nullptr) : A(a), B(b) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args) {
    if(args.size() == 1) {
      return std::shared_ptr<FieldGenerator>( new  FieldATan(args.front()));
    }else if(args.size() == 2) {
      return std::shared_ptr<FieldGenerator>( new  FieldATan(args.front(), args.back()));
    }
    throw ParseException("Incorrect number of arguments to atan function. Expecting 1 or 2, got %d", args.size());
  }
  BoutReal generate(double x, double y, double z, double t) {
    if(B == nullptr)
      return atan(A->generate(x,y,z,t));
    return atan2(A->generate(x,y,z,t), B->generate(x,y,z,t));
  }
private:
  std::shared_ptr<FieldGenerator> A, B;
};

/// Hyperbolic sine function
class FieldSinh : public FieldGenerator {
public:
  FieldSinh(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Hyperbolic cosine
class FieldCosh : public FieldGenerator {
public:
  FieldCosh(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Hyperbolic tangent
class FieldTanh : public FieldGenerator {
public:
  FieldTanh(std::shared_ptr<FieldGenerator> g=nullptr) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Gaussian distribution, taking mean and width arguments
class FieldGaussian : public FieldGenerator {
public:
  FieldGaussian(std::shared_ptr<FieldGenerator> xin, std::shared_ptr<FieldGenerator> sin) : X(xin), s(sin) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> X, s;
};

/// Absolute value
class FieldAbs : public FieldGenerator {
public:
  FieldAbs(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Square root function
class FieldSqrt : public FieldGenerator {
public:
  FieldSqrt(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Heaviside function, switches between 0 and 1
class FieldHeaviside : public FieldGenerator {
public:
  FieldHeaviside(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
  const std::string str() {return std::string("H(")+gen->str()+std::string(")");}
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Generator for the error function erf
class FieldErf : public FieldGenerator {
public:
  FieldErf(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  std::shared_ptr<FieldGenerator> gen;
};

/// Minimum
class FieldMin : public FieldGenerator {
public:
  FieldMin() {}
  FieldMin(const list<std::shared_ptr<FieldGenerator> > args) : input(args) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args) {
    if(args.size() == 0) {
      throw ParseException("min function must have some inputs");
    }
    return std::shared_ptr<FieldGenerator>( new  FieldMin(args));
  }
  BoutReal generate(double x, double y, double z, double t) {
    list<std::shared_ptr<FieldGenerator> >::iterator it=input.begin();
    BoutReal result = (*it)->generate(x,y,z,t);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(x,y,z,t);
      if(val < result)
        result = val;
    }
    return result;
  }
private:
  list<std::shared_ptr<FieldGenerator> > input;
};

/// Maximum
class FieldMax : public FieldGenerator {
public:
  FieldMax() {}
  FieldMax(const list<std::shared_ptr<FieldGenerator> > args) : input(args) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args) {
    if(args.size() == 0) {
      throw ParseException("max function must have some inputs");
    }
    return std::shared_ptr<FieldGenerator>( new  FieldMax(args));
  }
  BoutReal generate(double x, double y, double z, double t) {
    list<std::shared_ptr<FieldGenerator> >::iterator it=input.begin();
    BoutReal result = (*it)->generate(x,y,z,t);
    for(;it != input.end(); it++) {
      BoutReal val = (*it)->generate(x,y,z,t);
      if(val > result)
        result = val;
    }
    return result;
  }
private:
  list<std::shared_ptr<FieldGenerator> > input;
};

/// Generator to round to the nearest integer
class FieldRound : public FieldGenerator {
public:
  FieldRound(std::shared_ptr<FieldGenerator> g) : gen(g) {}

  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args) {
    if(args.size() != 1)
      throw BoutException("round function must have one input");
    return std::shared_ptr<FieldGenerator>( new  FieldRound(args.front()));
  }
  BoutReal generate(double x, double y, double z, double t) {
    BoutReal val = gen->generate(x,y,z,t);
    if(val > 0.0) {
      return static_cast<int>(val + 0.5);
    }
    return static_cast<int>(val - 0.5);
  }
private:
  std::shared_ptr<FieldGenerator> gen;
};

//////////////////////////////////////////////////////////
// Ballooning transform
// Use a truncated Ballooning transform to enforce periodicity
// in doubly periodic domains

#include <bout/mesh.hxx>

class FieldBallooning : public FieldGenerator {
public:
  FieldBallooning(Mesh *m, std::shared_ptr<FieldGenerator> a = nullptr, int n = 3) : mesh(m), arg(a), ball_n(n) {}
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  Mesh *mesh;
  std::shared_ptr<FieldGenerator> arg;
  int ball_n;   // How many times around in each direction
};

//////////////////////////////////////////////////////////
// Mix of mode numbers
// This is similar to BOUT initialisation option 3

class FieldMixmode : public FieldGenerator {
public:
  FieldMixmode(std::shared_ptr<FieldGenerator> a = nullptr, BoutReal seed = 0.5);
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  /// Generate a random number between 0 and 1 (exclusive)
  /// given an arbitrary seed value
  ///
  /// This PRNG has no memory, i.e. you need to call it 
  /// with a different seed each time.
  BoutReal genRand(BoutReal seed);

  std::shared_ptr<FieldGenerator> arg;
  BoutReal phase[14];
};

//////////////////////////////////////////////////////////
// TanhHat
class FieldTanhHat : public FieldGenerator {
public:
  // Constructor
  FieldTanhHat(std::shared_ptr<FieldGenerator> xin,
               std::shared_ptr<FieldGenerator> widthin,
               std::shared_ptr<FieldGenerator> centerin,
               std::shared_ptr<FieldGenerator> steepnessin)
        : X(xin), width(widthin), center(centerin), steepness(steepnessin) {};
  // Clone containing the list of arguments
  std::shared_ptr<FieldGenerator> clone(const list<std::shared_ptr<FieldGenerator> > args);
  BoutReal generate(double x, double y, double z, double t);
private:
  // The (x,y,z,t) field
  std::shared_ptr<FieldGenerator> X;
  // Other input arguments to the FieldTanhHat
  std::shared_ptr<FieldGenerator> width, center, steepness;
};

#endif // __FIELDGENERATORS_H__
