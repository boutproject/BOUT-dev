/*!
 * \file fieldgenerators.hxx
 *
 * These classes are used by FieldFactory
 */

#ifndef __FIELDGENERATORS_H__
#define __FIELDGENERATORS_H__

#include <boutexception.hxx>
#include <field_factory.hxx>
#include <unused.hxx>

#include <cmath>

//////////////////////////////////////////////////////////
// Generators from values

/// Creates a Field Generator using a pointer to value
/// WARNING: The value pointed to must remain in scope until this generator is finished
class FieldValuePtr : public FieldGenerator {
public:
  FieldValuePtr(BoutReal* val) : ptr(val) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> UNUSED(args)) override {
    return std::make_shared<FieldValuePtr>(ptr);
  }
  BoutReal generate(const bout::generator::Context&) override {
    return *ptr;
  }

private:
  BoutReal* ptr;
};

//////////////////////////////////////////////////////////
// Functions

/// Template class to define generators around a C function
using single_arg_op = BoutReal (*)(BoutReal);
template<single_arg_op Op>
class FieldGenOneArg : public FieldGenerator { ///< Template for single-argument function
public:
  FieldGenOneArg(FieldGeneratorPtr g, const std::string& name = "function") : gen(g), name(name) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() != 1) {
      throw ParseException("Incorrect number of arguments to {:s}. Expecting 1, got {:d}",
                           name, args.size());
    }
    return std::make_shared<FieldGenOneArg<Op>>(args.front(), name);
  }
  BoutReal generate(const bout::generator::Context& pos) override {
    return Op(gen->generate(pos));
  }
  std::string str() const override {
    return name + std::string("(") + gen->str() + std::string(")");
  }

private:
  FieldGeneratorPtr gen;
  std::string name; ///< A string describing the function, to be printed in error messages
};

/// Template for a FieldGenerator with two input arguments
using double_arg_op = BoutReal (*)(BoutReal, BoutReal);
template <double_arg_op Op>
class FieldGenTwoArg : public FieldGenerator { ///< Template for two-argument function
public:
  FieldGenTwoArg(FieldGeneratorPtr a, FieldGeneratorPtr b, const std::string& name = "function") : A(a), B(b), name(name) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() != 2) {
      throw ParseException("Incorrect number of arguments to {:s}. Expecting 2, got {:d}",
                           name, args.size());
    }
    return std::make_shared<FieldGenTwoArg<Op>>(args.front(), args.back(), name);
  }
  BoutReal generate(const bout::generator::Context& pos) override {
    return Op(A->generate(pos), B->generate(pos));
  }
  std::string str() const override {
    return name + std::string("(") + A->str() + "," + B->str() + std::string(")");
  }

private:
  FieldGeneratorPtr A, B;
  std::string name; ///< The name of the function, to be printed in error messages
};

/// Arc (Inverse) tangent. Either one or two argument versions
class FieldATan : public FieldGenerator {
public:
  FieldATan(FieldGeneratorPtr a, FieldGeneratorPtr b = nullptr) : A(a), B(b) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() == 1) {
      return std::make_shared<FieldATan>(args.front());
    } else if (args.size() == 2) {
      return std::make_shared<FieldATan>(args.front(), args.back());
    }
    throw ParseException(
        "Incorrect number of arguments to atan function. Expecting 1 or 2, got {:d}",
        args.size());
  }
  BoutReal generate(const bout::generator::Context& pos) override {
    if (B == nullptr)
      return atan(A->generate(pos));
    return atan2(A->generate(pos), B->generate(pos));
  }

private:
  FieldGeneratorPtr A, B;
};

/// Gaussian distribution, taking mean and width arguments
class FieldGaussian : public FieldGenerator {
public:
  FieldGaussian(FieldGeneratorPtr xin, FieldGeneratorPtr sin) : X(xin), s(sin) {}

  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override;
  BoutReal generate(const bout::generator::Context& pos) override;

private:
  FieldGeneratorPtr X, s;
};

/// Heaviside function, switches between 0 and 1
class FieldHeaviside : public FieldGenerator {
public:
  FieldHeaviside(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override;
  BoutReal generate(const bout::generator::Context& pos) override;
  std::string str() const override {
    return std::string("H(") + gen->str() + std::string(")");
  }

private:
  FieldGeneratorPtr gen;
};

/// Minimum
class FieldMin : public FieldGenerator {
public:
  FieldMin() = default;
  FieldMin(const std::list<FieldGeneratorPtr> args) : input(args) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() == 0) {
      throw ParseException("min function must have some inputs");
    }
    return std::make_shared<FieldMin>(args);
  }
  BoutReal generate(const bout::generator::Context& pos) override {
    auto it = input.begin();
    BoutReal result = (*it)->generate(pos);
    for (; it != input.end(); it++) {
      BoutReal val = (*it)->generate(pos);
      if (val < result)
        result = val;
    }
    return result;
  }

private:
  std::list<FieldGeneratorPtr> input;
};

/// Maximum
class FieldMax : public FieldGenerator {
public:
  FieldMax() = default;
  FieldMax(const std::list<FieldGeneratorPtr> args) : input(args) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() == 0) {
      throw ParseException("max function must have some inputs");
    }
    return std::make_shared<FieldMax>(args);
  }
  BoutReal generate(const bout::generator::Context& pos) override {
    auto it = input.begin();
    BoutReal result = (*it)->generate(pos);
    for (; it != input.end(); it++) {
      BoutReal val = (*it)->generate(pos);
      if (val > result)
        result = val;
    }
    return result;
  }

private:
  std::list<FieldGeneratorPtr> input;
};

/// Generator to round to the nearest integer
class FieldRound : public FieldGenerator {
public:
  FieldRound(FieldGeneratorPtr g) : gen(g) {}

  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() != 1) {
      throw ParseException("round function must have one input");
    }
    return std::make_shared<FieldRound>(args.front());
  }
  BoutReal generate(const bout::generator::Context& pos) override {
    BoutReal val = gen->generate(pos);
    if (val > 0.0) {
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
  FieldBallooning(Mesh* m, FieldGeneratorPtr a = nullptr, int n = 3)
      : mesh(m), arg(a), ball_n(n) {}
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override;
  BoutReal generate(const bout::generator::Context& pos) override;

private:
  Mesh* mesh;
  FieldGeneratorPtr arg;
  int ball_n; // How many times around in each direction
};

//////////////////////////////////////////////////////////
// Mix of mode numbers
// This is similar to BOUT initialisation option 3

class FieldMixmode : public FieldGenerator {
public:
  FieldMixmode(FieldGeneratorPtr a = nullptr, BoutReal seed = 0.5);
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override;
  BoutReal generate(const bout::generator::Context& pos) override;

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
  FieldTanhHat(FieldGeneratorPtr xin, FieldGeneratorPtr widthin,
               FieldGeneratorPtr centerin, FieldGeneratorPtr steepnessin)
      : X(xin), width(widthin), center(centerin), steepness(steepnessin) {};
  // Clone containing the list of arguments
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override;
  BoutReal generate(const bout::generator::Context& pos) override;

private:
  // The (x,y,z,t) field
  FieldGeneratorPtr X;
  // Other input arguments to the FieldTanhHat
  FieldGeneratorPtr width, center, steepness;
};

class FieldWhere : public FieldGenerator {
public:
  FieldWhere(FieldGeneratorPtr test, FieldGeneratorPtr gt0, FieldGeneratorPtr lt0)
      : test(test), gt0(gt0), lt0(lt0){};

  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> args) override {
    if (args.size() != 3) {
      throw ParseException("where expects 3 arguments (test, gt0, lt0) but got %lu",
                           static_cast<unsigned long>(args.size()));
    }
    auto arg_it = std::begin(args);
    auto first = *arg_it++;
    auto second = *arg_it++;
    auto third = *arg_it;

    return std::make_shared<FieldWhere>(first, second, third);
  }
  double generate(const bout::generator::Context& pos) override {
    if (test->generate(pos) > 0.0) {
      return gt0->generate(pos);
    }
    return lt0->generate(pos);
  }

  std::string str() const override {
    return std::string("where(") + test->str() + std::string(",") + gt0->str()
           + std::string(",") + lt0->str() + std::string(")");
  }

private:
  FieldGeneratorPtr test, gt0, lt0;
};


#endif // __FIELDGENERATORS_H__
