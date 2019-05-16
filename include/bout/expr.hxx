/**************************************************************************
 *
 * Operators, and support for template expressions
 * 
 * Originally based on article by Klaus Kreft & Angelika Langer
 * http://www.angelikalanger.com/Articles/Cuj/ExpressionTemplates/ExpressionTemplates.htm
 * 
 * Parts adapted from Blitz++ library
 *
 **************************************************************************/

#ifndef __EXPR_H__
#define __EXPR_H__

#warning expr.hxx is deprecated. Do not use!

#include <field3d.hxx>
#include <field2d.hxx>
#include <bout/mesh.hxx>

/// Literal class to capture BoutReal values in expressions
class Literal {
 public:
  /// Type of this expression
  using type = Literal;

  Literal(BoutReal v) : val(v) {}
  ~Literal() {}
  BoutReal operator()(int x, int y, int z) const {return val;}
private:
  const BoutReal val;
};

class Field3DExpr {
public:
  using type = Field3D;
  
  Field3DExpr(const Field3D &f) : data(&f(0,0,0)) {}
  const BoutReal& operator()(int x, int y, int z) const { return data[(x*bout::globals::mesh->LocalNy + y)*bout::globals::mesh->LocalNz + z]; }
private:
  const BoutReal *data;
};

class Field2DExpr {
public:
  using type = Field2D;
  
  Field2DExpr(const Field2D &f) : data(&f(0,0)) {}
  const BoutReal& operator()(int x, int y, int z) const { return data[x*bout::globals::mesh->LocalNy + y]; }
private:
  const BoutReal *data;
};

/// Expression traits, to convert doubles etc. to Literal

template <class ExprT>
struct exprTraits {
  using expr_type = ExprT;
};

template <>
struct exprTraits<double> {
  using expr_type = Literal;
};

template <>
struct exprTraits<float> {
  using expr_type = Literal;
};

template <>
struct exprTraits<int> {
  using expr_type = Literal;
};

///////////////////////////////////////////////
// asExpr: convert objects to expressions

template <typename T>
struct asExpr {
  using type = T;
  static const T& getExpr(const T& x) {return x;}
};

template <>
struct asExpr<int> {
  using type = Literal;
  static const Literal getExpr(const int& x) {return Literal(x);}
};

template <>
struct asExpr<double> {
  using type = Literal;
  static const Literal getExpr(const double& x) {return Literal(x);}
};

template <>
struct asExpr<float> {
  using type = Literal;
  static const Literal getExpr(const float& x) {return Literal(x);}
};

template <>
struct asExpr<Field3D> {
  using type = Field3DExpr;
  static const Field3DExpr getExpr(const Field3D& x) {return Field3DExpr(x);}
};

/////////////////////////////////////////////////////////////
// Type promotion. Work out the type of a calculation,
// based on the type of the arguments

template<typename Lhs, typename Rhs> // If in doubt, convert to Field3D
struct PromoteType {
  using type = Field3D;
};

/////////////////////////////////////////////////////////////
// Binary expressions

template <class ExprT1, class ExprT2, class BinOp> 
class BinaryExpr { 
public:
  BinaryExpr(const ExprT1 &e1, const ExprT2 &e2) 
    : _expr1(e1),_expr2(e2) {
  }
  
  // Work out the type of the inputs
  using ltype = typename exprTraits<ExprT1>::expr_type;
  using rtype = typename exprTraits<ExprT2>::expr_type;
  
  /// Type of the resulting expression
  using type = typename PromoteType<ltype, rtype>::type;
  
  BoutReal operator()(int x, int y, int z) const {
    return BinOp::apply((_expr1)(x,y,z),(_expr2)(x,y,z));
  }
  
private:
  ltype const _expr1;
  rtype const _expr2;
};

template<typename ExprT1, typename ExprT2, class name>
struct BinaryResult {
  using arg1 = typename asExpr<ExprT1>::type;
  using arg2 = typename asExpr<ExprT2>::type;
  using type = BinaryExpr<arg1, arg2,name>;
};

/// Binary operator classes

#define DEFINE_BINARY_OP(name,op)               \
struct name {                                   \
  template<typename T>                          \
  static inline T                               \
  apply(T a, T b)                               \
  { return a op b; }                            \
};

DEFINE_BINARY_OP(Add,+)
DEFINE_BINARY_OP(Subtract,-)
DEFINE_BINARY_OP(Multiply,*)
DEFINE_BINARY_OP(Divide,/)

struct Power {
  template<typename T> 
  static inline T apply(T a, T b) {
    return pow(a, b); 
  }
};

/// Define functions add, mul which use operator structs
#define DEFINE_OVERLOAD_FUNC(name, func)                                   \
  template  <typename ExprT1, typename ExprT2>                             \
  typename BinaryResult<ExprT1,ExprT2,name>::type                          \
  func(const ExprT1 &e1, const ExprT2 &e2) {                               \
    using type = typename BinaryResult<ExprT1,ExprT2,name>::type;          \
    return type(asExpr<ExprT1>::getExpr(e1), asExpr<ExprT2>::getExpr(e2)); \
  }

/// Addition of two Expressions
DEFINE_OVERLOAD_FUNC(Add, add);
/// Multiplication of two Expressions
DEFINE_OVERLOAD_FUNC(Multiply, mul);

/// A function to evaluate expressions
template<typename Expr>
const Field3D eval3D(Expr e) {
  Field3D result;
  result.allocate();
  for(int i=0;i<bout::globals::mesh->LocalNx;i++)
    for(int j=0;j<bout::globals::mesh->LocalNy;j++)
      for(int k=0;k<bout::globals::mesh->LocalNz;k++)
	result(i,j,k) = e(i,j,k);
  return result;
}

#endif // __EXPR_H__
