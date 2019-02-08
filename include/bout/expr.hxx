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
  typedef Literal type; ///< Type of this expression

  Literal(BoutReal v) : val(v) {}
  ~Literal() {}
  BoutReal operator()(int x, int y, int z) const {return val;}
private:
  const BoutReal val;
};

class Field3DExpr {
public:
  typedef Field3D type;
  
  Field3DExpr(const Field3D &f) : data(&f(0,0,0)) {}
  const BoutReal& operator()(int x, int y, int z) const { return data[(x*bout::globals::mesh->LocalNy + y)*bout::globals::mesh->LocalNz + z]; }
private:
  const BoutReal *data;
};

class Field2DExpr {
public:
  typedef Field2D type;
  
  Field2DExpr(const Field2D &f) : data(&f(0,0)) {}
  const BoutReal& operator()(int x, int y, int z) const { return data[x*bout::globals::mesh->LocalNy + y]; }
private:
  const BoutReal *data;
};

/// Expression traits, to convert doubles etc. to Literal

template <class ExprT>
struct exprTraits {
  typedef ExprT expr_type;
};

template <>
struct exprTraits<double> {
  typedef Literal expr_type;
};

template <>
struct exprTraits<float> {
  typedef Literal expr_type;
};

template <>
struct exprTraits<int> {
  typedef Literal expr_type;
};

///////////////////////////////////////////////
// asExpr: convert objects to expressions

template <typename T>
struct asExpr {
  typedef T type;
  static const T& getExpr(const T& x) {return x;}
};

template <>
struct asExpr<int> {
  typedef Literal type;
  static const Literal getExpr(const int& x) {return Literal(x);}
};

template <>
struct asExpr<double> {
  typedef Literal type;
  static const Literal getExpr(const double& x) {return Literal(x);}
};

template <>
struct asExpr<float> {
  typedef Literal type;
  static const Literal getExpr(const float& x) {return Literal(x);}
};

template <>
struct asExpr<Field3D> {
  typedef Field3DExpr type;
  static const Field3DExpr getExpr(const Field3D& x) {return Field3DExpr(x);}
};

/////////////////////////////////////////////////////////////
// Type promotion. Work out the type of a calculation,
// based on the type of the arguments

template<typename Lhs, typename Rhs> // If in doubt, convert to Field3D
struct PromoteType {
  typedef Field3D type;
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
  typedef typename exprTraits<ExprT1>::expr_type ltype;
  typedef typename exprTraits<ExprT2>::expr_type rtype;
  
  /// Type of the resulting expression
  typedef typename PromoteType<ltype, rtype>::type type;
  
  BoutReal operator()(int x, int y, int z) const {
    return BinOp::apply((_expr1)(x,y,z),(_expr2)(x,y,z));
  }
  
private:
  ltype const _expr1;
  rtype const _expr2;
};

template<typename ExprT1, typename ExprT2, class name>
struct BinaryResult {
  typedef typename asExpr<ExprT1>::type arg1;
  typedef typename asExpr<ExprT2>::type arg2;
  typedef BinaryExpr<arg1, arg2,name> type;
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
    typedef typename BinaryResult<ExprT1,ExprT2,name>::type type;          \
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
