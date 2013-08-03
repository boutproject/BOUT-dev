
#include "fieldgenerators.hxx"

#include <bout/constants.hxx>
#include <utils.hxx>

BoutReal FieldX::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return fieldmesh->GlobalX(x);
}

BoutReal FieldY::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return TWOPI*fieldmesh->GlobalY(y);
}

BoutReal FieldZ::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return TWOPI*((BoutReal) z) / ((BoutReal) (fieldmesh->ngz-1));
}

//////////////////////////////////////////////////////////

FieldBinary::~FieldBinary() {
  if(lhs)
    delete lhs;
  if(rhs)
    delete rhs;
}

FieldGenerator* FieldBinary::clone(const list<FieldGenerator*> args) {
  if(args.size() != 2)
    return NULL;
  
  return new FieldBinary(args.front(), args.back(), op);
}

BoutReal FieldBinary::generate(const Mesh *fieldmesh, int x, int y, int z) {
  BoutReal lval = lhs->generate(fieldmesh, x,y,z);
  BoutReal rval = rhs->generate(fieldmesh, x,y,z);
  switch(op) {
  case '+': return lval + rval;
  case '-': return lval - rval;
  case '*': return lval * rval;
  case '/': return lval / rval;
  case '^': return pow(lval, rval);
  }
  // Unknown operator. Throw an error?
  return 0.;
}

FieldGenerator* FieldSin::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to sin function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldSin(args.front());
}

BoutReal FieldSin::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return sin(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldCos::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to cos function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldCos(args.front());
}

BoutReal FieldCos::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return cos(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldSinh::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to sinh function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldSinh(args.front());
}

BoutReal FieldSinh::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return sinh(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldCosh::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to cosh function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldCosh(args.front());
}

BoutReal FieldCosh::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return cosh(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldTanh::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to tanh function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  return new FieldTanh(args.front());
}

BoutReal FieldTanh::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return tanh(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldGaussian::clone(const list<FieldGenerator*> args) {
  if((args.size() < 1) || (args.size() > 2)) {
    output << "FieldFactory error: Incorrect number of arguments to gaussian function. Expecting 1 or 2, got " << args.size() << endl;
    return NULL;
  }
  
  FieldGenerator *xin = args.front();
  FieldGenerator *sin;
  if(args.size() == 2) {
    sin = args.back(); // Optional second argument
  }else
    sin = new FieldValue(1.0);
  
  return new FieldGaussian(xin, sin);
}

BoutReal FieldGaussian::generate(const Mesh *fieldmesh, int x, int y, int z) {
  BoutReal sigma = s->generate(fieldmesh, x,y,z);
  return exp(-SQ(X->generate(fieldmesh, x,y,z)/sigma)/2.) / (sqrt(TWOPI) * sigma);
}

FieldGenerator* FieldAbs::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to abs function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldAbs(args.front());
}

BoutReal FieldAbs::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return fabs(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldSqrt::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to sqrt function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldSqrt(args.front());
}

BoutReal FieldSqrt::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return sqrt(gen->generate(fieldmesh, x,y,z));
}

FieldGenerator* FieldHeaviside::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    output << "FieldFactory error: Incorrect number of arguments to heaviside function. Expecting 1, got " << args.size() << endl;
    return NULL;
  }
  
  return new FieldHeaviside(args.front());
}

BoutReal FieldHeaviside::generate(const Mesh *fieldmesh, int x, int y, int z) {
  return (gen->generate(fieldmesh, x,y,z) > 0.0) ? 1.0 : 0.0;
}
