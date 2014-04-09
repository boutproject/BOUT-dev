
#include "fieldgenerators.hxx"

#include <bout/constants.hxx>
#include <utils.hxx>

//////////////////////////////////////////////////////////

FieldGenerator* FieldSin::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to sin function. Expecting 1, got %d", args.size());
  }
  
  return new FieldSin(args.front());
}

BoutReal FieldSin::generate(double x, double y, double z, double t) {
  return sin(gen->generate(x,y,z,t));
}

FieldGenerator* FieldCos::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to cos function. Expecting 1, got %d", args.size());
  }
  
  return new FieldCos(args.front());
}

BoutReal FieldCos::generate(double x, double y, double z, double t) {
  return cos(gen->generate(x,y,z,t));
}

FieldGenerator* FieldSinh::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to sinh function. Expecting 1, got %d", args.size());
  }
  
  return new FieldSinh(args.front());
}

BoutReal FieldSinh::generate(double x, double y, double z, double t) {
  return sinh(gen->generate(x,y,z,t));
}

FieldGenerator* FieldCosh::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to cosh function. Expecting 1, got %d", args.size());
  }
  
  return new FieldCosh(args.front());
}

BoutReal FieldCosh::generate(double x, double y, double z, double t) {
  return cosh(gen->generate(x,y,z,t));
}

FieldGenerator* FieldTanh::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to tanh function. Expecting 1, got ", args.size());
  }
  return new FieldTanh(args.front());
}

BoutReal FieldTanh::generate(double x, double y, double z, double t) {
  return tanh(gen->generate(x,y,z,t));
}

FieldGenerator* FieldGaussian::clone(const list<FieldGenerator*> args) {
  if((args.size() < 1) || (args.size() > 2)) {
    throw ParseException("Incorrect number of arguments to gaussian function. Expecting 1 or 2, got ", args.size());
  }
  
  FieldGenerator *xin = args.front();
  FieldGenerator *sin;
  if(args.size() == 2) {
    sin = args.back(); // Optional second argument
  }else
    sin = new FieldValue(1.0);
  
  return new FieldGaussian(xin, sin);
}

BoutReal FieldGaussian::generate(double x, double y, double z, double t) {
  BoutReal sigma = s->generate(x,y,z,t);
  return exp(-SQ(X->generate(x,y,z,t)/sigma)/2.) / (sqrt(TWOPI) * sigma);
}

FieldGenerator* FieldAbs::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to abs function. Expecting 1, got %d", args.size());
  }
  
  return new FieldAbs(args.front());
}

BoutReal FieldAbs::generate(double x, double y, double z, double t) {
  return fabs(gen->generate(x,y,z,t));
}

FieldGenerator* FieldSqrt::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to sqrt function. Expecting 1, got %d", args.size());
  }
  
  return new FieldSqrt(args.front());
}

BoutReal FieldSqrt::generate(double x, double y, double z, double t) {
  return sqrt(gen->generate(x,y,z,t));
}

FieldGenerator* FieldHeaviside::clone(const list<FieldGenerator*> args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to heaviside function. Expecting 1, got %d", args.size());
  }
  
  return new FieldHeaviside(args.front());
}

BoutReal FieldHeaviside::generate(double x, double y, double z, double t) {
  return (gen->generate(x,y,z,t) > 0.0) ? 1.0 : 0.0;
}
