
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

//////////////////////////////////////////////////////////
// Ballooning transform
// Use a truncated Ballooning transform to enforce periodicity in y and z

FieldGenerator* FieldBallooning::clone(const list<FieldGenerator*> args) {
  int n = ball_n;
  switch(args.size()) {
  case 2: {
    // Second optional argument is ball_n, an integer
    // This should probably warn if arg isn't constant
    n = ROUND( args.back()->generate(0,0,0,0) );
  } // Fall through
  case 1: {
    return new FieldBallooning(mesh, args.front(), n);
    break;
  }
  };

  throw ParseException("ballooning function must have one or two arguments");
}

BoutReal FieldBallooning::generate(double x, double y, double z, double t) {
  if(!mesh)
    throw BoutException("ballooning function needs a valid mesh");
  if(ball_n < 1)
    throw BoutException("ballooning function ball_n less than 1");

  BoutReal ts; // Twist-shift angle

  // Need to find the nearest flux surface (x index)
  // This assumes that mesh->GlobalX is linear in x index
  BoutReal dx = (mesh->GlobalX(mesh->xend) - mesh->GlobalX(mesh->xstart)) /
    (mesh->xend - mesh->xstart);
  int jx = ROUND((x - mesh->GlobalX(0)) / dx);

  if(mesh->periodicY(jx, ts)) {
    // Start with the value at this point
    BoutReal value = arg->generate(x,y,z,t);

    for(int i=1; i<= ball_n; i++) {
      // y - i * 2pi
      value += arg->generate(x,y - i*TWOPI,z + i*ts,t);

      // y + i * 2pi
      value += arg->generate(x,y + i*TWOPI,z - i*ts,t);
    }
    return value;
  }

  // Open surfaces. Not sure what to do, so set to zero
  return 0.0;
}

////////////////////////////////////////////////////////////////

FieldMixmode::FieldMixmode(FieldGenerator* a, BoutReal seed) : arg(a) {
  // Calculate the phases -PI to +PI
  // using genRand [0,1]

  for(int i=0;i<14;i++)
    phase[i] = PI * (2.*genRand(seed + i) - 1.);
}

FieldGenerator* FieldMixmode::clone(const list<FieldGenerator*> args) {
  BoutReal seed = 0.5;
  switch(args.size()) {
  case 2: {
    // Second optional argument is the seed, which should be a constant
    seed = args.back()->generate(0,0,0,0);
  } // Fall through
  case 1: {
    return new FieldMixmode(args.front(), seed);
  }
  };

  throw ParseException("mixmode function must have one or two arguments");
}

BoutReal FieldMixmode::generate(double x, double y, double z, double t) {
  BoutReal result = 0.0;

  // A mixture of mode numbers
  for(int i=0;i<14;i++) {
    // This produces a spectrum which is peaked around mode number 4
    result += ( 1./SQ(1. + abs(i - 4)) ) *
      cos(i * arg->generate(x,y,z,t) + phase[i]);
  }

  return result;
}

BoutReal FieldMixmode::genRand(BoutReal seed) {
  // Make sure seed is
  if(seed < 0.0)
    seed *= -1;

  // Round the seed to get the number of iterations
  int niter = 11 + (23 + ROUND(seed)) % 79;

  // Start x between 0 and 1
  const BoutReal A = 0.01, B = 1.23456789;
  BoutReal x = (A + fmod(seed,B)) / (B - 2.*A);

  // Iterate logistic map
  for(int i=0;i!=niter;++i)
    x = 3.99 * x * (1. - x);

  return x;
}
