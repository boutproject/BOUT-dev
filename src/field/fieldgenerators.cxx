
#include "fieldgenerators.hxx"

#include <bout/constants.hxx>
#include <utils.hxx>

//////////////////////////////////////////////////////////

FieldGeneratorPtr FieldSin::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to sin function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldSin>(args.front());
}

BoutReal FieldSin::generate(double x, double y, double z, double t) {
  return sin(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldCos::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to cos function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldCos>(args.front());
}

BoutReal FieldCos::generate(double x, double y, double z, double t) {
  return cos(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldSinh::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to sinh function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldSinh>(args.front());
}

BoutReal FieldSinh::generate(double x, double y, double z, double t) {
  return sinh(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldCosh::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to cosh function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldCosh>(args.front());
}

BoutReal FieldCosh::generate(double x, double y, double z, double t) {
  return cosh(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldTanh::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to tanh function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }
  return std::make_shared<FieldTanh>(args.front());
}

BoutReal FieldTanh::generate(double x, double y, double z, double t) {
  return tanh(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldGaussian::clone(const std::list<FieldGeneratorPtr> args) {
  if ((args.size() < 1) || (args.size() > 2)) {
    throw ParseException(
        "Incorrect number of arguments to gaussian function. Expecting 1 or 2, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  FieldGeneratorPtr xin = args.front();
  FieldGeneratorPtr sin;
  if(args.size() == 2) {
    sin = args.back(); // Optional second argument
  }else
    sin = std::make_shared<FieldValue>(1.0);

  return std::make_shared<FieldGaussian>(xin, sin);
}

BoutReal FieldGaussian::generate(double x, double y, double z, double t) {
  BoutReal sigma = s->generate(x,y,z,t);
  return exp(-SQ(X->generate(x,y,z,t)/sigma)/2.) / (sqrt(TWOPI) * sigma);
}

FieldGeneratorPtr FieldAbs::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to abs function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldAbs>(args.front());
}

BoutReal FieldAbs::generate(double x, double y, double z, double t) {
  return std::fabs(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldSqrt::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to sqrt function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldSqrt>(args.front());
}

BoutReal FieldSqrt::generate(double x, double y, double z, double t) {
  return sqrt(gen->generate(x, y, z, t));
}

FieldGeneratorPtr FieldHeaviside::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to heaviside function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldHeaviside>(args.front());
}

BoutReal FieldHeaviside::generate(double x, double y, double z, double t) {
  return (gen->generate(x, y, z, t) > 0.0) ? 1.0 : 0.0;
}

FieldGeneratorPtr FieldErf::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to erf function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldErf>(args.front());
}

BoutReal FieldErf::generate(double x, double y, double z, double t) {
  return erf(gen->generate(x,y,z,t));
}

//////////////////////////////////////////////////////////
// Ballooning transform
// Use a truncated Ballooning transform to enforce periodicity in y and z

FieldGeneratorPtr FieldBallooning::clone(const std::list<FieldGeneratorPtr> args) {
  int n = ball_n;
  switch(args.size()) {
  case 2: {
    // Second optional argument is ball_n, an integer
    // This should probably warn if arg isn't constant
    n = ROUND( args.back()->generate(0,0,0,0) );
  } // Fall through
  case 1: {
    return std::make_shared<FieldBallooning>(mesh, args.front(), n);
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
  Coordinates* coords = mesh->getCoordinates();

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
      value += arg->generate(x,y - i*TWOPI,z + i*ts*TWOPI/coords->zlength(),t);

      // y + i * 2pi
      value += arg->generate(x,y + i*TWOPI,z - i*ts*TWOPI/coords->zlength(),t);
    }
    return value;
  }

  // Open surfaces. Not sure what to do, so set to zero
  return 0.0;
}

////////////////////////////////////////////////////////////////

FieldMixmode::FieldMixmode(FieldGeneratorPtr a, BoutReal seed) : arg(std::move(a)) {
  // Calculate the phases -PI to +PI
  // using genRand [0,1]

  for(int i=0;i<14;i++)
    phase[i] = PI * (2.*genRand(seed + i) - 1.);
}

FieldGeneratorPtr FieldMixmode::clone(const std::list<FieldGeneratorPtr> args) {
  BoutReal seed = 0.5;
  switch(args.size()) {
  case 2: {
    // Second optional argument is the seed, which should be a constant
    seed = args.back()->generate(0,0,0,0);
  } // Fall through
  case 1: {
    return std::make_shared<FieldMixmode>(args.front(), seed);
  }
  };

  throw ParseException("mixmode function must have one or two arguments");
}

BoutReal FieldMixmode::generate(double x, double y, double z, double t) {
  BoutReal result = 0.0;

  // A mixture of mode numbers
  for(int i=0;i<14;i++) {
    // This produces a spectrum which is peaked around mode number 4
    result += ( 1./SQ(1. + std::abs(i - 4)) ) *
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

  // Start x between 0 and 1 (exclusive)
  const BoutReal A = 0.01, B = 1.23456789;
  BoutReal x = (A + fmod(seed,B)) / (B + 2.*A);

  // Iterate logistic map
  for(int i=0;i!=niter;++i)
    x = 3.99 * x * (1. - x);

  return x;
}

//////////////////////////////////////////////////////////
// TanhHat
FieldGeneratorPtr FieldTanhHat::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 4) {
    throw ParseException(
        "Incorrect number of arguments to TanhHat function. Expecting 4, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  // As lists are not meant to be indexed, we may use an iterator to get the
  // input arguments instead
  // Create the iterator
  auto it = args.begin();
  // Assign the input arguments to the input of the constructor and advance the
  // iterator
  FieldGeneratorPtr xin = *it;
  std::advance(it, 1);
  FieldGeneratorPtr widthin = *it;
  std::advance(it, 1);
  FieldGeneratorPtr centerin = *it;
  std::advance(it, 1);
  FieldGeneratorPtr steepnessin = *it;

  // Call the constructor
  return std::make_shared<FieldTanhHat>(xin, widthin, centerin, steepnessin);
}

BoutReal FieldTanhHat::generate(double x, double y, double z, double t) {
  // The following are constants
  BoutReal w = width    ->generate(0,0,0,0);
  BoutReal c = center   ->generate(0,0,0,0);
  BoutReal s = steepness->generate(0,0,0,0);
  return 0.5*(
                 tanh( s*(X->generate(x,y,z,t) - (c - 0.5*w)) )
               - tanh( s*(X->generate(x,y,z,t) - (c + 0.5*w)) )
             );
}
