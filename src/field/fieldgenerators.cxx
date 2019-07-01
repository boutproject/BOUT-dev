
#include "fieldgenerators.hxx"

#include <bout/constants.hxx>
#include <utils.hxx>

using bout::generator::Context;

//////////////////////////////////////////////////////////

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

BoutReal FieldGaussian::generate(const Context& pos) {
  BoutReal sigma = s->generate(pos);
  return exp(-SQ(X->generate(pos)/sigma)/2.) / (sqrt(TWOPI) * sigma);
}

FieldGeneratorPtr FieldHeaviside::clone(const std::list<FieldGeneratorPtr> args) {
  if (args.size() != 1) {
    throw ParseException(
        "Incorrect number of arguments to heaviside function. Expecting 1, got %lu",
        static_cast<unsigned long>(args.size()));
  }

  return std::make_shared<FieldHeaviside>(args.front());
}

BoutReal FieldHeaviside::generate(const Context& pos) {
  return (gen->generate(pos) > 0.0) ? 1.0 : 0.0;
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
    n = ROUND( args.back()->generate(Context()) );
  } // Fall through
  case 1: {
    return std::make_shared<FieldBallooning>(mesh, args.front(), n);
  }
  };

  throw ParseException("ballooning function must have one or two arguments");
}

BoutReal FieldBallooning::generate(const Context& pos) {
  Mesh *localmesh = pos.getMesh();
  if (!localmesh)
    throw BoutException("ballooning function needs a valid mesh");
  if (ball_n < 1)
    throw BoutException("ballooning function ball_n less than 1");

  BoutReal ts; // Twist-shift angle
  Coordinates* coords = localmesh->getCoordinates();

  // Need to find the nearest flux surface (x index)
  // This assumes that localmesh->GlobalX is linear in x index
  BoutReal dx = (localmesh->GlobalX(localmesh->xend) - localmesh->GlobalX(localmesh->xstart))
                / (localmesh->xend - localmesh->xstart);
  int jx = ROUND((pos.x() - localmesh->GlobalX(0)) / dx);

  if (localmesh->periodicY(jx, ts)) {
    // Start with the value at this point
    BoutReal value = arg->generate(pos);

    for (int i = 1; i <= ball_n; i++) {
      // y - i * 2pi
      value += arg->generate(Context(pos).set(
          "y", pos.y() - i * TWOPI,
          "z", pos.z() + i * ts * TWOPI / coords->zlength()));

      value += arg->generate(Context(pos).set(
          "y", pos.y() + i * TWOPI,
          "z", pos.z() - i * ts * TWOPI / coords->zlength()));
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

  for (int i = 0; i < 14; i++)
    phase[i] = PI * (2. * genRand(seed + i) - 1.);
}

FieldGeneratorPtr FieldMixmode::clone(const std::list<FieldGeneratorPtr> args) {
  BoutReal seed = 0.5;
  switch(args.size()) {
  case 2: {
    // Second optional argument is the seed, which should be a constant
    seed = args.back()->generate(Context());
  } // Fall through
  case 1: {
    return std::make_shared<FieldMixmode>(args.front(), seed);
  }
  };

  throw ParseException("mixmode function must have one or two arguments");
}

BoutReal FieldMixmode::generate(const Context& pos) {
  BoutReal result = 0.0;

  // A mixture of mode numbers
  for(int i=0;i<14;i++) {
    // This produces a spectrum which is peaked around mode number 4
    result += ( 1./SQ(1. + std::abs(i - 4)) ) *
      cos(i * arg->generate(pos) + phase[i]);
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

BoutReal FieldTanhHat::generate(const Context& pos) {
  // The following are constants
  BoutReal w = width    ->generate(Context());
  BoutReal c = center   ->generate(Context());
  BoutReal s = steepness->generate(Context());
  return 0.5*(
                 tanh( s*(X->generate(pos) - (c - 0.5*w)) )
               - tanh( s*(X->generate(pos) - (c + 0.5*w)) )
             );
}
