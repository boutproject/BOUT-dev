#include "fieldgenerators.hxx"

#include <bout/constants.hxx>
#include <utils.hxx>

//////////////////////////////////////////////////////////

std::shared_ptr<FieldGenerator> FieldSin::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to sin function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldSin(args.front()));
}

BoutReal FieldSin::generate(double x, double y, double z, double t) {
  return sin(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldCos::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to cos function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldCos(args.front()));
}

BoutReal FieldCos::generate(double x, double y, double z, double t) {
  return cos(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldSinh::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to sinh function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldSinh(args.front()));
}

BoutReal FieldSinh::generate(double x, double y, double z, double t) {
  return sinh(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldCosh::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to cosh function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldCosh(args.front()));
}

BoutReal FieldCosh::generate(double x, double y, double z, double t) {
  return cosh(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldTanh::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to tanh function. Expecting 1, got ", args.size());
  }
  return std::shared_ptr<FieldGenerator>( new FieldTanh(args.front()));
}

BoutReal FieldTanh::generate(double x, double y, double z, double t) {
  return tanh(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldGaussian::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if((args.size() < 1) || (args.size() > 2)) {
    throw ParseException("Incorrect number of arguments to gaussian function. Expecting 1 or 2, got ", args.size());
  }

  std::shared_ptr<FieldGenerator> xin = args.front();
  std::shared_ptr<FieldGenerator> sin;
  if(args.size() == 2) {
    sin = args.back(); // Optional second argument
  }else
    sin = std::shared_ptr<FieldGenerator>( new FieldValue(1.0));

  return std::shared_ptr<FieldGenerator>( new FieldGaussian(xin, sin));
}

BoutReal FieldGaussian::generate(double x, double y, double z, double t) {
  BoutReal sigma = s->generate(x,y,z,t);
  return exp(-SQ(X->generate(x,y,z,t)/sigma)/2.) / (sqrt(TWOPI) * sigma);
}

std::shared_ptr<FieldGenerator> FieldAbs::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to abs function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldAbs(args.front()));
}

BoutReal FieldAbs::generate(double x, double y, double z, double t) {
  return fabs(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldSqrt::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to sqrt function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldSqrt(args.front()));
}

BoutReal FieldSqrt::generate(double x, double y, double z, double t) {
  return sqrt(gen->generate(x,y,z,t));
}

std::shared_ptr<FieldGenerator> FieldHeaviside::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to heaviside function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldHeaviside(args.front()));
}

BoutReal FieldHeaviside::generate(double x, double y, double z, double t) {
  return (gen->generate(x,y,z,t) > 0.0) ? 1.0 : 0.0;
}

std::shared_ptr<FieldGenerator> FieldErf::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 1) {
    throw ParseException("Incorrect number of arguments to erf function. Expecting 1, got %d", args.size());
  }

  return std::shared_ptr<FieldGenerator>( new FieldErf(args.front()));
}

BoutReal FieldErf::generate(double x, double y, double z, double t) {
  return erf(gen->generate(x,y,z,t));
}

//////////////////////////////////////////////////////////
// Ballooning transform
// Use a truncated Ballooning transform to enforce periodicity in y and z

std::shared_ptr<FieldGenerator> FieldBallooning::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  int n = ball_n;
  switch(args.size()) {
  case 2: {
    // Second optional argument is ball_n, an integer
    // This should probably warn if arg isn't constant
    n = ROUND( args.back()->generate(0,0,0,0) );
  } // Fall through
  case 1: {
    return std::shared_ptr<FieldGenerator>( new FieldBallooning(mesh, args.front(), n));
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
      value += arg->generate(x,y - i*TWOPI,z + i*ts*TWOPI/mesh->coordinates()->zlength(),t);

      // y + i * 2pi
      value += arg->generate(x,y + i*TWOPI,z - i*ts*TWOPI/mesh->coordinates()->zlength(),t);
    }
    return value;
  }

  // Open surfaces. Not sure what to do, so set to zero
  return 0.0;
}

////////////////////////////////////////////////////////////////

FieldMixmode::FieldMixmode(std::shared_ptr<FieldGenerator> a, BoutReal seed)
    : arg(std::move(a)) {
  // Calculate the phases -PI to +PI
  // using genRand [0,1]

  for(int i=0;i<14;i++)
    phase[i] = PI * (2.*genRand(seed + i) - 1.);
}

std::shared_ptr<FieldGenerator> FieldMixmode::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  BoutReal seed = 0.5;
  switch(args.size()) {
  case 2: {
    // Second optional argument is the seed, which should be a constant
    seed = args.back()->generate(0,0,0,0);
  } // Fall through
  case 1: {
    return std::shared_ptr<FieldGenerator>( new FieldMixmode(args.front(), seed));
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
std::shared_ptr<FieldGenerator> FieldTanhHat::clone(const list<std::shared_ptr<FieldGenerator> > args) {
  if(args.size() != 4) {
    throw ParseException("Incorrect number of arguments to TanhHat function. Expecting 4, got %d", args.size());
  }

  // As lists are not meant to be indexed, we may use an iterator to get the
  // input arguments instead
  // Create the iterator
  list<std::shared_ptr<FieldGenerator> >::const_iterator it = args.begin();
  // Assign the input arguments to the input of the constructor and advance the
  // iterator
  std::shared_ptr<FieldGenerator> xin = *it;
  std::advance(it, 1);
  std::shared_ptr<FieldGenerator> widthin = *it;
  std::advance(it, 1);
  std::shared_ptr<FieldGenerator> centerin = *it;
  std::advance(it, 1);
  std::shared_ptr<FieldGenerator> steepnessin = *it;

  // Call the constructor
  return std::shared_ptr<FieldGenerator>( new FieldTanhHat(xin, widthin, centerin, steepnessin));
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

BoutReal FieldRealY::generate(double x, double y, double z, double t,
                              const DataIterator &i, Mesh *localmesh) {
  BoutReal *offset = doCache(localmesh);
  int gny = localmesh->GlobalNy - mesh->ystart * 2;
  int iy = y / 2 / PI * gny * 2 + 10.5;
  BoutReal res;
  Field2D &dy = localmesh->coordinates()->dy;
  res = offset[i.x];
  for (int y = 0; y < i.y; ++y) {
    res += dy(i.x, y);
  }
  if (iy % 2) { // centered
    res -= dy[i] * .5;
  } else {
    // res+=dy[i];
  }
  // printf("%g %8g %4d %4d %6g\n",offset[i.x],y/2/PI*gny*2,i.y,iy,res);
  return res;
}

BoutReal *FieldRealY::doCache(Mesh *localmesh) {
  auto it = cache.begin();
  for (const auto &cm : cached) {
    if (cm == localmesh) {
      return *it;
    }
    ++it;
  }
  cached.push_back(localmesh);
  BoutReal *offset = (BoutReal *)malloc(sizeof(BoutReal) * localmesh->LocalNx);
  BoutReal tmp[localmesh->LocalNx];
  BoutReal dummy1, dummy2;
  if (localmesh->firstY()) {
    for (int x = 0; x < localmesh->LocalNx; ++x) {
      offset[x] = 0;
      for (int y = localmesh->ystart - 1; y >= 0; --y) {
        offset[x] -= localmesh->coordinates()->dy(x, y);
      }
    }
  } else {
    localmesh->wait(localmesh->receiveFromProc(localmesh->getXProcIndex(),
                                               localmesh->getYProcIndex() - 1, offset,
                                               localmesh->LocalNx, 666));
    // printf("%d: receive %g\n",localmesh->getYProcIndex(),offset[localmesh->xstart]);
  }
  if (!localmesh->lastY()) {
    for (int x = 0; x < localmesh->LocalNx; ++x) {
      tmp[x] = offset[x];
      for (int y = 0; y < localmesh->LocalNy - localmesh->ystart * 2; ++y) {
        tmp[x] += localmesh->coordinates()->dy(x, y);
      }
    }
    // printf("%d: send %g\n",localmesh->getYProcIndex(),tmp[localmesh->xstart]);
    MPI_Request out =
        localmesh->sendToProc(localmesh->getXProcIndex(), localmesh->getYProcIndex() + 1,
                              tmp, localmesh->LocalNx, 666);
    MPI_Status stat;
    MPI_Wait(&out, &stat);
    // ASSERT2(stat == MPI_SUCCESS);
  }
  // dummy1=localmesh->getYProcIndex();
  // dummy2=-1;
  // if (! localmesh->firstY()){
  //   localmesh->sendToProc(localmesh->getXProcIndex(),
  //   localmesh->getYProcIndex()-1,&dummy1,1,667);
  // }
  // if (! localmesh->lastY()) {
  //   localmesh->wait(localmesh->receiveFromProc(localmesh->getXProcIndex(),
  //   localmesh->getYProcIndex()+1,&dummy2,1,667));
  // }
  // printf("%d: debug %g %g\n",localmesh->getYProcIndex(),dummy1,dummy2);

  cache.push_back(offset);
  return offset;
}
