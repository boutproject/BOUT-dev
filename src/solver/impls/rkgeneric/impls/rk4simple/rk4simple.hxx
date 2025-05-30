
class RK4SIMPLEScheme;

#ifndef BOUT_RK4SIMPLE_SCHEME_H
#define BOUT_RK4SIMPLE_SCHEME_H

#include <bout/rkscheme.hxx>
#include <bout/utils.hxx>

class RK4SIMPLEScheme : public RKScheme {
public:
  RK4SIMPLEScheme(Options* options);

  BoutReal setOutputStates(const Array<BoutReal>& start, BoutReal dt,
                           Array<BoutReal>& resultFollow);
};

namespace {
RegisterRKScheme<RK4SIMPLEScheme> registerrkscheme4simple(RKSCHEME_RK4);
}

#endif // BOUT_RK4SIMPLE_SCHEME_H
