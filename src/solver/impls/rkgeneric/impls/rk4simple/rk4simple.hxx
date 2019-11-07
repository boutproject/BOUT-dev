
class RK4SIMPLEScheme;

#ifndef __RK4SIMPLE_SCHEME_H__
#define __RK4SIMPLE_SCHEME_H__

#include <bout/rkscheme.hxx>
#include <utils.hxx>

#include "../../rkschemefactory.hxx"

class RK4SIMPLEScheme : public RKScheme{
public:
  RK4SIMPLEScheme(Options *options);

  BoutReal setOutputStates(const Array<BoutReal> &start,BoutReal dt, Array<BoutReal> &resultFollow);
};

namespace {
RegisterRKScheme<RK4SIMPLEScheme> registerrkscheme4simple(RKSCHEME_RK4);
}

#endif // __RK4SIMPLE_SCHEME_H__
