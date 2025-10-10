
class RKF34Scheme;

#ifndef BOUT_RKF34_SCHEME_H
#define BOUT_RKF34_SCHEME_H

#include <bout/rkscheme.hxx>
#include <bout/utils.hxx>

class RKF34Scheme : public RKScheme {
public:
  RKF34Scheme(Options* options);
};

namespace {
RegisterRKScheme<RKF34Scheme> registerrkschemef34(RKSCHEME_RKF34);
}

#endif // BOUT_RKF34_SCHEME_H
