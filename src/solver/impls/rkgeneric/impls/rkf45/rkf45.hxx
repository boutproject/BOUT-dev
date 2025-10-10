
class RKF45Scheme;

#ifndef BOUT_RKF45_SCHEME_H
#define BOUT_RKF45_SCHEME_H

#include <bout/rkscheme.hxx>
#include <bout/utils.hxx>

class RKF45Scheme : public RKScheme {
public:
  RKF45Scheme(Options* options);
};

namespace {
RegisterRKScheme<RKF45Scheme> registerrkschemef45(RKSCHEME_RKF45);
}

#endif // BOUT_RKF45_SCHEME_H
