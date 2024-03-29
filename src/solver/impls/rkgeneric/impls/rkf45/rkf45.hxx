
class RKF45Scheme;

#ifndef __RKF45_SCHEME_H__
#define __RKF45_SCHEME_H__

#include <bout/rkscheme.hxx>
#include <bout/utils.hxx>

class RKF45Scheme : public RKScheme {
public:
  RKF45Scheme(Options* options);
};

namespace {
RegisterRKScheme<RKF45Scheme> registerrkschemef45(RKSCHEME_RKF45);
}

#endif // __RKF45_SCHEME_H__
