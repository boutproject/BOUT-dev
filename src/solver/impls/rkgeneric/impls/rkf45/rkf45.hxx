
class RKF45Scheme;

#ifndef __RKF45_SCHEME_H__
#define __RKF45_SCHEME_H__

#include <bout/rkscheme.hxx>
#include <utils.hxx>

class RKF45Scheme : public RKScheme {
public:
  RKF45Scheme(Options* options);
};

#endif // __RKF45_SCHEME_H__
