
class RKF34Scheme;

#ifndef __RKF34_SCHEME_H__
#define __RKF34_SCHEME_H__

#include <bout/rkscheme.hxx>
#include <utils.hxx>

class RKF34Scheme : public RKScheme {
public:
  RKF34Scheme(Options* options);
};

#endif // __RKF34_SCHEME_H__
