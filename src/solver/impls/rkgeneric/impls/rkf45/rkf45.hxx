#ifndef __RKF45_SCHEME_H__
#define __RKF45_SCHEME_H__

#include "bout/rkscheme.hxx"

class Options;

class RKF45Scheme : public RKScheme{
 public:
  RKF45Scheme(Options *options);
  ~RKF45Scheme();
 private:

};

#endif // __RKF45_SCHEME_H__
