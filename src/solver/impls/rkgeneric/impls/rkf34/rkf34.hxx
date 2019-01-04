#ifndef __RKF34_SCHEME_H__
#define __RKF34_SCHEME_H__

#include "bout/rkscheme.hxx"
class Options;

class RKF34Scheme : public RKScheme{
 public:
  RKF34Scheme(Options *options);
  ~RKF34Scheme();
 private:

};

#endif // __RKF34_SCHEME_H__
