
class RK4SIMPLEScheme;

#ifndef __RK4SIMPLE_SCHEME_H__
#define __RK4SIMPLE_SCHEME_H__

#include <bout/rkscheme.hxx>
#include <utils.hxx>

class RK4SIMPLEScheme : public RKScheme{
 public:
  RK4SIMPLEScheme(Options *options);
  ~RK4SIMPLEScheme();
  
  BoutReal setOutputStates(const Array<BoutReal> &start,BoutReal dt, Array<BoutReal> &resultFollow);
 private:
  
};

#endif // __RK4SIMPLE_SCHEME_H__
