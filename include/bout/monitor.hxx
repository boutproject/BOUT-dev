#ifndef __MONITOR_H__
#define __MONITOR_H__

#include "bout_types.hxx"
#include "bout/assert.hxx"

#include <cmath>

class Solver;

/// Return true if either \p a is a multiple of \p b or vice-versa
///
/// Assumes both arguments are greater than zero
inline bool isMultiple(BoutReal a, BoutReal b) {
  ASSERT2(a > 0);
  ASSERT2(b > 0);
  return (std::fmod(a, b) == 0) || (std::fmod(b, a) == 0);
}

/// Monitor baseclass for the Solver
/// Can be called ether with a specified frequency, or with the
/// frequency of the BOUT++ output monitor.
class Monitor{
  friend class Solver; ///< needs access to timestep and freq
public:
  /// A timestep of -1 defaults to the the frequency of the BOUT++
  /// output monitor
  Monitor(BoutReal timestep_=-1):timestep(timestep_){};
  virtual ~Monitor(){};
  /// call is called by the solver after timestep_ has passed
  virtual int call(Solver * solver, BoutReal time, int iter, int nout)=0;
  /// cleanup is called when a clean shutdown is initiated
  virtual void cleanup(){};
  /// compare two monitors, check whether they are identicall
  bool operator==(const Monitor& rhs) const;
private:
  BoutReal timestep;
  int freq;
};

#endif // __MONITOR_H__
