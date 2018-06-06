#pragma once

#include "bout_types.hxx"

class Solver;
/// Check if two numbers are multiples of each other, with an accuracy
/// of 1e-12
bool isMultiple(BoutReal a, BoutReal b);

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

/// signature of legacy functions
typedef int (* MonitorFuncPointer )(Solver *solver, BoutReal simtime, int iter, int NOUT);

/// Class to wrap legacy monitor functions
class MonitorFunc: public Monitor{
public:
  MonitorFunc(MonitorFuncPointer pntr);
  virtual ~MonitorFunc(){};
  virtual int call(Solver * solver, BoutReal time, int iter, int nout) override;
private:
  MonitorFuncPointer callFunc;
};
