#pragma once

#include "bout_types.hxx"

class Solver;

bool isMultiple(BoutReal a, BoutReal b);

class Monitor{
  friend class Solver; // needs access to timestep and freq
public:
  Monitor(BoutReal timestep_=-1):timestep(timestep_){};
  virtual ~Monitor(){};
  virtual int call(Solver * solver, BoutReal time, int iter, int nout)=0;
  virtual void cleanup(){};
  bool operator==(const Monitor& rhs) const;
private:
  BoutReal timestep;
  int freq;
};

typedef int (* MonitorFuncPointer )(Solver *solver, BoutReal simtime, int iter, int NOUT);

class MonitorFunc: public Monitor{
public:
  MonitorFunc(MonitorFuncPointer pntr);
  virtual ~MonitorFunc(){};
  virtual int call(Solver * solver, BoutReal time, int iter, int nout);
private:
  MonitorFuncPointer callFunc;
};
