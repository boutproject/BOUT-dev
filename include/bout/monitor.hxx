#pragma once

#include "bout_types.hxx"

class Solver;

template <typename t>
t bout_abs(t a){
  return a>0? a:-1;
}

template <typename t>
bool isMultiple(t a, t b){
  t min = a>b?b:a;
  t max = a>b?a:b;
  int ratio = (max/min)+.5;
  t error = ratio*min - max;
  if (bout_abs(error/max) > 1e-5)
    return false;
  return true;
}


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
