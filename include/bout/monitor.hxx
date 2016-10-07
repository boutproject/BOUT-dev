#pragma once

template <typename t>
t abs(t a){
  return a>0? a:-1;
}
template <typename t>
bool isMultiple(t a, t b){
  t min = a>b?b:a;
  t max = a>b?a:b;
  int ratio = (max/min)+.5;
  t error = ratio*min - max;
  if (abs(error/max) > 1e-5)
    return false;
  return true;
}


class Monitor{
public:
  virtual int call(Solver * solver, BoutReal time, int iter, int nout)=0;
  bool operator==(const Monitor& rhs) const;
  BoutReal timestep;
  int freq;
};

typedef int (* MonitorFuncPointer )(Solver *solver, BoutReal simtime, int iter, int NOUT);

class MonitorFunc: public Monitor{
public:
  MonitorFunc(MonitorFuncPointer pntr);
  virtual int call(Solver * solver, BoutReal time, int iter, int nout);
private:
  MonitorFuncPointer callFunc;
};
