#include <bout/solver.hxx>
#include <bout/assert.hxx>
#include <cstdlib>
#include <cmath>

MonitorFunc::MonitorFunc(MonitorFuncPointer pntr){
  callFunc=pntr;
}


int MonitorFunc::call(Solver * solver, BoutReal time, int iter, int nout){
  return (callFunc)(solver,time,iter,nout);
}



bool isMultiple(BoutReal a, BoutReal b){
  ASSERT2(a>0);
  ASSERT2(b>0);
  auto min = a>b?b:a;
  auto max = a>b?a:b;
  auto ratio = std::round(max/min);
  auto error = ratio*min - max;
  if (std::abs(error/max) > 1e-12) {
    return false;
  }
  return true;
}
