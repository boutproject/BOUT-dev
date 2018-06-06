#include <bout/solver.hxx>

MonitorFunc::MonitorFunc(MonitorFuncPointer pntr){
  callFunc=pntr;
}


int MonitorFunc::call(Solver * solver, BoutReal time, int iter, int nout){
  return (callFunc)(solver,time,iter,nout);
}
