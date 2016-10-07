#include <bout/solver.hxx>

// int Monitor::call(Solver * solver, BoutReal time, int iter, int nout){
//   throw BoutException("Monitor::call - not implemented");
// }

// bool Monitor::operator==(const Monitor& rhs) const{
//   return this == &rhs;
// }


//MonitorFunc::MonitorFunc(int (&pntr)(Solver *solver, BoutReal simtime, int iter, int NOUT)){
MonitorFunc::MonitorFunc(MonitorFuncPointer pntr){
  callFunc=pntr;
}


int MonitorFunc::call(Solver * solver, BoutReal time, int iter, int nout){
  return (callFunc)(solver,time,iter,nout);
}
