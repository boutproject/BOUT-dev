// Integration check against a patched SUNDIALS 7.9.0 shared-library build.
// tests/nvtx3/nvToolsExt.h records NVTX calls without an NVIDIA dependency.
#include <bout/profiling.hxx>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cvode/cvode.h>
#include <cvode/cvode_ls.h>
#include <map>
#include <nvector/nvector_serial.h>
#include <string>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <vector>
std::vector<std::string> stack;
std::map<std::string, int> counts;
int nvtxRangePushEx(const nvtxEventAttributes_t* event) {
  assert(event->color == bout::profiling::depthColor(stack.size()));
  stack.emplace_back(event->message.ascii);
  ++counts[stack.back()];
  return stack.size() - 1;
}
int nvtxRangePop() {
  assert(!stack.empty());
  stack.pop_back();
  return stack.size();
}
extern "C" void bout_nvtx_range_push(const char* name) {
  bout::profiling::nvtxPushColor(
      name, bout::profiling::depthColor(bout::profiling::rangeDepth()));
}
extern "C" void bout_nvtx_range_pop() { bout::profiling::nvtxPop(); }
int rhs(sunrealtype, N_Vector y, N_Vector dy, void*) {
  bout::profiling::ProfileRange range("Model::rhs");
  NV_Ith_S(dy, 0) = -100 * NV_Ith_S(y, 0) + NV_Ith_S(y, 1);
  NV_Ith_S(dy, 1) = -NV_Ith_S(y, 1);
  return 0;
}
extern "C" int SUNBoutNVTXInstrumented(void);
int main() {
  assert(SUNBoutNVTXInstrumented());
  SUNContext ctx;
  assert(SUNContext_Create(SUN_COMM_NULL, &ctx) == 0);
  N_Vector y = N_VNew_Serial(2, ctx);
  N_VConst(1, y);
  void* cv = CVodeCreate(CV_BDF, ctx);
  assert(CVodeInit(cv, rhs, 0, y) == 0);
  assert(CVodeSStolerances(cv, 1e-6, 1e-10) == 0);
  SUNLinearSolver ls = SUNLinSol_SPGMR(y, SUN_PREC_NONE, 5, ctx);
  assert(CVodeSetLinearSolver(cv, ls, nullptr) == 0);
  sunrealtype t;
  {
    bout::profiling::ProfileRange run("CvodeSolver::run");
    bout::profiling::ProfileRange integrate("CvodeSolver::integrate");
    assert(CVode(cv, 1, y, &t, CV_NORMAL) == CV_SUCCESS);
    assert(std::abs(NV_Ith_S(y, 1) - std::exp(-1.)) < 1e-5);
  }
  assert(stack.empty());
  assert(bout::profiling::rangeDepth() == 0);
  for (const char* name : {"CVode", "CVStep", "StepAttempt", "CVNls", "NewtonSolve",
                           "NewtonIteration", "CVLsSolve", "SPGMRSolve", "GMRESIteration",
                           "CVLsATimes", "CVLsDQJtimes", "N_VLinearSum", "N_VScale",
                           "N_VDotProd", "N_VWrmsNorm", "N_VLinearCombination"}) {
    assert(counts[std::string("SUNDIALS::") + name] > 0);
    std::printf("%s: %d ranges\n", name, counts[std::string("SUNDIALS::") + name]);
  }
  // A failed integration must close all scopes, too.
  assert(CVodeSetMaxNumSteps(cv, 1) == 0);
  assert(CVode(cv, 100, y, &t, CV_NORMAL) == CV_TOO_MUCH_WORK);
  assert(stack.empty());
  assert(bout::profiling::rangeDepth() == 0);
  CVodeFree(&cv);
  SUNLinSolFree(ls);
  N_VDestroy(y);
  SUNContext_Free(&ctx);
  puts("Passed: real CVODE/Newton/GMRES ranges, depth colors, solution, and failure "
       "cleanup.");
}
