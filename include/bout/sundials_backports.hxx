// Backports for SUNDIALS compatibility between versions 3-6
//
// These are common backports shared between the CVode, ARKode, and IDA solvers
//
// Copyright 2022 Peter Hill, BOUT++ Team
// SPDX-License-Identifier: LGPLv3

#ifndef BOUT_SUNDIALS_BACKPORTS_H
#define BOUT_SUNDIALS_BACKPORTS_H

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_types.h>

#if SUNDIALS_VERSION_MAJOR >= 3
#include <sunlinsol/sunlinsol_spgmr.h>
#endif

#if SUNDIALS_VERSION_MAJOR >= 4
#include <sundials/sundials_nonlinearsolver.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#endif

#include "unused.hxx"

#if SUNDIALS_VERSION_MAJOR < 3
using SUNLinearSolver = int*;
inline void SUNLinSolFree(MAYBE_UNUSED(SUNLinearSolver solver)) {}
using sunindextype = long int;
#endif

#if SUNDIALS_VERSION_MAJOR < 4
using SUNNonlinearSolver = int*;
inline void SUNNonlinSolFree(MAYBE_UNUSED(SUNNonlinearSolver solver)) {}
#endif

#if SUNDIALS_VERSION_MAJOR < 6
namespace sundials {
struct Context {
  Context(MPI_Comm comm MAYBE_UNUSED() = MPI_COMM_NULL) {}
};
} // namespace sundials

using SUNContext = sundials::Context;

constexpr auto SUN_PREC_RIGHT = PREC_RIGHT;
constexpr auto SUN_PREC_LEFT = PREC_LEFT;
constexpr auto SUN_PREC_NONE = PREC_NONE;

inline N_Vector N_VNew_Parallel(MPI_Comm comm, sunindextype local_length,
                                sunindextype global_length,
                                MAYBE_UNUSED(SUNContext sunctx)) {
  return N_VNew_Parallel(comm, local_length, global_length);
}

#if SUNDIALS_VERSION_MAJOR >= 3
inline SUNLinearSolver SUNLinSol_SPGMR(N_Vector y, int pretype, int maxl,
                                       MAYBE_UNUSED(SUNContext sunctx)) {
#if SUNDIALS_VERSION_MAJOR == 3
  return SUNSPGMR(y, pretype, maxl);
#else
  return SUNLinSol_SPGMR(y, pretype, maxl);
#endif
}
#if SUNDIALS_VERSION_MAJOR >= 4
inline SUNNonlinearSolver SUNNonlinSol_FixedPoint(N_Vector y, int m,
                                                  MAYBE_UNUSED(SUNContext sunctx)) {
  return SUNNonlinSol_FixedPoint(y, m);
}

inline SUNNonlinearSolver SUNNonlinSol_Newton(N_Vector y,
                                              MAYBE_UNUSED(SUNContext sunctx)) {
  return SUNNonlinSol_Newton(y);
}
#endif // SUNDIALS_VERSION_MAJOR >= 4
#endif // SUNDIALS_VERSION_MAJOR >= 3
#endif // SUNDIALS_VERSION_MAJOR < 6

#endif // BOUT_SUNDIALS_BACKPORTS_H
