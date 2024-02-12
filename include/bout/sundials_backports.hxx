// Backports for SUNDIALS compatibility between versions 4-7
//
// These are common backports shared between the CVode, ARKode, and IDA solvers
//
// Copyright 2022 Peter Hill, BOUT++ Team
// SPDX-License-Identifier: LGPLv3

#ifndef BOUT_SUNDIALS_BACKPORTS_H
#define BOUT_SUNDIALS_BACKPORTS_H

#include <type_traits>

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#if SUNDIALS_VERSION_MAJOR >= 6
#include <sundials/sundials_context.hpp>
#endif

#include "bout/unused.hxx"

static_assert(std::is_same<BoutReal,
#if SUNDIALS_VERSION_MAJOR < 6
              realtype
#else
              sunrealtype
#endif
              >::value,
              "BOUT++ and SUNDIALS real types do not match");

#if SUNDIALS_VERSION_MAJOR < 6
constexpr auto SUN_PREC_RIGHT = PREC_RIGHT;
constexpr auto SUN_PREC_LEFT = PREC_LEFT;
constexpr auto SUN_PREC_NONE = PREC_NONE;

namespace sundials {
using Context = std::nullptr_t;
} // namespace sundials
#endif

inline sundials::Context createSUNContext(MAYBE_UNUSED(MPI_Comm& comm)) {
#if SUNDIALS_VERSION_MAJOR < 6
  return nullptr;
#elif SUNDIALS_VERSION_MAJOR < 7
  return sundials::Context(static_cast<void*>(&comm));
#else
  return sundials::Context(comm);
#endif
}

template<typename Func, typename... Args>
inline decltype(auto) callWithSUNContext(Func f,
                                        MAYBE_UNUSED(sundials::Context& ctx),
                                        Args&&... args) {
#if SUNDIALS_VERSION_MAJOR < 6
  return f(std::forward<Args>(args)...);
#else
  return f(std::forward<Args>(args)..., ctx);
#endif
}

#endif // BOUT_SUNDIALS_BACKPORTS_H
