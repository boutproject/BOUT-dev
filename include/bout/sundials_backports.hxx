// Backports for SUNDIALS compatibility between versions 4-7
//
// These are common backports shared between the CVode, ARKode, and IDA solvers
//
// Copyright 2022 Peter Hill, BOUT++ Team
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef BOUT_SUNDIALS_BACKPORTS_H
#define BOUT_SUNDIALS_BACKPORTS_H

#include "bout/bout_types.hxx"

#include <type_traits>

#include <nvector/nvector_parallel.h>
#include <sundials/sundials_config.h>
#include <sundials/sundials_iterative.h>
#include <sundials/sundials_nonlinearsolver.h>
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>

#ifndef SUNDIALS_VERSION
#error "Unable to determine SUNDIALS version"
#endif

#define SUNDIALS_VERSION_AT_LEAST(major, minor, patch) \
  (major < SUNDIALS_VERSION_MAJOR                      \
   || (major == SUNDIALS_VERSION_MAJOR                 \
       && (minor < SUNDIALS_VERSION_MINOR              \
           || (minor == SUNDIALS_VERSION_MINOR && patch <= SUNDIALS_VERSION_PATCH))))
#define SUNDIALS_VERSION_LESS_THAN(major, minor, patch) \
  (!SUNDIALS_VERSION_AT_LEAST(major, minor, patch))

#if SUNDIALS_VERSION_AT_LEAST(6, 0, 0)
#include <sundials/sundials_context.hpp>
#endif

#if SUNDIALS_VERSION_LESS_THAN(6, 0, 0)
using sundials_real_type = realtype;
#else
using sundials_real_type = sunrealtype;
#endif

static_assert(std::is_same_v<BoutReal, sundials_real_type>,
              "BOUT++ and SUNDIALS real types do not match");

#if SUNDIALS_VERSION_LESS_THAN(6, 0, 0)
constexpr auto SUN_PREC_RIGHT = PREC_RIGHT;
constexpr auto SUN_PREC_LEFT = PREC_LEFT;
constexpr auto SUN_PREC_NONE = PREC_NONE;

namespace sundials {
using Context = std::nullptr_t;
} // namespace sundials
#endif

inline sundials::Context createSUNContext([[maybe_unused]] MPI_Comm& comm) {
#if SUNDIALS_VERSION_LESS_THAN(6, 0, 0)
  return nullptr;
#elif SUNDIALS_VERSION_LESS_THAN(7, 0, 0)
  return sundials::Context(static_cast<void*>(&comm));
#else
  return sundials::Context(comm);
#endif
}

template <typename Func, typename... Args>
inline decltype(auto) callWithSUNContext(Func f, [[maybe_unused]] sundials::Context& ctx,
                                         Args&&... args) {
#if SUNDIALS_VERSION_LESS_THAN(6, 0, 0)
  return f(std::forward<Args>(args)...);
#else
  return f(std::forward<Args>(args)..., ctx);
#endif
}

#endif // BOUT_SUNDIALS_BACKPORTS_H
