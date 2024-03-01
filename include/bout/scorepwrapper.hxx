#ifndef BOUT_SCOREP_H
#define BOUT_SCOREP_H

#include "bout/build_config.hxx"

#include "bout/msg_stack.hxx"
#include <bout/bout_types.hxx>

#if BOUT_HAS_SCOREP
#include <scorep/SCOREP_User.h>
#endif

#ifndef SCOREPLVL
#define SCOREPLVL 0
#endif

/// Instrument a function with scorep
///
/// The scorep call is identical for all levels, so just define it here.
/// If we don't have scorep support then just define a null function
#if BOUT_HAS_SCOREP
#define SCOREP_BASE_CALL(...) \
  SCOREP_USER_REGION(__thefunc__, SCOREP_USER_REGION_TYPE_FUNCTION)
#else
#define SCOREP_BASE_CALL(...)
#endif

/// This is always defined
#if SCOREPLVL >= 0
#define SCOREP0(...) SCOREP_BASE_CALL(__VA_ARGS__)
#else
#define SCOREP0(...)
#endif

#if SCOREPLVL >= 1
#define SCOREP1(...) SCOREP_BASE_CALL(__VA_ARGS__)
#else
#define SCOREP1(...)
#endif

#if SCOREPLVL >= 2
#define SCOREP2(...) SCOREP_BASE_CALL(__VA_ARGS__)
#else
#define SCOREP2(...)
#endif

#if SCOREPLVL >= 3
#define SCOREP3(...) SCOREP_BASE_CALL(__VA_ARGS__)
#else
#define SCOREP3(...)
#endif

/// Instrument a region with scorep
#if BOUT_HAS_SCOREP
#define BOUT_SCOREP_REGION(...) \
  SCOREP_USER_REGION(__VA_ARGS__, SCOREP_USER_REGION_TYPE_COMMON)
#else
#define BOUT_SCOREP_REGION(...)
#endif

#endif
