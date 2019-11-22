#ifndef __BOUT_SCOREP_H__
#define __BOUT_SCOREP_H__

#include <bout_types.hxx>
#include "msg_stack.hxx"

#ifdef BOUT_HAS_SCOREP
#include <scorep/SCOREP_User.h>
#endif

#ifndef SCOREPLVL
#define SCOREPLVL 0
#endif

/// Instrument a region/function with scorep
///
/// The scorep call is identical for all levels, so just define it here.
/// If we don't have scorep support then just define a null function
#ifdef BOUT_HAS_SCOREP
#define SCOREP_BASE_CALL(...)						\
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

#endif
