#ifndef __BOUT_SCOREP_H__
#define __BOUT_SCOREP_H__

#ifdef BOUT_HAS_SCOREP
#include <scorep/SCOREP_User.h>
#endif

#ifndef SCOREPLVL
#define SCOREPLVL 0
#endif

//The __PRETTY_FUNCTION__ variable is defined by GCC (and some other families) but is not a part 
//of the standard. The __func__ variable *is* a part of the c++11 standard so we'd like to fall back
//to this if possible. However as these are variables/constants and not macros we can't just
//check if __PRETTY_FUNCITON__ is defined or not. Instead we need to say if we support this
//or not by defining HAS_PRETTY_FUNCTION (to be implemented in configure)

#ifdef HAS_PRETTY_FUNCTION
#define __thefunc__ __PRETTY_FUNCTION__ 
#else
#define __thefunc__ __func__
#endif

//The scorep call is identical for all levels, so just define it here.
//If we don't have scorep support then just define a null function
#ifdef BOUT_HAS_SCOREP
#define SCOREP_BASE_CALL(...)						\
  SCOREP_USER_REGION(__thefunc__, SCOREP_USER_REGION_TYPE_FUNCTION)
#else
#define SCOREP_BASE_CALL(...)
#endif


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
