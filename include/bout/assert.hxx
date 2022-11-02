/*!
 * Defines a macro ASSERT which throws a BoutException if a given
 * condition is false. Whether the assertion is tested depends on
 * the checking level, so assetions can be removed for optimised runs.
 * 
 * ASSERT<level> ( condition )
 *
 *   level     - An integer known at compile time.
 *               condition tested if level >= CHECK
 *
 *   condition - The expression to test
 * 
 * e.g. ASSERT2( condition ) will only test condition if CHECK >= 2
 * 
 */

#ifndef __BOUT_ASSERT_H__
#define __BOUT_ASSERT_H__

#include "../boutexception.hxx"

#ifndef CHECK
#define CHECKLEVEL 0
#else
#define CHECKLEVEL CHECK
#endif

#if CHECKLEVEL >= 0
#define ASSERT0(condition)     \
  if(!(condition)) {           \
    throw BoutException("Assertion failed in {:s}, line {:d}: {:s}", __FILE__, __LINE__, #condition);  \
  }
#else // CHECKLEVEL >= 0
#define ASSERT0(condition)
#endif

#if CHECKLEVEL >= 1
#define ASSERT1(condition)     \
  if(!(condition)) {           \
    throw BoutException("Assertion failed in {:s}, line {:d}: {:s}", __FILE__, __LINE__, #condition);  \
  }
#else // CHECKLEVEL >= 1
#define ASSERT1(condition)
#endif

#if CHECKLEVEL >= 2
#define ASSERT2(condition)     \
  if(!(condition)) {           \
    throw BoutException("Assertion failed in {:s}, line {:d}: {:s}", __FILE__, __LINE__, #condition);  \
  }
#else // CHECKLEVEL >= 2
#define ASSERT2(condition)
#endif

#if CHECKLEVEL >= 3
#define ASSERT3(condition)     \
  if(!(condition)) {           \
    throw BoutException("Assertion failed in {:s}, line {:d}: {:s}", __FILE__, __LINE__, #condition);  \
  }
#else // CHECKLEVEL >= 3
#define ASSERT3(condition)
#endif

#endif // __BOUT_ASSERT_H__
