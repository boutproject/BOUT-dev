#ifndef __BOUT_FORMAT_H__
#define __BOUT_FORMAT_H__

/// Tell GCC that a function has a printf-style like argument
/// The first argument is the position of format string, and the
/// second is the position of the first variadic argument
/// Note that it seems to start counting from 1, and also counts a
/// *this pointer, as the first argument, so often 2 would be the
/// first argument.
#if defined(__GNUC__)
# define BOUT_FORMAT_ARGS(i,j) __attribute__ ((format (printf, i, j)))
#else
# define BOUT_FORMAT_ARGS(i,j)
#endif

#endif //__BOUT_FORMAT_H__
