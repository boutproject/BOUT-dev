
#ifndef __MACRO_FOR_EACH_H__
#define __MACRO_FOR_EACH_H__

// Provides a macro MACRO_FOR_EACH which applies a
// macro to each argument in a VA_ARGS list
// 
// Based on this blog post by Daniel Hardman:
//    https://codecraft.co/2014/11/25/variadic-macros-tricks/
// This answer on StackOverflow:
//    https://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments/11763277
// See also: 
//    https://github.com/cormacc/va_args_iterators
//


/// _me_x set of macros expand a number of arguments without ';' between them
#define _me_1(_call, x) _call(x)
#define _me_2(_call, x, ...) _call(x) _me_1(_call, __VA_ARGS__)
#define _me_3(_call, x, ...) _call(x) _me_2(_call, __VA_ARGS__)
#define _me_4(_call, x, ...) _call(x) _me_3(_call, __VA_ARGS__)
#define _me_5(_call, x, ...) _call(x) _me_4(_call, __VA_ARGS__)
#define _me_6(_call, x, ...) _call(x) _me_5(_call, __VA_ARGS__)
#define _me_7(_call, x, ...) _call(x) _me_6(_call, __VA_ARGS__)
#define _me_8(_call, x, ...) _call(x) _me_7(_call, __VA_ARGS__)
#define _me_9(_call, x, ...) _call(x) _me_8(_call, __VA_ARGS__)
#define _me_10(_call, x, ...) _call(x) _me_9(_call, __VA_ARGS__)

/// _fe_x set of macros expand a number of arguments with ';' between them
#define _fe_1(_call, x) _call(x);
#define _fe_2(_call, x, ...) _call(x); _fe_1(_call, __VA_ARGS__)
#define _fe_3(_call, x, ...) _call(x); _fe_2(_call, __VA_ARGS__)
#define _fe_4(_call, x, ...) _call(x); _fe_3(_call, __VA_ARGS__)
#define _fe_5(_call, x, ...) _call(x); _fe_4(_call, __VA_ARGS__)
#define _fe_6(_call, x, ...) _call(x); _fe_5(_call, __VA_ARGS__)
#define _fe_7(_call, x, ...) _call(x); _fe_6(_call, __VA_ARGS__)
#define _fe_8(_call, x, ...) _call(x); _fe_7(_call, __VA_ARGS__)
#define _fe_9(_call, x, ...) _call(x); _fe_8(_call, __VA_ARGS__)
#define _fe_10(_call, x, ...) _call(x); _fe_9(_call, __VA_ARGS__)

/// When called with __VA_ARGS__ first, this evaluates to an argument which depends
/// on the length of __VA_ARGS__. This is used to find the appropriate macro to
/// begin the expansion.
#define _GET_FOR_EACH_EXPANSION(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, NAME, ...) NAME

/// Apply a macro (first argument) to each
/// of the following arguments.
/// Currently supports up to 10 arguments.
///
/// Example:
///
///    MACRO_FOR_EACH(test, a, b, c)
///
/// expands to
///
///    test(a) test(b) test(c)
///
/// Notes:
///
///   - No semicolon is inserted after each expansion
///   - No braces are put around the expansion. These
///     should usually be added in the top-level macro
///     to avoid surprising results.
/// 
#define MACRO_FOR_EACH(mac, ...)                                        \
  _GET_FOR_EACH_EXPANSION(__VA_ARGS__,                                  \
                          _me_10, _me_9, _me_8, _me_7, _me_6, _me_5,    \
                          _me_4, _me_3, _me_2, _me_1)                   \
  (mac, __VA_ARGS__)

/// Apply a function (first argument) to each
/// of the following arguments.
/// Currently supports up to 10 arguments.
///
/// Example:
///
///    MACRO_FOR_EACH_FN(test, a, b, c)
///
/// expands to
///
///    test(a); test(b); test(c);
///
/// Notes:
///
///   - A ; is inserted after each expansion
///   - No braces are put around the expansion. These
///     should usually be added in the top-level macro
///     to avoid surprising results.
///
#define MACRO_FOR_EACH_FN(fn, ...)                                      \
  _GET_FOR_EACH_EXPANSION(__VA_ARGS__,                                  \
                          _fe_10, _fe_9, _fe_8, _fe_7, _fe_6, _fe_5,    \
                          _fe_4, _fe_3, _fe_2, _fe_1)                   \
  (fn, __VA_ARGS__)

#endif
