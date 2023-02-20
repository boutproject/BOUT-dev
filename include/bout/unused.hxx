#ifndef __UNUSED_H__
#define __UNUSED_H__

/// Mark a function parameter as unused in the function body
///
/// For GCC, expands to
///
///     UNUSED_x __attribute__((unused))
///
/// Macro taken from http://stackoverflow.com/q/7090998/2043465
///
/// This will add the "unused" attribute to parameters in function
/// signatures, telling the compiler that we know the parameter isn't
/// used. This should cut down on false positives when using
/// -Wunused-parameters.
///
/// Additionally, this macro will also rename the
/// parameter so that if it is accidentally used, the compiler will
/// throw an error.
///
/// A better way to do this might be to detect how to silence the
/// warning in configure and use that in the macro instead.
///
/// Example
/// -------
///
///     void someFunction(int UNUSED(x)) {};
#if defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(_MSC_VER)
#define UNUSED(x) __pragma(warning(suppress : 4100)) UNUSED_ ## x
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#elif defined(__cplusplus)
# define UNUSED(x)
#else
# define UNUSED(x) x
#endif

/// Mark a function parameter as possibly unused in the function body
///
/// Unlike `UNUSED`, this has to go around the type as well:
///
///    MAYBE_UNUSED(int foo);
#ifdef __has_cpp_attribute
#if __has_cpp_attribute(maybe_unused)
# define MAYBE_UNUSED(x) [[maybe_unused]] x
#endif
#endif
#ifndef MAYBE_UNUSED
#if defined(__GNUC__)
# define MAYBE_UNUSED(x) [[gnu::unused]] x
#elif defined(_MSC_VER)
# define MAYBE_UNUSED(x) __pragma(warning(suppress : 4100)) x
#else
# define MAYBE_UNUSED(x) x
#endif
#endif

#endif //__UNUSED_H__
