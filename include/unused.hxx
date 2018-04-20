#ifndef __UNUSED_H__
#define __UNUSED_H__

/// Mark a function parameter as unused in the function body
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
#ifdef UNUSED
#elif defined(__GNUC__)
# define UNUSED(x) UNUSED_ ## x __attribute__((unused))
#elif defined(__LCLINT__)
# define UNUSED(x) /*@unused@*/ x
#elif defined(__cplusplus)
# define UNUSED(x)
#else
# define UNUSED(x) x
#endif

#endif //__UNUSED_H__
