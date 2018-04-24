#ifndef __DEPRECATED_H__
#define __DEPRECATED_H__

/// Mark functions for future removal
///
/// On gcc, expands to
///
///     func __attribute__ ((deprecated))
///
/// Example
/// -------
///
///     class SomeClass {
///      public:
///       DEPRECATED(int someFunction(const string &input));
///     }
#ifdef __GNUC__
#define DEPRECATED(func) __attribute__ ((deprecated)) func
#elif defined(_MSC_VER)
#define DEPRECATED(func) __declspec(deprecated) func
#else
//#pragma message("WARNING: You need to implement DEPRECATED for this compiler")
#define DEPRECATED(func) func
#endif

#endif // __DEPRECATED_H__
