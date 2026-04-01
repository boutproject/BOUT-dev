/// Support for i18n using GNU gettext

#ifndef BOUT_GETTEXT_H
#define BOUT_GETTEXT_H

#include "bout/build_defines.hxx"

#if BOUT_HAS_GETTEXT

#include <clocale> // IWYU pragma: keep

#include <libintl.h>

#define GETTEXT_PACKAGE "libbout"

// If we have C++23, we can get fmt to do compile-time checks of our format
// strings, _and_ have gettext do runtime replacement
#if __cpp_if_consteval >= 202106L
constexpr const char* dgettext_wrap(const char* __domainname, const char* __msgid) __THROW
    __attribute_format_arg__(2);

constexpr const char* dgettext_wrap(const char* __domainname, const char* __msgid) {
  if consteval {
    return __msgid;
  }
  return dgettext(__domainname, __msgid);
}

/// Gettext i18n macro for text containing fmt format specifiers
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define _f(string) dgettext_wrap(GETTEXT_PACKAGE, string)

#else
// We're pre-C++23, so all our i18n text must be fmt runtime formats
#include "fmt/base.h" // IWYU pragma: keep

/// Gettext i18n macro for text containing fmt format specifiers
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define _f(string) fmt::runtime(dgettext(GETTEXT_PACKAGE, string))
#endif

/// Gettext i18n macro for plain text that doesn't need formatting
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define _(string) dgettext(GETTEXT_PACKAGE, string)

#else // BOUT_HAS_GETTEXT

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define _f(string) string
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define _(string) string

#endif // BOUT_HAS_GETTEXT
#endif // BOUT_GETTEXT_H
