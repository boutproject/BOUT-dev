/// Support for i18n using GNU gettext

#ifndef BOUT_GETTEXT_H
#define BOUT_GETTEXT_H

#include "bout/build_defines.hxx"

#if BOUT_HAS_GETTEXT

#include <clocale>
#include <libintl.h>

#define GETTEXT_PACKAGE "libbout"

#define _(string) dgettext(GETTEXT_PACKAGE, string)

#else

#define _(string) string

#endif // BOUT_HAS_GETTEXT
#endif // BOUT_GETTEXT_H
