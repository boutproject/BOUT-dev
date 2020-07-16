/// Support for i18n using GNU gettext

#ifndef __BOUT_GETTEXT_H__
#define __BOUT_GETTEXT_H__

#include "bout/build_config.hxx"

#if BOUT_HAS_GETTEXT

#include <clocale>
#include <libintl.h>

#define GETTEXT_PACKAGE "libbout"

#define _(string) dgettext(GETTEXT_PACKAGE, string)

#else

#define _(string) string

#endif // BOUT_HAS_GETTEXT
#endif // __BOUT_GETTEXT_H__
