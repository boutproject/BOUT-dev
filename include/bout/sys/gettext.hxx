/// Support for i18n using GNU gettext

#ifndef __BOUT_GETTEXT_H__
#define __BOUT_GETTEXT_H__

#if BOUT_HAS_GETTEXT

#include <libintl.h>
#include <locale.h>

#define _(string) gettext(string)

#else

#define _(string) string

#endif // BOUT_HAS_GETTEXT
#endif // __BOUT_GETTEXT_H__
