/*!*************************************************************************
 * \file string_utils.hxx
 *
 * A mix of short utilities for strings
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#ifndef __STRING_UTILS_H__
#define __STRING_UTILS_H__

#include "bout_types.hxx"

#include <string>
#include <list>
#include <memory>
#include <sstream>


/// Allocate memory and copy string \p s
char *copy_string(const char *s);

/// Convert a value to a string by writing to a stringstream
template <class T> const std::string toString(const T &val) {
  std::stringstream ss;
  ss << val;
  return ss.str();
}

/// Convert a string to lower case
const std::string lowercase(const std::string &str);

/// Convert to lower case, except inside quotes (" or ')
const std::string lowercasequote(const std::string &str);

/// Convert a string to a BoutReal
///
/// Throws BoutException if can't be done
///
/// @param[in] s String containing digits
///
/// @return BoutReal representation of string
BoutReal stringToReal(const std::string &s);

/// Convert a string to an int
///
/// Throws BoutException if can't be done
///
/// @param[in] s String containing digits
///
/// @return Integer representation of string
int stringToInt(const std::string &s);

/// Split a string on a given delimiter
///
/// @param[in] s     The string to split (not modified by call)
/// @param[in] delim The delimiter to split on (single char)
/// @param[in, out] elems  A list to which the pieces will be appended using push_back
///
/// @return List of substrings
std::list<std::string> &strsplit(const std::string &s, char delim,
                                 std::list<std::string> &elems);

/// Split a string on a given delimiter
///
/// @param[in] s     The string to split (not modified by call)
/// @param[in] delim The delimiter to split on (single char)
///
/// @return List of substrings
std::list<std::string> strsplit(const std::string &s, char delim);

/// Strips leading and trailing spaces from a string
///
/// @param[in] s   The string to trim (not modified)
/// @param[in] c   Collection of characters to remove
///
/// @return Trimmed string
std::string trim(const std::string &s, const std::string &c=" \t\r");

/// Strips leading spaces from a string
///
/// @param[in] s   The string to trim (not modified)
/// @param[in] c   Collection of characters to remove
///
/// @return Trimmed string
std::string trimLeft(const std::string &s, const std::string &c=" \t");

/// Strips leading spaces from a string
///
/// @param[in] s   The string to trim (not modified)
/// @param[in] c   Collection of characters to remove
///
/// @return Trimmed string
std::string trimRight(const std::string &s, const std::string &c=" \t\r");

/// Strips the comments from a string
///
/// Removes anything after the first appearance of one
/// of the characters in \p c
///
/// @param[in] s   The string to trim (not modified)
/// @param[in] c   Collection of comment characters
///
/// @return Trimmed string
std::string trimComments(const std::string &s, const std::string &c="#;");

/// Format a string using C printf-style formatting
///
/// Can be used to replace C-style variadic functions that use the
/// bout_vsnprintf macro:
///
/// `include/datafile.hxx`:
///
///     bool openr(const char *filename, ...);
///
/// `src/fileio/datafile.cxx`:
///
///     bool Datafile::openr(const char *format, ...) {
///       if(format == (const char*) NULL)
///         throw BoutException("Datafile::open: No argument given for opening file!");
///       bout_vsnprintf(filename,filenamelen, format);
///
/// which can be replaced with:
///
/// `include/datafile.hxx`:
///
///     template <typename... Args>
///     bool openr(const std::string &filename, Args... args) {
///       return openr(string_format(filename, args...));
///     }
///     bool openr(const std::string &filename);
///
/// `filename` can then be used directly in `Datafile::openr`.
///
/// Warning: \p args should not contain anything of type std::string, as
/// the printf functions cannot handle them; if you need to call string_format
/// with a std::string, make sure to pass `foo.c_str()` instead
///
/// Taken from http://stackoverflow.com/a/26221725/2043465
///
/// @param format printf-style format string
/// @param args Arguments to format
///
/// @return Formatted string
template <typename... Args>
std::string string_format(const std::string &format, Args... args) {
  // Extra space for '\0'
  size_t size = snprintf(nullptr, 0, format.c_str(), args...) + 1;
  std::unique_ptr<char[]> buffer(new char[size]);
  snprintf(buffer.get(), size, format.c_str(), args...);
  // We don't want the '\0' inside
  return std::string(buffer.get(), buffer.get() + size - 1);
}

/// the bout_vsnprintf macro:
/// The first argument is an char * buffer of length len.
/// It needs to have been allocated with new[], as it may be
/// reallocated.
/// len: the length of said buffer. May be changed, mussn't be const.
/// fmt: the const char * descriping the format.
/// note that fmt should be the first argument of the function of type
/// const char * and has to be directly followed by the variable arguments.
#define bout_vsnprintf(buf, len, fmt)                                                    \
  {                                                                                      \
    va_list va;                                                                          \
    va_start(va, fmt);                                                                   \
    int _vsnprintflen = vsnprintf(buf, len, fmt, va);                                    \
    va_end(va);                                                                          \
    if (_vsnprintflen + 1 > len) {                                                       \
      _vsnprintflen += 1;                                                                \
      delete[] buf;                                                                      \
      buf = new char[_vsnprintflen];                                                     \
      len = _vsnprintflen;                                                               \
      va_start(va, fmt);                                                                 \
      vsnprintf(buf, len, fmt, va);                                                      \
      va_end(va);                                                                        \
    }                                                                                    \
  }

#endif // __STRING_UTILS_H__
