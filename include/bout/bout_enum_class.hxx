/**************************************************************************
 * Copyright 2019 B.D.Dudson, J.T.Omotani, P.Hill
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
 **************************************************************************/

#ifndef __BOUT_ENUM_CLASS_H__
#define __BOUT_ENUM_CLASS_H__

#include "bout/macro_for_each.hxx"
#include "boutexception.hxx"
#include "msg_stack.hxx"
#include "options.hxx"

#include <map>
#include <string>

/// Create some macro magic similar to bout/macro_for_each.hxx, but allowing for the enum
/// class name to be passed through to each _call
/// _ec_expand_x set of macros expand a number of arguments without ';' between them
#define _ec_expand_1(_call, enumname, x) _call(enumname, x)
#define _ec_expand_2(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_1(_call, enumname, __VA_ARGS__)
#define _ec_expand_3(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_2(_call, enumname, __VA_ARGS__)
#define _ec_expand_4(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_3(_call, enumname, __VA_ARGS__)
#define _ec_expand_5(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_4(_call, enumname, __VA_ARGS__)
#define _ec_expand_6(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_5(_call, enumname, __VA_ARGS__)
#define _ec_expand_7(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_6(_call, enumname, __VA_ARGS__)
#define _ec_expand_8(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_7(_call, enumname, __VA_ARGS__)
#define _ec_expand_9(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_8(_call, enumname, __VA_ARGS__)
#define _ec_expand_10(_call, enumname, x, ...) \
  _call(enumname, x) _ec_expand_9(_call, enumname, __VA_ARGS__)

#define BOUT_ENUM_CLASS_MAP_ARGS(mac, enumname, ...)                                    \
  BOUT_EXPAND(_GET_FOR_EACH_EXPANSION(                                                  \
    __VA_ARGS__, _ec_expand_10, _ec_expand_9, _ec_expand_8, _ec_expand_7, _ec_expand_6, \
    _ec_expand_5, _ec_expand_4, _ec_expand_3, _ec_expand_2, _ec_expand_1)               \
  (mac, enumname, __VA_ARGS__))

#define BOUT_ENUM_CLASS_STR(enumname, val) {enumname::val, lowercase(#val)},
#define BOUT_STR_ENUM_CLASS(enumname, val) {lowercase(#val), enumname::val},

#define BOUT_MAKE_FROMSTRING_NAME(enumname) enumname ## FromString

/// Create an enum class with toString and <enum name>FromString functions, and an
/// Options::as<enum> overload to read the enum
#define BOUT_ENUM_CLASS(enumname, ...)                                      \
enum class enumname { __VA_ARGS__ };                                        \
                                                                            \
inline std::string toString(enumname e) {                                   \
  AUTO_TRACE();                                                             \
  const static std::map<enumname, std::string> toString_map = {             \
    BOUT_ENUM_CLASS_MAP_ARGS(BOUT_ENUM_CLASS_STR, enumname, __VA_ARGS__)    \
  };                                                                        \
  auto found = toString_map.find(e);                                        \
  if (found == toString_map.end()) {                                        \
    throw BoutException("Did not find enum {:d}", static_cast<int>(e));     \
  }                                                                         \
  return found->second;                                                     \
}                                                                           \
                                                                            \
inline enumname BOUT_MAKE_FROMSTRING_NAME(enumname)(const std::string& s) { \
  AUTO_TRACE();                                                             \
  const static std::map<std::string, enumname> fromString_map = {           \
    BOUT_ENUM_CLASS_MAP_ARGS(BOUT_STR_ENUM_CLASS, enumname, __VA_ARGS__)    \
  };                                                                        \
  auto found = fromString_map.find(s);                                      \
  if (found == fromString_map.end()) {                                      \
    throw BoutException("Did not find enum {:s}", s);                       \
  }                                                                         \
  return found->second;                                                     \
}                                                                           \
                                                                            \
template <> inline enumname Options::as<enumname>(const enumname&) const {  \
  return BOUT_MAKE_FROMSTRING_NAME(enumname)(this->as<std::string>());      \
}                                                                           \
                                                                            \
inline std::ostream& operator<<(std::ostream& out, const enumname& e) {     \
  return out << toString(e);                                                \
}

#endif // __BOUT_ENUM_CLASS_H__
