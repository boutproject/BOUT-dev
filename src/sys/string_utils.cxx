/**************************************************************************
 * String utility routines
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

#include "boutexception.hxx"
#include "string_utils.hxx"

#include <algorithm>
#include <string.h>
#include <stdio.h>

// Allocate memory for a copy of given string
char *copy_string(const char *s) {
  char *s2;
  int n;

  if (s == NULL)
    return NULL;

  n = strlen(s);
  s2 = (char *)malloc(n + 1);
  strcpy(s2, s);
  return s2;
}

/// Convert a string to lower case
const std::string lowercase(const std::string &str) {
  std::string strlow(str);

  std::transform(strlow.begin(), strlow.end(), strlow.begin(), ::tolower);
  return strlow;
}

/// Convert to lowercase, except for inside strings
const std::string lowercasequote(const std::string &str) {
  std::string strlow(str);

  bool quote = false, dquote = false;
  for (std::string::size_type i = 0; i < strlow.length(); i++) {
    if (strlow[i] == '\'') {
      quote ^= true;
    } else if (strlow[i] == '"') {
      dquote ^= true;
    } else if ((!quote) && (!dquote)) {
      strlow[i] = tolower(strlow[i]);
    }
  }
  return strlow;
}

BoutReal stringToReal(const std::string &s) {
  std::stringstream ss(s);
  BoutReal val;
  if (!(ss >> val)) {
    throw BoutException("Could not convert string '%s' to BoutReal\n", s.c_str());
  }
  return val;
}

int stringToInt(const std::string &s) {
  std::stringstream ss(s);
  int val;
  if (!(ss >> val)) {
    throw BoutException("Could not convert string '%s' to int\n", s.c_str());
  }
  return val;
}

std::list<std::string> &strsplit(const std::string &s, char delim,
                                 std::list<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::list<std::string> strsplit(const std::string &s, char delim) {
  std::list<std::string> elems;
  return strsplit(s, delim, elems);
}

// Strips leading and trailing spaces from a string
std::string trim(const std::string &s, const std::string &c) {
  return trimLeft(trimRight(s, c), c);
}

std::string trimRight(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.erase(s.find_last_not_of(c) + 1);
}

std::string trimLeft(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.erase(0, s.find_first_not_of(c));
}

// Strips the comments from a string
// This is the compliment of trimLeft
std::string trimComments(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.substr(0, s.find_first_of(c));
}
