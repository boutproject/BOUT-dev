/**************************************************************************
 * Memory allocation routines
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

#include <utils.hxx>
#include <cstring>
#include <cstdlib>
#include <cctype>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cmath>
#include <ctime>
#include <iomanip>

#include "fmt/chrono.h"

/**************************************************************************
 * String routines
 **************************************************************************/

// Allocate memory for a copy of given string
char* copy_string(const char* s) {

  if (s == nullptr)
    return nullptr;

  const auto n = strlen(s);
  auto s2 = static_cast<char*>(malloc(n + 1));
  strcpy(s2, s);
  return s2;
}

// Convert a string to lower case
const std::string lowercase(const std::string &str) {
  std::string strlow(str);

  std::transform(strlow.begin(), strlow.end(), strlow.begin(), ::tolower);
  return strlow;
}

// Convert a string to upper case
const std::string uppercase(const std::string& str) {
  std::string strup(str);

  std::transform(strup.begin(), strup.end(), strup.begin(), ::toupper);
  return strup;
}

// Convert to lowercase, except for inside strings
const std::string lowercasequote(const std::string &str) {
  std::string strlow(str);

  bool quote = false, dquote = false;
  for (char &i : strlow) {
    if (i == '\'') {
      quote ^= true;
    } else if (i == '"') {
      dquote ^= true;
    } else if ((!quote) && (!dquote)) {
      i = static_cast<char>(tolower(i));
    }
  }
  return strlow;
}

BoutReal stringToReal(const std::string &s) {
  std::stringstream ss(s);
  BoutReal val;
  if(!(ss >> val)) {
    throw BoutException("Could not convert string '{:s}' to BoutReal\n", s);
  }
  return val;
}

int stringToInt(const std::string &s) {
  std::stringstream ss(s);
  int val;
  if(!(ss >> val)) {
    throw BoutException("Could not convert string '{:s}' to int\n", s);
  }
  return val;
}

std::list<std::string> &strsplit(const std::string &s, char delim, std::list<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
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
  return str.erase(s.find_last_not_of(c)+1);
}

std::string trimLeft(const std::string &s, const std::string &c) {
  std::string str(s);
  return str.erase(0, s.find_first_not_of(c));
}

// Strips the comments from a string
// This is the compliment of trimLeft
std::string trimComments(const std::string &s, const std::string &c) {
  return s.substr(0, s.find_first_of(c));
}

std::string toString(const time_t& time) {
  return fmt::format("{:%c}", *std::localtime(&time));
}

std::string::size_type editDistance(const std::string& str1, const std::string& str2) {

  using str_size_t = std::string::size_type;

  const auto str1_size = str1.size() + 1;
  const auto str2_size = str2.size() + 1;

  auto distance = Matrix<str_size_t>(str1_size, str2_size);

  // Initialise zeroth column and row with string index
  for (str_size_t i = 0; i < str1_size; ++i) {
    distance(i, 0) = i;
  }
  for (str_size_t j = 0; j < str2_size; ++j) {
    distance(0, j) = j;
  }

  // Wikipedia uses 1-indexing for the input strings, but 0-indexing
  // for the `d` matrix, so the input strings have an additional `-1`
  // when indexing them
  for (str_size_t i = 1; i < str1_size; ++i) {
    for (str_size_t j = 1; j < str2_size; ++j) {
      const str_size_t cost = (str1[i - 1] == str2[j - 1]) ? 0 : 1;

      distance(i, j) = std::min({
          distance(i - 1, j) + 1,       // deletion
          distance(i, j - 1) + 1,       // insertion
          distance(i - 1, j - 1) + cost // substitution
      });

      if (i > 1 and j > 1 and (str1[i - 1] == str2[j - 2])
          and (str1[i - 2] == str2[j - 1])) {
        // transposition
        distance(i, j) = std::min(distance(i, j), distance(i - 2, j - 2) + 1);
      }
    }
  }

  return distance(str1_size - 1, str2_size - 1);
}
