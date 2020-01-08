/**************************************************************************
 * Output, for printing messages/errors etc.
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

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <output.hxx>
#include <utils.hxx>

void Output::enable() {
  add(std::cout);
  enabled = true;
}

void Output::disable() {
  remove(std::cout);
  enabled = false;
}

int Output::open(const std::string& filename) {
  if (filename.empty()) {
    return 1;
  }

  close();

  file.open(filename);

  if (!file.is_open()) {
    fmt::print(stderr, "Could not open log file '{}'\n", filename);
    return 1;
  }

  add(file);

  return 0;
}

void Output::close() {
  if (!file.is_open()) {
    return;
  }

  remove(file);
  file.close();
}

void Output::write(const std::string& message) {
  multioutbuf_init::buf()->sputn(message.c_str(), message.length());
}

void Output::print(const std::string& message) {
  if (!enabled) {
    return; // Only output if to screen
  }
  std::cout << message;
  std::cout.flush();
}

Output *Output::getInstance() {
  static Output instance;
  return &instance;
}

void ConditionalOutput::write(const std::string& message) {
  if (enabled) {
    ASSERT1(base != nullptr);
    base->write(message);
  }
}

void ConditionalOutput::print(const std::string& message) {
  if (enabled) {
    ASSERT1(base != nullptr);
    base->print(message);
  }
}

#ifdef DEBUG_ENABLED
ConditionalOutput output_debug(Output::getInstance());
#else
DummyOutput output_debug;
#endif
ConditionalOutput output_warn(Output::getInstance());
ConditionalOutput output_info(Output::getInstance());
ConditionalOutput output_progress(Output::getInstance());
ConditionalOutput output_error(Output::getInstance());
ConditionalOutput output_verbose(Output::getInstance(), false);
ConditionalOutput output(Output::getInstance());

#undef bout_vsnprint_pre
