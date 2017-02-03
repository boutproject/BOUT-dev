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

#include <stdarg.h>
#include <stdio.h>
#include <string.h>
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

int Output::open(const char* fname, ...) {
  
  if(fname == (const char*) NULL)
    return 1;

  bout_vsnprintf(buffer, buffer_len, fname);

  close();

  file.open(buffer);

  if(!file.is_open()) {
    fprintf(stderr, "Could not open log file '%s'\n", buffer);
    return 1;
  }

  add(file);

  return 0;
}

void Output::close() {
  if(!file.is_open())
    return;
  
  remove(file);
  file.close();
}

void Output::write(const char* string, ...) {

  if(string == (const char*) NULL)
    return;
  
  bout_vsnprintf(buffer, buffer_len, string);

  multioutbuf_init::buf()->sputn(buffer, strlen(buffer));
}

void Output::print(const char* string, ...) {
  va_list ap;  // List of arguments

  if(!enabled)
    return; // Only output if to screen

  if(string == (const char*) NULL)
    return;
  
  va_start(ap, string);
    vprintf(string, ap);
  va_end(ap);
  
  fflush(stdout);
}


Output* Output::instance = NULL;

Output* Output::getInstance() {
  if(instance == NULL) {
    // Create the instance
    instance = new Output();
  }
  return instance;
}

void Output::cleanup() {
  if(instance == NULL)
    return;
  
  delete instance;
  instance = NULL;
}
