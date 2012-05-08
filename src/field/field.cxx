/**************************************************************************
 * Base class for fields
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

//#include <globals.hxx>

#include <stdarg.h>

#include <field.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <boutexception.hxx>

Field::Field() {
#ifdef CHECK
  bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true;
#endif
}

/////////////////// PROTECTED ////////////////////

char err_buffer[512];

// Report an error occurring
void Field::error(const char *s, ...) const {
  va_list ap;  // List of arguments

  if(s == (const char*) NULL) {
    output.write("Unspecified error in field\n");
  }else {
  
    va_start(ap, s);
      vsprintf(err_buffer, s, ap);
      va_end(ap);

#ifdef TRACK
      output.write("Error in '%s': %s", name.c_str(), err_buffer);
#else
      output.write("Error in field: %s", err_buffer);
#endif
  }
  
  throw BoutException("Error in field: %s", err_buffer);
}

