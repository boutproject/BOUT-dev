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
#include <utils.hxx>
#include <bout/mesh.hxx>

Field::Field() : fieldmesh(nullptr), is_const(false){
#if CHECK > 0
  bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true;
#endif
}

Field::Field(Mesh * msh) : fieldmesh(msh), is_const(false){
  if (fieldmesh ==nullptr){
    fieldmesh=mesh;
  }
#if CHECK > 0
  bndry_xin = bndry_xout = bndry_yup = bndry_ydown = true;
#endif
}

int Field::getNx() const{
  return getMesh()->LocalNx;
};

int Field::getNy() const{
  return getMesh()->LocalNy;
};

int Field::getNz() const{
  return getMesh()->LocalNz;
};

/////////////////// PROTECTED ////////////////////


// Report an error occurring
void Field::error(const char *s, ...) const {
  int buf_len=512;
  char * err_buffer=new char[buf_len];

  if(s == (const char*) NULL) {
    output_error.write("Unspecified error in field\n");
  }else {
  
    bout_vsnprintf(err_buffer,buf_len, s);

#ifdef TRACK
      output_error.write("Error in '%s': %s", name.c_str(), err_buffer);
#else
      output_error.write("Error in field: %s", err_buffer);
#endif
  }
  std::string msg="Error in field: ";
  msg+=err_buffer;
  delete[] err_buffer;
  throw BoutException(msg);
}

void Field::makeAutoConstant(){
  BoutReal val = (*this)[Indices{0,0,0}];
  for (auto i : *this){
    if (val != (*this)[i]){
      is_const=false;
      output_info.write("difference!");
      return;
    }
  }
  output_info.write("field is const\n");
  is_const=true;
}
