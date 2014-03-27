/**************************************************************************
 * Generate a field with specified values, mainly for creating
 * initial perturbations
 * 
 * 
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


class FieldFactory;

#ifndef __FIELD_FACTORY_H__
#define __FIELD_FACTORY_H__

#include "bout/mesh.hxx"

#include "bout/sys/expressionparser.hxx"

#include "field2d.hxx"
#include "field3d.hxx"
#include "options.hxx"

#include <string>
#include <map>
#include <list>

// Utility routines to create generators from values

FieldGenerator* generator(BoutReal value);
FieldGenerator* generator(BoutReal *ptr);

//////////////////////////////////////////////////////////
// Create a tree of generators from an input string

class FieldFactory : public ExpressionParser {
public:
  FieldFactory(Mesh *m);
  ~FieldFactory();
  
  const Field2D create2D(const std::string &value, Options *opt = NULL);
  const Field3D create3D(const std::string &value, Options *opt = NULL);
  
protected:
  // These functions called by the parser
  FieldGenerator* resolve(std::string &name);
  
private:
  Mesh *fieldmesh;
  
  Options *options;
  std::list<std::string> lookup; // Names currently being parsed
  
  // Cache parsed strings
  std::map<std::string, FieldGenerator* > cache;
  
  FieldGenerator* parse(const std::string &input, Options *opt=NULL);
};

#endif // __FIELD_FACTORY_H__
