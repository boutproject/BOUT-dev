/**************************************************************************
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

#include <field_factory.hxx>

#include <cmath>

#include <output.hxx>
#include <bout/constants.hxx>
#include <utils.hxx>

#include "bout/constants.hxx"

#include "fieldgenerators.hxx"

FieldGenerator* generator(BoutReal value) {
  return new FieldValue(value);
}

FieldGenerator* generator(BoutReal *ptr) {
  return new FieldValuePtr(ptr);
}

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory(Mesh *m) : fieldmesh(m) {
  
  // Useful values
  addGenerator("pi", new FieldValue(PI));

  // Some standard functions
  addGenerator("sin", new FieldSin(NULL));
  addGenerator("cos", new FieldCos(NULL));
  addGenerator("tan", new FieldGenOneArg<tan>(NULL));
  
  addGenerator("acos", new FieldGenOneArg<acos>(NULL));
  addGenerator("asin", new FieldGenOneArg<asin>(NULL));
  addGenerator("atan", new FieldATan(NULL));
  
  addGenerator("sinh", new FieldSinh(NULL));
  addGenerator("cosh", new FieldCosh(NULL));
  addGenerator("tanh", new FieldTanh());

  addGenerator("exp", new FieldGenOneArg<exp>(NULL));
  addGenerator("log", new FieldGenOneArg<log>(NULL));
  addGenerator("gauss", new FieldGaussian(NULL, NULL));
  addGenerator("abs", new FieldAbs(NULL));
  addGenerator("sqrt", new FieldSqrt(NULL));
  addGenerator("h", new FieldHeaviside(NULL));

  addGenerator("min", new FieldMin());
  addGenerator("max", new FieldMax());
  
  addGenerator("power", new FieldGenTwoArg<pow>(NULL,NULL));
}

FieldFactory::~FieldFactory() {
  
}

const Field2D FieldFactory::create2D(const string &value, Options *opt) {
  Field2D result = 0.;

  FieldGenerator* gen = parse(value, opt);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 2D field from '"
           << value
           << "'" << endl;
    return result;
  }
  
  for(int x=0;x<fieldmesh->ngx;x++)
    for(int y=0;y<fieldmesh->ngy;y++)
      result[x][y] = gen->generate(fieldmesh->GlobalX(x),
                                   TWOPI*fieldmesh->GlobalY(y),
                                   0.0,  // Z
                                   0.0); // T
  
  // Don't delete the generator, as will be cached

  return result;
}

const Field3D FieldFactory::create3D(const string &value, Options *opt) {
  Field3D result = 0.;

  FieldGenerator* gen = parse(value, opt);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 3D field from '"
           << value
           << "'" << endl;
    return result;
  }
  
  for(int x=0;x<fieldmesh->ngx;x++)
    for(int y=0;y<fieldmesh->ngy;y++)
      for(int z=0;z<fieldmesh->ngz;z++)
        result[x][y][z] = gen->generate(fieldmesh->GlobalX(x),
                                        TWOPI*fieldmesh->GlobalY(y),
                                        TWOPI*((BoutReal) z) / ((BoutReal) (fieldmesh->ngz-1)),  // Z
                                        0.0); // T
  // Don't delete generator

  return result;
}

FieldGenerator* FieldFactory::resolve(string &name) {
  if(options) {
    // Look up in options
    
    // Check if already looking up this symbol
    for(list<string>::const_iterator it=lookup.begin(); it != lookup.end(); it++)
      if( name.compare(*it) == 0 ) {
        // Name matches, so already looking up
        output << "ExpressionParser lookup stack:\n";
        for(list<string>::const_iterator it=lookup.begin(); it != lookup.end(); it++) {
          output << *it << " -> ";
        }
        output << name << endl;
        throw BoutException("ExpressionParser: Infinite recursion in parsing '%s'", name.c_str());
      }
    
    // Syntax for sections?
    if(options->isSet(name)) {
      // Add to lookup list
      lookup.push_back(name);
      
      // Get the string from options
      string val;
      options->get(name, val, "");
      
      // Parse
      FieldGenerator *g = parse(val, options);
      
      // Remove from lookup list
      lookup.pop_back();
      
      return g;
    }
  }
  output << "ExpressionParser error: Can't find generator '" << name << "'" << endl;
  return NULL;
}

FieldGenerator* FieldFactory::parse(const string &input, Options *opt) {

  // Check if in the cache
  
  string key = input;
  if(opt)
    key += opt->str(); // Include options in key
  
  map<string, FieldGenerator*>::iterator it = cache.find(key);
  if(it != cache.end()) {
    // Found in cache
    output << "Found '" << key << "' in cache\n";
    return it->second;
  }

  // Save the current options
  Options *oldoptions = options;

  // Store the options tree for token lookups
  options = opt;
  
  // Parse
  FieldGenerator *expr = parseString(input);
  
  output << "Adding '" << key << "' to cache\n";
  // Add to cache
  cache[key] = expr;

  // Restore the old options
  options = oldoptions;
  
  return expr;
}
