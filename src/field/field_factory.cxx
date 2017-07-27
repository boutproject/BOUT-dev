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
#include <globals.hxx>

#include <field_factory.hxx>

#include <cmath>

#include <output.hxx>
#include <bout/constants.hxx>
#include <utils.hxx>

#include "bout/constants.hxx"

#include "fieldgenerators.hxx"

/// Helper function to create a FieldValue generator from a BoutReal
FieldGenerator* generator(BoutReal value) {
  return new FieldValue(value);
}

/// Helper function to create a FieldValuePtr from a pointer to BoutReal
FieldGenerator* generator(BoutReal *ptr) {
  return new FieldValuePtr(ptr);
}

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory(Mesh *m, Options *opt) : fieldmesh(m), options(opt) {

  if(options == NULL)
    options = Options::getRoot();

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
  addGenerator("erf", new FieldErf(NULL));

  addGenerator("min", new FieldMin());
  addGenerator("max", new FieldMax());

  addGenerator("power", new FieldGenTwoArg<pow>(NULL,NULL));

  addGenerator("round", new FieldRound(NULL));
  
  // Ballooning transform
  addGenerator("ballooning", new FieldBallooning(fieldmesh));

  // Mixmode function
  addGenerator("mixmode", new FieldMixmode());

  // TanhHat function
  addGenerator("tanhhat", new FieldTanhHat(NULL, NULL, NULL, NULL));
}

FieldFactory::~FieldFactory() {

}

const Field2D FieldFactory::create2D(const string &value, Options *opt, Mesh *m, CELL_LOC loc, BoutReal t) {
  Field2D result = 0.;

  if(mesh->StaggerGrids == false){
    loc = CELL_CENTRE ;
  }
  result.setLocation(loc);

  if(m == NULL)
    m = fieldmesh;
  if(m == NULL)
    throw BoutException("Not a valid mesh");

  FieldGenerator* gen = parse(value, opt);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 2D field from '"
           << value
           << "'" << endl;
    return result;
  }

  switch(loc)  {
  case CELL_XLOW: {
    for(auto i : result) {
      BoutReal xpos = 0.5*(m->GlobalX(i.x-1) + m->GlobalX(i.x));
      result[i] = gen->generate(xpos,
                                TWOPI*m->GlobalY(i.y),
                                0.0,  // Z
                                t); // T
    }
    break;
  }
  case CELL_YLOW: {
    for(auto i : result) {
      BoutReal ypos = TWOPI*0.5*(m->GlobalY(i.y-1) + m->GlobalY(i.y));
      result[i] = gen->generate(m->GlobalX(i.x),
                                ypos,
                                0.0,  // Z
                                t); // T
    }
    break;
  }
  default: {// CELL_CENTRE or CELL_ZLOW
    for(auto i : result) {
      result[i] = gen->generate(m->GlobalX(i.x),
                                TWOPI*m->GlobalY(i.y),
                                0.0,  // Z
                                t); // T
    }
  }
  };

  // Don't delete the generator, as will be cached

  return result;
}

const Field3D FieldFactory::create3D(const string &value, Options *opt, Mesh *m, CELL_LOC loc, BoutReal t) {
  
  if(m == NULL)
    m = fieldmesh;
  if(m == NULL)
    throw BoutException("Not a valid mesh");
  
  // Create a Field3D over mesh "m"
  Field3D result(m);
  
  // Ensure that data is allocated and unique
  result.allocate();
  
  if(mesh->StaggerGrids == false){
    loc = CELL_CENTRE ;
  }
  result.setLocation(loc);

  // Parse expression to create a tree of generators
  FieldGenerator* gen = parse(value, opt);
  if(!gen) {
    throw BoutException("FieldFactory error: Couldn't create 3D field from '%s'", value.c_str());
  }

  switch(loc)  {
  case CELL_XLOW: {
    for(auto i : result) {
      BoutReal xpos = 0.5*(m->GlobalX(i.x-1) + m->GlobalX(i.x));
      result[i] = gen->generate(xpos,
                                TWOPI*m->GlobalY(i.y),
                                TWOPI*((BoutReal) i.z) / ((BoutReal) (m->LocalNz)),  // Z
                                t); // T
    }
    break;
  }
  case CELL_YLOW: {
    for(auto i : result) {
      BoutReal ypos = TWOPI*0.5*(m->GlobalY(i.y-1) + m->GlobalY(i.y));
      result[i] = gen->generate(m->GlobalX(i.x),
                                ypos,
                                TWOPI*((BoutReal) i.z) / ((BoutReal) (m->LocalNz)),  // Z
                                t); // T
    }
    break;
  }
  case CELL_ZLOW: {
    for(auto i : result) {
      result[i] = gen->generate(m->GlobalX(i.x),
                                TWOPI*m->GlobalY(i.y),
                                TWOPI*(((BoutReal) i.z) - 0.5) / ((BoutReal) (m->LocalNz)),  // Z
                                t); // T
    }
    break;
  }
  default: {// CELL_CENTRE
    for(auto i : result) {
      result[i] = gen->generate(m->GlobalX(i.x),
                                TWOPI*m->GlobalY(i.y),
                                TWOPI*((BoutReal) i.z) / ((BoutReal) (m->LocalNz)),  // Z
                                t); // T
    }
  }
  };

  // Don't delete generator
  
  // Transform from field aligned coordinates, to be compatible with
  // older BOUT++ inputs. This is not a particularly "nice" solution.
  try {
    result = m->fromFieldAligned(result);
  }catch(BoutException &e) {
    // might fail if not possible to shift coordinates
    // e.g. FCI
  }

  return result;
}

Options* FieldFactory::findOption(Options *opt, const string &name, string &val) {
  // Find an Options object which contains the given name

  Options *result = opt;

  // Check if name contains a section separator ':'
  size_t pos = name.find(':');
  if(pos == string::npos) {
    // No separator. Try this section, and then go through parents

    while(!result->isSet(name)) {
      result = result->getParent();
      if(result == NULL)
        throw ParseException("Cannot find variable '%s'", name.c_str());
    }
    result->get(name, val, "");

  }else {
    // Go to the root, and go up through sections
    result = Options::getRoot();

    size_t lastpos = 0;
    while(pos != string::npos) {
      string sectionname = name.substr(lastpos,pos);
      if( sectionname.length() > 0 ) {
        result = result->getSection(sectionname);
      }
      lastpos = pos+1;
      pos = name.find(':', lastpos);
    }
    // Now look for the name in this section

    string varname = name.substr(lastpos);

    if(!result->isSet(varname)) {
      // Not in this section
      throw ParseException("Cannot find variable '%s'", name.c_str());
    }

    result->get(varname, val, "");
  }

  return result;
}

FieldGenerator* FieldFactory::resolve(string &name) {
  if(options) {
    // Check if in cache
    string key;
    if(name.find(':') != string::npos) {
      // Already has section
      key = name;
    }else {
      key = options->str();
      if(key.length() > 0)
        key += string(":");
      key += name;
    }

    map<string, FieldGenerator*>::iterator it = cache.find(key);
    if(it != cache.end()) {
      // Found in cache
      return it->second;
    }

    // Look up in options

    // Check if already looking up this symbol
    for(const auto& it : lookup) {
      if( key.compare(it) == 0 ) {
        // Name matches, so already looking up
        output << "ExpressionParser lookup stack:\n";
        for(const auto& it : lookup) {
          output << it << " -> ";
        }
        output << name << endl;
        throw BoutException("ExpressionParser: Infinite recursion in parsing '%s'", name.c_str());
      }
    }
    
    // Find the option, including traversing sections.
    // Throws exception if not found
    string value;
    Options *section = findOption(options, name, value);

    // Add to lookup list
    lookup.push_back(key);

    // Parse
    FieldGenerator *g = parse(value, section);

    // Cache
    cache[key] = g;

    // Remove from lookup list
    lookup.pop_back();

    return g;
  }
  output << "ExpressionParser error: Can't find generator '" << name << "'" << endl;
  return NULL;
}

FieldGenerator* FieldFactory::parse(const string &input, Options *opt) {

  //output.write("FieldFactory::parse('%s')", input.c_str());

  // Check if in the cache

  string key = string("#") + input;
  if(opt)
    key = opt->str()+key; // Include options context in key

  map<string, FieldGenerator*>::iterator it = cache.find(key);
  if(it != cache.end()) {
    // Found in cache
    return it->second;
  }

  // Save the current options
  Options *oldoptions = options;

  // Store the options tree for token lookups
  if(opt)
    options = opt;

  // Parse
  FieldGenerator *expr = parseString(input);

  // Add to cache
  cache[key] = expr;

  // Restore the old options
  options = oldoptions;

  return expr;
}

FieldFactory* FieldFactory::get() {
  static FieldFactory instance(NULL, Options::getRoot());

  return &instance;
}
