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
std::shared_ptr<FieldGenerator> generator(BoutReal value) {
  return std::shared_ptr<FieldGenerator>( new FieldValue(value));
}

/// Helper function to create a FieldValuePtr from a pointer to BoutReal
std::shared_ptr<FieldGenerator> generator(BoutReal *ptr) {
  return std::shared_ptr<FieldGenerator>( new FieldValuePtr(ptr));
}

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory(Mesh * localmesh, Options *opt) : fieldmesh(localmesh), options(opt) {

  if(options == NULL)
    options = Options::getRoot();

  // Useful values
  addGenerator("pi", std::shared_ptr<FieldGenerator>( new FieldValue(PI)));

  // Some standard functions
  addGenerator("sin", std::shared_ptr<FieldGenerator>( new FieldSin(NULL)));
  addGenerator("cos", std::shared_ptr<FieldGenerator>( new FieldCos(NULL)));
  addGenerator("tan", std::shared_ptr<FieldGenerator>( new FieldGenOneArg<tan>(NULL)));

  addGenerator("acos", std::shared_ptr<FieldGenerator>( new FieldGenOneArg<acos>(NULL)));
  addGenerator("asin", std::shared_ptr<FieldGenerator>( new FieldGenOneArg<asin>(NULL)));
  addGenerator("atan", std::shared_ptr<FieldGenerator>( new FieldATan(NULL)));

  addGenerator("sinh", std::shared_ptr<FieldGenerator>( new FieldSinh(NULL)));
  addGenerator("cosh", std::shared_ptr<FieldGenerator>( new FieldCosh(NULL)));
  addGenerator("tanh", std::shared_ptr<FieldGenerator>( new FieldTanh()));

  addGenerator("exp", std::shared_ptr<FieldGenerator>( new FieldGenOneArg<exp>(NULL)));
  addGenerator("log", std::shared_ptr<FieldGenerator>( new FieldGenOneArg<log>(NULL)));
  addGenerator("gauss", std::shared_ptr<FieldGenerator>( new FieldGaussian(NULL, NULL)));
  addGenerator("abs", std::shared_ptr<FieldGenerator>( new FieldAbs(NULL)));
  addGenerator("sqrt", std::shared_ptr<FieldGenerator>( new FieldSqrt(NULL)));
  addGenerator("h", std::shared_ptr<FieldGenerator>( new FieldHeaviside(NULL)));
  addGenerator("erf", std::shared_ptr<FieldGenerator>( new FieldErf(NULL)));

  addGenerator("min", std::shared_ptr<FieldGenerator>( new FieldMin()));
  addGenerator("max", std::shared_ptr<FieldGenerator>( new FieldMax()));

  addGenerator("power", std::shared_ptr<FieldGenerator>( new FieldGenTwoArg<pow>(NULL,NULL)));

  addGenerator("round", std::shared_ptr<FieldGenerator>( new FieldRound(NULL)));
  
  // Ballooning transform
  addGenerator("ballooning", std::shared_ptr<FieldGenerator>( new FieldBallooning(fieldmesh)));

  // Mixmode function
  addGenerator("mixmode", std::shared_ptr<FieldGenerator>( new FieldMixmode()));

  // TanhHat function
  addGenerator("tanhhat", std::shared_ptr<FieldGenerator>( new FieldTanhHat(NULL, NULL, NULL, NULL)));

  // Real X and Real Y
  addGenerator("realx", std::shared_ptr<FieldGenerator>( new FieldRealX()));

  addGenerator("realy", std::shared_ptr<FieldGenerator>( new FieldRealY()));
}

FieldFactory::~FieldFactory() {
}

const Field2D FieldFactory::create2D(const string &value, Options *opt,
                                     Mesh *localmesh, CELL_LOC loc,
                                     BoutReal t) {

  if(localmesh == nullptr)
    localmesh = fieldmesh;
  if(localmesh == nullptr)
    throw BoutException("Not a valid mesh");

  Field2D result(0.,localmesh);

  if(localmesh->StaggerGrids == false){
    loc = CELL_CENTRE ;
  }
  result.setLocation(loc);

  std::shared_ptr<FieldGenerator> gen = parse(value, opt);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 2D field from '"
           << value
           << "'" << endl;
    return result;
  }

  switch(loc)  {
  case CELL_XLOW: {
    for(auto i : result) {
      BoutReal xpos = 0.5*(localmesh->GlobalX(i.x-1) + localmesh->GlobalX(i.x));
      result[i] = gen->generate(xpos,
                                TWOPI*localmesh->GlobalY(i.y),
                                0.0,  // Z
                                t, // T
                                i,localmesh);
    }
    break;
  }
  case CELL_YLOW: {
    for(auto i : result) {
      BoutReal ypos = TWOPI*0.5*(localmesh->GlobalY(i.y-1) + localmesh->GlobalY(i.y));
      result[i] = gen->generate(localmesh->GlobalX(i.x),
                                ypos,
                                0.0,  // Z
                                t, // T
                                i,localmesh);
    }
    break;
  }
  default: {// CELL_CENTRE or CELL_ZLOW
    for(auto i : result) {
      result[i] = gen->generate(localmesh->GlobalX(i.x),
                                TWOPI*localmesh->GlobalY(i.y),
                                0.0,  // Z
                                t, // T
                                i,localmesh);
    }
  }
  };

  // Don't delete the generator, as will be cached

  return result;
}

const Field3D FieldFactory::create3D(const string &value, Options *opt,
                                     Mesh *localmesh, CELL_LOC loc,
                                     BoutReal t) {

  if(localmesh == nullptr)
    localmesh = fieldmesh;
  if(localmesh == nullptr)
    throw BoutException("Not a valid mesh");

  // Create a Field3D over mesh "localmesh"
  Field3D result(localmesh);
  
  // Ensure that data is allocated and unique
  result.allocate();

  result.setLocation(loc);

  // Parse expression to create a tree of generators
  std::shared_ptr<FieldGenerator> gen = parse(value, opt);
  if(!gen) {
    throw BoutException("FieldFactory error: Couldn't create 3D field from '%s'", value.c_str());
  }

  switch(loc)  {
  case CELL_XLOW: {
    for(auto i : result) {
      BoutReal xpos = 0.5*(localmesh->GlobalX(i.x-1) + localmesh->GlobalX(i.x));
      result[i] = gen->generate(xpos,
                                TWOPI*localmesh->GlobalY(i.y),
                                TWOPI*static_cast<BoutReal>(i.z) / static_cast<BoutReal>(localmesh->LocalNz),  // Z
                                t, // T
                                i,localmesh);
    }
    break;
  }
  case CELL_YLOW: {
    for(auto i : result) {
      BoutReal ypos = TWOPI*0.5*(localmesh->GlobalY(i.y-1) + localmesh->GlobalY(i.y));
      result[i] = gen->generate(localmesh->GlobalX(i.x),
                                ypos,
                                TWOPI*static_cast<BoutReal>(i.z) / static_cast<BoutReal>(localmesh->LocalNz),  // Z
                                t, // T
                                i,localmesh);
    }
    break;
  }
  case CELL_ZLOW: {
    for(auto i : result) {
      result[i] = gen->generate(localmesh->GlobalX(i.x),
                                TWOPI*localmesh->GlobalY(i.y),
                                TWOPI*(static_cast<BoutReal>(i.z) - 0.5) / static_cast<BoutReal>(localmesh->LocalNz),  // Z
                                t, // T
                                i,localmesh);
    }
    break;
  }
  default: {// CELL_CENTRE
    for(auto i : result) {
      result[i] = gen->generate(localmesh->GlobalX(i.x),
                                TWOPI*localmesh->GlobalY(i.y),
                                TWOPI*static_cast<BoutReal>(i.z) / static_cast<BoutReal>(localmesh->LocalNz),  // Z
                                t, // T
                                i,localmesh);
    }
  }
  };

  // Don't delete generator
  
  if (localmesh->canToFromFieldAligned()){ // Ask wheter it is possible
    // Transform from field aligned coordinates, to be compatible with
    // older BOUT++ inputs. This is not a particularly "nice" solution.
    result = localmesh->fromFieldAligned(result);
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

std::shared_ptr<FieldGenerator> FieldFactory::resolve(string &name) {
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

    auto cached_value = cache.find(key);
    if (cached_value != cache.end()) {
      // Found in cache
      return cached_value->second;
    }

    // Look up in options

    // Check if already looking up this symbol
    for (const auto &lookup_value : lookup) {
      if (key.compare(lookup_value) == 0) {
        // Name matches, so already looking up
        output << "ExpressionParser lookup stack:\n";
        for (const auto &stack_value : lookup) {
          output << stack_value << " -> ";
        }
        output << name << endl;
        throw BoutException("ExpressionParser: Infinite recursion in parsing '%s'",
                            name.c_str());
      }
    }

    // Find the option, including traversing sections.
    // Throws exception if not found
    string value;
    Options *section = findOption(options, name, value);

    // Add to lookup list
    lookup.push_back(key);

    // Parse
    std::shared_ptr<FieldGenerator> g = parse(value, section);

    // Cache
    cache[key] = g;

    // Remove from lookup list
    lookup.pop_back();

    return g;
  }
  output << "ExpressionParser error: Can't find generator '" << name << "'" << endl;
  return NULL;
}

std::shared_ptr<FieldGenerator> FieldFactory::parse(const string &input, Options *opt) {

  // Check if in the cache
  string key = string("#") + input;
  if(opt)
    key = opt->str()+key; // Include options context in key

  auto it = cache.find(key);
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
  std::shared_ptr<FieldGenerator> expr = parseString(input);

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

void FieldFactory::cleanCache() {
  cache.clear();
}
