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
#include "field_factory.hxx"

#include <cmath>
#include <cstddef>
#include <list>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "bout/constants.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include "bout/sys/expressionparser.hxx"

#include "bout_types.hxx"
#include "boutexception.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "fieldgenerators.hxx"
#include "options.hxx"
#include "output.hxx"

/// Helper function to create a FieldValue generator from a BoutReal
FieldGeneratorPtr generator(BoutReal value) {
  return std::make_shared<FieldValue>(value);
}

/// Helper function to create a FieldValuePtr from a pointer to BoutReal
FieldGeneratorPtr generator(BoutReal *ptr) {
  return std::make_shared<FieldValuePtr>(ptr);
}

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory(Mesh * localmesh, Options *opt) : fieldmesh(localmesh), options(opt) {

  if (options == nullptr)
    options = Options::getRoot();

  // Useful values
  addGenerator("pi", std::make_shared<FieldValue>(PI));
  addGenerator("π", std::make_shared<FieldValue>(PI));
  
  // Some standard functions
  addGenerator("sin", std::make_shared<FieldSin>(nullptr));
  addGenerator("cos", std::make_shared<FieldCos>(nullptr));
  addGenerator("tan", std::make_shared<FieldGenOneArg<tan>>(nullptr));

  addGenerator("acos", std::make_shared<FieldGenOneArg<acos>>(nullptr));
  addGenerator("asin", std::make_shared<FieldGenOneArg<asin>>(nullptr));
  addGenerator("atan", std::make_shared<FieldATan>(nullptr));

  addGenerator("sinh", std::make_shared<FieldSinh>(nullptr));
  addGenerator("cosh", std::make_shared<FieldCosh>(nullptr));
  addGenerator("tanh", std::make_shared<FieldTanh>());

  addGenerator("exp", std::make_shared<FieldGenOneArg<exp>>(nullptr));
  addGenerator("log", std::make_shared<FieldGenOneArg<log>>(nullptr));
  addGenerator("gauss", std::make_shared<FieldGaussian>(nullptr, nullptr));
  addGenerator("abs", std::make_shared<FieldAbs>(nullptr));
  addGenerator("sqrt", std::make_shared<FieldSqrt>(nullptr));
  addGenerator("h", std::make_shared<FieldHeaviside>(nullptr));
  addGenerator("erf", std::make_shared<FieldErf>(nullptr));
  addGenerator("fmod", std::make_shared<FieldGenTwoArg<fmod>>(nullptr, nullptr));

  addGenerator("min", std::make_shared<FieldMin>());
  addGenerator("max", std::make_shared<FieldMax>());

  addGenerator("power", std::make_shared<FieldGenTwoArg<pow>>(nullptr, nullptr));

  addGenerator("round", std::make_shared<FieldRound>(nullptr));

  // Ballooning transform
  addGenerator("ballooning", std::make_shared<FieldBallooning>(fieldmesh));

  // Mixmode function
  addGenerator("mixmode", std::make_shared<FieldMixmode>());

  // TanhHat function
  addGenerator("tanhhat",
               std::make_shared<FieldTanhHat>(nullptr, nullptr, nullptr, nullptr));
}

FieldFactory::~FieldFactory() {

}

const Field2D FieldFactory::create2D(const std::string &value, const Options *opt,
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

  FieldGeneratorPtr gen = parse(value, opt);
  if(!gen) {
    output << "FieldFactory error: Couldn't create 2D field from '"
           << value
           << "'" << endl;
    return result;
  }

  switch(loc)  {
  case CELL_XLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal xpos = 0.5 * (localmesh->GlobalX(i.x() - 1) + localmesh->GlobalX(i.x()));
      result[i] = gen->generate(xpos, TWOPI * localmesh->GlobalY(i.y()),
                                0.0, // Z
                                t);  // T
    }
    break;
  }
  case CELL_YLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal ypos =
          TWOPI * 0.5 * (localmesh->GlobalY(i.y() - 1) + localmesh->GlobalY(i.y()));
      result[i] = gen->generate(localmesh->GlobalX(i.x()), ypos,
                                0.0, // Z
                                t);  // T
    }
    break;
  }
  default: {// CELL_CENTRE or CELL_ZLOW
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), TWOPI * localmesh->GlobalY(i.y()),
                        0.0, // Z
                        t);  // T
    }
  }
  };

  // Don't delete the generator, as will be cached

  return result;
}

const Field3D FieldFactory::create3D(const std::string &value, const Options *opt,
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
  FieldGeneratorPtr gen = parse(value, opt);
  if(!gen) {
    throw BoutException("FieldFactory error: Couldn't create 3D field from '%s'", value.c_str());
  }

  switch(loc)  {
  case CELL_XLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal xpos = 0.5 * (localmesh->GlobalX(i.x() - 1) + localmesh->GlobalX(i.x()));
      result[i] = gen->generate(xpos, TWOPI * localmesh->GlobalY(i.y()),
                                TWOPI * static_cast<BoutReal>(i.z()) /
                                    static_cast<BoutReal>(localmesh->LocalNz), // Z
                                t);                                            // T
    }
    break;
  }
  case CELL_YLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal ypos =
          TWOPI * 0.5 * (localmesh->GlobalY(i.y() - 1) + localmesh->GlobalY(i.y()));
      result[i] = gen->generate(localmesh->GlobalX(i.x()), ypos,
                                TWOPI * static_cast<BoutReal>(i.z()) /
                                    static_cast<BoutReal>(localmesh->LocalNz), // Z
                                t);                                            // T
    }
    break;
  }
  case CELL_ZLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), TWOPI * localmesh->GlobalY(i.y()),
                        TWOPI * (static_cast<BoutReal>(i.z()) - 0.5) /
                            static_cast<BoutReal>(localmesh->LocalNz), // Z
                        t);                                            // T
    }
    break;
  }
  default: {// CELL_CENTRE
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), TWOPI * localmesh->GlobalY(i.y()),
                        TWOPI * static_cast<BoutReal>(i.z()) /
                            static_cast<BoutReal>(localmesh->LocalNz), // Z
                        t);                                            // T
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

const Options* FieldFactory::findOption(const Options *opt, const std::string &name, std::string &val) {
  // Find an Options object which contains the given name

  const Options *result = opt;

  // Check if name contains a section separator ':'
  size_t pos = name.find(':');
  if(pos == std::string::npos) {
    // No separator. Try this section, and then go through parents

    while(!result->isSet(name)) {
      result = result->getParent();
      if (result == nullptr)
        throw ParseException("Cannot find variable '%s'", name.c_str());
    }
    result->get(name, val, "");

  }else {
    // Go to the root, and go up through sections
    result = Options::getRoot();

    size_t lastpos = 0;
    while(pos != std::string::npos) {
      std::string sectionname = name.substr(lastpos,pos);
      if( sectionname.length() > 0 ) {
        result = result->getSection(sectionname);
      }
      lastpos = pos+1;
      pos = name.find(':', lastpos);
    }
    // Now look for the name in this section

    std::string varname = name.substr(lastpos);

    if(!result->isSet(varname)) {
      // Not in this section
      throw ParseException("Cannot find variable '%s'", name.c_str());
    }

    result->get(varname, val, "");
  }

  return result;
}

FieldGeneratorPtr FieldFactory::resolve(std::string &name) {
  if (options) {
    // Check if in cache
    std::string key;
    if(name.find(':') != std::string::npos) {
      // Already has section
      key = name;
    }else {
      key = options->str();
      if(key.length() > 0)
        key += ":";
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
    std::string value;
    const Options *section = findOption(options, name, value);

    // Add to lookup list
    lookup.push_back(key);

    // Parse
    FieldGeneratorPtr g = parse(value, section);

    // Cache
    cache[key] = g;

    // Remove from lookup list
    lookup.pop_back();

    return g;
  }
  output << "ExpressionParser error: Can't find generator '" << name << "'" << endl;
  return nullptr;
}

FieldGeneratorPtr FieldFactory::parse(const std::string &input, const Options *opt) {

  // Check if in the cache
  std::string key = "#" + input;
  if (opt)
    key = opt->str() + key; // Include options context in key

  auto it = cache.find(key);
  if (it != cache.end()) {
    // Found in cache
    return it->second;
  }

  // Save the current options
  const Options *oldoptions = options;

  // Store the options tree for token lookups
  if (opt)
    options = opt;

  // Parse
  FieldGeneratorPtr expr = parseString(input);

  // Add to cache
  cache[key] = expr;

  // Restore the old options
  options = oldoptions;

  return expr;
}

FieldFactory* FieldFactory::get() {
  static FieldFactory instance(nullptr, Options::getRoot());

  return &instance;
}

void FieldFactory::cleanCache() {
  cache.clear();
}
