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

#include <bout/constants.hxx>
#include <output.hxx>
#include <utils.hxx>

#include "bout/constants.hxx"

#include "fieldgenerators.hxx"

/// Helper function to create a FieldValue generator from a BoutReal
FieldGeneratorPtr generator(BoutReal value) {
  return std::make_shared<FieldValue>(value);
}

/// Helper function to create a FieldValuePtr from a pointer to BoutReal
FieldGeneratorPtr generator(BoutReal* ptr) {
  return std::make_shared<FieldValuePtr>(ptr);
}

//////////////////////////////////////////////////////////
// FieldFactory public functions

FieldFactory::FieldFactory(Mesh* localmesh, Options* opt)
    : fieldmesh(localmesh == nullptr ? bout::globals::mesh : localmesh),
      options(opt == nullptr ? Options::getRoot() : opt) {

  // Set options
  // Note: don't use 'options' here because 'options' is a 'const Options*'
  // pointer, so this would fail if the "input" section is not present.
  Options& nonconst_options{opt == nullptr ? Options::root() : *opt};
  transform_from_field_aligned
    = nonconst_options["input"]["transform_from_field_aligned"].withDefault(true);

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

Field2D FieldFactory::create2D(const std::string& value, const Options* opt,
                               Mesh* localmesh, CELL_LOC loc, BoutReal t) const {
  return create2D(parse(value, opt), localmesh, loc, t);
}

Field2D FieldFactory::create2D(FieldGeneratorPtr gen, Mesh* localmesh, CELL_LOC loc,
                               BoutReal t) const {
  AUTO_TRACE();

  if (localmesh == nullptr) {
    if (fieldmesh == nullptr) {
      throw BoutException("FieldFactory not created with mesh and no mesh passed in");
    }
    localmesh = fieldmesh;
  }

  if (!gen) {
    throw BoutException("Couldn't create 2D field from null generator");
  }

  Field2D result{localmesh};

  result.allocate();
  result.setLocation(loc);

  constexpr BoutReal z_position{0.0};

  switch(loc)  {
  case CELL_XLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal xpos = 0.5 * (localmesh->GlobalX(i.x() - 1) + localmesh->GlobalX(i.x()));
      result[i] = gen->generate(xpos, TWOPI * localmesh->GlobalY(i.y()), z_position, t);
    }
    break;
  }
  case CELL_YLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal ypos =
          TWOPI * 0.5 * (localmesh->GlobalY(i.y() - 1) + localmesh->GlobalY(i.y()));
      result[i] = gen->generate(localmesh->GlobalX(i.x()), ypos, z_position, t);
    }
    break;
  }
  default: { // CELL_CENTRE or CELL_ZLOW
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] = gen->generate(localmesh->GlobalX(i.x()),
                                TWOPI * localmesh->GlobalY(i.y()), z_position, t);
    }
  }
  };

  return result;
}

Field3D FieldFactory::create3D(const std::string& value, const Options* opt,
                               Mesh* localmesh, CELL_LOC loc, BoutReal t) const {
  return create3D(parse(value, opt), localmesh, loc, t);
}

Field3D FieldFactory::create3D(FieldGeneratorPtr gen, Mesh* localmesh, CELL_LOC loc,
                               BoutReal t) const {
  AUTO_TRACE();

  if (localmesh == nullptr) {
    if (fieldmesh == nullptr) {
      throw BoutException("FieldFactory not created with mesh and no mesh passed in");
    }
    localmesh = fieldmesh;
  }

  if (!gen) {
    throw BoutException("Couldn't create 3D field from null generator");
  }

  Field3D result(localmesh);
  result.allocate();
  result.setLocation(loc);

  switch (loc) {
  case CELL_XLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal xpos = 0.5 * (localmesh->GlobalX(i.x() - 1) + localmesh->GlobalX(i.x()));
      result[i] = gen->generate(xpos, TWOPI * localmesh->GlobalY(i.y()),
                                TWOPI * static_cast<BoutReal>(i.z())
                                    / static_cast<BoutReal>(localmesh->LocalNz),
                                t);
    }
    break;
  }
  case CELL_YLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal ypos =
          TWOPI * 0.5 * (localmesh->GlobalY(i.y() - 1) + localmesh->GlobalY(i.y()));
      result[i] = gen->generate(localmesh->GlobalX(i.x()), ypos,
                                TWOPI * static_cast<BoutReal>(i.z())
                                    / static_cast<BoutReal>(localmesh->LocalNz),
                                t);
    }
    break;
  }
  case CELL_ZLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), TWOPI * localmesh->GlobalY(i.y()),
                        TWOPI * (static_cast<BoutReal>(i.z()) - 0.5)
                            / static_cast<BoutReal>(localmesh->LocalNz),
                        t);
    }
    break;
  }
  default: { // CELL_CENTRE
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), TWOPI * localmesh->GlobalY(i.y()),
                        TWOPI * static_cast<BoutReal>(i.z())
                            / static_cast<BoutReal>(localmesh->LocalNz),
                        t);
    }
  }
  };

  if (transform_from_field_aligned) {
    auto coords = result.getCoordinates();
    if (coords == nullptr) {
      throw BoutException("Unable to transform result: Mesh does not have Coordinates set");
    }
    if (coords->getParallelTransform().canToFromFieldAligned()) {
      // Transform from field aligned coordinates, to be compatible with
      // older BOUT++ inputs. This is not a particularly "nice" solution.
      result = fromFieldAligned(result, "RGN_ALL");
    }
  }

  return result;
}

FieldPerp FieldFactory::createPerp(const std::string& value, const Options* opt,
    Mesh* localmesh, CELL_LOC loc, BoutReal t) const {
  return createPerp(parse(value, opt), localmesh, loc, t);
}

FieldPerp FieldFactory::createPerp(FieldGeneratorPtr gen, Mesh* localmesh, CELL_LOC loc,
                                   BoutReal t) const {
  AUTO_TRACE();

  if (localmesh == nullptr) {
    if (fieldmesh == nullptr) {
      throw BoutException("FieldFactory not created with mesh and no mesh passed in");
    }
    localmesh = fieldmesh;
  }

  if (!gen) {
    throw BoutException("Couldn't create FieldPerp from null generator");
  }

  FieldPerp result(localmesh);
  result.allocate();
  result.setLocation(loc);

  switch (loc) {
  case CELL_XLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      BoutReal xpos = 0.5 * (localmesh->GlobalX(i.x() - 1) + localmesh->GlobalX(i.x()));
      result[i] = gen->generate(xpos, 0.,
                                TWOPI * static_cast<BoutReal>(i.z())
                                    / static_cast<BoutReal>(localmesh->LocalNz),
                                t);
    }
    break;
  }
  case CELL_YLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] = gen->generate(localmesh->GlobalX(i.x()), 0.,
                                TWOPI * static_cast<BoutReal>(i.z())
                                    / static_cast<BoutReal>(localmesh->LocalNz),
                                t);
    }
    break;
  }
  case CELL_ZLOW: {
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), 0.,
                        TWOPI * (static_cast<BoutReal>(i.z()) - 0.5)
                            / static_cast<BoutReal>(localmesh->LocalNz),
                        t);
    }
    break;
  }
  default: { // CELL_CENTRE
    BOUT_FOR(i, result.getRegion("RGN_ALL")) {
      result[i] =
          gen->generate(localmesh->GlobalX(i.x()), 0.,
                        TWOPI * static_cast<BoutReal>(i.z())
                            / static_cast<BoutReal>(localmesh->LocalNz),
                        t);
    }
  }
  };

  if (transform_from_field_aligned) {
    auto coords = result.getCoordinates();
    if (coords == nullptr) {
      throw BoutException("Unable to transform result: Mesh does not have Coordinates set");
    }
    if (coords->getParallelTransform().canToFromFieldAligned()) {
      // Transform from field aligned coordinates, to be compatible with
      // older BOUT++ inputs. This is not a particularly "nice" solution.
      result = fromFieldAligned(result, "RGN_ALL");
    }
  }

  return result;
}

const Options* FieldFactory::findOption(const Options* opt, const std::string& name,
                                        std::string& val) const {
  const Options* result = opt;

  // Check if name contains a section separator ':'
  size_t pos = name.find(':');
  if (pos == std::string::npos) {
    // No separator. Try this section, and then go through parents

    while (!result->isSet(name)) {
      result = result->getParent();
      if (result == nullptr)
        throw ParseException("Cannot find variable '%s'", name.c_str());
    }
    result->get(name, val, "");

  } else {
    // Go to the root, and go up through sections
    result = Options::getRoot();

    size_t lastpos = 0;
    while (pos != std::string::npos) {
      std::string sectionname = name.substr(lastpos, pos);
      if (sectionname.length() > 0) {
        result = result->getSection(sectionname);
      }
      lastpos = pos + 1;
      pos = name.find(':', lastpos);
    }
    // Now look for the name in this section

    std::string varname = name.substr(lastpos);

    if (!result->isSet(varname)) {
      // Not in this section
      throw ParseException("Cannot find variable '%s'", name.c_str());
    }

    result->get(varname, val, "");
  }

  return result;
}

FieldGeneratorPtr FieldFactory::resolve(std::string& name) const {
  if (options != nullptr) {
    // Check if in cache
    std::string key;
    if (name.find(':') != std::string::npos) {
      // Already has section
      key = name;
    } else {
      key = options->str();
      if (key.length() > 0) {
        key += ":";
      }
      key += name;
    }

    auto cached_value = cache.find(key);
    if (cached_value != cache.end()) {
      // Found in cache
      return cached_value->second;
    }

    // Look up in options

    // Check if already looking up this symbol
    for (const auto& lookup_value : lookup) {
      if (key.compare(lookup_value) == 0) {
        // Name matches, so already looking up
        output_error << "ExpressionParser lookup stack:\n";
        for (const auto& stack_value : lookup) {
          output_error << stack_value << " -> ";
        }
        output_error << name << endl;
        throw BoutException("ExpressionParser: Infinite recursion in parsing '%s'",
                            name.c_str());
      }
    }

    // Find the option, including traversing sections.
    // Throws exception if not found
    std::string value;
    const Options* section = findOption(options, name, value);

    lookup.push_back(key);

    FieldGeneratorPtr g = parse(value, section);

    cache[key] = g;

    lookup.pop_back();

    return g;
  }
  output << "ExpressionParser error: Can't find generator '" << name << "'" << endl;
  return nullptr;
}

FieldGeneratorPtr FieldFactory::parse(const std::string& input, const Options* opt) const {

  // Check if in the cache
  std::string key = "#" + input;
  if (opt != nullptr) {
    key = opt->str() + key; // Include options context in key
  }

  auto it = cache.find(key);
  if (it != cache.end()) {
    return it->second;
  }

  // Save the current options
  const Options* oldoptions = options;

  // Store the options tree for token lookups
  if (opt != nullptr) {
    options = opt;
  }

  FieldGeneratorPtr expr = parseString(input);

  cache[key] = expr;

  options = oldoptions;

  return expr;
}

FieldFactory* FieldFactory::get() {
  static FieldFactory instance(nullptr, Options::getRoot());

  return &instance;
}

void FieldFactory::cleanCache() { cache.clear(); }
