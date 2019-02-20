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

#include "unused.hxx"

#include <string>
#include <map>
#include <list>

// Utility routines to create generators from values

FieldGeneratorPtr generator(BoutReal value);
FieldGeneratorPtr generator(BoutReal *ptr);

//////////////////////////////////////////////////////////
// Create a tree of generators from an input string

class FieldFactory : public ExpressionParser {
public:
  FieldFactory(Mesh* mesh = nullptr, Options* opt = nullptr);
  ~FieldFactory() override = default;

  /// Create a 2D field by parsing a string and evaluating the expression
  /// using the given options \p opt, over Mesh \p m at time \p t.
  /// The resulting field is at cell location \p loc.
  Field2D create2D(const std::string& value, const Options* opt = nullptr,
                   Mesh* m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0);

  /// Create a 3D field by parsing a string and evaluating the expression
  /// using the given options \p opt, over Mesh \p m at time \p t.
  /// The resulting field is at cell location \p loc.
  Field3D create3D(const std::string& value, const Options* opt = nullptr,
                   Mesh* m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0);

  /// Parse a string into a tree of generators
  FieldGeneratorPtr parse(const std::string& input, const Options* opt = nullptr);

  /// Create a 2D field from a generator, over a given mesh
  /// at a given cell location and time.
  Field2D create2D(FieldGeneratorPtr generator, Mesh* m = nullptr,
                   CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0);

  /// Create a 3D field from a generator, over a given mesh
  /// at a given cell location and time.
  Field3D create3D(FieldGeneratorPtr generator, Mesh* m = nullptr,
                   CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0);

  /// Get the Singleton object
  static FieldFactory* get();

  /// clean the cache of parsed strings
  void cleanCache();
protected:
  /// These functions called by the parser to resolve unknown symbols.
  /// This is used to enable options to be referred to in expressions.
  FieldGeneratorPtr resolve(std::string &name) override;

private:
  Mesh *fieldmesh;        ///< The default mesh for create functions.
  const Options *options; ///< Set in parse() and used in resolve()

  std::list<std::string> lookup; // Names currently being parsed
  
  /// Cache parsed strings so repeated evaluations
  /// don't result in allocating more generators.
  std::map<std::string, FieldGeneratorPtr > cache;
  
  const Options* findOption(const Options *opt, const std::string &name, std::string &val);
};

//////////////////////////////////////////////////////////
// Generator from function

class FieldFunction : public FieldGenerator {
public:
  FieldFunction() = delete;
  FieldFunction(FuncPtr userfunc) : func(userfunc) {}
  double generate(double x, double y, double z, double t) override {
    return func(t, x, y, z);
  }
private:
  FuncPtr func;
};

//////////////////////////////////////////////////////////
// Null generator 

class FieldNull : public FieldGenerator {
public:
  double generate(double UNUSED(x), double UNUSED(y), double UNUSED(z),
                  double UNUSED(t)) override {
    return 0.0;
  }
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> UNUSED(args)) override {
    return get();
  }
  /// Singeton
  static FieldGeneratorPtr get() {
    static FieldGeneratorPtr instance = nullptr;

    if(!instance)
      instance = std::make_shared<FieldNull>();
    return instance;
  }
private:
};

#endif // __FIELD_FACTORY_H__
