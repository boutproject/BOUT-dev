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

#include <list>
#include <map>
#include <string>

// Utility routines to create generators from values

FieldGeneratorPtr generator(BoutReal value);
FieldGeneratorPtr generator(BoutReal* ptr);

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
                   Mesh* m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0) const;

  /// Create a 3D field by parsing a string and evaluating the expression
  /// using the given options \p opt, over Mesh \p m at time \p t.
  /// The resulting field is at cell location \p loc.
  Field3D create3D(const std::string& value, const Options* opt = nullptr,
                   Mesh* m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0) const;

  /// Create a FieldPerp by parsing a string and evaluating the expression
  /// using the given options \p opt, over Mesh \p m at time \p t.
  /// The resulting field is at cell location \p loc.
  FieldPerp createPerp(const std::string& value, const Options* opt = nullptr,
      Mesh* m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0) const;

  /// Parse a string into a tree of generators
  FieldGeneratorPtr parse(const std::string& input, const Options* opt = nullptr) const;

  /// Create a 2D field from a generator, over a given mesh
  /// at a given cell location and time.
  Field2D create2D(FieldGeneratorPtr generator, Mesh* m = nullptr,
                   CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0) const;

  /// Create a 3D field from a generator, over a given mesh
  /// at a given cell location and time.
  Field3D create3D(FieldGeneratorPtr generator, Mesh* m = nullptr,
                   CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0) const;

  /// Create a FieldPerp from a generator, over a given mesh
  /// at a given cell location and time.
  FieldPerp createPerp(FieldGeneratorPtr generator, Mesh* m = nullptr,
      CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0) const;

  /// Get the Singleton object
  static FieldFactory* get();

  /// clean the cache of parsed strings
  void cleanCache();

protected:
  /// These functions called by the parser to resolve unknown symbols.
  /// This is used to enable options to be referred to in expressions.
  FieldGeneratorPtr resolve(std::string& name) const override;

private:
  /// The default mesh for create functions.
  Mesh* fieldmesh;

  /// Should we transform input from field-aligned coordinates (if possible)?
  bool transform_from_field_aligned{true};

  /// The default options used in resolve(), can be *temporarily*
  /// overridden in parse()/create2D()/create3D()
  mutable const Options* options;

  /// Names currently being parsed
  mutable std::list<std::string> lookup;

  /// Cache parsed strings so repeated evaluations
  /// don't result in allocating more generators.
  mutable std::map<std::string, FieldGeneratorPtr> cache;

  /// Find an Options object which contains the given \p name
  const Options* findOption(const Options* opt, const std::string& name,
                            std::string& val) const;
};

//////////////////////////////////////////////////////////
// Generator from function

class FieldFunction : public FieldGenerator {
public:
  FieldFunction() = delete;
  FieldFunction(FuncPtr userfunc) : func(userfunc) {}
  BoutReal generate(const bout::generator::Context& pos) override {
    return func(pos.t(), pos.x(), pos.y(), pos.z());
  }

private:
  FuncPtr func;
};

//////////////////////////////////////////////////////////
// Null generator

class FieldNull : public FieldGenerator {
public:
  FieldNull() = default;
  BoutReal generate(const bout::generator::Context&) override {
    return 0.0;
  }
  FieldGeneratorPtr clone(const std::list<FieldGeneratorPtr> UNUSED(args)) override {
    return get();
  }
  /// Singeton
  static FieldGeneratorPtr get() {
    static FieldGeneratorPtr instance = std::make_shared<FieldNull>();
    return instance;
  }
};

#endif // __FIELD_FACTORY_H__
