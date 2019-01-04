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

#ifndef __FIELD_FACTORY_H__
#define __FIELD_FACTORY_H__

#include <string>
#include <map>
#include <list>

#include "bout/sys/expressionparser.hxx"
#include "bout_types.hxx"

class Field2D;
class Field3D;
class Mesh;
class Options;

// Utility routines to create generators from values

FieldGeneratorPtr generator(BoutReal value);
FieldGeneratorPtr generator(BoutReal *ptr);

//////////////////////////////////////////////////////////
// Create a tree of generators from an input string

class FieldFactory : public ExpressionParser {
public:
  FieldFactory(Mesh *m, Options *opt = nullptr);
  ~FieldFactory() override;

  const Field2D create2D(const std::string &value, const Options *opt = nullptr,
                         Mesh *m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0);
  const Field3D create3D(const std::string &value, const Options *opt = nullptr,
                         Mesh *m = nullptr, CELL_LOC loc = CELL_CENTRE, BoutReal t = 0.0);

  // Parse a string into a tree of generators
  FieldGeneratorPtr parse(const std::string &input, const Options *opt = nullptr);

  // Singleton object
  static FieldFactory *get();

  /// clean the cache of parsed strings
  void cleanCache();
protected:
  // These functions called by the parser
  FieldGeneratorPtr resolve(std::string &name) override;

private:
  Mesh *fieldmesh;  
  const Options *options;

  std::list<std::string> lookup; // Names currently being parsed
  
  // Cache parsed strings
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
