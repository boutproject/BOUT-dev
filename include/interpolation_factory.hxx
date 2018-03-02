#ifndef __INTERP_FACTORY_H__
#define __INTERP_FACTORY_H__

#include "interpolation.hxx"
#include "options.hxx"

#include <map>
#include <string>

using std::string;

class Mesh;

class InterpolationFactory {
public:
  /// Callback function definition for creating Interpolation objects
  typedef Interpolation* (*CreateInterpCallback)(Mesh*);
private:
  /// Prevent instantiation of this class
  InterpolationFactory();
  /// The only instance of this class (singleton)
  static InterpolationFactory* instance;

  /// Database of available interpolations
  std::map<string, CreateInterpCallback> interp_map;
  /// Look up interpolations in database
  CreateInterpCallback findInterpolation(const string &name);
public:
  ~InterpolationFactory() {};

  /// Return a pointer to the only instance
  static InterpolationFactory* getInstance();

  static void cleanup();

  inline string getDefaultInterpType();

  /// Create an interpolation object
  Interpolation *create(Mesh *mesh) {
    return create(nullptr, mesh);
  }
  Interpolation *create(Options *options) {
    return create(options, nullptr);
  }
  Interpolation *create(Options *options = nullptr, Mesh *mesh = nullptr);
  Interpolation *create(const string &name, Options *options = nullptr,
                        Mesh *mesh = nullptr);

  /// Add available interpolations to database
  void add(CreateInterpCallback interp, const string &name);
};

#endif //__INTERP_FACTORY_H__
