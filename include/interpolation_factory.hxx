#ifndef __INTERP_FACTORY_H__
#define __INTERP_FACTORY_H__

#include "interpolation.hxx"
#include "options.hxx"

#include <map>
#include <string>

using std::string;

class InterpolationFactory {
public:
  /// Callback function definition for creating Interpolation objects
  typedef Interpolation* (*CreateInterpCallback)();
private:
  /// Prevent instantiation of this class
  InterpolationFactory();
  /// The only instance of this class (singleton)
  static InterpolationFactory* instance;

  /// Database of available interpolations
  std::map<string, CreateInterpCallback> interp_map;
  /// Look up interpolations in database
  Interpolation* findInterpolation(const string &name);
public:
  ~InterpolationFactory() {};

  /// Return a pointer to the only instance
  static InterpolationFactory* getInstance();

  static void cleanup();

  inline string getDefaultInterpType();

  /// Create an interpolation object
  Interpolation* create(Options *options = nullptr);
  Interpolation* create(const string &name, Options *options = nullptr);

  /// Add available interpolations to database
  void add(CreateInterpCallback interp, const string &name);
};

#endif //__INTERP_FACTORY_H__
