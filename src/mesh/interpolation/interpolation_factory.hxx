#ifndef __INTERP_FACTORY_H__
#define __INTERP_FACTORY_H__

#include "interpolation.hxx"
#include "options.hxx"

#include <map>
#include <string>

using std::string;

class InterpolationFactory {
  /// Prevent instantiation of this class
  InterpolationFactory();
  /// The only instance of this class (singleton)
  static InterpolationFactory* instance;

  /// Database of available interpolations
  std::map<string, Interpolation*> interp_map;
  /// Look up interpolations in database
  Interpolation* findInterpolation(const string &name);
public:
  ~InterpolationFactory();

  /// Return a pointer to the only instance
  static InterpolationFactory* getInstance();

  /// Create an interpolation object
  Interpolation* createInterpolation(Options *opts = nullptr);
  Interpolation* createInterpolation(const string &name, Options *opts = nullptr);

  /// Add available interpolations to database
  void add(Interpolation* interp, const string &name);
};

#endif //__INTERP_FACTORY_H__
