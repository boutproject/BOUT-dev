#ifndef __INTERP_FACTORY_H__
#define __INTERP_FACTORY_H__

#include "interpolation.hxx"
#include "options.hxx"

#include <map>
#include <string>

class Mesh;

class InterpolationFactory {
public:
  /// Callback function definition for creating Interpolation objects
  using CreateInterpCallback = Interpolation* (*)(Mesh*);

private:
  /// Add the available interpolation methods to the internal map
  ///
  /// Private default constructor to prevent instantiation of this
  /// class
  InterpolationFactory();
  /// The only instance of this class (singleton)
  static InterpolationFactory* instance;

  /// Database of available interpolations
  std::map<std::string, CreateInterpCallback> interp_map;

  /// Find an interpolation method in the list of available methods
  ///
  /// @param name Name of the interpolation method
  ///
  /// @return A pointer to the Interpolation object in the map
  CreateInterpCallback findInterpolation(const std::string& name);

public:
  ~InterpolationFactory() = default;
  ;

  /// Create or get the singleton instance of the factory
  static InterpolationFactory* getInstance();

  /// Destroy the singleton instance
  static void cleanup();

  /// A string representing the default interpolation type
  inline std::string getDefaultInterpType() { return "hermitespline"; }

  /// Create an interpolation object
  Interpolation* create(Mesh* mesh) { return create(nullptr, mesh); }
  Interpolation* create(Options* options = nullptr, Mesh* mesh = nullptr);

  /// Create an Interpolation object
  ///
  /// @param name    The name of the interpolation method
  /// @param options An Options object (e.g. an input file)
  /// @param mesh    A Mesh object to construct the interpolation on
  ///
  /// @return A new copy of an Interpolation object
  Interpolation* create(const std::string& name, Options* options = nullptr,
                        Mesh* mesh = nullptr);

  /// Add available interpolations to database
  void add(CreateInterpCallback interp, const std::string& name);
};

#endif //__INTERP_FACTORY_H__
