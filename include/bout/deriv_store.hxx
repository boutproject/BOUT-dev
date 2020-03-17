/*!************************************************************************
 * \file deriv_store.hxx
 *
 * Definition of derivative methods storage class
 *
 **************************************************************************
 * Copyright 2018
 *    D.Dickinson, P.Hill, B.Dudson
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

#ifndef __DERIV_STORE_HXX__
#define __DERIV_STORE_HXX__

#include <functional>
#include <map>
#include <set>
#include <unordered_map>

#include <bout/scorepwrapper.hxx>

#include <bout_types.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <options.hxx>

/// Here we have a templated singleton that is used to store DerivativeFunctions
/// for all types of derivatives. It is templated on the FieldType (2D or 3D) as
/// the function interfaces also depend on this. It provides public routines for
/// registering and fetching derivative methods defined by a string key (e.g. "C2")
/// a DIRECTION (e.g. DIRECTION::X) and a STAGGER (e.g. STAGGER::None). There is
/// one routine for each class of derivative (standard, standard2nd, standard4th,
/// upwind and flux).
template <typename FieldType>
struct DerivativeStore {
  using standardFunc = std::function<void(const FieldType&, FieldType&, const std::string&)>;
  using flowFunc =
      std::function<void(const FieldType&, const FieldType&, FieldType&, const std::string&)>;
  using upwindFunc = flowFunc;
  using fluxFunc = flowFunc;

#ifdef USE_ORDERED_MAP_FOR_DERIVATIVE_STORE
  template <typename K, typename V>
  using storageType = std::map<K, V>;
#else
  template <typename K, typename V>
  using storageType = std::unordered_map<K, V>;
#endif

  // No copy constructor allowed
  DerivativeStore(const DerivativeStore& junk) = delete;

  // Singleton method
  static DerivativeStore& getInstance() {
    static DerivativeStore instance;
    return instance;
  }

  /// Report if store has any registered methods
  bool isEmpty() const {
    AUTO_TRACE();
    return registeredMethods.empty();
  };

  /// Report if store has any registered methods for specific type determined by key
  bool isEmpty(std::size_t key) const {
    AUTO_TRACE();
    return registeredMethods.count(key) == 0;
  }

  /// Report if store has any registered methods for specific type
  bool isEmpty(DERIV derivType, DIRECTION direction,
               STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();

    // Get the key
    auto key = getKey(direction, stagger, toString(derivType));
    return isEmpty(key);
  }

  /// Returns a vector of all registered method names for the
  /// specified derivative type, direction and stagger.
  std::set<std::string> getAvailableMethods(DERIV derivType, DIRECTION direction,
                                            STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();

    // Get the key
    auto key = getKey(direction, stagger, toString(derivType));
    if (isEmpty(key)) {
      return std::set<std::string>{};
    } else {
      return registeredMethods.at(key);
    }
  };

  /// Outputs a list of all registered method names for the
  /// specified derivative type, direction and stagger.
  void listAvailableMethods(DERIV derivType, DIRECTION direction,
                            STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();

    // Introductory information
    output_info << "Available methods for derivative type '";
    output_info << toString(derivType);
    output_info << "' in direction ";
    output_info << toString(direction);
    output_info << " ( Staggering : ";
    output_info << toString(stagger) << " ) : \n";

    for (const auto& i : getAvailableMethods(derivType, direction, stagger)) {
      output_info << "\t" << i << "\n";
    };
  };

  /// Register a function with standardFunc interface. Which map is used
  /// depends on the derivType input.
  void registerDerivative(standardFunc func, DERIV derivType, DIRECTION direction,
                          STAGGER stagger, std::string methodName) {
    AUTO_TRACE();
    const auto key = getKey(direction, stagger, methodName);

    switch (derivType) {
    case (DERIV::Standard):
      if (standard.count(key) != 0) {
        throw BoutException("Trying to override standard derivative : "
                            "direction %s, stagger %s, key %s",
                            toString(direction).c_str(),
                            toString(stagger).c_str(), methodName.c_str());
      }
      standard[key] = func;
      break;
    case (DERIV::StandardSecond):
      if (standardSecond.count(key) != 0) {
        throw BoutException("Trying to override standardSecond derivative : "
                            "direction %s, stagger %s, key %s",
                            toString(direction).c_str(),
                            toString(stagger).c_str(), methodName.c_str());
      }
      standardSecond[key] = func;
      break;
    case (DERIV::StandardFourth):
      if (standardFourth.count(key) != 0) {
        throw BoutException("Trying to override standardFourth derivative : "
                            "direction %s, stagger %s, key %s",
                            toString(direction).c_str(),
                            toString(stagger).c_str(), methodName.c_str());
      }
      standardFourth[key] = func;
      break;
    default:
      throw BoutException("Invalid function signature in registerDerivative : Function "
                          "signature 'standard' but derivative type %s passed",
                          toString(derivType).c_str());
    };

    // Register this method name in lookup of known methods
    registeredMethods[getKey(direction, stagger, toString(derivType))].insert(
        methodName);

  };

  /// Register a function with upwindFunc/fluxFunc interface. Which map is used
  /// depends on the derivType input.
  void registerDerivative(upwindFunc func, DERIV derivType, DIRECTION direction,
                          STAGGER stagger, std::string methodName) {
    AUTO_TRACE();
    const auto key = getKey(direction, stagger, methodName);

    switch (derivType) {
    case (DERIV::Upwind):
      if (upwind.count(key) != 0) {
        throw BoutException("Trying to override upwind derivative : "
                            "direction %s, stagger %s, key %s",
                            toString(direction).c_str(),
                            toString(stagger).c_str(), methodName.c_str());
      }
      upwind[key] = func;
      break;
    case (DERIV::Flux):
      if (flux.count(key) != 0) {
        throw BoutException("Trying to override flux derivative : "
                            "direction %s, stagger %s, key %s",
                            toString(direction).c_str(),
                            toString(stagger).c_str(), methodName.c_str());
      }
      flux[key] = func;
      break;
    default:
      throw BoutException("Invalid function signature in registerDerivative : Function "
                          "signature 'upwind/flux' but derivative type %s passed",
                          toString(derivType).c_str());
    };

    // Register this method name in lookup of known methods
    registeredMethods[getKey(direction, stagger, toString(derivType))].insert(
        methodName);
  };

  /// Templated versions of the above registration routines.
  template <typename Direction, typename Stagger, typename Method>
  void registerDerivative(standardFunc func, Direction direction, Stagger stagger,
                          Method method) {
    AUTO_TRACE();
    registerDerivative(func, method.meta.derivType, direction.lookup(), stagger.lookup(),
                       method.meta.key);
  };
  template <typename Direction, typename Stagger, typename Method>
  void registerDerivative(upwindFunc func, Direction direction, Stagger stagger,
                          Method method) {
    AUTO_TRACE();
    registerDerivative(func, method.meta.derivType, direction.lookup(), stagger.lookup(),
                       method.meta.key);
  };

  /// Routines to return a specific differential operator. Note we
  /// have to have a separate routine for different methods as they
  /// have different return types. As such we choose to use a
  /// different name for each of the method-classes so everything is
  /// consistently treated
  standardFunc getStandardDerivative(std::string name, DIRECTION direction,
                                     STAGGER stagger = STAGGER::None,
                                     DERIV derivType = DERIV::Standard) const {

    AUTO_TRACE();
    const auto realName = nameLookup(
        name, defaultMethods.at(getKey(direction, stagger, toString(derivType))));
    const auto key = getKey(direction, stagger, realName);

    const storageType<std::size_t, standardFunc>* theMap = nullptr;

    if (derivType == DERIV::Standard) {
      theMap = &standard;
    } else if (derivType == DERIV::StandardSecond) {
      theMap = &standardSecond;
    } else if (derivType == DERIV::StandardFourth) {
      theMap = &standardFourth;
    } else {
      throw BoutException("getStandardDerivative only works for derivType in {Standard, "
                          "StandardSecond, StandardFourth} but receieved %s",
                          toString(derivType).c_str());
    };

    const auto resultOfFind = theMap->find(key);
    if (resultOfFind != theMap->end())
      return resultOfFind->second;

    throw BoutException(
        "Couldn't find requested method %s in map for standard derivative of type %s.",
        getMethodName(realName, direction, stagger).c_str(),
        toString(derivType).c_str());
  };

  standardFunc getStandard2ndDerivative(std::string name, DIRECTION direction,
                                        STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();
    return getStandardDerivative(name, direction, stagger, DERIV::StandardSecond);
  };

  standardFunc getStandard4thDerivative(std::string name, DIRECTION direction,
                                        STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();
    return getStandardDerivative(name, direction, stagger, DERIV::StandardFourth);
  };

  flowFunc getFlowDerivative(std::string name, DIRECTION direction,
                             STAGGER stagger = STAGGER::None,
                             DERIV derivType = DERIV::Upwind) const {
    AUTO_TRACE();
    const auto realName = nameLookup(
        name, defaultMethods.at(getKey(direction, stagger, toString(derivType))));
    const auto key = getKey(direction, stagger, realName);

    const storageType<std::size_t, flowFunc>* theMap = nullptr;

    if (derivType == DERIV::Upwind) {
      theMap = &upwind;
    } else if (derivType == DERIV::Flux) {
      theMap = &flux;
    } else {
      throw BoutException(
          "getFlowDerivative only works for derivType in {Upwind, Flux} but receieved %s",
          toString(derivType).c_str());
    };

    const auto resultOfFind = theMap->find(key);
    if (resultOfFind != theMap->end())
      return resultOfFind->second;

    throw BoutException(
        "Couldn't find requested method %s in map for standard flow of type %s.",
        getMethodName(realName, direction, stagger).c_str(),
        toString(derivType).c_str());
  }

  upwindFunc getUpwindDerivative(std::string name, DIRECTION direction,
                                 STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();
    return getFlowDerivative(name, direction, stagger, DERIV::Upwind);
  };

  fluxFunc getFluxDerivative(std::string name, DIRECTION direction,
                             STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();
    return getFlowDerivative(name, direction, stagger, DERIV::Flux);
  };

  void initialise(Options* options) {
    AUTO_TRACE();

    // To replicate the existing behaviour we first search for a section called
    //"dd?" and if the option isn't in there we search a section called "diff"
    auto backupSection = options->getSection("diff");

    std::map<DIRECTION, std::string> directions = {{DIRECTION::X, "ddx"},
                                                   {DIRECTION::Y, "ddy"},
                                                   {DIRECTION::YOrthogonal, "ddy"},
                                                   {DIRECTION::Z, "ddz"}};

    std::map<DERIV, std::string> derivTypes = {{DERIV::Standard, "First"},
                                               {DERIV::StandardSecond, "Second"},
                                               {DERIV::StandardFourth, "Fourth"},
                                               {DERIV::Upwind, "Upwind"},
                                               {DERIV::Flux, "Flux"}};

    for (const auto& direction : directions) {
      for (const auto& deriv : derivTypes) {
        // This corresponds to the key in the input file
        auto derivName = deriv.second;

        const auto theDirection = direction.first;
        const auto theDerivType = deriv.first;
        const auto theDerivTypeString = toString(theDerivType);

        // Note both staggered and unstaggered have the same fallback default currently
        std::string theDefault{
            defaultMethods[getKey(theDirection, STAGGER::None, theDerivTypeString)]};

        //-------------------------------------------------------------
        // Unstaggered
        //-------------------------------------------------------------

        // The direction specific section to consider
        auto specificSection = options->getSection(direction.second);

        // Find the appropriate value for theDefault either from
        // the input file or if not found then use the value in
        // initialDefaultMethods
        // If neither branch matched (i.e. not in options) then we use
        // the hard-coded default already in the defaultMethods section
        if (specificSection->isSet(derivName)) {
          specificSection->get(derivName, theDefault, "");
        } else if (backupSection->isSet(derivName)) {
          backupSection->get(derivName, theDefault, "");
        }

        // Now we have the default method we should store it in defaultMethods
        theDefault = uppercase(theDefault);
        defaultMethods[getKey(theDirection, STAGGER::None, theDerivTypeString)] =
            theDefault;
        output_verbose << "The default method for derivative type " << theDerivTypeString
                       << " in direction " << toString(theDirection) << " is "
                       << theDefault << "\n";

        //-------------------------------------------------------------
        // Staggered
        //-------------------------------------------------------------

        // The direction specific section to consider with staggering
        specificSection = options->getSection(direction.second + "stag");

        // Find the appropriate value for theDefault either from
        // the input file. Note if the specific section for staggering isn't
        // found then we leave the default as for the non-staggered version
        if (specificSection->isSet(derivName)) {
          specificSection->get(derivName, theDefault, "");
        }

        // Now we have the default method we should store it in defaultMethods
        theDefault = uppercase(theDefault);
        defaultMethods[getKey(theDirection, STAGGER::L2C, theDerivTypeString)] =
            theDefault;
        defaultMethods[getKey(theDirection, STAGGER::C2L, theDerivTypeString)] =
            theDefault;
        output_verbose << "The default method for staggered derivative type "
                       << theDerivTypeString << " in direction "
                       << toString(theDirection) << " is " << theDefault << "\n";
      }
    }
  }

  /// Provide a method to override/force a specific default method
  void forceDefaultMethod(std::string methodName, DERIV deriv, DIRECTION direction,
                          STAGGER stagger = STAGGER::None) {
    const auto key = getKey(direction, stagger, toString(deriv));
    defaultMethods[key] = uppercase(methodName);
  }

  /// Empty all member storage
  void clear() {
    defaultMethods.clear();
    standard.clear();
    standardSecond.clear();
    standardFourth.clear();
    upwind.clear();
    flux.clear();
    registeredMethods.clear();
  }

  /// Reset to initial state
  void reset() {
    clear();

    setDefaults();
  }

private:
  // Make empty constructor private so we can't make instances outside
  // of the struct
  DerivativeStore() {
    // Ensure the default methods are set on construction
    // This populates the defaultMethods map
    setDefaults();
  }

  storageType<std::size_t, standardFunc> standard;
  storageType<std::size_t, standardFunc> standardSecond;
  storageType<std::size_t, standardFunc> standardFourth;
  storageType<std::size_t, upwindFunc> upwind;
  storageType<std::size_t, fluxFunc> flux;

  storageType<std::size_t, std::set<std::string>> registeredMethods;

  /// The following stores what actual method to use when DIFF_DEFAULT
  /// is passed. The key is determined using the getKey routine here,
  /// where the name we pass is determined by the type of method (standard,
  /// upwind etc.). Note for now we'll always use STAGGER::None as we
  /// currently assume the default method is independent of staggering --
  /// it might be useful to relax this assumption!
  storageType<std::size_t, std::string> defaultMethods;

  void setDefaults() {
    std::map<DERIV, std::string> initialDefaultMethods = {{DERIV::Standard, "C2"},
                                                          {DERIV::StandardSecond, "C2"},
                                                          {DERIV::StandardFourth, "C2"},
                                                          {DERIV::Upwind, "U1"},
                                                          {DERIV::Flux, "U1"}};

    std::map<DIRECTION, std::string> directions = {{DIRECTION::X, "ddx"},
                                                   {DIRECTION::Y, "ddy"},
                                                   {DIRECTION::YOrthogonal, "ddy"},
                                                   {DIRECTION::Z, "ddz"}};

    std::map<DERIV, std::string> derivTypes = {{DERIV::Standard, "First"},
                                               {DERIV::StandardSecond, "Second"},
                                               {DERIV::StandardFourth, "Fourth"},
                                               {DERIV::Upwind, "Upwind"},
                                               {DERIV::Flux, "Flux"}};

    for (const auto& direction : directions) {
      for (const auto& deriv : derivTypes) {
        const auto theDirection = direction.first;
        const auto theDerivTypeString = toString(deriv.first);
        const std::string theDefault{uppercase(initialDefaultMethods[deriv.first])};

        //-------------------------------------------------------------
        // Unstaggered and Staggered -- both get same default currently
        //-------------------------------------------------------------

        // Now we have the default method we should store it in defaultMethods
        defaultMethods[getKey(theDirection, STAGGER::None, theDerivTypeString)] =
            theDefault;
        defaultMethods[getKey(theDirection, STAGGER::L2C, theDerivTypeString)] =
            theDefault;
        defaultMethods[getKey(theDirection, STAGGER::C2L, theDerivTypeString)] =
            theDefault;
      }
    }
  };

  std::string getMethodName(std::string name, DIRECTION direction,
                            STAGGER stagger = STAGGER::None) const {
    AUTO_TRACE();
    return name + " (" + toString(direction) + ", " + toString(stagger)
           + ")";
  };

  std::string nameLookup(const std::string name, const std::string defaultName) const {
    return name != toString(DIFF_DEFAULT) ? name : defaultName;
  }

  /// Provides a routine to produce a unique key given information
  /// about the specific type required. This is templated so requires
  /// compile-time information. Need to also supply a non-templated
  /// version to account for run-time choices Note : We could include
  /// the derivType in the key -- this would allow us to store all
  /// methods with the same function interface in the same map, which
  /// might be nice.
  std::size_t getKey(DIRECTION direction, STAGGER stagger, std::string key) const {
    AUTO_TRACE();
    // Note this key is indepedent of the field type (and hence the
    // key is the same for 3D/2D fields) as we have to use different
    // maps to store the different field types as the signature is
    // different.
    std::size_t result;
    result = std::hash<std::string>{}(toString(direction));
    result = result ^ std::hash<std::string>{}(toString(stagger));
    result = result ^ std::hash<std::string>{}(key);
    return result;
  }

  /// Provides a routine to produce a unique key given information
  /// about the specific type required. This is templated so requires
  /// compile-time information. Makes use of a non-templated version
  /// that can be used to account for run-time choices
  template <typename Direction, typename Stagger, typename Method>
  std::size_t getKey() const {
    AUTO_TRACE();
    // Note this key is indepedent of the field type (and hence the
    // key is the same for 3D/2D fields) as we have to use different
    // maps to store the different field types as the signature is
    // different.
    return getKey(Direction{}.lookup(), Stagger{}.lookup(), Method{}.meta.key);
  }
};

#endif
