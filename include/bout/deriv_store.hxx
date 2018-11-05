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
template <typename FieldType> struct DerivativeStore {
  using standardFunc = std::function<void(const FieldType &, FieldType &, const REGION)>;
  using upwindFunc = std::function<void(const FieldType &, const FieldType &, FieldType &,
					const REGION)>;
  using fluxFunc = std::function<void(const FieldType &, const FieldType &, FieldType &,
				      const REGION)>;
  
#ifdef USE_ORDERED_MAP_FOR_DERIVATIVE_STORE
  template<typename K, typename V>
  using storageType = std::map<K, V>;
#else
  template<typename K, typename V>  
  using storageType = std::unordered_map<K, V>;
#endif

  // Singleton method
  static DerivativeStore& getInstance() {
    static DerivativeStore instance;
    return instance;
  }

  /// Returns a vector of all registered method names for the
  /// specified derivative type, direction and stagger.
  std::vector<std::string> getAvailableMethods(DERIV derivType, DIRECTION direction,
					       STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    // Get the key
    auto key = getKey(direction, stagger, DERIV_STRING(derivType));
    return registeredMethods.at(key);
  };

  /// Outputs a list of all registered method names for the
  /// specified derivative type, direction and stagger.
  void listAvailableMethods(DERIV derivType, DIRECTION direction,
                            STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);

    // Introductory information
    output_info << "Available methods for derivative type '";
    output_info << DERIV_STRING(derivType);
    output_info << "' in direction ";
    output_info << DIRECTION_STRING(direction);
    output_info << " ( Staggering : ";
    output_info << STAGGER_STRING(stagger) << " ) : \n";

    for (const auto &i :
         getAvailableMethods(derivType, direction, stagger)) {
      output_info << "\t" << i << "\n";
    };
  };

  /// Register a function with standardFunc interface. Which map is used
  /// depends on the derivType input.
  void registerDerivative(standardFunc func, DERIV derivType,
                          DIRECTION direction, STAGGER stagger,
                          std::string methodName) {
    TRACE("%s", __thefunc__);
    const auto key = getKey(direction, stagger, methodName);

    // Register this method name in lookup of known methods
    registeredMethods[getKey(direction, stagger, DERIV_STRING(derivType))]
      .push_back(methodName);

    switch (derivType) {
    case (DERIV::Standard):
      standard[key] = func;
      return;
    case (DERIV::StandardSecond):
      standardSecond[key] = func;
      return;
    case (DERIV::StandardFourth):
      standardFourth[key] = func;
      return;
    default:
      throw BoutException(
          "Invalid function signature in registerDerivative : Function "
          "signature 'standard' but derivative type %s passed",
          DERIV_STRING(derivType).c_str());
    };
    return;
  };

  /// Register a function with upwindFunc/fluxFunc interface. Which map is used
  /// depends on the derivType input.
  void registerDerivative(upwindFunc func, DERIV derivType, DIRECTION direction,
                          STAGGER stagger, std::string methodName) {
    TRACE("%s", __thefunc__);
    const auto key = getKey(direction, stagger, methodName);

    // Register this method name in lookup of known methods
    registeredMethods[getKey(direction, stagger, DERIV_STRING(derivType))]
        .push_back(methodName);

    switch (derivType) {
    case (DERIV::Upwind):
      upwind[key] = func;
      return;
    case (DERIV::Flux):
      flux[key] = func;
      return;
    default:
      throw BoutException(
          "Invalid function signature in registerDerivative : Function "
          "signature 'upwind/flux' but derivative type %s passed",
          DERIV_STRING(derivType).c_str());
    };
    return;
  };

  /// Templated versions of the above registration routines.
  template <typename Direction, typename Stagger, typename Method>
  void registerDerivative(standardFunc func) {
    TRACE("%s", __thefunc__);
    registerDerivative(func, Method{}.meta.derivType, Direction{}.lookup(),
                       Stagger{}.lookup(), Method{}.meta.key);
  };
  template <typename Direction, typename Stagger, typename Method>
  void registerDerivative(upwindFunc func) {
    TRACE("%s", __thefunc__);
    registerDerivative(func, Method{}.meta.derivType, Direction{}.lookup(),
                       Stagger{}.lookup(), Method{}.meta.key);
  };

  /// Routines to return a specific differential operator. Note we have to have
  /// a separate
  /// routine for different
  /// methods as they have different return types. As such we choose to use a
  /// different
  /// name for each of the
  /// method-classes so everything is consistently treated
  standardFunc getStandardDerivative(std::string name, DIRECTION direction,
                                     STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    const auto realName =
        nameLookup(name,
                   defaultMethods.at(getKey(direction, stagger,
                                            DERIV_STRING(DERIV::Standard))));
    const auto key = getKey(direction, stagger, realName);
    const auto resultOfFind = standard.find(key);
    if (resultOfFind != standard.end())
      return resultOfFind->second;
    throw BoutException(
        "Couldn't find requested method %s in map for standard derivative.",
        getMethodName(realName, direction, stagger).c_str());
  };

  standardFunc getStandard2ndDerivative(std::string name, DIRECTION direction,
                                        STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    const auto realName = nameLookup(
        name,
        defaultMethods.at(
            getKey(direction, stagger, DERIV_STRING(DERIV::StandardSecond))));
    const auto key = getKey(direction, stagger, realName);
    const auto resultOfFind = standardSecond.find(key);
    if (resultOfFind != standardSecond.end())
      return resultOfFind->second;
    throw BoutException("Couldn't find requested method %s in map for "
                        "standardSecond derivative.",
                        getMethodName(realName, direction, stagger).c_str());
  };

  standardFunc getStandard4thDerivative(std::string name, DIRECTION direction,
                                        STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    const auto realName = nameLookup(
        name,
        defaultMethods.at(
            getKey(direction, stagger, DERIV_STRING(DERIV::StandardFourth))));
    const auto key = getKey(direction, stagger, realName);
    const auto resultOfFind = standardFourth.find(key);
    if (resultOfFind != standardFourth.end())
      return resultOfFind->second;
    throw BoutException("Couldn't find requested method %s in map for "
                        "standardFourth derivative.",
                        getMethodName(realName, direction, stagger).c_str());
  };
  
  upwindFunc getUpwindDerivative(std::string name, DIRECTION direction,
                                 STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    const auto realName =
        nameLookup(name,
                   defaultMethods.at(getKey(direction, stagger,
                                            DERIV_STRING(DERIV::Upwind))));
    const auto key = getKey(direction, stagger, realName);
    const auto resultOfFind = upwind.find(key);
    if (resultOfFind != upwind.end())
      return resultOfFind->second;
    throw BoutException(
        "Couldn't find requested method %s in map for upwind derivative.",
        getMethodName(realName, direction, stagger).c_str());
  };
  fluxFunc getFluxDerivative(std::string name, DIRECTION direction,
                             STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    const auto realName =
        nameLookup(name,
                   defaultMethods.at(
                       getKey(direction, stagger, DERIV_STRING(DERIV::Flux))));
    const auto key = getKey(direction, stagger, realName);
    const auto resultOfFind = flux.find(key);
    if (resultOfFind != flux.end())
      return resultOfFind->second;
    throw BoutException(
        "Couldn't find requested method %s in map for flux derivative.",
        getMethodName(realName, direction, stagger).c_str());
  };

  void initialise(Options *options) {
    TRACE("%s", __thefunc__);

    // To replicate the existing behaviour we first search for a section called
    //"dd?" and if the option isn't in there we search a section called "diff"
    auto backupSection = options->getSection("diff");

    std::map<DERIV, std::string> initialDefaultMethods = {
        {DERIV::Standard, "C2"},
        {DERIV::StandardSecond, "C2"},
        {DERIV::StandardFourth, "C2"},
        {DERIV::Upwind, "U1"},
        {DERIV::Flux, "U1"}};

    std::map<DIRECTION, std::string> directions = {
        {DIRECTION::X, "ddx"},
        {DIRECTION::Y, "ddy"},
        {DIRECTION::YOrthogonal, "ddy"},
        {DIRECTION::Z, "ddz"}};

    std::map<DERIV, std::string> derivTypes = {
        {DERIV::Standard, "First"},
        {DERIV::StandardSecond, "Second"},
        {DERIV::StandardFourth, "Fourth"},
        {DERIV::Upwind, "Upwind"},
        {DERIV::Flux, "Flux"}};

    for (const auto &direction : directions) {
      for (const auto &deriv : derivTypes) {
        std::string theDefault;

        // This corresponds to the key in the input file
        auto derivName = deriv.second;

        const auto theDirection = direction.first;
        const auto theDerivType = deriv.first;
        const auto theDerivTypeString = DERIV_STRING(theDerivType);

        //-------------------------------------------------------------
        // Unstaggered
        //-------------------------------------------------------------

        // The direction specific section to consider
        auto specificSection = options->getSection(direction.second);

        // Find the appropriate value for theDefault either from
        // the input file or if not found then use the value in
        // initialDefaultMethods
        if (specificSection->isSet(derivName)) {
          specificSection->get(derivName, theDefault, "");
        } else if (backupSection->isSet(derivName)) {
          backupSection->get(derivName, theDefault, "");
        } else {
          theDefault = initialDefaultMethods[theDerivType];
        }

        // Now we have the default method we should store it in defaultMethods
        theDefault = uppercase(theDefault);
        defaultMethods[getKey(theDirection, STAGGER::None,
                              theDerivTypeString)] = theDefault;
        output_info << "The default method for derivative type "
                    << theDerivTypeString;
        output_info << " in direction " << DIRECTION_STRING(theDirection);
        output_info << " is " << theDefault << "\n";

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
        output_info << "The default method for staggered derivative type "
                    << theDerivTypeString;
        output_info << " in direction " << DIRECTION_STRING(theDirection);
        output_info << " is " << theDefault << "\n";
      }
    }
  }

  /// The following stores what actual method to use when DIFF_DEFAULT
  /// is passed. The key is determined using the getKey routine here,
  /// where the name we pass is determined by the type of method (standard,
  /// upwind etc.). Note for now we'll always use STAGGER::None as we
  /// currently assume the default method is independent of staggering --
  /// it might be useful to relax this assumption!
  storageType<std::size_t, std::string> defaultMethods;

private:
  // Make the constructor private so we can't make instances outside
  // of the struct
  DerivativeStore(){};
  // No copy constructor allowed
  DerivativeStore(const DerivativeStore &junk) = delete;

  storageType<std::size_t, standardFunc> standard;
  storageType<std::size_t, standardFunc> standardSecond;
  storageType<std::size_t, standardFunc> standardFourth;
  storageType<std::size_t, upwindFunc> upwind;
  storageType<std::size_t, fluxFunc> flux;

  storageType<std::size_t, std::vector<std::string>> registeredMethods;

  std::string getMethodName(std::string name, DIRECTION direction,
                            STAGGER stagger = STAGGER::None) const {
    TRACE("%s", __thefunc__);
    return name + " (" + DIRECTION_STRING(direction) + ", " +
           STAGGER_STRING(stagger) + ")";
  };

  std::string nameLookup(const std::string name,
                         const std::string defaultName) const {
    return name != DIFF_METHOD_STRING(DIFF_DEFAULT) ? name : defaultName;
  }

  /// Provides a routine to produce a unique key given information about the
  /// specific type
  /// required. This is templated so requires compile-time information. Need to
  /// also
  /// supply
  /// a non-templated version to account for run-time choices
  /// Note : We could include the derivType in the key -- this would allow us to
  /// store
  /// all methods with the same function interface in the same map, which might
  /// be nice.
  std::size_t getKey(DIRECTION direction, STAGGER stagger,
                     std::string key) const {
    TRACE("%s", __thefunc__);
    // Note this key is indepedent of the field type (and hence the key is the
    // same for
    // 3D/2D
    // fields) as we have to use different maps to store the different field
    // types as the
    // signature is different.
    std::size_t result;
    result = std::hash<std::string>{}(DIRECTION_STRING(direction));
    result = result ^ std::hash<std::string>{}(STAGGER_STRING(stagger));
    result = result ^ std::hash<std::string>{}(key);
    return result;
  }

  /// Provides a routine to produce a unique key given information about the
  /// specific type
  /// required. This is templated so requires compile-time information. Makes
  /// use of
  /// a non-templated version that can be used to account for run-time choices
  template <typename Direction, typename Stagger, typename Method>
  std::size_t getKey() const {
    TRACE("%s", __thefunc__);
    // Note this key is indepedent of the field type (and hence the key is the
    // same for
    // 3D/2D
    // fields) as we have to use different maps to store the different field
    // types as the
    // signature is different.
    return getKey(Direction{}.lookup(), Stagger{}.lookup(), Method{}.meta.key);
  }
};

#endif
