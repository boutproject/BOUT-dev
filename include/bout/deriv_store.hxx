#ifndef __DERIV_STORE_HXX__
#define __DERIV_STORE_HXX__

#include <functional>
#include <map>

#include <bout/scorepwrapper.hxx>
#include <bout_types.hxx>
#include <boutexception.hxx>
#include <msg_stack.hxx>

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

  // Singleton method
  static DerivativeStore &getInstance() {
    static DerivativeStore instance;
    return instance;
  }

  /// Register a function with standardFunc interface. Which map is used
  /// depends on the derivType input.
  void registerDerivative(const standardFunc func, const DERIV derivType,
			  const DIRECTION direction, const STAGGER stagger,
			  const std::string methodName) {
    TRACE("%s", __thefunc__);
    const auto key = getKey(direction, stagger, methodName);

    switch (derivType) {
    case (DERIV::Standard):
      getInstance().standard[key] = func;
      return;
    case (DERIV::StandardSecond):
      getInstance().standardSecond[key] = func;
      return;
    case (DERIV::StandardFourth):
      getInstance().standardFourth[key] = func;
      return;
    default:
      throw BoutException("Invalid function signature in registerDerivative : Function "
			  "signature 'standard' but derivative type %s passed",
			  DERIV_STRING(derivType).c_str());
    };
    return;
  };

  /// Register a function with upwindFunc/fluxFunc interface. Which map is used
  /// depends on the derivType input.
  void registerDerivative(const upwindFunc func, const DERIV derivType,
			  const DIRECTION direction, const STAGGER stagger,
			  const std::string methodName) {
    TRACE("%s", __thefunc__);
    const auto key = getKey(direction, stagger, methodName);
    switch (derivType) {
    case (DERIV::Upwind):
      getInstance().upwind[key] = func;
      return;
    case (DERIV::Flux):
      getInstance().flux[key] = func;
      return;
    default:
      throw BoutException("Invalid function signature in registerDerivative : Function "
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

  std::string getMethodName(std::string name, DIRECTION direction,
			    STAGGER stagger = STAGGER::None) {
    TRACE("%s", __thefunc__);
    return name + " (" + DIRECTION_STRING(direction) + ", " + STAGGER_STRING(stagger) +
	   ")";
  };

  /// Routines to return a specific differential operator. Note we have to have a separate
  /// routine for different
  /// methods as they have different return types. As such we choose to use a different
  /// name for each of the
  /// method-classes so everything is consistently treated
  standardFunc getStandardDerivative(std::string name, DIRECTION direction,
				     STAGGER stagger = STAGGER::None) {
    TRACE("%s", __thefunc__);
    const auto instance = getInstance();
    const auto key = getKey(direction, stagger, name);
    const auto resultOfFind = instance.standard.find(key);
    if (resultOfFind != instance.standard.end())
      return resultOfFind->second;
    throw BoutException(
	"Couldn't find requested method %s in map for standard derivative.",
	getMethodName(name, direction, stagger).c_str());
  };

  standardFunc getStandard2ndDerivative(std::string name, DIRECTION direction,
					STAGGER stagger = STAGGER::None) {
    TRACE("%s", __thefunc__);
    const auto instance = getInstance();
    const auto key = getKey(direction, stagger, name);
    const auto resultOfFind = instance.standardSecond.find(key);
    if (resultOfFind != instance.standardSecond.end())
      return resultOfFind->second;
    throw BoutException(
	"Couldn't find requested method %s in map for standardSecond derivative.",
	getMethodName(name, direction, stagger).c_str());
  };

  standardFunc getStandard4thDerivative(std::string name, DIRECTION direction,
					STAGGER stagger = STAGGER::None) {
    TRACE("%s", __thefunc__);
    const auto instance = getInstance();
    const auto key = getKey(direction, stagger, name);
    const auto resultOfFind = instance.standardFourth.find(key);
    if (resultOfFind != instance.standardFourth.end())
      return resultOfFind->second;
    throw BoutException(
	"Couldn't find requested method %s in map for standardFourth derivative.",
	getMethodName(name, direction, stagger).c_str());
  };
  upwindFunc getUpwindDerivative(std::string name, DIRECTION direction,
				 STAGGER stagger = STAGGER::None) {
    TRACE("%s", __thefunc__);
    const auto instance = getInstance();
    const auto key = getKey(direction, stagger, name);
    const auto resultOfFind = instance.upwind.find(key);
    if (resultOfFind != instance.upwind.end())
      return resultOfFind->second;
    throw BoutException("Couldn't find requested method %s in map for upwind derivative.",
			getMethodName(name, direction, stagger).c_str());
  };
  fluxFunc getFluxDerivative(std::string name, DIRECTION direction,
			     STAGGER stagger = STAGGER::None) {
    TRACE("%s", __thefunc__);
    const auto instance = getInstance();
    const auto key = getKey(direction, stagger, name);
    const auto resultOfFind = instance.flux.find(key);
    if (resultOfFind != instance.flux.end())
      return resultOfFind->second;
    throw BoutException("Couldn't find requested method %s in map for flux derivative.",
			getMethodName(name, direction, stagger).c_str());
  };

private:
  std::map<std::size_t, standardFunc> standard;
  std::map<std::size_t, standardFunc> standardSecond;
  std::map<std::size_t, standardFunc> standardFourth;
  std::map<std::size_t, upwindFunc> upwind;
  std::map<std::size_t, fluxFunc> flux;

  /// Provides a routine to produce a unique key given information about the specific type
  /// required. This is templated so requires compile-time information. Need to also
  /// supply
  /// a non-templated version to account for run-time choices
  /// Note : We could include the derivType in the key -- this would allow us to store
  /// all methods with the same function interface in the same map, which might be nice.
  std::size_t getKey(DIRECTION direction, STAGGER stagger, std::string key) {
    TRACE("%s", __thefunc__);
    // Note this key is indepedent of the field type (and hence the key is the same for
    // 3D/2D
    // fields) as we have to use different maps to store the different field types as the
    // signature is different.
    std::size_t result;
    result = std::hash<std::string>{}(DIRECTION_STRING(direction));
    result = result ^ std::hash<std::string>{}(STAGGER_STRING(stagger));
    result = result ^ std::hash<std::string>{}(key);
    return result;
  }

  /// Provides a routine to produce a unique key given information about the specific type
  /// required. This is templated so requires compile-time information. Makes use of
  /// a non-templated version that can be used to account for run-time choices
  template <typename Direction, typename Stagger, typename Method> std::size_t getKey() {
    TRACE("%s", __thefunc__);
    // Note this key is indepedent of the field type (and hence the key is the same for
    // 3D/2D
    // fields) as we have to use different maps to store the different field types as the
    // signature is different.
    return getKey(Direction{}.lookup(), Stagger{}.lookup(), Method{}.meta.key);
  }
};

#endif
