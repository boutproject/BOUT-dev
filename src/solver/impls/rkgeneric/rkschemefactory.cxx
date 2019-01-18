#include "rkschemefactory.hxx"

#include "impls/rkf45/rkf45.hxx"
#include "impls/cashkarp/cashkarp.hxx"
#include "impls/rk4simple/rk4simple.hxx"
#include "impls/rkf34/rkf34.hxx"

#include <boutexception.hxx>

RKSchemeFactory *RKSchemeFactory::instance = nullptr;

RKSchemeFactory* RKSchemeFactory::getInstance() {
  if (instance == nullptr) {
    // Create the singleton object
    instance = new RKSchemeFactory();
  }
  return instance;
}

inline RKSchemeType RKSchemeFactory::getDefaultRKSchemeType() {
  RKSchemeType type = nullptr;
  type = RKSCHEME_RKF45;
  return type;
}

RKScheme* RKSchemeFactory::createRKScheme(Options *options) {
  RKSchemeType type = getDefaultRKSchemeType();

  if (options == nullptr)
    options = Options::getRoot()->getSection("solver");
  
  std::string scheme;
  options->get("scheme", scheme, "");

  if(!scheme.empty()) type = scheme.c_str();

  return createRKScheme(type, options);
}

RKScheme* RKSchemeFactory::createRKScheme(RKSchemeType &type, Options *options) {
  if (options == nullptr)
    options = Options::getRoot()->getSection("solver");
  
  if(!strcasecmp(type, RKSCHEME_RKF45)) {
    return new RKF45Scheme(options);
  }else if(!strcasecmp(type, RKSCHEME_CASHKARP)) {
    return new CASHKARPScheme(options);
  }else if(!strcasecmp(type, RKSCHEME_RK4)) {
    return new RK4SIMPLEScheme(options);
  }else if(!strcasecmp(type, RKSCHEME_RKF34)) {
    return new RKF34Scheme(options);
  };

  // Need to throw an error saying 'Supplied option "type"' was not found
  throw BoutException("No such scheme exists in this build, type: %s", type);
}
