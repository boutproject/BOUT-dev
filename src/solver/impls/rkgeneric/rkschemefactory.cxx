#include "RKSchemefactory.hxx"

#include "impls/rkf45/rkf45.hxx"
#include "impls/cashkarp/cashkarp.hxx"

#include <boutexception.hxx>

RKSchemefactory* RKSchemefactory::instance = NULL;

RKSchemefactory* RKSchemefactory::getInstance() {
  if(instance == NULL) {
    // Create the singleton object
    instance = new RKSchemefactory();
  }
  return instance;
}

inline RKSchemeType RKSchemefactory::getDefaultRKSchemeType() {
  RKSchemeType type = NULL;
  type = RKSCHEME_RKF45;
  return type;
}

RKScheme* RKSchemefactory::createRKScheme(Options *options) {
  RKSchemeType type = getDefaultRKSchemeType();

  if(options == NULL) 
    options = Options::getRoot()->getSection("solver");
  
  string scheme;
  options->get("scheme", scheme, "");

  if(!scheme.empty()) type = scheme.c_str();

  return createRKScheme(type, options);
}

RKScheme* RKSchemefactory::createRKScheme(RKSchemeType &type, Options *options) {
  if(options == NULL)
    options = Options::getRoot()->getSection("solver");
  
  if(!strcasecmp(type, RKSCHEME_RKF45)) {
    return new RKF45Scheme(options);
  } else if(!strcasecmp(type, RKSCHEME_CASHKARP)) {
    return new CashKarpScheme(options);
  };

  // Need to throw an error saying 'Supplied option "type"' was not found
  throw BoutException("No such scheme exists in this build, type: %s", type);
}
