
class RKSchemeFactory;

#ifndef __RKSCHEME_FACTORY_H__
#define __RKSCHEME_FACTORY_H__

#include <bout/rkscheme.hxx>

#include <string.h>

class RKSchemeFactory {
 public:
  /// Return a pointer to the only instance
  static RKSchemeFactory* getInstance();
  
  RKSchemeType getDefaultRKSchemeType();
  
  RKScheme* createRKScheme(Options *opts = NULL);
  RKScheme* createRKScheme(RKSchemeType &, Options *opts = NULL);
  
 private:
  RKSchemeFactory() {} // Prevent instantiation of this class
  static RKSchemeFactory* instance; ///< The only instance of this class (Singleton)

};

#endif // __RKSCHEME_FACTORY_H__

