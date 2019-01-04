#ifndef __RKSCHEME_FACTORY_H__
#define __RKSCHEME_FACTORY_H__

#include "bout/rkscheme.hxx"

class Options;

class RKSchemeFactory {
 public:
  /// Return a pointer to the only instance
  static RKSchemeFactory* getInstance();
  
  RKSchemeType getDefaultRKSchemeType();

  RKScheme *createRKScheme(Options *opts = nullptr);
  RKScheme *createRKScheme(RKSchemeType &, Options *opts = nullptr);

private:
  RKSchemeFactory() {} // Prevent instantiation of this class
  static RKSchemeFactory* instance; ///< The only instance of this class (Singleton)

};

#endif // __RKSCHEME_FACTORY_H__

