class ParDerivFactory;

#ifndef __PARDERIV_FACTORY_H__
#define __PARDERIV_FACTORY_H__

#include <bout/invert_parderiv.hxx>
#include <bout/options.hxx>

class ParDerivFactory {
 public:
  /// Return a pointer to the only instance
  static ParDerivFactory* getInstance();
  
  InvertPar* createInvertPar();
  InvertPar* createInvertPar(const char* type, Options *opt = NULL);
  InvertPar* createInvertPar(Options *opts);
 private:
  ParDerivFactory() {} // Prevent instantiation of this class
  static ParDerivFactory* instance; ///< The only instance of this class (Singleton)
};

#endif // __PARDERIV_FACTORY_H__
