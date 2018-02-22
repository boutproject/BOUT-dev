
class LaplaceFactory;

#ifndef __LAPLACE_FACTORY_H__
#define __LAPLACE_FACTORY_H__

#include <bout/invert_laplace.hxx>

class LaplaceFactory {
 public:
  /// Return a pointer to the only instance
  static LaplaceFactory* getInstance();
  
  Laplacian* createLaplacian(Options *options = NULL);
  
private:
  LaplaceFactory() {} // Prevent instantiation of this class
  static LaplaceFactory* instance; ///< The only instance of this class (Singleton)
};

#endif // __LAPLACE_FACTORY_H__

