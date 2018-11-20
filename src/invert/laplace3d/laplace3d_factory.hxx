
class Laplace3DFactory;

#ifndef __LAPLACE3D_FACTORY_H__
#define __LAPLACE3D_FACTORY_H__

#include <bout/invert/laplace3d.hxx>
#include <bout/options.hxx>

class Laplace3DFactory {
 public:
  /// Return a pointer to the only instance
  static Laplace3DFactory* getInstance();
  
  Laplace3D* createLaplace3D(Options *options = nullptr);
  
private:
  Laplace3DFactory() {} // Prevent instantiation of this class
  static Laplace3DFactory* instance; ///< The only instance of this class (Singleton)
};

#endif // __LAPLACE3D_FACTORY_H__

