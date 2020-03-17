
class LaplaceFactory;

#ifndef __LAPLACE_FACTORY_H__
#define __LAPLACE_FACTORY_H__

#include <invert_laplace.hxx>

class LaplaceFactory {
 public:
  /// Return a pointer to the only instance
  static LaplaceFactory* getInstance();

  Laplacian *createLaplacian(Options *options = nullptr, const CELL_LOC loc = CELL_CENTRE,
      Mesh *mesh_in = nullptr);

private:
  LaplaceFactory() {} // Prevent instantiation of this class
  static LaplaceFactory* instance; ///< The only instance of this class (Singleton)
};

#endif // __LAPLACE_FACTORY_H__

