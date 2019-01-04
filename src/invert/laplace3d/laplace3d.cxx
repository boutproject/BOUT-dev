
#include "bout/invert/laplace3d.hxx"

#include "laplace3d_factory.hxx"

Laplace3D* Laplace3D::create(Options *opts) {
  return Laplace3DFactory::getInstance()->createLaplace3D(opts);
}
