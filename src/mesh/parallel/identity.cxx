#include "bout/paralleltransform.hxx"
#include "bout/mesh.hxx"

void ParallelTransformIdentity::calcYUpDown(Field3D& f) {
  f.splitYupYdown();

  for (int i = 0; i < f.getMesh()->ystart; ++i) {
    f.yup(i) = f;
    f.ydown(i) = f;
  }
}
