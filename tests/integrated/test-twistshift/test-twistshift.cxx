#include "bout.hxx"
#include "field_factory.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D test = FieldFactory::get()->create3D("test");

  Field3D test_aligned = toFieldAligned(test);

  using bout::globals::mesh;

  // zero guard cells to check that communication is doing something
  for (int x=0; x<mesh->LocalNx; x++) {
    for (int z=0; z<mesh->LocalNz; z++) {
      for (int y=0; y<mesh->ystart; y++) {
        test_aligned(x, y, z) = 0.;
      }
      for (int y=mesh->yend+1; y<mesh->LocalNy; y++) {
        test_aligned(x, y, z) = 0.;
      }
    }
  }

  mesh->communicate(test_aligned);

  Field3D result = fromFieldAligned(test_aligned);

  SAVE_ONCE(test, test_aligned, result);

  bout::globals::dump.write();

  BoutFinalise();
}
