#include "bout.hxx"
#include "derivs.hxx"
#include "field_factory.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  using bout::globals::mesh;

  Field3D input{FieldFactory::get()->create3D("input", Options::getRoot(), mesh)};
  Field3D solution{FieldFactory::get()->create3D("solution", Options::getRoot(), mesh)};

  // Communicate to calculate parallel transform
  mesh->communicate(input);

  Field3D result{Grad_par(input)};
  Field3D error{result - solution};
  BoutReal l_2{sqrt(mean(SQ(error), true, "RGN_NOBNDRY"))};
  BoutReal l_inf{max(abs(error), true, "RGN_NOBNDRY")};

  SAVE_ONCE6(input, solution, result, error, l_2, l_inf);

  for (int slice = 1; slice < mesh->ystart; ++slice) {
    SAVE_ONCE2(input.ynext(-slice), input.ynext(slice));
  }

  bout::globals::dump.write();

  BoutFinalise();
}
