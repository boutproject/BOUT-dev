#include "bout.hxx"
#include "derivs.hxx"
#include "field_factory.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  using bout::globals::mesh;

  Field3D input{FieldFactory::get()->create3D("input_field", Options::getRoot(), mesh)};
  Field3D solution{FieldFactory::get()->create3D("solution", Options::getRoot(), mesh)};

  // Communicate to calculate parallel transform
  mesh->communicate(input);

  Field3D result{Grad_par(input)};
  Field3D error{result - solution};

  Options dump;

  dump["l_2"] = sqrt(mean(SQ(error), true, "RGN_NOBNDRY"));
  dump["l_inf"] = max(abs(error), true, "RGN_NOBNDRY");

  dump["result"] = result;
  dump["error"] = error;
  dump["input"] = input;
  dump["solution"] = solution;

  for (int slice = 1; slice < mesh->ystart; ++slice) {
    dump[fmt::format("input.ynext(-{})", slice)] = input.ynext(-slice);
    dump[fmt::format("input.ynext({})", slice)] = input.ynext(slice);
  }

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
}
