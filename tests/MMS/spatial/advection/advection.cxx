#include "bout/bout.hxx"
#include "bout/derivs.hxx"
#include "bout/field_factory.hxx"

using bout::globals::mesh;

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field3D g{FieldFactory::get()->create3D("g", Options::getRoot(), mesh)};
  Field3D f{FieldFactory::get()->create3D("f", Options::getRoot(), mesh)};
  Field3D solution{FieldFactory::get()->create3D("solution", Options::getRoot(), mesh)};

  mesh->communicate(f, g);

  auto method = static_cast<BRACKET_METHOD>(Options::root()["method"].as<int>());

  Field3D result{bracket(g, f, method)};
  Field3D error{result - solution};
  BoutReal l_2{sqrt(mean(SQ(error), true, "RGN_NOBNDRY"))};
  BoutReal l_inf{max(abs(error), true, "RGN_NOBNDRY")};

  Options dump;
  dump["f"] = f;
  dump["g"] = g;
  dump["solution"] = solution;
  dump["result"] = result;
  dump["error"] = error;
  dump["l_2"] = l_2;
  dump["l_inf"] = l_inf;
  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
}
