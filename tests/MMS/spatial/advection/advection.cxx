#include "bout.hxx"
#include "derivs.hxx"
#include "field_factory.hxx"

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

  SAVE_ONCE2(f, g);
  SAVE_ONCE5(solution, result, error, l_2, l_inf);

  dump.write();

  BoutFinalise();
}
