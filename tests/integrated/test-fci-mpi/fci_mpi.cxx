#include "bout/bout.hxx"
#include "bout/derivs.hxx"
#include "bout/field_factory.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  {
    using bout::globals::mesh;
    Options* options = Options::getRoot();
    int i = 0;
    const std::string default_str{"not_set"};
    Options dump;
    while (true) {
      std::string temp_str;
      options->get(fmt::format("input_{:d}:function", i), temp_str, default_str);
      if (temp_str == default_str) {
        break;
      }
      Field3D input{FieldFactory::get()->create3D(fmt::format("input_{:d}:function", i),
                                                  Options::getRoot(), mesh)};
      // options->get(fmt::format("input_{:d}:boundary_perp", i), temp_str, s"free_o3");
      mesh->communicate(input);
      input.applyParallelBoundary("parallel_neumann");
      for (int slice = -mesh->ystart; slice <= mesh->ystart; ++slice) {
        if (slice != 0) {
          Field3D tmp{0.};
          BOUT_FOR(i, tmp.getRegion("RGN_NOBNDRY")) {
            tmp[i] = input.ynext(slice)[i.yp(slice)];
          }
          dump[fmt::format("output_{:d}_{:+d}", i, slice)] = tmp;
        }
      }
      ++i;
    }
    bout::writeDefaultOutputFile(dump);
  }
  BoutFinalise();
}
