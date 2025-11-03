#include "fmt/format.h"
#include "bout/bout.hxx"
#include "bout/field_factory.hxx"

namespace {
auto fci_mpi_test(int num, Options& dump) {
  using bout::globals::mesh;
  Field3D input{FieldFactory::get()->create3D(fmt::format("input_{:d}:function", num),
                                              Options::getRoot(), mesh)};
  mesh->communicate(input);

  input.applyParallelBoundary("parallel_neumann_o2");

  for (int slice = -mesh->ystart; slice <= mesh->ystart; ++slice) {
    if (slice == 0) {
      continue;
    }
    Field3D tmp{0.};
    BOUT_FOR(i, tmp.getRegion("RGN_NOBNDRY")) {
      tmp[i] = input.ynext(slice)[i.yp(slice)];
    }
    dump[fmt::format("output_{:d}_{:+d}", num, slice)] = tmp;
  }
}
} // namespace

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  Options dump;

  for (auto num : {0, 1, 2, 3}) {
    fci_mpi_test(num, dump);
  }

  bout::writeDefaultOutputFile(dump);
  BoutFinalise();
}
