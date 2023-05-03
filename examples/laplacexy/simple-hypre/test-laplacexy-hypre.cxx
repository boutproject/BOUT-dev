#include <bout/bout.hxx>
#include <bout/field_factory.hxx>
#include <bout/invert/laplacexy2_hypre.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  {
    /// Create a LaplaceXY object
    LaplaceXY2Hypre laplacexy(bout::globals::mesh);

    /// Generate rhs function
    Field2D rhs = FieldFactory::get()->create2D("laplacexy:rhs", Options::getRoot(),
                                                bout::globals::mesh);

    /// Solution
    Field2D solution = 0.0;

    solution = laplacexy.solve(rhs, solution);

    Options dump;
    dump["rhs"] = rhs;
    dump["x"] = solution;
    bout::writeDefaultOutputFile(dump);
  }
  BoutFinalise();
#if BOUT_USE_CUDA
  cudaDeviceReset();
#endif
  return 0;
}
