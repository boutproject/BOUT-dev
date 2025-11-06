#include "bout/bout.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx"
#include "bout/output.hxx"
#include "bout/parallel_boundary_region.hxx"
#include <fmt/format.h>
#include <vector>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  using bout::globals::mesh;

  std::vector<Field3D> fields(static_cast<int>(BoundaryParType::SIZE), Field3D{0.0});

  Options dump;
  for (int i = 0; i < fields.size(); i++) {
    fields[i].allocate();
    const auto boundary = static_cast<BoundaryParType>(i);
    const auto boundary_name = toString(boundary);
    mesh->communicate(fields[i]);
    for (const auto& bndry_par : mesh->getBoundariesPar(boundary)) {
      output.write("{:s} region\n", boundary_name);
      for (const auto& pnt : *bndry_par) {
        fields[i][pnt.ind()] += 1;
        output.write("{:s} increment\n", boundary_name);
      }
    }
    output.write("{:s} done\n", boundary_name);

    dump[fmt::format("field_{:s}", boundary_name)] = fields[i];
  }

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
}
