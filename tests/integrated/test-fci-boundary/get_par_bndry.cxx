#include "bout.hxx"
#include "derivs.hxx"
#include "field_factory.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  using bout::globals::mesh;

  std::vector<Field3D> fields;
  fields.resize(static_cast<int>(BoundaryParType::SIZE));
  for (int i=0; i< fields.size(); i++){
    fields[i] = Field3D{0.0};
    mesh->communicate(fields[i]);
    for (const auto &bndry_par : mesh->getBoundariesPar(static_cast<BoundaryParType>(i))) {
      output.write("{:s} region\n", toString(static_cast<BoundaryParType>(i)));
      for (bndry_par->first(); !bndry_par->isDone(); bndry_par->next()) {
	fields[i](bndry_par->x, bndry_par->y, bndry_par->z) += 1;
	output.write("{:s} increment\n", toString(static_cast<BoundaryParType>(i)));
      }
    }
    output.write("{:s} done\n", toString(static_cast<BoundaryParType>(i)));
    bout::globals::dump.addOnce(fields[i], fmt::format("field_{:s}", toString(static_cast<BoundaryParType>(i))));
  }
  
  bout::globals::dump.write();

  BoutFinalise();
}
