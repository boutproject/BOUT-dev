#include "bout/bout.hxx"
#include "bout/difops.hxx"
#include "bout/field.hxx"
#include "bout/field3d.hxx"
#include "bout/field_factory.hxx"
#include "bout/globals.hxx"
#include "bout/options.hxx"
#include "bout/petsc_operators.hxx"

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  const PetscOperators ops(bout::globals::mesh);

  // Parallel operators
  const auto parallel = ops.getParallel();

  // Test field
  Field3D f = FieldFactory::get()->create3D("sin(y) + cos(z)");
  bout::globals::mesh->communicate(f);
  f.applyParallelBoundary("parallel_dirichlet_o2");

  Field3D f_neumann = FieldFactory::get()->create3D("sin(y) + cos(z)");
  bout::globals::mesh->communicate(f_neumann);
  f_neumann.applyParallelBoundary("parallel_neumann_o1");

  auto forward_op = parallel.Restrict_minus * parallel.Forward;

  Options dump;
  dump["f"] = f;
  dump["dV"] = parallel.dV;

  // Test1 : (parallel.Restrict_minus * parallel.Inject_plus) is the identity in the interior (not X boundaries)

  dump["grad_par_op"] = parallel.Grad_par(f);
  Field3D forward{0.0};
  BOUT_FOR(i, forward.getRegion("RGN_NOBNDRY")) { forward[i] = f.yup()[i.yp()]; }

  dump["forward"] = forward;

  dump["forward_op"] = forward_op(f);

  dump["grad_par_yud"] = Grad_par(f);

  // Test2 : Sum of dV * Div_par over domain is zero

  dump["div_par_op"] = parallel.Div_par(f);

  auto* coords = bout::globals::mesh->getCoordinates();
  bout::globals::mesh->communicate(coords->Bxy);
  coords->Bxy.applyParallelBoundary("parallel_neumann_o1");
  dump["div_par_yud"] = Div_par(f);

  dump["div_par_grad_par_op"] = parallel.Div_par_Grad_par(f_neumann);
  dump["div_par_grad_par_yud"] = Div_par_K_Grad_par(1.0, f_neumann);

  bout::writeDefaultOutputFile(dump);

  BoutFinalise();
  return 0;
}
