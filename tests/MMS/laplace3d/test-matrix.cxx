#include <bout.hxx>

#include <invert_laplace.hxx>
#include <field_factory.hxx>
#include <bout/constants.hxx>
#include <derivs.hxx>
#include "../../../src/invert/laplace/impls/petscamg/petscamg.hxx"

Field3D this_Grad_perp2(const Field3D &f);
Field3D this_Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g);

int main(int argc, char **argv) {
  int init_err = BoutInitialise(argc, argv);
  if (init_err < 0) {
    return 0;
  } else if (init_err > 0) {
    return init_err;
  }

  // Create a Laplacian inversion solver
  LaplacePetscAmg *lap = (LaplacePetscAmg*) Laplacian::create();

  FieldFactory fact(mesh);

  std::shared_ptr<FieldGenerator> gen = fact.parse("input");
  output << "GEN = " << gen->str() << endl;

  Field3D rhs = fact.create3D("input");

  Field3D x = fact.create3D("solution");

  Field3D bout_rhs = this_Grad_perp2(x);

  Field3D petsc_rhs(mesh);
  petsc_rhs.allocate();
  for (int j=mesh->ystart; j<=mesh->yend; j++) {
    FieldPerp xslice = sliceXZ(x, j);
    FieldPerp result = lap->multiplyAx(xslice);
    petsc_rhs = result;
  }

  SAVE_ONCE4(x, rhs, bout_rhs, petsc_rhs);
  dump.write();

  BoutFinalise();
  return 0;
}

// Need custom version of perpendicular Laplacian operator, which is the inverse of the multigrid solver
// Delp2 uses FFT z-derivatives and Laplace includes y-derivatives, so can't use those
// The function is a copy of Laplace() with the y-derivatives deleted
Field3D this_Grad_perp2(const Field3D &f) {
  Field3D result = mesh->coordinates()->G1 * ::DDX(f) +  mesh->coordinates()->G3 * ::DDZ(f) +
                   mesh->coordinates()->g11 * ::D2DX2(f) + mesh->coordinates()->g33 * ::D2DZ2(f) +
                   2.0 * mesh->coordinates()->g13 * ::D2DXDZ(f);

  return result;
}

Field3D this_Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g) {
  Field3D result = mesh->coordinates()->g11 * ::DDX(f) * ::DDX(g) + mesh->coordinates()->g33 * ::DDZ(f) * ::DDZ(g)
                   + mesh->coordinates()->g13 * (DDX(f)*DDZ(g) + DDZ(f)*DDX(g));
  
  return result;
}
