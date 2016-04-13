#include <bout.hxx>

#include <bout/paralleltransform.hxx>

#include <derivs.hxx>

const Field3D DDY_yud(const Field3D &f) {
  Field3D result;
  result.allocate();

  for(int i=0;i<mesh->ngx;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->ngz-1;k++)
    result(i,j,k) = 0.5*(f.yup()(i,j+1,k) - f.ydown()(i,j-1,k));

  return result;
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  ShiftedMetric s(*mesh);

  // Read variable from mesh
  Field3D var;
  mesh->get(var, "var");

  // var now field aligned
  Field3D ddy = DDY(var);

  // Shift into X-Z orthogonal
  //Field3D varxz = var.shiftZ(true);
  Field3D varxz = mesh->fromFieldAligned(var);

  // Calculate yup and ydown
  s.calcYUpDown(varxz);

  Field3D ddy2 = DDY_yud(varxz);

  // Shift back to field aligned
  ddy2 = ddy2.shiftZ(false);

  SAVE_ONCE2(ddy, ddy2);
  dump.write();

  BoutFinalise();

  return 0;
}
