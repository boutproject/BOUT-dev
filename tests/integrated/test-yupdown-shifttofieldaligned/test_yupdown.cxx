#include <bout.hxx>

#include <bout/paralleltransform.hxx>

#include <derivs.hxx>

// Y derivative assuming field is aligned in Y
const Field3D DDY_aligned(const Field3D &f) {
  Field3D result;
  result.allocate();
  result = 0.0;
  result.setCoordinateSystem(f.getCoordinateSystem());

  for(int i=0;i<mesh->LocalNx;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++)
      for(int k=0;k<mesh->LocalNz;k++)
	result(i,j,k) = 0.5*(f(i,j+1,k) - f(i,j-1,k));

  return result;
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  ShiftedMetric s(*mesh);

  // Read variable from mesh
  Field3D var;
  mesh->get(var, "var");

  // Var starts in orthogonal X-Z coordinates

  // Calculate d/dy using ShiftToFieldAligned
  Field3D ddy = DDY(var);

  // Change into field-aligned coordinates
  Field3D var_aligned = s.toFieldAligned(var);

  // var now field aligned
  Field3D ddy2 = DDY_aligned(var_aligned);

  // Shift back to orthogonal X-Z coordinates
  ddy2 = s.fromFieldAligned(ddy2);

  SAVE_ONCE2(ddy, ddy2);
  dump.write();

  BoutFinalise();

  return 0;
}
