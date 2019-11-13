#include <bout.hxx>

#include <bout/paralleltransform.hxx>
#include <derivs.hxx>

// Y derivative using yup() and ydown() fields
const Field3D DDY_yud(const Field3D& f) {
  Field3D result = emptyFrom(f);

  for (int i = 0; i < mesh->LocalNx; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++)
      for (int k = 0; k < mesh->LocalNz; k++)
        result(i, j, k) = 0.5 * (f.yup()(i, j + 1, k) - f.ydown()(i, j - 1, k));

  return result;
}

// Y derivative assuming field is aligned in Y
const Field3D DDY_aligned(const Field3D& f) {
  Field3D result = emptyFrom(f);

  for (int i = 0; i < mesh->LocalNx; i++)
    for (int j = mesh->ystart; j <= mesh->yend; j++)
      for (int k = 0; k < mesh->LocalNz; k++)
        result(i, j, k) = 0.5 * (f(i, j + 1, k) - f(i, j - 1, k));

  return result;
}

int main(int argc, char** argv) {

  BoutInitialise(argc, argv);

  Field2D zShift;
  mesh->get(zShift, "zShift");

  ShiftedMetric s(*mesh, CELL_CENTRE, zShift, mesh->getCoordinates()->zlength());

  // Read variable from mesh
  Field3D var;
  mesh->get(var, "var");

  Field3D var2 = copy(var);

  // Var starts in orthogonal X-Z coordinates

  // Calculate yup and ydown
  s.calcParallelSlices(var);

  // Calculate d/dy using yup() and ydown() fields
  Field3D ddy = DDY(var);

  // Calculate d/dy by transform to field-aligned coordinates
  // (var2 has no yup/ydown fields)
  Field3D ddy2 = DDY(var2);

  // Change into field-aligned coordinates
  Field3D var_aligned = toFieldAligned(var);

  // var now field aligned
  Field3D ddy_check = DDY_aligned(var_aligned);

  // Shift back to orthogonal X-Z coordinates
  ddy_check = fromFieldAligned(ddy_check);

  SAVE_ONCE3(ddy, ddy2, ddy_check);
  dump.write();

  BoutFinalise();

  return 0;
}
