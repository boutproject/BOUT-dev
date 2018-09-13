#include <bout.hxx>
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

  // Read variable from mesh
  Field3D var;
  mesh->get(var, "var");

  Field3D var2 = copy(var);

  // Var starts in orthogonal X-Z coordinates

  // Calculate yup and ydown
  mesh->communicate(var);

  // Calculate d/dy ysing yup() and ydown() fields
  Field3D ddy = DDY_yud(var);

  // Change into field-aligned coordinates
  Field3D var_aligned = mesh->toFieldAligned(var);
  var_aligned.applyBoundary("neumann");
  mesh->communicate(var_aligned);

  // var now field aligned
  Field3D ddy2 = DDY_aligned(var_aligned);
  mesh->communicate(ddy2);

  // Shift back to orthogonal X-Z coordinates
  ddy_check = fromFieldAligned(ddy_check);

  SAVE_ONCE3(ddy, ddy2, ddy_check);
  dump.write();

  BoutFinalise();

  return 0;
}
