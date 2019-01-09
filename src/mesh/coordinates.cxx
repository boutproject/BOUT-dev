/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include <bout/assert.hxx>
#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <msg_stack.hxx>
#include <output.hxx>
#include <utils.hxx>

#include <derivs.hxx>
#include <fft.hxx>
#include <interpolation.hxx>

#include <globals.hxx>

Coordinates::Coordinates(Mesh *mesh)
    : dx(1, mesh), dy(1, mesh), dz(1), d1_dx(mesh), d1_dy(mesh), J(1, mesh), Bxy(1, mesh),
      // Identity metric tensor
      g11(1, mesh), g22(1, mesh), g33(1, mesh), g12(0, mesh), g13(0, mesh), g23(0, mesh),
      g_11(1, mesh), g_22(1, mesh), g_33(1, mesh), g_12(0, mesh), g_13(0, mesh),
      g_23(0, mesh), G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh),
      G1_23(mesh), G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh),
      G2_23(mesh), G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh),
      G3_23(mesh), G1(mesh), G2(mesh), G3(mesh), ShiftTorsion(mesh),
      IntShiftTorsion(mesh), localmesh(mesh), location(CELL_CENTRE) {

  if (mesh->get(dx, "dx", 1.0, false)) {
    output_warn.write("\tWARNING: differencing quantity 'dx' not found. Set to 1.0\n");
  }
  mesh->communicateXZ(dx);

  if (mesh->get(dy, "dy", 1.0, false)) {
    output_warn.write("\tWARNING: differencing quantity 'dy' not found. Set to 1.0\n");
  }
  mesh->communicateXZ(dy);

  nz = mesh->LocalNz;

  if (mesh->get(dz, "dz")) {
    // Couldn't read dz from input
    int zperiod;
    BoutReal ZMIN, ZMAX;
    Options *options = Options::getRoot();
    if (options->isSet("zperiod")) {
      OPTION(options, zperiod, 1);
      ZMIN = 0.0;
      ZMAX = 1.0 / static_cast<BoutReal>(zperiod);
    } else {
      OPTION(options, ZMIN, 0.0);
      OPTION(options, ZMAX, 1.0);

      zperiod = ROUND(1.0 / (ZMAX - ZMIN));
    }

    dz = (ZMAX - ZMIN) * TWOPI / nz;
  }

  // Diagonal components of metric tensor g^{ij} (default to 1)
  mesh->get(g11, "g11", 1.0, false);
  mesh->get(g22, "g22", 1.0, false);
  mesh->get(g33, "g33", 1.0, false);

  // Off-diagonal elements. Default to 0
  mesh->get(g12, "g12", 0.0, false);
  mesh->get(g13, "g13", 0.0, false);
  mesh->get(g23, "g23", 0.0, false);

  mesh->communicateXZ(g11, g22, g33, g12, g13, g23);

  // Check input metrics
  if ((!finite(g11)) || (!finite(g22)) || (!finite(g33))) {
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if ((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if ((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }

  /// Find covariant metric components
  // Check if any of the components are present
  if (mesh->sourceHasVar("g_11") or mesh->sourceHasVar("g_22") or
      mesh->sourceHasVar("g_33") or mesh->sourceHasVar("g_12") or
      mesh->sourceHasVar("g_13") or mesh->sourceHasVar("g_23")) {
    // Check that all components are present
    if (mesh->sourceHasVar("g_11") and mesh->sourceHasVar("g_22") and
        mesh->sourceHasVar("g_33") and mesh->sourceHasVar("g_12") and
        mesh->sourceHasVar("g_13") and mesh->sourceHasVar("g_23")) {
      mesh->get(g_11, "g_11", 1.0, false);
      mesh->get(g_22, "g_22", 1.0, false);
      mesh->get(g_33, "g_33", 1.0, false);
      mesh->get(g_12, "g_12", 0.0, false);
      mesh->get(g_13, "g_13", 0.0, false);
      mesh->get(g_23, "g_23", 0.0, false);

      mesh->communicateXZ(g_11, g_22, g_33, g_12, g_13, g_23);

      output_warn.write("\tWARNING! Covariant components of metric tensor set manually. "
                        "Contravariant components NOT recalculated\n");

    } else {
      output_warn.write("Not all covariant components of metric tensor found. "
                        "Calculating all from the contravariant tensor\n");
      /// Calculate contravariant metric components if not found
      if (calcCovariant()) {
        throw BoutException("Error in calcCovariant call");
      }
    }
  } else {
    /// Calculate contravariant metric components if not found
    if (calcCovariant()) {
      throw BoutException("Error in calcCovariant call");
    }
  }

  /// Calculate Jacobian and Bxy
  if (jacobian())
    throw BoutException("Error in jacobian call");

  // Attempt to read J from the grid file
  auto Jcalc = J;
  if (mesh->get(J, "J", 0.0, false)) {
    output_warn.write("\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    J = Jcalc;
  } else {
    // Compare calculated and loaded values
    output_warn.write("\tMaximum difference in J is %e\n", max(abs(J - Jcalc)));

    mesh->communicateXZ(J);

    // Re-evaluate Bxy using new J
    Bxy = sqrt(g_22) / J;
  }

  // Attempt to read Bxy from the grid file
  auto Bcalc = Bxy;
  if (mesh->get(Bxy, "Bxy", 0.0, false)) {
    output_warn.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from "
                      "metric tensor\n");
    Bxy = Bcalc;
  } else {
    mesh->communicateXZ(Bxy);

    output_warn.write("\tMaximum difference in Bxy is %e\n", max(abs(Bxy - Bcalc)));
    // Check Bxy
    if (!finite(Bxy)) {
      throw BoutException("\tERROR: Bxy not finite everywhere!\n");
    }
  }

  if (mesh->get(ShiftTorsion, "ShiftTorsion", 0.0, false)) {
    output_warn.write("\tWARNING: No Torsion specified for zShift. Derivatives may not be correct\n");
  }
  mesh->communicateXZ(ShiftTorsion);

  //////////////////////////////////////////////////////

  if (mesh->IncIntShear) {
    if (mesh->get(IntShiftTorsion, "IntShiftTorsion", 0.0, false)) {
      output_warn.write("\tWARNING: No Integrated torsion specified\n");
    }
    mesh->communicateXZ(IntShiftTorsion);
  }
}

// use anonymous namespace so this utility function is not available outside this file
namespace {
  /// Interpolate a Field2D to a new CELL_LOC with interp_to.
  /// Communicates to set internal guard cells.
  /// Boundary guard cells are set equal to the nearest grid point (equivalent to
  /// 2nd order accurate Neumann boundary condition).
  /// Corner guard cells are set to BoutNaN
Coordinates::metric_field_type
interpolateAndNeumann(const Coordinates::metric_field_type& f, CELL_LOC location) {
#ifndef COORDINATES_USE_3D
  Mesh* localmesh = f.getMesh();
  auto result = interp_to(f, location, RGN_NOBNDRY);
  localmesh->communicate(result);

  // Copy nearest value into boundaries so that differential geometry terms can
  // be interpolated if necessary
  // Note: cannot use applyBoundary("neumann") here because applyBoundary()
  // would try to create a new Coordinates object since we have not finished
  // initializing yet, leading to an infinite recursion
  for (auto bndry : localmesh->getBoundaries()) {
    if (bndry->bx != 0) {
      // If bx!=0 we are on an x-boundary, inner if bx>0 and outer if bx<0
      for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
        for (int i = 0; i < localmesh->xstart; i++)
          result(bndry->x + i * bndry->bx, bndry->y) =
              result(bndry->x + (i - 1) * bndry->bx, bndry->y - bndry->by);
      }
    }
    if (bndry->by != 0) {
      // If by!=0 we are on a y-boundary, upper if by>0 and lower if by<0
      for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
        for (int i = 0; i < localmesh->ystart; i++)
          result(bndry->x, bndry->y + i * bndry->by) =
              result(bndry->x - bndry->bx, bndry->y + (i - 1) * bndry->by);
      }
    }
  }

  // Set corner guard cells
  for (int i = 0; i < localmesh->xstart; i++) {
    for (int j = 0; j < localmesh->ystart; j++) {
      result(i, j) = BoutNaN;
      result(i, localmesh->LocalNy - 1 - j) = BoutNaN;
      result(localmesh->LocalNx - 1 - i, j) = BoutNaN;
      result(localmesh->LocalNx - 1 - i, localmesh->LocalNy - 1 - j) = BoutNaN;
    }
  }

  return result;
#else
  throw BoutException(
      "Staggered coordinates locations not currently supported with 3D metrics.");
#endif
  }
}

Coordinates::Coordinates(Mesh *mesh, const CELL_LOC loc, const Coordinates* coords_in)
    : dx(1, mesh), dy(1, mesh), dz(1), d1_dx(mesh), d1_dy(mesh), J(1, mesh), Bxy(1, mesh),
      // Identity metric tensor
      g11(1, mesh), g22(1, mesh), g33(1, mesh), g12(0, mesh), g13(0, mesh), g23(0, mesh),
      g_11(1, mesh), g_22(1, mesh), g_33(1, mesh), g_12(0, mesh), g_13(0, mesh),
      g_23(0, mesh), G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh),
      G1_23(mesh), G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh),
      G2_23(mesh), G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh),
      G3_23(mesh), G1(mesh), G2(mesh), G3(mesh), ShiftTorsion(mesh),
      IntShiftTorsion(mesh), localmesh(mesh), location(loc) {

  dx = interpolateAndNeumann(coords_in->dx, location);
  dy = interpolateAndNeumann(coords_in->dy, location);

  nz = mesh->LocalNz;

  dz = coords_in->dz;

  // Diagonal components of metric tensor g^{ij}
  g11 = interpolateAndNeumann(coords_in->g11, location);
  g22 = interpolateAndNeumann(coords_in->g22, location);
  g33 = interpolateAndNeumann(coords_in->g33, location);

  // Off-diagonal elements.
  g12 = interpolateAndNeumann(coords_in->g12, location);
  g13 = interpolateAndNeumann(coords_in->g13, location);
  g23 = interpolateAndNeumann(coords_in->g23, location);

  // Check input metrics
  if ((!finite(g11, RGN_NOBNDRY)) || (!finite(g22, RGN_NOBNDRY)) || (!finite(g33, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Interpolated diagonal metrics are not finite!\n");
  }
  if ((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Interpolated diagonal metrics are negative!\n");
  }
  if ((!finite(g12, RGN_NOBNDRY)) || (!finite(g13, RGN_NOBNDRY)) || (!finite(g23, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Interpolated off-diagonal metrics are not finite!\n");
  }

  /// Always calculate contravariant metric components so that they are
  /// consistent with the interpolated covariant components
  if (calcCovariant()) {
    throw BoutException("Error in calcCovariant call");
  }

  /// Calculate Jacobian and Bxy
  if (jacobian())
    throw BoutException("Error in jacobian call");

  ShiftTorsion = interpolateAndNeumann(coords_in->ShiftTorsion, location);

  //////////////////////////////////////////////////////

  if (mesh->IncIntShear) {
    IntShiftTorsion = interpolateAndNeumann(coords_in->IntShiftTorsion, location);
  }
}

void Coordinates::outputVars(Datafile &file) {
  file.add(dx, "dx", false);
  file.add(dy, "dy", false);
  file.add(dz, "dz", false);

  file.add(g11, "g11", false);
  file.add(g22, "g22", false);
  file.add(g33, "g33", false);
  file.add(g12, "g12", false);
  file.add(g13, "g13", false);
  file.add(g23, "g23", false);

  file.add(g_11, "g_11", false);
  file.add(g_22, "g_22", false);
  file.add(g_33, "g_33", false);
  file.add(g_12, "g_12", false);
  file.add(g_13, "g_13", false);
  file.add(g_23, "g_23", false);

  file.add(J, "J", false);
}

int Coordinates::geometry() {
  TRACE("Coordinates::geometry");
  localmesh->communicate(dx, dy, g11, g22, g33, g12, g13, g23);
  localmesh->communicate(g_11, g_22, g_33, g_12, g_13, g_23, J, Bxy);

  output_progress.write("Calculating differential geometry terms\n");

  if (min(abs(dx)) < 1e-8)
    throw BoutException("dx magnitude less than 1e-8");

  if (min(abs(dy)) < 1e-8)
    throw BoutException("dy magnitude less than 1e-8");

  if (fabs(dz) < 1e-8)
    throw BoutException("dz magnitude less than 1e-8");

  // Check input metrics
  if ((!finite(g11, RGN_NOBNDRY)) || (!finite(g22, RGN_NOBNDRY)) || (!finite(g33, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if ((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if ((!finite(g12, RGN_NOBNDRY)) || (!finite(g13, RGN_NOBNDRY)) || (!finite(g23, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }

  if ((!finite(g_11, RGN_NOBNDRY)) || (!finite(g_22, RGN_NOBNDRY)) || (!finite(g_33, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Diagonal g_ij metrics are not finite!\n");
  }
  if ((min(g_11) <= 0.0) || (min(g_22) <= 0.0) || (min(g_33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal g_ij metrics are negative!\n");
  }
  if ((!finite(g_12, RGN_NOBNDRY)) || (!finite(g_13, RGN_NOBNDRY)) || (!finite(g_23, RGN_NOBNDRY))) {
    throw BoutException("\tERROR: Off-diagonal g_ij metrics are not finite!\n");
  }

  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11 = 0.5 * g11 * DDX(g_11) + g12 * (DDX(g_12) - 0.5 * DDY(g_11)) +
          g13 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G1_22 = g11 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g12 * DDY(g_22) +
          g13 * (DDY(g_23) - 0.5 * DDZ(g_22));
  G1_33 = g11 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g12 * (DDZ(g_23) - 0.5 * DDY(g_33)) +
          0.5 * g13 * DDZ(g_33);
  G1_12 = 0.5 * g11 * DDY(g_11) + 0.5 * g12 * DDX(g_22) +
          0.5 * g13 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G1_13 = 0.5 * g11 * DDZ(g_11) + 0.5 * g12 * (DDZ(g_12) + DDX(g_23) - DDY(g_13)) +
          0.5 * g13 * DDX(g_33);
  G1_23 = 0.5 * g11 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) +
          0.5 * g12 * (DDZ(g_22) + DDY(g_23) - DDY(g_23))
          // + 0.5 *g13*(DDZ(g_32) + DDY(g_33) - DDZ(g_23));
          // which equals
          + 0.5 * g13 * DDY(g_33);

  G2_11 = 0.5 * g12 * DDX(g_11) + g22 * (DDX(g_12) - 0.5 * DDY(g_11)) +
          g23 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G2_22 = g12 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g22 * DDY(g_22) +
          g23 * (DDY(g23) - 0.5 * DDZ(g_22));
  G2_33 = g12 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g22 * (DDZ(g_23) - 0.5 * DDY(g_33)) +
          0.5 * g23 * DDZ(g_33);
  G2_12 = 0.5 * g12 * DDY(g_11) + 0.5 * g22 * DDX(g_22) +
          0.5 * g23 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G2_13 =
      // 0.5 *g21*(DDZ(g_11) + DDX(g_13) - DDX(g_13))
      // which equals
      0.5 * g12 * (DDZ(g_11) + DDX(g_13) - DDX(g_13))
      // + 0.5 *g22*(DDZ(g_21) + DDX(g_23) - DDY(g_13))
      // which equals
      + 0.5 * g22 * (DDZ(g_12) + DDX(g_23) - DDY(g_13))
      // + 0.5 *g23*(DDZ(g_31) + DDX(g_33) - DDZ(g_13));
      // which equals
      + 0.5 * g23 * DDX(g_33);
  G2_23 = 0.5 * g12 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) + 0.5 * g22 * DDZ(g_22) +
          0.5 * g23 * DDY(g_33);

  G3_11 = 0.5 * g13 * DDX(g_11) + g23 * (DDX(g_12) - 0.5 * DDY(g_11)) +
          g33 * (DDX(g_13) - 0.5 * DDZ(g_11));
  G3_22 = g13 * (DDY(g_12) - 0.5 * DDX(g_22)) + 0.5 * g23 * DDY(g_22) +
          g33 * (DDY(g_23) - 0.5 * DDZ(g_22));
  G3_33 = g13 * (DDZ(g_13) - 0.5 * DDX(g_33)) + g23 * (DDZ(g_23) - 0.5 * DDY(g_33)) +
          0.5 * g33 * DDZ(g_33);
  G3_12 =
      // 0.5 *g31*(DDY(g_11) + DDX(g_12) - DDX(g_12))
      // which equals to
      0.5 * g13 * DDY(g_11)
      // + 0.5 *g32*(DDY(g_21) + DDX(g_22) - DDY(g_12))
      // which equals to
      + 0.5 * g23 * DDX(g_22)
      //+ 0.5 *g33*(DDY(g_31) + DDX(g_32) - DDZ(g_12));
      // which equals to
      + 0.5 * g33 * (DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G3_13 = 0.5 * g13 * DDZ(g_11) + 0.5 * g23 * (DDZ(g_12) + DDX(g_23) - DDY(g_13)) +
          0.5 * g33 * DDX(g_33);
  G3_23 = 0.5 * g13 * (DDZ(g_12) + DDY(g_13) - DDX(g_23)) + 0.5 * g23 * DDZ(g_22) +
          0.5 * g33 * DDY(g_33);

  auto tmp = J * g12;
  localmesh->communicate(tmp);
  G1 = (DDX(J * g11) + DDY(tmp) + DDZ(J * g13)) / J;
  tmp = J * g22;
  localmesh->communicate(tmp);
  G2 = (DDX(J * g12) + DDY(tmp) + DDZ(J * g23)) / J;
  tmp = J * g23;
  localmesh->communicate(tmp);
  G3 = (DDX(J * g13) + DDY(tmp) + DDZ(J * g33)) / J;

  // Communicate christoffel symbol terms
  output_progress.write("\tCommunicating connection terms\n");

  FieldGroup com;

  com.add(G1_11);
  com.add(G1_22);
  com.add(G1_33);
  com.add(G1_12);
  com.add(G1_13);
  com.add(G1_23);

  com.add(G2_11);
  com.add(G2_22);
  com.add(G2_33);
  com.add(G2_12);
  com.add(G2_13);
  com.add(G2_23);

  com.add(G3_11);
  com.add(G3_22);
  com.add(G3_33);
  com.add(G3_12);
  com.add(G3_13);
  com.add(G3_23);

  com.add(G1);
  com.add(G2);
  com.add(G3);

  localmesh->communicateXZ(com);

  //////////////////////////////////////////////////////
  /// Non-uniform meshes. Need to use DDX, DDY

  OPTION(Options::getRoot(), non_uniform, true);

  Coordinates::metric_field_type d2x, d2y; // d^2 x / d i^2
  // Read correction for non-uniform meshes
  if (localmesh->get(d2x, "d2x", 0.0, false)) {
    output_warn.write(
        "\tWARNING: differencing quantity 'd2x' not found. Calculating from dx\n");
    d1_dx = bout::derivatives::index::DDX(1. / dx); // d/di(1/dx)
  } else {
    localmesh->communicateXZ(d2x);
    d1_dx = -d2x / (dx * dx);
  }

  if (localmesh->get(d2y, "d2y", 0.0, false)) {
    output_warn.write(
        "\tWARNING: differencing quantity 'd2y' not found. Calculating from dy\n");
    auto tmp2 = 1. / dy;
    localmesh->communicate(tmp2);
    d1_dy = bout::derivatives::index::DDY(tmp2); // d/di(1/dy)
  } else {
    localmesh->communicateXZ(d2y);
    d1_dy = -d2y / (dy * dy);
  }
  localmesh->communicateXZ(d1_dx, d1_dy);

  return 0;
}

int Coordinates::calcCovariant() {
  TRACE("Coordinates::calcCovariant");

  // Make sure metric elements are allocated
  g_11.allocate();
  g_22.allocate();
  g_33.allocate();
  g_12.allocate();
  g_13.allocate();
  g_23.allocate();

  g_11.setLocation(location);
  g_22.setLocation(location);
  g_33.setLocation(location);
  g_12.setLocation(location);
  g_13.setLocation(location);
  g_23.setLocation(location);

  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  for (int jx = 0; jx < localmesh->LocalNx; jx++) {
    for (int jy = 0; jy < localmesh->LocalNy; jy++) {
#ifndef COORDINATES_USE_3D
      {
        int jz = 0;
#else
      for (int jz = 0; jz < localmesh->LocalNz; jz++) { // Inefficient for 2D metric type
#endif
        // set elements of g
        a(0, 0) = g11(jx, jy, jz);
        a(1, 1) = g22(jx, jy, jz);
        a(2, 2) = g33(jx, jy, jz);

        a(0, 1) = a(1, 0) = g12(jx, jy, jz);
        a(1, 2) = a(2, 1) = g23(jx, jy, jz);
        a(0, 2) = a(2, 0) = g13(jx, jy, jz);

        // invert
        if (invert3x3(a)) {
          output_error.write("\tERROR: metric tensor is singular at (%d, %d, %d)\n", jx,
                             jy, jz);
          return 1;
        }

        // put elements into g_{ij}
        g_11(jx, jy, jz) = a(0, 0);
        g_22(jx, jy, jz) = a(1, 1);
        g_33(jx, jy, jz) = a(2, 2);

        g_12(jx, jy, jz) = a(0, 1);
        g_13(jx, jy, jz) = a(0, 2);
        g_23(jx, jy, jz) = a(1, 2);
      }
    }
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tLocal maximum error in diagonal inversion is %e\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tLocal maximum error in off-diagonal inversion is %e\n", maxerr);

  return 0;
}

int Coordinates::calcContravariant() {
  TRACE("Coordinates::calcContravariant");

  // Make sure metric elements are allocated
  g11.allocate();
  g22.allocate();
  g33.allocate();
  g12.allocate();
  g13.allocate();
  g23.allocate();

  // Perform inversion of g_{ij} to get g^{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  auto a = Matrix<BoutReal>(3, 3);

  for (int jx = 0; jx < localmesh->LocalNx; jx++) {
    for (int jy = 0; jy < localmesh->LocalNy; jy++) {
#ifndef COORDINATES_USE_3D
      {
        int jz = 0;
#else
      for (int jz = 0; jz < localmesh->LocalNz; jz++) { // Inefficient for 2D metric type
#endif
        // set elements of g
        a(0, 0) = g_11(jx, jy, jz);
        a(1, 1) = g_22(jx, jy, jz);
        a(2, 2) = g_33(jx, jy, jz);

        a(0, 1) = a(1, 0) = g_12(jx, jy, jz);
        a(1, 2) = a(2, 1) = g_23(jx, jy, jz);
        a(0, 2) = a(2, 0) = g_13(jx, jy, jz);

        // invert
        if (invert3x3(a)) {
          output_error.write("\tERROR: metric tensor is singular at (%d, %d, %d)\n", jx,
                             jy, jz);
          return 1;
        }

        // put elements into g_{ij}
        g11(jx, jy, jz) = a(0, 0);
        g22(jx, jy, jz) = a(1, 1);
        g33(jx, jy, jz) = a(2, 2);

        g12(jx, jy, jz) = a(0, 1);
        g13(jx, jy, jz) = a(0, 2);
        g23(jx, jy, jz) = a(1, 2);
      }
    }
  }

  BoutReal maxerr;
  maxerr = BOUTMAX(max(abs((g_11 * g11 + g_12 * g12 + g_13 * g13) - 1)),
                   max(abs((g_12 * g12 + g_22 * g22 + g_23 * g23) - 1)),
                   max(abs((g_13 * g13 + g_23 * g23 + g_33 * g33) - 1)));

  output_info.write("\tMaximum error in diagonal inversion is %e\n", maxerr);

  maxerr = BOUTMAX(max(abs(g_11 * g12 + g_12 * g22 + g_13 * g23)),
                   max(abs(g_11 * g13 + g_12 * g23 + g_13 * g33)),
                   max(abs(g_12 * g13 + g_22 * g23 + g_23 * g33)));

  output_info.write("\tMaximum error in off-diagonal inversion is %e\n", maxerr);
  return 0;
}

int Coordinates::jacobian() {
  TRACE("Coordinates::jacobian");
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)

  auto g = g11 * g22 * g33 + 2.0 * g12 * g13 * g23 - g11 * g23 * g23 - g22 * g13 * g13
           - g33 * g12 * g12;

  // Check that g is positive
  if (min(g) < 0.0) {
    throw BoutException("The determinant of g^ij is somewhere less than 0.0");
  }
  J = 1. / sqrt(g);

  // Check jacobian
  if (!finite(J, RGN_NOBNDRY)) {
    throw BoutException("\tERROR: Jacobian not finite everywhere!\n");
  }
  if (min(abs(J)) < 1.0e-10) {
    throw BoutException("\tERROR: Jacobian becomes very small\n");
  }

  if (min(g_22) < 0.0) {
    throw BoutException("g_22 is somewhere less than 0.0");
  }
  Bxy = sqrt(g_22) / J;

  return 0;
}

/*******************************************************************************
 * Operators
 *
 *******************************************************************************/

const Coordinates::metric_field_type Coordinates::DDX(const Field2D& f, CELL_LOC loc,
                                                      const std::string& method,
                                                      REGION region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  return bout::derivatives::index::DDX(f, loc, method, region) / dx;
}
const Field3D Coordinates::DDX(const Field3D& f, CELL_LOC outloc,
                               const std::string& method, REGION region) {

  auto result = bout::derivatives::index::DDX(f, outloc, method, region);
  result /= dx;

  if (f.getMesh()->IncIntShear) {
    // Using BOUT-06 style shifting
    result += IntShiftTorsion * DDZ(f, outloc, method, region);
  }

  return result;
};

const Coordinates::metric_field_type Coordinates::DDY(const Field2D& f, CELL_LOC loc,
                                                      const std::string& method,
                                                      REGION region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  return bout::derivatives::index::DDY(f, loc, method, region) / dy;
}
const Field3D Coordinates::DDY(const Field3D& f, CELL_LOC outloc,
                               const std::string& method, REGION region) {
  return bout::derivatives::index::DDY(f, outloc, method, region) / dy;
};

const Coordinates::metric_field_type Coordinates::DDZ(MAYBE_UNUSED(const Field2D& f),
                                                      CELL_LOC loc,
                                                      const std::string& UNUSED(method),
                                                      REGION UNUSED(region)) {
  ASSERT1(location == loc || loc == CELL_DEFAULT);
  ASSERT1(f.getMesh() == localmesh);
  auto result = Coordinates::metric_field_type(0.0, localmesh);
  result.setLocation(location);
  return result;
}
const Field3D Coordinates::DDZ(const Field3D& f, CELL_LOC outloc,
                               const std::string& method, REGION region) {
  return bout::derivatives::index::DDZ(f, outloc, method, region) / dz;
};

/////////////////////////////////////////////////////////
// Parallel gradient

const Coordinates::metric_field_type
Coordinates::Grad_par(const Field2D& var, MAYBE_UNUSED(CELL_LOC outloc),
                      const std::string& UNUSED(method)) {
  TRACE("Coordinates::Grad_par( Field2D )");
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == var.getLocation()));

  return DDY(var) / sqrt(g_22);
}

const Field3D Coordinates::Grad_par(const Field3D &var, CELL_LOC outloc,
                                    const std::string &method) {
  TRACE("Coordinates::Grad_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  return ::DDY(var, outloc, method) / sqrt(g_22);
}

/////////////////////////////////////////////////////////
// Vpar_Grad_par
// vparallel times the parallel derivative along unperturbed B-field

const Coordinates::metric_field_type
Coordinates::Vpar_Grad_par(const Field2D& v, const Field2D& f,
                           MAYBE_UNUSED(CELL_LOC outloc),
                           const std::string& UNUSED(method)) {
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()));
  return VDDY(v, f) / sqrt(g_22);
}

const Field3D Coordinates::Vpar_Grad_par(const Field3D &v, const Field3D &f, CELL_LOC outloc,
                                         const std::string &method) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return VDDY(v, f, outloc, method) / sqrt(g_22);
}

/////////////////////////////////////////////////////////
// Parallel divergence

const Coordinates::metric_field_type
Coordinates::Div_par(const Field2D& f, CELL_LOC outloc, const std::string& method) {
  TRACE("Coordinates::Div_par( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy;

  return Bxy * Grad_par(f / Bxy_floc, outloc, method);
}

const Field3D Coordinates::Div_par(const Field3D &f, CELL_LOC outloc,
                                   const std::string &method) {
  TRACE("Coordinates::Div_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  
  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy;

  if (!f.hasYupYdown()) {
    // No yup/ydown fields. The Grad_par operator will
    // shift to field aligned coordinates
    return Bxy * Grad_par(f / Bxy_floc, outloc, method);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  if (&f.yup() == &f) {
    // Identity, yup and ydown point to same field
    f_B.mergeYupYdown();
  } else {
    // Distinct fields
    f_B.splitYupYdown();
    f_B.yup() = f.yup() / Bxy_floc;
    f_B.ydown() = f.ydown() / Bxy_floc;
  }
  return Bxy * Grad_par(f_B, outloc, method);
}

/////////////////////////////////////////////////////////
// second parallel derivative (b dot Grad)(b dot Grad)
// Note: For parallel Laplacian use Laplace_par

const Coordinates::metric_field_type
Coordinates::Grad2_par2(const Field2D& f, CELL_LOC outloc, const std::string& method) {
  TRACE("Coordinates::Grad2_par2( Field2D )");
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()));

  auto invSg = 1.0 / sqrt(g_22);
  localmesh->communicate(invSg);
  auto result = DDY(invSg, outloc, method) * DDY(f, outloc, method) * invSg
                + D2DY2(f, outloc, method) / g_22;

  return result;
}

const Field3D Coordinates::Grad2_par2(const Field3D &f, CELL_LOC outloc, const std::string &method) {
  TRACE("Coordinates::Grad2_par2( Field3D )");
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  ASSERT1(location == outloc);

  auto sg = sqrt(g_22);
  auto invSg = 1.0 / sg;
  localmesh->communicate(invSg);
  sg = DDY(invSg, outloc, method) * invSg;

  Field3D result(localmesh), r2(localmesh);

  result = ::DDY(f, outloc, method);

  r2 = D2DY2(f, outloc, method) / g_22;

  result = sg * result + r2;

  ASSERT2(result.getLocation() == outloc);

  return result;
}

/////////////////////////////////////////////////////////
// perpendicular Laplacian operator

#include <invert_laplace.hxx> // Delp2 uses same coefficients as inversion code

const Coordinates::metric_field_type Coordinates::Delp2(const Field2D& f,
                                                        CELL_LOC outloc) {
  TRACE("Coordinates::Delp2( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  auto result = G1 * DDX(f, outloc) + g11 * D2DX2(f, outloc);

  return result;
}

const Field3D Coordinates::Delp2(const Field3D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Delp2( Field3D )");
#ifndef COORDINATES_USE_3D
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  ASSERT1(location == outloc);

  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
    // copy mesh, location, etc
    return f*0;
  }
  ASSERT2(localmesh->xstart > 0); // Need at least one guard cell

  ASSERT2(f.getLocation() == outloc);

  Field3D result(localmesh);
  result.allocate();
  result.setLocation(f.getLocation());

  int ncz = localmesh->LocalNz;

  // Allocate memory
  auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
  auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);

  // Loop over all y indices
  for (int jy = 0; jy < localmesh->LocalNy; jy++) {

    // Take forward FFT

    for (int jx = 0; jx < localmesh->LocalNx; jx++)
      rfft(&f(jx, jy, 0), ncz, &ft(jx, 0));

    // Loop over kz
    for (int jz = 0; jz <= ncz / 2; jz++) {
      dcomplex a, b, c;

      // No smoothing in the x direction
      for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
        // Perform x derivative

        laplace_tridag_coefs(jx, jy, jz, a, b, c, nullptr, nullptr, outloc);

        delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
      }
    }

    // Reverse FFT
    for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {

      irfft(&delft(jx, 0), ncz, &result(jx, jy, 0));
    }

    // Boundaries
    for (int jz = 0; jz < ncz; jz++) {
      for (int jx = 0; jx < localmesh->xstart; jx++) {
        result(jx, jy, jz) = 0.0;
      }
      for (int jx = localmesh->xend + 1; jx < localmesh->LocalNx; jx++) {
        result(jx, jy, jz) = 0.0;
      }
    }
  }

  ASSERT2(result.getLocation() == f.getLocation());

  return result;
#else
  throw BoutException("Delp2(Field3D) currently relies on the tridiagonal code that is "
                      "not supported with 3D metrics.");
#endif
}

const FieldPerp Coordinates::Delp2(const FieldPerp &f, CELL_LOC outloc) {
  TRACE("Coordinates::Delp2( FieldPerp )");

  if (outloc == CELL_DEFAULT) outloc = f.getLocation();

  ASSERT1(location == outloc);
  ASSERT2(f.getLocation() == outloc);

  FieldPerp result(localmesh);
  result.allocate();
  result.setLocation(outloc);

  int jy = f.getIndex();
  result.setIndex(jy);

  int ncz = localmesh->LocalNz;

  // Allocate memory
  auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
  auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);

  // Take forward FFT
  for (int jx = 0; jx < localmesh->LocalNx; jx++)
    rfft(&f(jx, 0), ncz, &ft(jx, 0));

  // Loop over kz
  for (int jz = 0; jz <= ncz / 2; jz++) {

    // No smoothing in the x direction
    for (int jx = 2; jx < (localmesh->LocalNx - 2); jx++) {
      // Perform x derivative

      dcomplex a, b, c;
      laplace_tridag_coefs(jx, jy, jz, a, b, c);

      delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
    }
  }

  // Reverse FFT
  for (int jx = 1; jx < (localmesh->LocalNx - 1); jx++) {
    irfft(&delft(jx, 0), ncz, &result(jx, 0));
  }

  // Boundaries
  for (int jz = 0; jz < ncz; jz++) {
    result(0, jz) = 0.0;
    result(localmesh->LocalNx - 1, jz) = 0.0;
  }

  return result;
}

const Coordinates::metric_field_type Coordinates::Laplace_par(const Field2D& f,
                                                              CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return D2DY2(f, outloc) / g_22 + DDY(J / g_22, outloc) * DDY(f, outloc) / J;
}

const Field3D Coordinates::Laplace_par(const Field3D &f, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);
  return D2DY2(f, outloc) / g_22 + DDY(J / g_22, outloc) * ::DDY(f, outloc) / J;
}

// Full Laplacian operator on scalar field

const Coordinates::metric_field_type Coordinates::Laplace(const Field2D& f,
                                                          CELL_LOC outloc) {
  TRACE("Coordinates::Laplace( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  auto result = G1 * DDX(f, outloc) + G2 * DDY(f, outloc) + g11 * D2DX2(f, outloc)
                + g22 * D2DY2(f, outloc) + 2.0 * g12 * D2DXDY(f, outloc);

  return result;
}

const Field3D Coordinates::Laplace(const Field3D &f, CELL_LOC outloc) {
  TRACE("Coordinates::Laplace( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT);

  Field3D result = G1 * ::DDX(f, outloc) + G2 * ::DDY(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * D2DX2(f, outloc) +
                   g22 * D2DY2(f, outloc) + g33 * D2DZ2(f, outloc) +
                   2.0 * (g12 * D2DXDY(f, outloc) + g13 * D2DXDZ(f, outloc) + g23 * D2DYDZ(f, outloc));

  ASSERT2(result.getLocation() == f.getLocation());

  return result;
}

 bool Coordinates::is3D() {
#ifdef COORDINATES_USE_3D
  return true;
#else
  return false;
#endif
}
