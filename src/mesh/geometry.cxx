/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

//#include <bout/assert.hxx>
//#include <bout/constants.hxx>
#include <bout/coordinates.hxx>
#include <bout/geometry.hxx>
//#include <bout/msg_stack.hxx>
//#include <bout/sys/timer.hxx>
//#include <bout/utils.hxx>
//
//#include <bout/derivs.hxx>
//#include <bout/fft.hxx>
//#include <bout/interpolation.hxx>
//
//#include <bout/globals.hxx>
//
//#include "parallel/fci.hxx"
//#include "parallel/shiftedmetricinterp.hxx"

//// use anonymous namespace so this utility function is not available outside this file
//namespace {

//template <typename T, typename... Ts>
//// Use sendY()/sendX() and wait() instead of Mesh::communicate() to ensure we
//// don't try to calculate parallel slices as Coordinates are not constructed yet
//void communicate(T& t, Ts... ts) {
//  coordinates.communicate(t, ts);
//}

///// Interpolate a Field2D to a new CELL_LOC with interp_to.
///// Communicates to set internal guard cells.
///// Boundary guard cells are set by extrapolating from the grid, like
///// 'free_o3' boundary conditions
///// Corner guard cells are set to BoutNaN
//const Field2D interpolateAndExtrapolate(const Field2D& f, CELL_LOC location,
//                                        bool extrapolate_x, bool extrapolate_y,
//                                        bool no_extra_interpolate,
//                                        ParallelTransform* UNUSED(pt) = nullptr,
//                                        const std::string& region = "RGN_NOBNDRY") {
//
//  Mesh* localmesh = f.getMesh();
//  Field2D result = interp_to(f, location, region);
//  // Ensure result's data is unique. Otherwise result might be a duplicate of
//  // f (if no interpolation is needed, e.g. if interpolation is in the
//  // z-direction); then f would be communicated. Since this function is used
//  // on geometrical quantities that might not be periodic in y even on closed
//  // field lines (due to dependence on integrated shear), we don't want to
//  // communicate f. We will sort out result's boundary guard cells below, but
//  // not f's so we don't want to change f.
//  result.allocate();
//  communicate(result);
//
//  // Extrapolate into boundaries (if requested) so that differential calculateGeometry
//  // terms can be interpolated if necessary
//  // Note: cannot use applyBoundary("free_o3") here because applyBoundary()
//  // would try to create a new Coordinates object since we have not finished
//  // initializing yet, leading to an infinite recursion.
//  // Also, here we interpolate for the boundary points at xstart/ystart and
//  // (xend+1)/(yend+1) instead of extrapolating.
//  for (auto& bndry : localmesh->getBoundaries()) {
//    if ((extrapolate_x and bndry->bx != 0) or (extrapolate_y and bndry->by != 0)) {
//      int extrap_start = 0;
//      if (not no_extra_interpolate) {
//        // Can use no_extra_interpolate argument to skip the extra interpolation when we
//        // want to extrapolate the Christoffel symbol terms which come from derivatives so
//        // don't have the extra point set already
//        if ((location == CELL_XLOW) && (bndry->bx > 0)) {
//          extrap_start = 1;
//        } else if ((location == CELL_YLOW) && (bndry->by > 0)) {
//          extrap_start = 1;
//        }
//      }
//      for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
//        // interpolate extra boundary point that is missed by interp_to, if
//        // necessary.
//        // Only interpolate this point if we are actually changing location. E.g.
//        // when we use this function to extrapolate J and Bxy on staggered grids,
//        // this point should already be set correctly because the metric
//        // components have been interpolated to here.
//        if (extrap_start > 0 and f.getLocation() != location) {
//          ASSERT1(bndry->bx == 0 or localmesh->xstart > 1)
//          ASSERT1(bndry->by == 0 or localmesh->ystart > 1)
//          // note that either bx or by is >0 here
//          result(bndry->x, bndry->y) =
//              (9.
//                   * (f(bndry->x - bndry->bx, bndry->y - bndry->by)
//                      + f(bndry->x, bndry->y))
//               - f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by)
//               - f(bndry->x + bndry->bx, bndry->y + bndry->by))
//              / 16.;
//        }
//
//        // set boundary guard cells
//        if ((bndry->bx != 0 && localmesh->GlobalNx - 2 * bndry->width >= 3)
//            || (bndry->by != 0
//                && localmesh->GlobalNy - localmesh->numberOfYBoundaries() * bndry->width
//                       >= 3)) {
//          if (bndry->bx != 0 && localmesh->LocalNx == 1 && bndry->width == 1) {
//            throw BoutException(
//                "Not enough points in the x-direction on this "
//                "processor for extrapolation needed to use staggered grids. "
//                "Increase number of x-guard cells MXG or decrease number of "
//                "processors in the x-direction NXPE.");
//          }
//          if (bndry->by != 0 && localmesh->LocalNy == 1 && bndry->width == 1) {
//            throw BoutException(
//                "Not enough points in the y-direction on this "
//                "processor for extrapolation needed to use staggered grids. "
//                "Increase number of y-guard cells MYG or decrease number of "
//                "processors in the y-direction NYPE.");
//          }
//          // extrapolate into boundary guard cells if there are enough grid points
//          for (int i = extrap_start; i < bndry->width; i++) {
//            int xi = bndry->x + i * bndry->bx;
//            int yi = bndry->y + i * bndry->by;
//            result(xi, yi) = 3.0 * result(xi - bndry->bx, yi - bndry->by)
//                             - 3.0 * result(xi - 2 * bndry->bx, yi - 2 * bndry->by)
//                             + result(xi - 3 * bndry->bx, yi - 3 * bndry->by);
//          }
//        } else {
//          // not enough grid points to extrapolate, set equal to last grid point
//          for (int i = extrap_start; i < bndry->width; i++) {
//            result(bndry->x + i * bndry->bx, bndry->y + i * bndry->by) =
//                result(bndry->x - bndry->bx, bndry->y - bndry->by);
//          }
//        }
//      }
//    }
//  }
//#if CHECK > 0
//  if (not(
//          // if include_corner_cells=true, then we extrapolate valid data into the
//          // corner cells if they are not already filled
//          localmesh->include_corner_cells
//
//          // if we are not extrapolating at all, the corner cells should contain valid
//          // data
//          or (not extrapolate_x and not extrapolate_y))) {
//    // Invalidate corner guard cells
//    for (int i = 0; i < localmesh->xstart; i++) {
//      for (int j = 0; j < localmesh->ystart; j++) {
//        result(i, j) = BoutNaN;
//        result(i, localmesh->LocalNy - 1 - j) = BoutNaN;
//        result(localmesh->LocalNx - 1 - i, j) = BoutNaN;
//        result(localmesh->LocalNx - 1 - i, localmesh->LocalNy - 1 - j) = BoutNaN;
//      }
//    }
//  }
//#endif
//
//  return result;
//}
/*

#if BOUT_USE_METRIC_3D
Field3D interpolateAndExtrapolate(const Field3D& f_, CELL_LOC location,
                                  bool extrapolate_x, bool extrapolate_y,
                                  bool no_extra_interpolate, ParallelTransform* pt_) {

  Mesh* localmesh = f_.getMesh();
  Field3D result;
  Field3D f = f_;
  ParallelTransform* pt_f;
  if (f.getCoordinates() == nullptr) {
    // if input f is member of the Coordinates we are currently constructing, it will not
    // have Coordinates and needs to use the passed-in ParallelTransform
    pt_f = pt_;
  } else {
    // if input f is from Coordinates at a different location, it will have its own
    // Coordinates, and we should use its ParallelTransform
    pt_f = &f.getCoordinates()->getParallelTransform();
  }
  if (f.getDirectionY() != YDirectionType::Standard) {
    if (pt_f->canToFromFieldAligned()) {
      f = pt_f->fromFieldAligned(f);
    } else {
      f.setDirectionY(YDirectionType::Standard);
    }
  }
  if (location == CELL_YLOW and f.getLocation() != CELL_YLOW) {
    auto f_aligned = pt_f->toFieldAligned(f, "RGN_NOX");
    result = interp_to(f_aligned, location, "RGN_NOBNDRY");
    ParallelTransform* pt_result;
    if (result.getCoordinates() == nullptr) {
      pt_result = pt_;
    } else {
      pt_result = &result.getCoordinates()->getParallelTransform();
    }
    result = pt_result->fromFieldAligned(result, "RGN_NOBNDRY");
  } else {
    result = interp_to(f, location, "RGN_NOBNDRY");
  }
  // Ensure result's data is unique. Otherwise result might be a duplicate of
  // f (if no interpolation is needed, e.g. if interpolation is in the
  // z-direction); then f would be communicated. Since this function is used
  // on geometrical quantities that might not be periodic in y even on closed
  // field lines (due to dependence on integrated shear), we don't want to
  // communicate f. We will sort out result's boundary guard cells below, but
  // not f's so we don't want to change f.
  result.allocate();
  communicate(result);

  // Extrapolate into boundaries (if requested) so that differential geometry
  // terms can be interpolated if necessary
  // Note: cannot use applyBoundary("free_o3") here because applyBoundary()
  // would try to create a new Coordinates object since we have not finished
  // initializing yet, leading to an infinite recursion.
  // Also, here we interpolate for the boundary points at xstart/ystart and
  // (xend+1)/(yend+1) instead of extrapolating.
  for (auto& bndry : localmesh->getBoundaries()) {
    if ((extrapolate_x and bndry->bx != 0) or (extrapolate_y and bndry->by != 0)) {
      int extrap_start = 0;
      if (not no_extra_interpolate) {
        // Can use no_extra_interpolate argument to skip the extra interpolation when we
        // want to extrapolate the Christoffel symbol terms which come from derivatives so
        // don't have the extra point set already
        if ((location == CELL_XLOW) && (bndry->bx > 0)) {
          extrap_start = 1;
        } else if ((location == CELL_YLOW) && (bndry->by > 0)) {
          extrap_start = 1;
        }
      }
      for (bndry->first(); !bndry->isDone(); bndry->next1d()) {
        // interpolate extra boundary point that is missed by interp_to, if
        // necessary.
        // Only interpolate this point if we are actually changing location. E.g.
        // when we use this function to extrapolate J and Bxy on staggered grids,
        // this point should already be set correctly because the metric
        // components have been interpolated to here.
        if (extrap_start > 0 and f.getLocation() != location) {
          ASSERT1(bndry->bx == 0 or localmesh->xstart > 1);
          ASSERT1(bndry->by == 0 or localmesh->ystart > 1);
          // note that either bx or by is >0 here
          for (int zi = 0; zi < localmesh->LocalNz; ++zi) {
            result(bndry->x, bndry->y, zi) =
                (9.
                     * (f(bndry->x - bndry->bx, bndry->y - bndry->by, zi)
                        + f(bndry->x, bndry->y, zi))
                 - f(bndry->x - 2 * bndry->bx, bndry->y - 2 * bndry->by, zi)
                 - f(bndry->x + bndry->bx, bndry->y + bndry->by, zi))
                / 16.;
          }
        }
        // set boundary guard cells
        if ((bndry->bx != 0 && localmesh->GlobalNx - 2 * bndry->width >= 3)
            || (bndry->by != 0 && localmesh->GlobalNy - 2 * bndry->width >= 3)) {
          if (bndry->bx != 0 && localmesh->LocalNx == 1 && bndry->width == 1) {
            throw BoutException(
                "Not enough points in the x-direction on this "
                "processor for extrapolation needed to use staggered grids. "
                "Increase number of x-guard cells MXG or decrease number of "
                "processors in the x-direction NXPE.");
          }
          if (bndry->by != 0 && localmesh->LocalNy == 1 && bndry->width == 1) {
            throw BoutException(
                "Not enough points in the y-direction on this "
                "processor for extrapolation needed to use staggered grids. "
                "Increase number of y-guard cells MYG or decrease number of "
                "processors in the y-direction NYPE.");
          }
          // extrapolate into boundary guard cells if there are enough grid points
          for (int i = extrap_start; i < bndry->width; i++) {
            int xi = bndry->x + i * bndry->bx;
            int yi = bndry->y + i * bndry->by;
            for (int zi = 0; zi < localmesh->LocalNz; ++zi) {
              result(xi, yi, zi) =
                  3.0 * result(xi - bndry->bx, yi - bndry->by, zi)
                  - 3.0 * result(xi - 2 * bndry->bx, yi - 2 * bndry->by, zi)
                  + result(xi - 3 * bndry->bx, yi - 3 * bndry->by, zi);
            }
          }
        } else {
          // not enough grid points to extrapolate, set equal to last grid point
          for (int i = extrap_start; i < bndry->width; i++) {
            for (int zi = 0; zi < localmesh->LocalNz; ++zi) {
              result(bndry->x + i * bndry->bx, bndry->y + i * bndry->by, zi) =
                  result(bndry->x - bndry->bx, bndry->y - bndry->by, zi);
            }
          }
        }
      }
    }
  }
#if CHECK > 0
  if (not(
          // if include_corner_cells=true, then we extrapolate valid data into the
          // corner cells if they are not already filled
          localmesh->include_corner_cells

          // if we are not extrapolating at all, the corner cells should contain valid
          // data
          or (not extrapolate_x and not extrapolate_y))) {
    // Invalidate corner guard cells
    for (int i = 0; i < localmesh->xstart; i++) {
      for (int j = 0; j < localmesh->ystart; j++) {
        for (int k = 0; k < localmesh->LocalNz; ++k) {
          result(i, j, k) = BoutNaN;
          result(i, localmesh->LocalNy - 1 - j, k) = BoutNaN;
          result(localmesh->LocalNx - 1 - i, j, k) = BoutNaN;
          result(localmesh->LocalNx - 1 - i, localmesh->LocalNy - 1 - j, k) = BoutNaN;
        }
      }
    }
  }
#endif // CHECK > 0

  return result;
}
#endif // BOUT_USE_METRIC_3D
*/
///*
//
//// If the CELL_CENTRE variable was read, the staggered version is required to
//// also exist for consistency
//void checkStaggeredGet(Mesh* mesh, const std::string& name, const std::string& suffix) {
//  if (mesh->sourceHasVar(name) != mesh->sourceHasVar(name + suffix)) {
//    throw BoutException("Attempting to read staggered fields from grid, but " + name
//                        + " is not present in both CELL_CENTRE and staggered versions.");
//  }
//}
//
//// convenience function for repeated code
//int getAtLoc(Mesh* mesh, Coordinates::FieldMetric& var, const std::string& name,
//             const std::string& suffix, CELL_LOC location, BoutReal default_value = 0.) {
//
//  checkStaggeredGet(mesh, name, suffix);
//  return mesh->get(var, name + suffix, default_value, false, location);
//}
//
//auto getAtLoc(Mesh* mesh, const std::string& name, const std::string& suffix,
//              CELL_LOC location, BoutReal default_value = 0.) {
//
//  checkStaggeredGet(mesh, name, suffix);
//  return mesh->get(name + suffix, default_value, false, location);
//}
//
//std::string getLocationSuffix(CELL_LOC location) {
//  switch (location) {
//  case CELL_CENTRE: {
//    return "";
//  }
//  case CELL_XLOW: {
//    return "_xlow";
//  }
//  case CELL_YLOW: {
//    return "_ylow";
//  }
//  case CELL_ZLOW: {
//    // in 2D metric, same as CELL_CENTRE
//    return bout::build::use_metric_3d ? "_zlow" : "";
//  }
//  default: {
//    throw BoutException(
//        "Incorrect location passed to "
//        "Coordinates(Mesh*,const CELL_LOC,const Coordinates*) constructor.");
//  }
//  }
//}
//
//} // anonymous namespace
//*/
//
//Coordinates::FieldMetric Coordinates::getAtLocOrUnaligned(Mesh* mesh,
//                                                          const std::string& name,
//                                                          BoutReal default_value,
//                                                          const std::string& suffix,
//                                                          CELL_LOC cell_location) {
//  if (cell_location == CELL_CENTRE) {
//    return getUnaligned(name, default_value);
//  }
//  // grid data source has staggered fields, so read instead of interpolating
//  // Diagonal components of metric tensor g^{ij} (default to 1)
//  return getAtLoc(mesh, name, suffix, cell_location, default_value);
//}
//
//Coordinates::FieldMetric Coordinates::getUnaligned(const std::string& name,
//                                                   BoutReal default_value) {
//
//  auto field = localmesh->get(name, default_value, false);
//  if (field.getDirectionY() == YDirectionType::Aligned
//      and transform->canToFromFieldAligned()) {
//    return transform->fromFieldAligned(field);
//  } else {
//    field.setDirectionY(YDirectionType::Standard);
//    return field;
//  }
//}
//
//Coordinates::FieldMetric Coordinates::getUnalignedAtLocationAndFillGuards(
//    Mesh* mesh, const std::string& name, BoutReal default_value,
//    const std::string& suffix, CELL_LOC cell_location, bool extrapolate_x,
//    bool extrapolate_y, bool no_extra_interpolate,
//    ParallelTransform* pParallelTransform) {
//
//  auto field = getAtLocOrUnaligned(mesh, name, default_value, suffix, cell_location);
//  if (suffix == "") {
//    no_extra_interpolate = false;
//    pParallelTransform = transform.get();
//  }
//  return interpolateAndExtrapolate(field, cell_location, extrapolate_x, extrapolate_y,
//                                   no_extra_interpolate, pParallelTransform);
//}

//Geometry::Geometry() {}

Geometry::Geometry(
    //    Coordinates& coordinates,
    FieldMetric J, FieldMetric Bxy, FieldMetric g11, FieldMetric g22, FieldMetric g33,
    FieldMetric g12, FieldMetric g13, FieldMetric g23, FieldMetric g_11, FieldMetric g_22,
    FieldMetric g_33, FieldMetric g_12, FieldMetric g_13, FieldMetric g_23)
    : //      coordinates(coordinates),
      contravariantMetricTensor(g11, g22, g33, g12, g13, g23),
      covariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23), this_J(J),
      this_Bxy(Bxy) {}

Geometry::Geometry(Mesh* mesh, const CELL_LOC cell_location)
    //bool force_interpolate_from_centre)
    : G1_11(mesh), G1_22(mesh), G1_33(mesh), G1_12(mesh), G1_13(mesh), G1_23(mesh),
      G2_11(mesh), G2_22(mesh), G2_33(mesh), G2_12(mesh), G2_13(mesh), G2_23(mesh),
      G3_11(mesh), G3_22(mesh), G3_33(mesh), G3_12(mesh), G3_13(mesh), G3_23(mesh),
      G1(mesh), G2(mesh), G3(mesh), contravariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      covariantMetricTensor(1., 1., 1., 0, 0, 0, mesh),
      // Identity metric tensor
      this_J(1., mesh), this_Bxy(1., mesh) {}

//
//void Coordinates::outputVars(Options& output_options) {
//  Timer time("io");
//  const std::string loc_string =
//      (location == CELL_CENTRE) ? "" : "_" + toString(location);
//
//  output_options["dx" + loc_string].force(dx, "Coordinates");
//  output_options["dy" + loc_string].force(dy, "Coordinates");
//  output_options["dz" + loc_string].force(dz, "Coordinates");
//
//  output_options["g11" + loc_string].force(contravariantMetricTensor.Getg11(),
//                                           "Coordinates");
//  output_options["g22" + loc_string].force(contravariantMetricTensor.Getg22(),
//                                           "Coordinates");
//  output_options["g33" + loc_string].force(contravariantMetricTensor.Getg33(),
//                                           "Coordinates");
//  output_options["g12" + loc_string].force(contravariantMetricTensor.Getg12(),
//                                           "Coordinates");
//  output_options["g13" + loc_string].force(contravariantMetricTensor.Getg13(),
//                                           "Coordinates");
//  output_options["g23" + loc_string].force(contravariantMetricTensor.Getg23(),
//                                           "Coordinates");
//
//  output_options["g_11" + loc_string].force(covariantMetricTensor.Getg11(),
//                                            "Coordinates");
//  output_options["g_22" + loc_string].force(covariantMetricTensor.Getg22(),
//                                            "Coordinates");
//  output_options["g_33" + loc_string].force(covariantMetricTensor.Getg33(),
//                                            "Coordinates");
//  output_options["g_12" + loc_string].force(covariantMetricTensor.Getg12(),
//                                            "Coordinates");
//  output_options["g_13" + loc_string].force(covariantMetricTensor.Getg13(),
//                                            "Coordinates");
//  output_options["g_23" + loc_string].force(covariantMetricTensor.Getg23(),
//                                            "Coordinates");
//
//  output_options["J" + loc_string].force(this_J, "Coordinates");
//  output_options["Bxy" + loc_string].force(this_Bxy, "Coordinates");
//
//  output_options["G1" + loc_string].force(G1, "Coordinates");
//  output_options["G2" + loc_string].force(G2, "Coordinates");
//  output_options["G3" + loc_string].force(G3, "Coordinates");
//
//  getParallelTransform().outputVars(output_options);
//}

//const Field2D& Coordinates::zlength() const {
//  BOUT_OMP(critical)
//  if (not zlength_cache) {
//    zlength_cache = std::make_unique<Field2D>(0., localmesh);
//
//#if BOUT_USE_METRIC_3D
//    BOUT_FOR_SERIAL(i, dz.getRegion("RGN_ALL")) { (*zlength_cache)[i] += dz[i]; }
//#else
//    (*zlength_cache) = dz * nz;
//#endif
//  }
//
//  return *zlength_cache;
//}

//int Geometry::calculateGeometry(bool recalculate_staggered, bool force_interpolate_from_centre) {
//  TRACE("Geometry::calculateGeometry");
//
//  communicate(dx, dy, dz, contravariantMetricTensor.Getg11(),
//              contravariantMetricTensor.Getg22(), contravariantMetricTensor.Getg33(),
//              contravariantMetricTensor.Getg12(), contravariantMetricTensor.Getg13(),
//              contravariantMetricTensor.Getg23(), covariantMetricTensor.Getg11(),
//              covariantMetricTensor.Getg22(), covariantMetricTensor.Getg33(),
//              covariantMetricTensor.Getg12(), covariantMetricTensor.Getg13(),
//              covariantMetricTensor.Getg23(), this_J, this_Bxy);
//
//  output_progress.write("Calculating differential calculateGeometry terms\n");
//
//  if (min(abs(dx)) < 1e-8) {
//    throw BoutException("dx magnitude less than 1e-8");
//  }
//
//  if (min(abs(dy)) < 1e-8) {
//    throw BoutException("dy magnitude less than 1e-8");
//  }
//
//  if (min(abs(dz)) < 1e-8) {
//    throw BoutException("dz magnitude less than 1e-8");
//  }
//
//  // Check input metrics
//  checkContravariant();
//  checkCovariant();
//  CalculateChristoffelSymbols();
//
//  auto tmp = this_J * contravariantMetricTensor.Getg12();
//  communicate(tmp);
//  G1 = (DDX(this_J * contravariantMetricTensor.Getg11()) + DDY(tmp)
//        + DDZ(this_J * contravariantMetricTensor.Getg13()))
//       / this_J;
//  tmp = this_J * contravariantMetricTensor.Getg22();
//  communicate(tmp);
//  G2 = (DDX(this_J * contravariantMetricTensor.Getg12()) + DDY(tmp)
//        + DDZ(this_J * contravariantMetricTensor.Getg23()))
//       / this_J;
//  tmp = this_J * contravariantMetricTensor.Getg23();
//  communicate(tmp);
//  G3 = (DDX(this_J * contravariantMetricTensor.Getg13()) + DDY(tmp)
//        + DDZ(this_J * contravariantMetricTensor.Getg33()))
//       / this_J;
//
//  // Communicate christoffel symbol terms
//  output_progress.write("\tCommunicating connection terms\n");
//
//  communicate(G1_11, G1_22, G1_33, G1_12, G1_13, G1_23, G2_11, G2_22, G2_33, G2_12, G2_13,
//              G2_23, G3_11, G3_22, G3_33, G3_12, G3_13, G3_23, G1, G2, G3);
//
//  // Set boundary guard cells of Christoffel symbol terms
//  // Ideally, when location is staggered, we would set the upper/outer boundary point
//  // correctly rather than by extrapolating here: e.g. if location==CELL_YLOW and we are
//  // at the upper y-boundary the x- and z-derivatives at yend+1 at the boundary can be
//  // calculated because the guard cells are available, while the y-derivative could be
//  // calculated from the CELL_CENTRE metric components (which have guard cells available
//  // past the boundary location). This would avoid the problem that the y-boundary on the
//  // CELL_YLOW grid is at a 'guard cell' location (yend+1).
//  // However, the above would require lots of special handling, so just extrapolate for
//  // now.
//  G1_11 = interpolateAndExtrapolate(G1_11, location, true, true, true, transform.get());
//  G1_22 = interpolateAndExtrapolate(G1_22, location, true, true, true, transform.get());
//  G1_33 = interpolateAndExtrapolate(G1_33, location, true, true, true, transform.get());
//  G1_12 = interpolateAndExtrapolate(G1_12, location, true, true, true, transform.get());
//  G1_13 = interpolateAndExtrapolate(G1_13, location, true, true, true, transform.get());
//  G1_23 = interpolateAndExtrapolate(G1_23, location, true, true, true, transform.get());
//
//  G2_11 = interpolateAndExtrapolate(G2_11, location, true, true, true, transform.get());
//  G2_22 = interpolateAndExtrapolate(G2_22, location, true, true, true, transform.get());
//  G2_33 = interpolateAndExtrapolate(G2_33, location, true, true, true, transform.get());
//  G2_12 = interpolateAndExtrapolate(G2_12, location, true, true, true, transform.get());
//  G2_13 = interpolateAndExtrapolate(G2_13, location, true, true, true, transform.get());
//  G2_23 = interpolateAndExtrapolate(G2_23, location, true, true, true, transform.get());
//
//  G3_11 = interpolateAndExtrapolate(G3_11, location, true, true, true, transform.get());
//  G3_22 = interpolateAndExtrapolate(G3_22, location, true, true, true, transform.get());
//  G3_33 = interpolateAndExtrapolate(G3_33, location, true, true, true, transform.get());
//  G3_12 = interpolateAndExtrapolate(G3_12, location, true, true, true, transform.get());
//  G3_13 = interpolateAndExtrapolate(G3_13, location, true, true, true, transform.get());
//  G3_23 = interpolateAndExtrapolate(G3_23, location, true, true, true, transform.get());
//
//  G1 = interpolateAndExtrapolate(G1, location, true, true, true, transform.get());
//  G2 = interpolateAndExtrapolate(G2, location, true, true, true, transform.get());
//  G3 = interpolateAndExtrapolate(G3, location, true, true, true, transform.get());
//  //
//  //  //////////////////////////////////////////////////////
//  //  /// Non-uniform meshes. Need to use DDX, DDY
//  //
//  //  OPTION(Options::getRoot(), non_uniform, true);
//  //
//  //  Coordinates::FieldMetric d2x(localmesh), d2y(localmesh),
//  //      d2z(localmesh); // d^2 x / d i^2
//  //
//  //  // Read correction for non-uniform meshes
//  //  std::string suffix = getLocationSuffix(location);
//  //  if (location == CELL_CENTRE
//  //      or (!force_interpolate_from_centre and localmesh->sourceHasVar("dx" + suffix))) {
//  //    bool extrapolate_x = not localmesh->sourceHasXBoundaryGuards();
//  //    bool extrapolate_y = not localmesh->sourceHasYBoundaryGuards();
//  //
//  //    if (localmesh->get(d2x, "d2x" + suffix, 0.0, false, location)) {
//  //      output_warn.write("\tWARNING: differencing quantity 'd2x' not found. "
//  //                        "Calculating from dx\n");
//  //      d1_dx = bout::derivatives::index::DDX(1. / dx); // d/di(1/dx)
//  //
//  //      communicate(d1_dx);
//  //      d1_dx =
//  //          interpolateAndExtrapolate(d1_dx, location, true, true, true, transform.get());
//  //    } else {
//  //      d2x.setLocation(location);
//  //      // set boundary cells if necessary
//  //      d2x = interpolateAndExtrapolate(d2x, location, extrapolate_x, extrapolate_y, false,
//  //                                      transform.get());
//  //
//  //      d1_dx = -d2x / (dx * dx);
//  //    }
//  //
//  //    if (localmesh->get(d2y, "d2y" + suffix, 0.0, false, location)) {
//  //      output_warn.write("\tWARNING: differencing quantity 'd2y' not found. "
//  //                        "Calculating from dy\n");
//  //      d1_dy = DDY(1. / dy); // d/di(1/dy)
//  //
//  //      communicate(d1_dy);
//  //      d1_dy =
//  //          interpolateAndExtrapolate(d1_dy, location, true, true, true, transform.get());
//  //    } else {
//  //      d2y.setLocation(location);
//  //      // set boundary cells if necessary
//  //      d2y = interpolateAndExtrapolate(d2y, location, extrapolate_x, extrapolate_y, false,
//  //                                      transform.get());
//  //
//  //      d1_dy = -d2y / (dy * dy);
//  //    }
//  //
//  //#if BOUT_USE_METRIC_3D
//  //    if (localmesh->get(d2z, "d2z" + suffix, 0.0, false)) {
//  //      output_warn.write("\tWARNING: differencing quantity 'd2z' not found. "
//  //                        "Calculating from dz\n");
//  //      d1_dz = bout::derivatives::index::DDZ(1. / dz);
//  //      communicate(d1_dz);
//  //      d1_dz =
//  //          interpolateAndExtrapolate(d1_dz, location, true, true, true, transform.get());
//  //    } else {
//  //      d2z.setLocation(location);
//  //      // set boundary cells if necessary
//  //      d2z = interpolateAndExtrapolate(d2z, location, extrapolate_x, extrapolate_y, false,
//  //                                      transform.get());
//  //
//  //      d1_dz = -d2z / (dz * dz);
//  //    }
//  //#else
//  //    d1_dz = 0;
//  //#endif
//  //  } else {
//  //    if (localmesh->get(d2x, "d2x", 0.0, false)) {
//  //      output_warn.write("\tWARNING: differencing quantity 'd2x' not found. "
//  //                        "Calculating from dx\n");
//  //      d1_dx = bout::derivatives::index::DDX(1. / dx); // d/di(1/dx)
//  //
//  //      communicate(d1_dx);
//  //      d1_dx =
//  //          interpolateAndExtrapolate(d1_dx, location, true, true, true, transform.get());
//  //    } else {
//  //      // Shift d2x to our location
//  //      d2x = interpolateAndExtrapolate(d2x, location, true, true, false, transform.get());
//  //
//  //      d1_dx = -d2x / (dx * dx);
//  //    }
//  //
//  //    if (localmesh->get(d2y, "d2y", 0.0, false)) {
//  //      output_warn.write("\tWARNING: differencing quantity 'd2y' not found. "
//  //                        "Calculating from dy\n");
//  //      d1_dy = DDY(1. / dy); // d/di(1/dy)
//  //
//  //      communicate(d1_dy);
//  //      d1_dy =
//  //          interpolateAndExtrapolate(d1_dy, location, true, true, true, transform.get());
//  //    } else {
//  //      // Shift d2y to our location
//  //      d2y = interpolateAndExtrapolate(d2y, location, true, true, false, transform.get());
//  //
//  //      d1_dy = -d2y / (dy * dy);
//  //    }
//  //
//  //#if BOUT_USE_METRIC_3D
//  //    if (localmesh->get(d2z, "d2z", 0.0, false)) {
//  //      output_warn.write("\tWARNING: differencing quantity 'd2z' not found. "
//  //                        "Calculating from dz\n");
//  //      d1_dz = bout::derivatives::index::DDZ(1. / dz);
//  //
//  //      communicate(d1_dz);
//  //      d1_dz =
//  //          interpolateAndExtrapolate(d1_dz, location, true, true, true, transform.get());
//  //    } else {
//  //      // Shift d2z to our location
//  //      d2z = interpolateAndExtrapolate(d2z, location, true, true, false, transform.get());
//  //
//  //      d1_dz = -d2z / (dz * dz);
//  //    }
//  //#else
//  //    d1_dz = 0;
//  //#endif
//  //  }
//  //  communicate(d1_dx, d1_dy, d1_dz);
//  //
//  //  if (location == CELL_CENTRE && recalculate_staggered) {
//  //    // Re-calculate interpolated Coordinates at staggered locations
//  //    localmesh->recalculateStaggeredCoordinates();
//  //  }
//  //
//  //  // Invalidate and recalculate cached variables
//  //  zlength_cache.reset();
//  //  Grad2_par2_DDY_invSgCache.clear();
//  //  invSgCache.reset();
//  //
//  //  return 0;
//  //}
//
void Geometry::CalculateChristoffelSymbols() {
  // Calculate Christoffel symbol terms (18 independent values)
  // Note: This calculation is completely general: metric
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11 = 0.5 * contravariantMetricTensor.Getg11() * DDX(covariantMetricTensor.Getg11())
          + contravariantMetricTensor.Getg12()
                * (DDX(covariantMetricTensor.Getg12())
                   - 0.5 * DDY(covariantMetricTensor.Getg11()))
          + contravariantMetricTensor.Getg13()
                * (DDX(covariantMetricTensor.Getg13())
                   - 0.5 * DDZ(covariantMetricTensor.Getg11()));
  G1_22 = contravariantMetricTensor.Getg11()
              * (DDY(covariantMetricTensor.Getg12())
                 - 0.5 * DDX(covariantMetricTensor.Getg22()))
          + 0.5 * contravariantMetricTensor.Getg12() * DDY(covariantMetricTensor.Getg22())
          + contravariantMetricTensor.Getg13()
                * (DDY(covariantMetricTensor.Getg23())
                   - 0.5 * DDZ(covariantMetricTensor.Getg22()));
  G1_33 =
      contravariantMetricTensor.Getg11()
          * (DDZ(covariantMetricTensor.Getg13())
             - 0.5 * DDX(covariantMetricTensor.Getg33()))
      + contravariantMetricTensor.Getg12()
            * (DDZ(covariantMetricTensor.Getg23())
               - 0.5 * DDY(covariantMetricTensor.Getg33()))
      + 0.5 * contravariantMetricTensor.Getg13() * DDZ(covariantMetricTensor.Getg33());
  G1_12 =
      0.5 * contravariantMetricTensor.Getg11() * DDY(covariantMetricTensor.Getg11())
      + 0.5 * contravariantMetricTensor.Getg12() * DDX(covariantMetricTensor.Getg22())
      + 0.5 * contravariantMetricTensor.Getg13()
            * (DDY(covariantMetricTensor.Getg13()) + DDX(covariantMetricTensor.Getg23())
               - DDZ(covariantMetricTensor.Getg12()));
  G1_13 =
      0.5 * contravariantMetricTensor.Getg11() * DDZ(covariantMetricTensor.Getg11())
      + 0.5 * contravariantMetricTensor.Getg12()
            * (DDZ(covariantMetricTensor.Getg12()) + DDX(covariantMetricTensor.Getg23())
               - DDY(covariantMetricTensor.Getg13()))
      + 0.5 * contravariantMetricTensor.Getg13() * DDX(covariantMetricTensor.Getg33());
  G1_23 =
      0.5 * contravariantMetricTensor.Getg11()
          * (DDZ(covariantMetricTensor.Getg12()) + DDY(covariantMetricTensor.Getg13())
             - DDX(covariantMetricTensor.Getg23()))
      + 0.5 * contravariantMetricTensor.Getg12()
            * (DDZ(covariantMetricTensor.Getg22()) + DDY(covariantMetricTensor.Getg23())
               - DDY(covariantMetricTensor.Getg23()))
      // + 0.5 *g13*(DDZ(g_32) + DDY(g_33) - DDZ(g_23));
      // which equals
      + 0.5 * contravariantMetricTensor.Getg13() * DDY(covariantMetricTensor.Getg33());

  G2_11 = 0.5 * contravariantMetricTensor.Getg12() * DDX(covariantMetricTensor.Getg11())
          + contravariantMetricTensor.Getg22()
                * (DDX(covariantMetricTensor.Getg12())
                   - 0.5 * DDY(covariantMetricTensor.Getg11()))
          + contravariantMetricTensor.Getg23()
                * (DDX(covariantMetricTensor.Getg13())
                   - 0.5 * DDZ(covariantMetricTensor.Getg11()));
  G2_22 = contravariantMetricTensor.Getg12()
              * (DDY(covariantMetricTensor.Getg12())
                 - 0.5 * DDX(covariantMetricTensor.Getg22()))
          + 0.5 * contravariantMetricTensor.Getg22() * DDY(covariantMetricTensor.Getg22())
          + contravariantMetricTensor.Getg23()
                * (DDY(contravariantMetricTensor.Getg23())
                   - 0.5 * DDZ(covariantMetricTensor.Getg22()));
  G2_33 =
      contravariantMetricTensor.Getg12()
          * (DDZ(covariantMetricTensor.Getg13())
             - 0.5 * DDX(covariantMetricTensor.Getg33()))
      + contravariantMetricTensor.Getg22()
            * (DDZ(covariantMetricTensor.Getg23())
               - 0.5 * DDY(covariantMetricTensor.Getg33()))
      + 0.5 * contravariantMetricTensor.Getg23() * DDZ(covariantMetricTensor.Getg33());
  G2_12 =
      0.5 * contravariantMetricTensor.Getg12() * DDY(covariantMetricTensor.Getg11())
      + 0.5 * contravariantMetricTensor.Getg22() * DDX(covariantMetricTensor.Getg22())
      + 0.5 * contravariantMetricTensor.Getg23()
            * (DDY(covariantMetricTensor.Getg13()) + DDX(covariantMetricTensor.Getg23())
               - DDZ(covariantMetricTensor.Getg12()));
  G2_13 =
      // 0.5 *g21*(DDZ(covariantMetricTensor.Getg11()) + DDX(covariantMetricTensor.Getg13()) - DDX(covariantMetricTensor.Getg13()))
      // which equals
      0.5 * contravariantMetricTensor.Getg12()
          * (DDZ(covariantMetricTensor.Getg11()) + DDX(covariantMetricTensor.Getg13())
             - DDX(covariantMetricTensor.Getg13()))
      // + 0.5 *g22*(DDZ(covariantMetricTensor.Getg21()) + DDX(covariantMetricTensor.Getg23()) - DDY(covariantMetricTensor.Getg13()))
      // which equals
      + 0.5 * contravariantMetricTensor.Getg22()
            * (DDZ(covariantMetricTensor.Getg12()) + DDX(covariantMetricTensor.Getg23())
               - DDY(covariantMetricTensor.Getg13()))
      // + 0.5 *g23*(DDZ(covariantMetricTensor.Getg31()) + DDX(covariantMetricTensor.Getg33()) - DDZ(g_13));
      // which equals
      + 0.5 * contravariantMetricTensor.Getg23() * DDX(covariantMetricTensor.Getg33());
  G2_23 =
      0.5 * contravariantMetricTensor.Getg12()
          * (DDZ(covariantMetricTensor.Getg12()) + DDY(covariantMetricTensor.Getg13())
             - DDX(covariantMetricTensor.Getg23()))
      + 0.5 * contravariantMetricTensor.Getg22() * DDZ(covariantMetricTensor.Getg22())
      + 0.5 * contravariantMetricTensor.Getg23() * DDY(covariantMetricTensor.Getg33());

  G3_11 = 0.5 * contravariantMetricTensor.Getg13() * DDX(covariantMetricTensor.Getg11())
          + contravariantMetricTensor.Getg23()
                * (DDX(covariantMetricTensor.Getg12())
                   - 0.5 * DDY(covariantMetricTensor.Getg11()))
          + contravariantMetricTensor.Getg33()
                * (DDX(covariantMetricTensor.Getg13())
                   - 0.5 * DDZ(covariantMetricTensor.Getg11()));
  G3_22 = contravariantMetricTensor.Getg13()
              * (DDY(covariantMetricTensor.Getg12())
                 - 0.5 * DDX(covariantMetricTensor.Getg22()))
          + 0.5 * contravariantMetricTensor.Getg23() * DDY(covariantMetricTensor.Getg22())
          + contravariantMetricTensor.Getg33()
                * (DDY(covariantMetricTensor.Getg23())
                   - 0.5 * DDZ(covariantMetricTensor.Getg22()));
  G3_33 =
      contravariantMetricTensor.Getg13()
          * (DDZ(covariantMetricTensor.Getg13())
             - 0.5 * DDX(covariantMetricTensor.Getg33()))
      + contravariantMetricTensor.Getg23()
            * (DDZ(covariantMetricTensor.Getg23())
               - 0.5 * DDY(covariantMetricTensor.Getg33()))
      + 0.5 * contravariantMetricTensor.Getg33() * DDZ(covariantMetricTensor.Getg33());
  G3_12 =
      // 0.5 *g31*(DDY(covariantMetricTensor.Getg11()) + DDX(covariantMetricTensor.Getg12()) - DDX(covariantMetricTensor.Getg12()))
      // which equals to
      0.5 * contravariantMetricTensor.Getg13() * DDY(covariantMetricTensor.Getg11())
      // + 0.5 *g32*(DDY(covariantMetricTensor.Getg21()) + DDX(covariantMetricTensor.Getg22()) - DDY(covariantMetricTensor.Getg12()))
      // which equals to
      + 0.5 * contravariantMetricTensor.Getg23() * DDX(covariantMetricTensor.Getg22())
      //+ 0.5 *g33*(DDY(covariantMetricTensor.Getg31()) + DDX(covariantMetricTensor.Getg32()) - DDZ(covariantMetricTensor.Getg12()));
      // which equals to
      + 0.5 * contravariantMetricTensor.Getg33()
            * (DDY(covariantMetricTensor.Getg13()) + DDX(covariantMetricTensor.Getg23())
               - DDZ(covariantMetricTensor.Getg12()));
  G3_13 =
      0.5 * contravariantMetricTensor.Getg13() * DDZ(covariantMetricTensor.Getg11())
      + 0.5 * contravariantMetricTensor.Getg23()
            * (DDZ(covariantMetricTensor.Getg12()) + DDX(covariantMetricTensor.Getg23())
               - DDY(covariantMetricTensor.Getg13()))
      + 0.5 * contravariantMetricTensor.Getg33() * DDX(covariantMetricTensor.Getg33());
  G3_23 =
      0.5 * contravariantMetricTensor.Getg13()
          * (DDZ(covariantMetricTensor.Getg12()) + DDY(covariantMetricTensor.Getg13())
             - DDX(covariantMetricTensor.Getg23()))
      + 0.5 * contravariantMetricTensor.Getg23() * DDZ(covariantMetricTensor.Getg22())
      + 0.5 * contravariantMetricTensor.Getg33() * DDY(covariantMetricTensor.Getg33());
}

void Geometry::calcCovariant(CELL_LOC cell_location, const std::string& region) {
  TRACE("Geometry::calcCovariant");
  covariantMetricTensor.setMetricTensor(
      contravariantMetricTensor.oppositeRepresentation(cell_location, region));
}

void Geometry::calcContravariant(CELL_LOC cell_location, const std::string& region) {
  TRACE("Geometry::calcContravariant");
  contravariantMetricTensor.setMetricTensor(
      covariantMetricTensor.oppositeRepresentation(cell_location, region));
}

//void Geometry::jacobian(bool extrapolate_x, bool extrapolate_y) {
//  TRACE("Geometry::jacobian");
//  try {
//
//    const auto j = recalculateJacobian(extrapolate_x, extrapolate_y);
//    // More robust to extrapolate derived quantities directly, rather than
//    // deriving from extrapolated covariant metric components
//    this_J = interpolateAndExtrapolate(j, extrapolate_x, extrapolate_y, false);
//
//    const auto Bxy = recalculateBxy(extrapolate_x, extrapolate_y);
////    CELL_LOC location, ParallelTransform* pParallelTransform
//    this_Bxy = interpolateAndExtrapolate(Bxy, location, extrapolate_x, extrapolate_y, false,
//                                     transform.get());
//  } catch (BoutException&) {
//    output_error.write("\tError in jacobian call\n");
//    throw;
//  }
//}

MetricTensor::FieldMetric Geometry::recalculateJacobian() {

  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  auto g = contravariantMetricTensor.Getg11() * contravariantMetricTensor.Getg22()
               * contravariantMetricTensor.Getg33()
           + 2.0 * contravariantMetricTensor.Getg12() * contravariantMetricTensor.Getg13()
                 * contravariantMetricTensor.Getg23()
           - contravariantMetricTensor.Getg11() * contravariantMetricTensor.Getg23()
                 * contravariantMetricTensor.Getg23()
           - contravariantMetricTensor.Getg22() * contravariantMetricTensor.Getg13()
                 * contravariantMetricTensor.Getg13()
           - contravariantMetricTensor.Getg33() * contravariantMetricTensor.Getg12()
                 * contravariantMetricTensor.Getg12();

  // Check that g is positive
  bout::checkPositive(g, "The determinant of g^ij", "RGN_NOBNDRY");

  return 1. / sqrt(g);
}

MetricTensor::FieldMetric Geometry::recalculateBxy() {

  return sqrt(covariantMetricTensor.Getg22()) / this_J;
}

//  namespace {
//  // Utility function for fixing up guard cells of zShift
//  void fixZShiftGuards(Field2D& zShift) {
//    auto localmesh = zShift.getMesh();
//
//    // extrapolate into boundary guard cells if necessary
//    zShift = interpolateAndExtrapolate(zShift, zShift.getLocation(),
//                                       not localmesh->sourceHasXBoundaryGuards(),
//                                       not localmesh->sourceHasYBoundaryGuards(), false);
//
//    // make sure zShift has been communicated
//    communicate(zShift);
//
//    // Correct guard cells for discontinuity of zShift at poloidal branch cut
//    for (int x = 0; x < localmesh->LocalNx; x++) {
//      const auto lower = localmesh->hasBranchCutLower(x);
//      if (lower.first) {
//        for (int y = 0; y < localmesh->ystart; y++) {
//          zShift(x, y) -= lower.second;
//        }
//      }
//      const auto upper = localmesh->hasBranchCutUpper(x);
//      if (upper.first) {
//        for (int y = localmesh->yend + 1; y < localmesh->LocalNy; y++) {
//          zShift(x, y) += upper.second;
//        }
//      }
//    }
//  }
//  } // namespace
//
//  //void Coordinates::setParallelTransform(Options* options) {
//  //
//  //  auto ptoptions = options->getSection("paralleltransform");
//  //
//  //  std::string ptstr;
//  //  ptoptions->get("type", ptstr, "identity");
//  //
//  //  // Convert to lower case for comparison
//  //  ptstr = lowercase(ptstr);
//  //
//  //  if (ptstr == "identity") {
//  //    // Identity method i.e. no transform needed
//  //    transform =
//  //        bout::utils::make_unique<ParallelTransformIdentity>(*localmesh, ptoptions);
//  //
//  //  } else if (ptstr == "shifted" or ptstr == "shiftedinterp") {
//  //    // Shifted metric method
//  //
//  //    Field2D zShift{localmesh};
//  //
//  //    // Read the zShift angle from the mesh
//  //    std::string suffix = getLocationSuffix(location);
//  //    if (localmesh->sourceHasVar("dx" + suffix)) {
//  //      // Grid file has variables at this location, so should be able to read
//  //      checkStaggeredGet(localmesh, "zShift", suffix);
//  //      if (localmesh->get(zShift, "zShift" + suffix, 0.0, false, location)) {
//  //        // No zShift variable. Try qinty in BOUT grid files
//  //        if (localmesh->get(zShift, "qinty" + suffix, 0.0, false, location)) {
//  //          // Failed to find either variable, cannot use ShiftedMetric
//  //          throw BoutException("Could not read zShift" + suffix + " from grid file");
//  //        }
//  //      }
//  //    } else {
//  //      if (location == CELL_YLOW and bout::build::use_metric_3d) {
//  //        throw BoutException("Cannot interpolate zShift to construct ShiftedMetric when "
//  //                            "using 3d metrics. You must provide zShift_ylow in the grid "
//  //                            "file.");
//  //      }
//  //      Field2D zShift_centre;
//  //      if (localmesh->get(zShift_centre, "zShift", 0.0, false)) {
//  //        // No zShift variable. Try qinty in BOUT grid files
//  //        if (localmesh->get(zShift_centre, "qinty", 0.0, false)) {
//  //          // Failed to find either variable, cannot use ShiftedMetric
//  //          throw BoutException("Could not read zShift from grid file");
//  //        }
//  //      }
//  //
//  //      fixZShiftGuards(zShift_centre);
//  //
//  //      zShift = interpolateAndExtrapolate(zShift_centre, location, true, true, false,
//  //                                         transform.get());
//  //    }
//  //
//  //    fixZShiftGuards(zShift);
//  //
//  //    if (ptstr == "shifted") {
//  //      transform = bout::utils::make_unique<ShiftedMetric>(*localmesh, location, zShift,
//  //                                                          getUniform(zlength()));
//  //    } else if (ptstr == "shiftedinterp") {
//  //      transform = bout::utils::make_unique<ShiftedMetricInterp>(
//  //          *localmesh, location, zShift, getUniform(zlength()));
//  //    }
//  //
//  //  } else if (ptstr == "fci") {
//  //
//  //    if (location != CELL_CENTRE) {
//  //      throw BoutException("FCITransform is not available on staggered grids.");
//  //    }
//  //
//  //    // Flux Coordinate Independent method
//  //    const bool fci_zperiodic = (*ptoptions)["z_periodic"].withDefault(true);
//  //    transform =
//  //        bout::utils::make_unique<FCITransform>(*localmesh, dy, fci_zperiodic, ptoptions);
//  //
//  //  } else {
//  //    throw BoutException(_("Unrecognised paralleltransform option.\n"
//  //                          "Valid choices are 'identity', 'shifted', 'fci'"));
//  //  }
//  //}
//
//  ///*******************************************************************************
//  // * Operators
//  // *
//  // *******************************************************************************/
//  //
//  //Coordinates::FieldMetric Coordinates::DDX(const Field2D& f, CELL_LOC loc,
//  //                                          const std::string& method,
//  //                                          const std::string& region) {
//  //  ASSERT1(location == loc || loc == CELL_DEFAULT)
//  //  return bout::derivatives::index::DDX(f, loc, method, region) / dx;
//  //}
//  //Field3D Coordinates::DDX(const Field3D& f, CELL_LOC outloc, const std::string& method,
//  //                         const std::string& region) {
//  //
//  //  auto result = bout::derivatives::index::DDX(f, outloc, method, region);
//  //  result /= dx;
//  //
//  //  if (f.getMesh()->IncIntShear) {
//  //    // Using BOUT-06 style shifting
//  //    result += IntShiftTorsion * DDZ(f, outloc, method, region);
//  //  }
//  //
//  //  return result;
//  //}
//  //
//  //Coordinates::FieldMetric Coordinates::DDY(const Field2D& f, CELL_LOC loc,
//  //                                          const std::string& method,
//  //                                          const std::string& region) const {
//  //  ASSERT1(location == loc || loc == CELL_DEFAULT)
//  //  return bout::derivatives::index::DDY(f, loc, method, region) / dy;
//  //}
//  //
//  //Field3D Coordinates::DDY(const Field3D& f, CELL_LOC outloc, const std::string& method,
//  //                         const std::string& region) const {
//  //#if BOUT_USE_METRIC_3D
//  //  if (!f.hasParallelSlices() and !transform->canToFromFieldAligned()) {
//  //    Field3D f_parallel = f;
//  //    transform->calcParallelSlices(f_parallel);
//  //    f_parallel.applyParallelBoundary("parallel_neumann");
//  //    return bout::derivatives::index::DDY(f_parallel, outloc, method, region);
//  //  }
//  //#endif
//  //  return bout::derivatives::index::DDY(f, outloc, method, region) / dy;
//  //}
//  //
//  //Coordinates::FieldMetric Coordinates::DDZ(const Field2D& f, CELL_LOC loc,
//  //                                          const std::string& UNUSED(method),
//  //                                          const std::string& UNUSED(region)) {
//  //  ASSERT1(location == loc || loc == CELL_DEFAULT)
//  //  ASSERT1(f.getMesh() == localmesh)
//  //  if (loc == CELL_DEFAULT) {
//  //    loc = f.getLocation();
//  //  }
//  //  return zeroFrom(f).setLocation(loc);
//  //}
//  //Field3D Coordinates::DDZ(const Field3D& f, CELL_LOC outloc, const std::string& method,
//  //                         const std::string& region) {
//  //  return bout::derivatives::index::DDZ(f, outloc, method, region) / dz;
//  //}
//
//  ///////////////////////////////////////////////////////////
//  //// Parallel gradient
//  //
//  //Coordinates::FieldMetric Coordinates::Grad_par(const Field2D& var,
//  //                                               MAYBE_UNUSED(CELL_LOC outloc),
//  //                                               const std::string& UNUSED(method)) {
//  //  TRACE("Coordinates::Grad_par( Field2D )");
//  //  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == var.getLocation()))
//  //
//  //  return DDY(var) * invSg();
//  //}
//  //
//  //Field3D Coordinates::Grad_par(const Field3D& var, CELL_LOC outloc,
//  //                              const std::string& method) {
//  //  TRACE("Coordinates::Grad_par( Field3D )");
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  return ::DDY(var, outloc, method) * invSg();
//  //}
//  //
//  ///////////////////////////////////////////////////////////
//  //// Vpar_Grad_par
//  //// vparallel times the parallel derivative along unperturbed B-field
//  //
//  //Coordinates::FieldMetric Coordinates::Vpar_Grad_par(const Field2D& v, const Field2D& f,
//  //                                                    MAYBE_UNUSED(CELL_LOC outloc),
//  //                                                    const std::string& UNUSED(method)) {
//  //  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()))
//  //
//  //  return VDDY(v, f) * invSg();
//  //}
//  //
//  //Field3D Coordinates::Vpar_Grad_par(const Field3D& v, const Field3D& f, CELL_LOC outloc,
//  //                                   const std::string& method) {
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  return VDDY(v, f, outloc, method) * invSg();
//  //}
//  //
//  ///////////////////////////////////////////////////////////
//  //// Parallel divergence
//  //
//  //Coordinates::FieldMetric Coordinates::Div_par(const Field2D& f, CELL_LOC outloc,
//  //                                              const std::string& method) {
//  //  TRACE("Coordinates::Div_par( Field2D )");
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  // Need Bxy at location of f, which might be different from location of this
//  //  // Coordinates object
//  //  auto Bxy_floc = f.getCoordinates()->Bxy();
//  //
//  //  return this_Bxy * Grad_par(f / Bxy_floc, outloc, method);
//  //}
//  //
//  //Field3D Coordinates::Div_par(const Field3D& f, CELL_LOC outloc,
//  //                             const std::string& method) {
//  //  TRACE("Coordinates::Div_par( Field3D )");
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  // Need Bxy at location of f, which might be different from location of this
//  //  // Coordinates object
//  //  auto Bxy_floc = f.getCoordinates()->Bxy();
//  //
//  //  if (!f.hasParallelSlices()) {
//  //    // No yup/ydown fields. The Grad_par operator will
//  //    // shift to field aligned coordinates
//  //    return this_Bxy * Grad_par(f / Bxy_floc, outloc, method);
//  //  }
//  //
//  //  // Need to modify yup and ydown fields
//  //  Field3D f_B = f / Bxy_floc;
//  //  f_B.splitParallelSlices();
//  //  for (int i = 0; i < f.getMesh()->ystart; ++i) {
//  //    f_B.yup(i) = f.yup(i) / Bxy_floc.yup(i);
//  //    f_B.ydown(i) = f.ydown(i) / Bxy_floc.ydown(i);
//  //  }
//  //  return this_Bxy * Grad_par(f_B, outloc, method);
//  //}
//  //
//  ///////////////////////////////////////////////////////////
//  //// second parallel derivative (b dot Grad)(b dot Grad)
//  //// Note: For parallel Laplacian use Laplace_par
//  //
//  //Coordinates::FieldMetric Coordinates::Grad2_par2(const Field2D& f, CELL_LOC outloc,
//  //                                                 const std::string& method) {
//  //  TRACE("Coordinates::Grad2_par2( Field2D )");
//  //  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()))
//  //
//  //  auto result = Grad2_par2_DDY_invSg(outloc, method) * DDY(f, outloc, method)
//  //                + D2DY2(f, outloc, method) / covariantMetricTensor.Getg22();
//  //
//  //  return result;
//  //}
//  //
//  //Field3D Coordinates::Grad2_par2(const Field3D& f, CELL_LOC outloc,
//  //                                const std::string& method) {
//  //  TRACE("Coordinates::Grad2_par2( Field3D )");
//  //  if (outloc == CELL_DEFAULT) {
//  //    outloc = f.getLocation();
//  //  }
//  //  ASSERT1(location == outloc)
//  //
//  //  Field3D result = ::DDY(f, outloc, method);
//  //
//  //  Field3D r2 = D2DY2(f, outloc, method) / covariantMetricTensor.Getg22();
//  //
//  //  result = Grad2_par2_DDY_invSg(outloc, method) * result + r2;
//  //
//  //  ASSERT2(result.getLocation() == outloc)
//  //
//  //  return result;
//  //}
//  //
//  ///////////////////////////////////////////////////////////
//  //// perpendicular Laplacian operator
//  //
//  //#include <bout/invert_laplace.hxx> // Delp2 uses same coefficients as inversion code
//  //
//  //Coordinates::FieldMetric Coordinates::Delp2(const Field2D& f, CELL_LOC outloc,
//  //                                            bool UNUSED(useFFT)) {
//  //  TRACE("Coordinates::Delp2( Field2D )");
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  auto result =
//  //      G1 * DDX(f, outloc) + contravariantMetricTensor.Getg11() * D2DX2(f, outloc);
//  //
//  //  return result;
//  //}
//  //
//  //Field3D Coordinates::Delp2(const Field3D& f, CELL_LOC outloc, bool useFFT) {
//  //  TRACE("Coordinates::Delp2( Field3D )");
//  //
//  //  if (outloc == CELL_DEFAULT) {
//  //    outloc = f.getLocation();
//  //  }
//  //
//  //  ASSERT1(location == outloc)
//  //  ASSERT1(f.getLocation() == outloc)
//  //
//  //  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
//  //    // copy mesh, location, etc
//  //    return f * 0;
//  //  }
//  //  ASSERT2(localmesh->xstart > 0) // Need at least one guard cell
//  //
//  //  Field3D result{emptyFrom(f).setLocation(outloc)};
//  //
//  //  if (useFFT and not bout::build::use_metric_3d) {
//  //    int ncz = localmesh->LocalNz;
//  //
//  //    // Allocate memory
//  //    auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
//  //    auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
//  //
//  //    // Loop over y indices
//  //    // Note: should not include y-guard or y-boundary points here as that would
//  //    // use values from corner cells in dx, which may not be initialised.
//  //    for (int jy = localmesh->ystart; jy <= localmesh->yend; jy++) {
//  //
//  //      // Take forward FFT
//  //
//  //      for (int jx = 0; jx < localmesh->LocalNx; jx++) {
//  //        rfft(&f(jx, jy, 0), ncz, &ft(jx, 0));
//  //      }
//  //
//  //      // Loop over kz
//  //      for (int jz = 0; jz <= ncz / 2; jz++) {
//  //
//  //        // No smoothing in the x direction
//  //        for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
//  //          // Perform x derivative
//  //
//  //          dcomplex a, b, c;
//  //          laplace_tridag_coefs(jx, jy, jz, a, b, c, nullptr, nullptr, outloc);
//  //
//  //          delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
//  //        }
//  //      }
//  //
//  //      // Reverse FFT
//  //      for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
//  //
//  //        irfft(&delft(jx, 0), ncz, &result(jx, jy, 0));
//  //      }
//  //    }
//  //  } else {
//  //    result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc)
//  //             + contravariantMetricTensor.Getg11() * ::D2DX2(f, outloc)
//  //             + contravariantMetricTensor.Getg33() * ::D2DZ2(f, outloc)
//  //             + 2 * contravariantMetricTensor.Getg13() * ::D2DXDZ(f, outloc);
//  //  }
//  //
//  //  ASSERT2(result.getLocation() == outloc)
//  //
//  //  return result;
//  //}
//  //
//  //FieldPerp Coordinates::Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
//  //  TRACE("Coordinates::Delp2( FieldPerp )");
//  //
//  //  if (outloc == CELL_DEFAULT) {
//  //    outloc = f.getLocation();
//  //  }
//  //
//  //  ASSERT1(location == outloc)
//  //  ASSERT1(f.getLocation() == outloc)
//  //
//  //  if (localmesh->GlobalNx == 1 && localmesh->GlobalNz == 1) {
//  //    // copy mesh, location, etc
//  //    return f * 0;
//  //  }
//  //  ASSERT2(localmesh->xstart > 0) // Need at least one guard cell
//  //
//  //  FieldPerp result{emptyFrom(f).setLocation(outloc)};
//  //
//  //  int jy = f.getIndex();
//  //  result.setIndex(jy);
//  //
//  //  if (useFFT) {
//  //    int ncz = localmesh->LocalNz;
//  //
//  //    // Allocate memory
//  //    auto ft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
//  //    auto delft = Matrix<dcomplex>(localmesh->LocalNx, ncz / 2 + 1);
//  //
//  //    // Take forward FFT
//  //    for (int jx = 0; jx < localmesh->LocalNx; jx++) {
//  //      rfft(&f(jx, 0), ncz, &ft(jx, 0));
//  //    }
//  //
//  //    // Loop over kz
//  //    for (int jz = 0; jz <= ncz / 2; jz++) {
//  //
//  //      // No smoothing in the x direction
//  //      for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
//  //        // Perform x derivative
//  //
//  //        dcomplex a, b, c;
//  //        laplace_tridag_coefs(jx, jy, jz, a, b, c);
//  //
//  //        delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
//  //      }
//  //    }
//  //
//  //    // Reverse FFT
//  //    for (int jx = localmesh->xstart; jx <= localmesh->xend; jx++) {
//  //      irfft(&delft(jx, 0), ncz, &result(jx, 0));
//  //    }
//  //
//  //  } else {
//  //    throw BoutException("Non-fourier Delp2 not currently implented for FieldPerp.");
//  //    // Would be the following but don't have standard derivative operators for FieldPerps
//  //    // yet
//  //    // result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * ::D2DX2(f, outloc)
//  //    //          + g33 * ::D2DZ2(f, outloc) + 2 * g13 * ::D2DXDZ(f, outloc);
//  //  }
//  //
//  //  return result;
//  //}
//  //
//  //Coordinates::FieldMetric Coordinates::Laplace_par(const Field2D& f, CELL_LOC outloc) {
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  return D2DY2(f, outloc) / covariantMetricTensor.Getg22()
//  //         + DDY(this_J / covariantMetricTensor.Getg22(), outloc) * DDY(f, outloc) / this_J;
//  //}
//  //
//  //Field3D Coordinates::Laplace_par(const Field3D& f, CELL_LOC outloc) {
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //  return D2DY2(f, outloc) / covariantMetricTensor.Getg22()
//  //         + DDY(this_J / covariantMetricTensor.Getg22(), outloc) * ::DDY(f, outloc)
//  //               / this_J;
//  //}
//  //
//  //// Full Laplacian operator on scalar field
//  //
//  //Coordinates::FieldMetric Coordinates::Laplace(const Field2D& f, CELL_LOC outloc,
//  //                                              const std::string& dfdy_boundary_conditions,
//  //                                              const std::string& dfdy_dy_region) {
//  //  TRACE("Coordinates::Laplace( Field2D )");
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  auto result = G1 * DDX(f, outloc) + G2 * DDY(f, outloc)
//  //                + contravariantMetricTensor.Getg11() * D2DX2(f, outloc)
//  //                + contravariantMetricTensor.Getg22() * D2DY2(f, outloc)
//  //                + 2.0 * contravariantMetricTensor.Getg12()
//  //                      * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY",
//  //                               dfdy_boundary_conditions, dfdy_dy_region);
//  //
//  //  return result;
//  //}
//  //
//  //Field3D Coordinates::Laplace(const Field3D& f, CELL_LOC outloc,
//  //                             const std::string& dfdy_boundary_conditions,
//  //                             const std::string& dfdy_dy_region) {
//  //  TRACE("Coordinates::Laplace( Field3D )");
//  //  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//  //
//  //  Field3D result = G1 * ::DDX(f, outloc) + G2 * ::DDY(f, outloc) + G3 * ::DDZ(f, outloc)
//  //                   + contravariantMetricTensor.Getg11() * D2DX2(f, outloc)
//  //                   + contravariantMetricTensor.Getg22() * D2DY2(f, outloc)
//  //                   + contravariantMetricTensor.Getg33() * D2DZ2(f, outloc)
//  //                   + 2.0
//  //                         * (contravariantMetricTensor.Getg12()
//  //                                * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY",
//  //                                         dfdy_boundary_conditions, dfdy_dy_region)
//  //                            + contravariantMetricTensor.Getg13() * D2DXDZ(f, outloc)
//  //                            + contravariantMetricTensor.Getg23() * D2DYDZ(f, outloc));
//  //
//  //  return result;
//  //}
//  //
//  //// Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
//  //// solver
//  //Field2D Coordinates::Laplace_perpXY(MAYBE_UNUSED(const Field2D& A),
//  //                                    MAYBE_UNUSED(const Field2D& f)) {
//  //  TRACE("Coordinates::Laplace_perpXY( Field2D )");
//  //#if not(BOUT_USE_METRIC_3D)
//  //  Field2D result;
//  //  result.allocate();
//  //  for (auto i : result.getRegion(RGN_NOBNDRY)) {
//  //    result[i] = 0.;
//  //
//  //    // outer x boundary
//  //    const auto outer_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xp()]); };
//  //    const BoutReal outer_x_A = outer_x_avg(A);
//  //    const BoutReal outer_x_J = outer_x_avg(this_J);
//  //    const BoutReal outer_x_g11 = outer_x_avg(contravariantMetricTensor.Getg11());
//  //    const BoutReal outer_x_dx = outer_x_avg(dx);
//  //    const BoutReal outer_x_value =
//  //        outer_x_A * outer_x_J * outer_x_g11 / (this_J[i] * outer_x_dx * dx[i]);
//  //    result[i] += outer_x_value * (f[i.xp()] - f[i]);
//  //
//  //    // inner x boundary
//  //    const auto inner_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xm()]); };
//  //    const BoutReal inner_x_A = inner_x_avg(A);
//  //    const BoutReal inner_x_J = inner_x_avg(this_J);
//  //    const BoutReal inner_x_g11 = inner_x_avg(contravariantMetricTensor.Getg11());
//  //    const BoutReal inner_x_dx = inner_x_avg(dx);
//  //    const BoutReal inner_x_value =
//  //        inner_x_A * inner_x_J * inner_x_g11 / (this_J[i] * inner_x_dx * dx[i]);
//  //    result[i] += inner_x_value * (f[i.xm()] - f[i]);
//  //
//  //    // upper y boundary
//  //    const auto upper_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.yp()]); };
//  //    const BoutReal upper_y_A = upper_y_avg(A);
//  //    const BoutReal upper_y_J = upper_y_avg(this_J);
//  //    const BoutReal upper_y_g_22 = upper_y_avg(covariantMetricTensor.Getg22());
//  //    const BoutReal upper_y_g23 = upper_y_avg(covariantMetricTensor.Getg23());
//  //    const BoutReal upper_y_g_23 = upper_y_avg(covariantMetricTensor.Getg23());
//  //    const BoutReal upper_y_dy = upper_y_avg(dy);
//  //    const BoutReal upper_y_value = -upper_y_A * upper_y_J * upper_y_g23 * upper_y_g_23
//  //                                   / (upper_y_g_22 * this_J[i] * upper_y_dy * dy[i]);
//  //    result[i] += upper_y_value * (f[i.yp()] - f[i]);
//  //
//  //    // lower y boundary
//  //    const auto lower_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.ym()]); };
//  //    const BoutReal lower_y_A = lower_y_avg(A);
//  //    const BoutReal lower_y_J = lower_y_avg(this_J);
//  //    const BoutReal lower_y_g_22 = lower_y_avg(covariantMetricTensor.Getg22());
//  //    const BoutReal lower_y_g23 = lower_y_avg(contravariantMetricTensor.Getg23());
//  //    const BoutReal lower_y_g_23 = lower_y_avg(covariantMetricTensor.Getg23());
//  //    const BoutReal lower_y_dy = lower_y_avg(dy);
//  //    const BoutReal lower_y_value = -lower_y_A * lower_y_J * lower_y_g23 * lower_y_g_23
//  //                                   / (lower_y_g_22 * this_J[i] * lower_y_dy * dy[i]);
//  //    result[i] += lower_y_value * (f[i.ym()] - f[i]);
//  //  }
//  //
//  //  return result;
//  //#else
//  //  throw BoutException("Coordinates::Laplace_perpXY for 3D metric not implemented");
//  //#endif
//  //}
//  //
//  //const Coordinates::FieldMetric& Coordinates::invSg() const {
//  //  if (invSgCache == nullptr) {
//  //    auto ptr = std::make_unique<FieldMetric>();
//  //    (*ptr) = 1.0 / sqrt(covariantMetricTensor.Getg22());
//  //    invSgCache = std::move(ptr);
//  //  }
//  //  return *invSgCache;
//  //}
//  //
//  //const Coordinates::FieldMetric&
//  //Coordinates::Grad2_par2_DDY_invSg(CELL_LOC outloc, const std::string& method) const {
//  //  if (auto search = Grad2_par2_DDY_invSgCache.find(method);
//  //      search != Grad2_par2_DDY_invSgCache.end()) {
//  //    return *search->second;
//  //  }
//  //  invSg();
//  //
//  //  // Communicate to get parallel slices
//  //  localmesh->communicate(*invSgCache);
//  //  invSgCache->applyParallelBoundary("parallel_neumann");
//  //
//  //  // cache
//  //  auto ptr = std::make_unique<FieldMetric>();
//  //  *ptr = DDY(*invSgCache, outloc, method) * invSg();
//  //  Grad2_par2_DDY_invSgCache[method] = std::move(ptr);
//  //  return *Grad2_par2_DDY_invSgCache[method];
//  //}

void Geometry::checkCovariant(int ystart) { covariantMetricTensor.check(ystart); }

void Geometry::checkContravariant(int ystart) { contravariantMetricTensor.check(ystart); }

void Geometry::setContravariantMetricTensor(MetricTensor metric_tensor,
                                            CELL_LOC cell_location,
                                            const std::string& region) {
  contravariantMetricTensor.setMetricTensor(metric_tensor);
  covariantMetricTensor.setMetricTensor(
      contravariantMetricTensor.oppositeRepresentation(cell_location, region));
}

void Geometry::setCovariantMetricTensor(MetricTensor metric_tensor,
                                        CELL_LOC cell_location,
                                        const std::string& region) {
  covariantMetricTensor.setMetricTensor(metric_tensor);
  contravariantMetricTensor.setMetricTensor(
      covariantMetricTensor.oppositeRepresentation(cell_location, region));
}

const MetricTensor::FieldMetric& Geometry::g_11() const {
  return covariantMetricTensor.Getg11();
}
const MetricTensor::FieldMetric& Geometry::g_22() const {
  return covariantMetricTensor.Getg22();
}
const MetricTensor::FieldMetric& Geometry::g_33() const {
  return covariantMetricTensor.Getg33();
}
const MetricTensor::FieldMetric& Geometry::g_12() const {
  return covariantMetricTensor.Getg12();
}
const MetricTensor::FieldMetric& Geometry::g_13() const {
  return covariantMetricTensor.Getg13();
}
const MetricTensor::FieldMetric& Geometry::g_23() const {
  return covariantMetricTensor.Getg23();
}

const MetricTensor::FieldMetric& Geometry::g11() const {
  return contravariantMetricTensor.Getg11();
}
const MetricTensor::FieldMetric& Geometry::g22() const {
  return contravariantMetricTensor.Getg22();
}
const MetricTensor::FieldMetric& Geometry::g33() const {
  return contravariantMetricTensor.Getg33();
}
const MetricTensor::FieldMetric& Geometry::g12() const {
  return contravariantMetricTensor.Getg12();
}
const MetricTensor::FieldMetric& Geometry::g13() const {
  return contravariantMetricTensor.Getg13();
}
const MetricTensor::FieldMetric& Geometry::g23() const {
  return contravariantMetricTensor.Getg23();
}

const Geometry::FieldMetric& Geometry::J() const { return this_J; }

const Geometry::FieldMetric& Geometry::Bxy() const { return this_Bxy; }

void Geometry::setJ(FieldMetric J) {
  //TODO: Calculate J and check value is close
  this_J = J;
}

void Geometry::setJ(BoutReal value, int x, int y) {
  //TODO: Calculate Bxy and check value is close
  this_J(x, y) = value;
}

void Geometry::setBxy(FieldMetric Bxy) {
  //TODO: Calculate Bxy and check value is close
  this_Bxy = Bxy;
}

const MetricTensor& Geometry::getContravariantMetricTensor() const {
  return contravariantMetricTensor;
}

void Geometry::applyToContravariantMetricTensor(
    std::function<const FieldMetric(const FieldMetric)> function) {
  contravariantMetricTensor.map(function);
}

void Geometry::applyToCovariantMetricTensor(
    std::function<const FieldMetric(const FieldMetric)> function) {
  covariantMetricTensor.map(function);
}

//MetricTensor& Geometry::getCovariantMetricTensor() { return covariantMetricTensor; }
