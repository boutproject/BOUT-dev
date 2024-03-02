
#include "bout/differential_operators.hxx"
//#include "bout/field2d.hxx"
//#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
//#include "bout/metricTensor.hxx"
//#include "bout/index_derivs_interface.hxx"
//#include "bout/paralleltransform.hxx"
#include <bout/derivs.hxx>

DifferentialOperators::DifferentialOperators(const DifferentialOperators& operators) {}

DifferentialOperators::DifferentialOperators(Mesh* mesh, FieldMetric intShiftTorsion,
                                             const CELL_LOC location, FieldMetric& dx,
                                             FieldMetric& dy, FieldMetric& dz)
    : mesh(mesh), intShiftTorsion(intShiftTorsion), location(location), dx(dx), dy(dy),
      dz(dz) {}

MetricTensor::FieldMetric DifferentialOperators::DDX(const Field2D& f, CELL_LOC loc,
                                                     const std::string& method,
                                                     const std::string& region) {
  ASSERT1(location == loc || loc == CELL_DEFAULT)
  return bout::derivatives::index::DDX(f, loc, method, region) / dx;
}

Field3D DifferentialOperators::DDX(const Field3D& f, CELL_LOC outloc,
                                   const std::string& method, const std::string& region) {
  auto result = bout::derivatives::index::DDX(f, outloc, method, region);
  result /= dx;

  if (f.getMesh()->IncIntShear) {
    // Using BOUT-06 style shifting
    result += intShiftTorsion * DDZ(f, outloc, method, region);
  }

  return result;
}

FieldMetric DifferentialOperators::DDY(const Field2D& f, CELL_LOC loc,
                                       const std::string& method,
                                       const std::string& region) const {
  ASSERT1(location == loc || loc == CELL_DEFAULT)
  return bout::derivatives::index::DDY(f, loc, method, region) / dy;
}

Field3D DifferentialOperators::DDY(const Field3D& f, CELL_LOC outloc,
                                   const std::string& method,
                                   const std::string& region) const {
#if BOUT_USE_METRIC_3D
  if (!f.hasParallelSlices() and !transform->canToFromFieldAligned()) {
    Field3D f_parallel = f;
    transform->calcParallelSlices(f_parallel);
    f_parallel.applyParallelBoundary("parallel_neumann");
    return bout::derivatives::index::DDY(f_parallel, outloc, method, region);
  }
#endif
  return bout::derivatives::index::DDY(f, outloc, method, region) / dy;
}

FieldMetric DifferentialOperators::DDZ(const Field2D& f, CELL_LOC loc,
                                       const std::string& UNUSED(method),
                                       const std::string& UNUSED(region)) {
  ASSERT1(location == loc || loc == CELL_DEFAULT)
  ASSERT1(f.getMesh() == mesh)
  if (loc == CELL_DEFAULT) {
    loc = f.getLocation();
  }
  return zeroFrom(f).setLocation(loc);
}
Field3D DifferentialOperators::DDZ(const Field3D& f, CELL_LOC outloc,
                                   const std::string& method, const std::string& region) {
  return bout::derivatives::index::DDZ(f, outloc, method, region) / dz;
}

/////////////////////////////////////////////////////////
// Parallel gradient

FieldMetric DifferentialOperators::Grad_par(const Field2D& var,
                                            const MetricTensor& covariantMetricTensor,
                                            [[maybe_unused]] CELL_LOC outloc,
                                            const std::string& UNUSED(method)) {
  TRACE("DifferentialOperators::Grad_par( Field2D )");

  return DDY(var) * invSg(covariantMetricTensor);
}

Field3D DifferentialOperators::Grad_par(const Field3D& var,
                                        const MetricTensor& covariantMetricTensor,
                                        CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Grad_par( Field3D )");

  return DDY(var, outloc, method) * invSg(covariantMetricTensor);
}

/////////////////////////////////////////////////////////
// Vpar_Grad_par
// vparallel times the parallel derivative along unperturbed B-field

Field2D DifferentialOperators::Vpar_Grad_par(const Field2D& v, const Field2D& f,
                                             const MetricTensor& covariantMetricTensor,
                                             [[maybe_unused]] CELL_LOC outloc,
                                             const std::string& UNUSED(method)) {
  return VDDY(v, f) * invSg(covariantMetricTensor);
}

Field3D DifferentialOperators::Vpar_Grad_par(const Field3D& v, const Field3D& f,
                                             const MetricTensor& covariantMetricTensor,
                                             CELL_LOC outloc, const std::string& method) {

  return VDDY(v, f, outloc, method) * invSg(covariantMetricTensor);
}

/////////////////////////////////////////////////////////
// Parallel divergence

FieldMetric DifferentialOperators::Div_par(const Field2D& f, const Field2D& Bxy,
                                           const MetricTensor& covariantMetricTensor,
                                           CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Div_par( Field2D )");

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy();

  return Bxy * Grad_par(f / Bxy_floc, covariantMetricTensor, outloc, method);
}

Field3D DifferentialOperators::Div_par(const Field3D& f, const Field2D& Bxy,
                                       const MetricTensor& covariantMetricTensor,
                                       CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Div_par( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy();

  if (!f.hasParallelSlices()) {
    // No yup/ydown fields. The Grad_par operator will
    // shift to field aligned coordinates
    return Bxy * Grad_par(f / Bxy_floc, covariantMetricTensor, outloc, method);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  f_B.splitParallelSlices();
  for (int i = 0; i < f.getMesh()->ystart; ++i) {
    f_B.yup(i) = f.yup(i) / Bxy_floc.yup(i);
    f_B.ydown(i) = f.ydown(i) / Bxy_floc.ydown(i);
  }
  return Bxy * Grad_par(f_B, covariantMetricTensor, outloc, method);
}

/////////////////////////////////////////////////////////
// second parallel derivative (b dot Grad)(b dot Grad)
// Note: For parallel Laplacian use Laplace_par

FieldMetric DifferentialOperators::Grad2_par2(const Field2D& f,
                                              MetricTensor& covariantMetricTensor,
                                              CELL_LOC outloc,
                                              const std::string& method) {
  TRACE("DifferentialOperators::Grad2_par2( Field2D )");
  ASSERT1(location == outloc || (outloc == CELL_DEFAULT && location == f.getLocation()))

  auto result =
      Grad2_par2_DDY_invSg(covariantMetricTensor, outloc, method) * DDY(f, outloc, method)
      + D2DY2(f, outloc, method) / covariantMetricTensor.Getg22();

  return result;
}

Field3D DifferentialOperators::Grad2_par2(const Field3D& f,
                                          MetricTensor& covariantMetricTensor,
                                          CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Grad2_par2( Field3D )");
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }
  ASSERT1(location == outloc)

  Field3D result = ::DDY(f, outloc, method);

  Field3D r2 = D2DY2(f, outloc, method) / covariantMetricTensor.Getg22();

  result = Grad2_par2_DDY_invSg(covariantMetricTensor, outloc, method) * result + r2;

  ASSERT2(result.getLocation() == outloc)

  return result;
}

///////////////////////////////////////////////////////////
//// perpendicular Laplacian operator
//
//#include <bout/invert_laplace.hxx> // Delp2 uses same coefficients as inversion code
//
//FieldMetric DifferentialOperators::Delp2(const Field2D& f, const Field2D& g11,
//                                         const Field2D& G1, CELL_LOC outloc,
//                                         bool UNUSED(useFFT)) {
//  TRACE("DifferentialOperators::Delp2( Field2D )");
//  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
//
//  auto result = G1 * DDX(f, outloc) + g11 * D2DX2(f, outloc);
//
//  return result;
//}

//Field3D DifferentialOperators::Delp2(const Field3D& f,
//                                     MetricTensor& covariantMetricTensor,
//                                     const Field3D& G1, const Field3D& G3,
//                                     CELL_LOC outloc, bool useFFT) {
//  TRACE("DifferentialOperators::Delp2( Field3D )");
//
//  if (outloc == CELL_DEFAULT) {
//    outloc = f.getLocation();
//  }
//
//  ASSERT1(location == outloc)
//  ASSERT1(f.getLocation() == outloc)
//
//  if (mesh->GlobalNx == 1 && mesh->GlobalNz == 1) {
//    // copy mesh, location, etc
//    return f * 0;
//  }
//  ASSERT2(mesh->xstart > 0) // Need at least one guard cell
//
//  Field3D result{emptyFrom(f).setLocation(outloc)};
//
//  if (useFFT and not bout::build::use_metric_3d) {
//    int ncz = mesh->LocalNz;
//
//    // Allocate memory
//    auto ft = Matrix<dcomplex>(mesh->LocalNx, ncz / 2 + 1);
//    auto delft = Matrix<dcomplex>(mesh->LocalNx, ncz / 2 + 1);
//
//    // Loop over y indices
//    // Note: should not include y-guard or y-boundary points here as that would
//    // use values from corner cells in dx, which may not be initialised.
//    for (int jy = mesh->ystart; jy <= mesh->yend; jy++) {
//
//      // Take forward FFT
//
//      for (int jx = 0; jx < mesh->LocalNx; jx++) {
//        rfft(&f(jx, jy, 0), ncz, &ft(jx, 0));
//      }
//
//      // Loop over kz
//      for (int jz = 0; jz <= ncz / 2; jz++) {
//
//        // No smoothing in the x direction
//        for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
//          // Perform x derivative
//
//          dcomplex a, b, c;
//          laplace_tridag_coefs(jx, jy, jz, a, b, c, nullptr, nullptr, outloc);
//
//          delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
//        }
//      }
//
//      // Reverse FFT
//      for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
//
//        irfft(&delft(jx, 0), ncz, &result(jx, jy, 0));
//      }
//    }
//  } else {
//    result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc)
//             + covariantMetricTensor.Getg11() * ::D2DX2(f, outloc)
//             + covariantMetricTensor.Getg33() * ::D2DZ2(f, outloc)
//             + 2 * covariantMetricTensor.Getg13() * ::D2DXDZ(f, outloc);
//  }
//
//  ASSERT2(result.getLocation() == outloc)
//
//  return result;
//}

//FieldPerp DifferentialOperators::Delp2(const FieldPerp& f, CELL_LOC outloc, bool useFFT) {
//  TRACE("DifferentialOperators::Delp2( FieldPerp )");
//
//  if (outloc == CELL_DEFAULT) {
//    outloc = f.getLocation();
//  }
//
//  ASSERT1(location == outloc)
//  ASSERT1(f.getLocation() == outloc)
//
//  if (mesh->GlobalNx == 1 && mesh->GlobalNz == 1) {
//    // copy mesh, location, etc
//    return f * 0;
//  }
//  ASSERT2(mesh->xstart > 0) // Need at least one guard cell
//
//  FieldPerp result{emptyFrom(f).setLocation(outloc)};
//
//  int jy = f.getIndex();
//  result.setIndex(jy);
//
//  if (useFFT) {
//    int ncz = mesh->LocalNz;
//
//    // Allocate memory
//    auto ft = Matrix<dcomplex>(mesh->LocalNx, ncz / 2 + 1);
//    auto delft = Matrix<dcomplex>(mesh->LocalNx, ncz / 2 + 1);
//
//    // Take forward FFT
//    for (int jx = 0; jx < mesh->LocalNx; jx++) {
//      rfft(&f(jx, 0), ncz, &ft(jx, 0));
//    }
//
//    // Loop over kz
//    for (int jz = 0; jz <= ncz / 2; jz++) {
//
//      // No smoothing in the x direction
//      for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
//        // Perform x derivative
//
//        dcomplex a, b, c;
//        laplace_tridag_coefs(jx, jy, jz, a, b, c);
//
//        delft(jx, jz) = a * ft(jx - 1, jz) + b * ft(jx, jz) + c * ft(jx + 1, jz);
//      }
//    }
//
//    // Reverse FFT
//    for (int jx = mesh->xstart; jx <= mesh->xend; jx++) {
//      irfft(&delft(jx, 0), ncz, &result(jx, 0));
//    }
//
//  } else {
//    throw BoutException("Non-fourier Delp2 not currently implented for FieldPerp.");
//    // Would be the following but don't have standard derivative operators for FieldPerps
//    // yet
//    // result = G1 * ::DDX(f, outloc) + G3 * ::DDZ(f, outloc) + g11 * ::D2DX2(f, outloc)
//    //          + g33 * ::D2DZ2(f, outloc) + 2 * g13 * ::D2DXDZ(f, outloc);
//  }
//
//  return result;
//}

FieldMetric DifferentialOperators::Laplace_par(const Field2D& f, const Field2D& g22,
                                               const Field2D& J, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  return D2DY2(f, outloc) / g22 + DDY(J / g22, outloc) * DDY(f, outloc) / J;
}

Field3D DifferentialOperators::Laplace_par(const Field3D& f, const Field3D& g22,
                                           const Field3D& J, CELL_LOC outloc) {
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)
  return D2DY2(f, outloc) / g22 + DDY(J / g22, outloc) * ::DDY(f, outloc) / J;
}

// Full Laplacian operator on scalar field

FieldMetric DifferentialOperators::Laplace(const Field2D& f,
                                           MetricTensor& covariantMetricTensor,
                                           const Field2D& G1, const Field2D& G2,
                                           const Field2D& G3, CELL_LOC outloc,
                                           const std::string& dfdy_boundary_conditions,
                                           const std::string& dfdy_dy_region) {
  TRACE("DifferentialOperators::Laplace( Field2D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  auto result = G1 * DDX(f, outloc) + G2 * DDY(f, outloc)
                + covariantMetricTensor.Getg11() * D2DX2(f, outloc)
                + covariantMetricTensor.Getg22() * D2DY2(f, outloc)
                + 2.0 * covariantMetricTensor.Getg12()
                      * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY",
                               dfdy_boundary_conditions, dfdy_dy_region);

  return result;
}

Field3D DifferentialOperators::Laplace(const Field3D& f,
                                       MetricTensor& covariantMetricTensor,
                                       const Field3D& G1, const Field3D& G2,
                                       const Field3D& G3, CELL_LOC outloc,
                                       const std::string& dfdy_boundary_conditions,
                                       const std::string& dfdy_dy_region) {
  TRACE("DifferentialOperators::Laplace( Field3D )");
  ASSERT1(location == outloc || outloc == CELL_DEFAULT)

  Field3D result = G1 * ::DDX(f, outloc) + G2 * ::DDY(f, outloc) + G3 * ::DDZ(f, outloc)
                   + covariantMetricTensor.Getg11() * D2DX2(f, outloc)
                   + covariantMetricTensor.Getg22() * D2DY2(f, outloc)
                   + covariantMetricTensor.Getg33() * D2DZ2(f, outloc)
                   + 2.0
                         * (covariantMetricTensor.Getg12()
                                * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY",
                                         dfdy_boundary_conditions, dfdy_dy_region)
                            + covariantMetricTensor.Getg13() * D2DXDZ(f, outloc)
                            + covariantMetricTensor.Getg23() * D2DYDZ(f, outloc));

  return result;
}

// Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
// solver
Field2D DifferentialOperators::Laplace_perpXY([[maybe_unused]] const Field2D& A,
                                              [[maybe_unused]] const Field2D& f), 
                                              MetricTensor& covariantMetricTensor,
                                              const Field2D& J) {
  TRACE("DifferentialOperators::Laplace_perpXY( Field2D )");
#if not(BOUT_USE_METRIC_3D)
  Field2D result;
  result.allocate();
  for (auto i : result.getRegion(RGN_NOBNDRY)) {
    result[i] = 0.;

    // outer x boundary
    const auto outer_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xp()]); };
    const BoutReal outer_x_A = outer_x_avg(A);
    const BoutReal outer_x_J = outer_x_avg(J);
    const BoutReal outer_x_g11 = outer_x_avg(covariantMetricTensor.Getg11());
    const BoutReal outer_x_dx = outer_x_avg(dx);
    const BoutReal outer_x_value =
        outer_x_A * outer_x_J * outer_x_g11 / (J[i] * outer_x_dx * dx[i]);
    result[i] += outer_x_value * (f[i.xp()] - f[i]);

    // inner x boundary
    const auto inner_x_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.xm()]); };
    const BoutReal inner_x_A = inner_x_avg(A);
    const BoutReal inner_x_J = inner_x_avg(J);
    const BoutReal inner_x_g11 = inner_x_avg(covariantMetricTensor.Getg11());
    const BoutReal inner_x_dx = inner_x_avg(dx);
    const BoutReal inner_x_value =
        inner_x_A * inner_x_J * inner_x_g11 / (J[i] * inner_x_dx * dx[i]);
    result[i] += inner_x_value * (f[i.xm()] - f[i]);

    // upper y boundary
    const auto upper_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.yp()]); };
    const BoutReal upper_y_A = upper_y_avg(A);
    const BoutReal upper_y_J = upper_y_avg(J);
    const BoutReal upper_y_g_22 = upper_y_avg(covariantMetricTensor.Getg22());
    const BoutReal upper_y_g23 = upper_y_avg(covariantMetricTensor.Getg23());
    const BoutReal upper_y_g_23 = upper_y_avg(covariantMetricTensor.Getg23());
    const BoutReal upper_y_dy = upper_y_avg(dy);
    const BoutReal upper_y_value = -upper_y_A * upper_y_J * upper_y_g23 * upper_y_g_23
                                   / (upper_y_g_22 * J[i] * upper_y_dy * dy[i]);
    result[i] += upper_y_value * (f[i.yp()] - f[i]);

    // lower y boundary
    const auto lower_y_avg = [&i](const auto& f) { return 0.5 * (f[i] + f[i.ym()]); };
    const BoutReal lower_y_A = lower_y_avg(A);
    const BoutReal lower_y_J = lower_y_avg(J);
    const BoutReal lower_y_g_22 = lower_y_avg(covariantMetricTensor.Getg22());
    const BoutReal lower_y_g23 = lower_y_avg(covariantMetricTensor.Getg23());
    const BoutReal lower_y_g_23 = lower_y_avg(covariantMetricTensor.Getg23());
    const BoutReal lower_y_dy = lower_y_avg(dy);
    const BoutReal lower_y_value = -lower_y_A * lower_y_J * lower_y_g23 * lower_y_g_23
                                   / (lower_y_g_22 * J[i] * lower_y_dy * dy[i]);
    result[i] += lower_y_value * (f[i.ym()] - f[i]);
  }

  return result;
#else
  throw BoutException(
      "DifferentialOperators::Laplace_perpXY for 3D metric not implemented");
#endif
}

FieldMetric&
DifferentialOperators::invSg(const MetricTensor& covariantMetricTensor) const {
  if (invSgCache == nullptr) {
    auto ptr = std::make_unique<FieldMetric>();
    (*ptr) = 1.0 / sqrt(covariantMetricTensor.Getg22());
    invSgCache = std::move(ptr);
  }
  return *invSgCache;
}

const FieldMetric&
DifferentialOperators::Grad2_par2_DDY_invSg(const MetricTensor& covariantMetricTensor,
                                            CELL_LOC outloc,
                                            const std::string& method) const {
  if (auto search = Grad2_par2_DDY_invSgCache.find(method);
      search != Grad2_par2_DDY_invSgCache.end()) {
    return *search->second;
  }
  invSg(covariantMetricTensor);

  // Communicate to get parallel slices
  mesh->communicate(*invSgCache);
  invSgCache->applyParallelBoundary("parallel_neumann");

  // cache
  auto ptr = std::make_unique<FieldMetric>();
  *ptr = DDY(*invSgCache, outloc, method) * invSg(covariantMetricTensor);
  Grad2_par2_DDY_invSgCache[method] = std::move(ptr);
  return *Grad2_par2_DDY_invSgCache[method];
}

void DifferentialOperators::invalidateAndRecalculateCachedVariables() {
  Grad2_par2_DDY_invSgCache.clear();
  invSgCache.reset();
}
