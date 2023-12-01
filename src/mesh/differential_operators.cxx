
#include "bout/differential_operators.hxx"
#include "bout/mesh.hxx"
#include <bout/derivs.hxx>

DifferentialOperators::DifferentialOperators(Mesh* mesh) : mesh(mesh) {}

Field2D DifferentialOperators::DDX(const Field2D& f, const Field2D& dx, CELL_LOC loc,
                                   const std::string& method,
                                   const std::string& region) const {
  return bout::derivatives::index::DDX(f, loc, method, region) / dx;
}

Field3D DifferentialOperators::DDX(const Field3D& f, const Field3D& dx, const Field3D& dz,
                                   const FieldMetric& intShiftTorsion, CELL_LOC outloc,
                                   const std::string& method,
                                   const std::string& region) const {
  auto result = bout::derivatives::index::DDX(f, outloc, method, region);
  result /= dx;

  if (f.getMesh()->IncIntShear) {
    // Using BOUT-06 style shifting
    result += intShiftTorsion * DDZ(f, dz, outloc, method, region);
  }

  return result;
}

Field2D DifferentialOperators::DDY(const Field2D& f, const Field2D& dy, CELL_LOC loc,
                                   const std::string& method,
                                   const std::string& region) const {
  return bout::derivatives::index::DDY(f, loc, method, region) / dy;
}

Field3D DifferentialOperators::DDY(const Field3D& f, const Field3D& dy, CELL_LOC outloc,
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

Field2D DifferentialOperators::DDZ(const Field2D& f, CELL_LOC loc,
                                   const std::string& UNUSED(method),
                                   const std::string& UNUSED(region)) const {
  if (loc == CELL_DEFAULT) {
    loc = f.getLocation();
  }
  return zeroFrom(f).setLocation(loc);
}

Field3D DifferentialOperators::DDZ(const Field3D& f, const Field3D& dz, CELL_LOC outloc,
                                   const std::string& method,
                                   const std::string& region) const {
  return bout::derivatives::index::DDZ(f, outloc, method, region) / dz;
}

/////////////////////////////////////////////////////////
// Parallel gradient

Field2D DifferentialOperators::Grad_par(const Field2D& var, const Field2D& dy,
                                        const MetricTensor& covariantMetricTensor,
                                        [[maybe_unused]] CELL_LOC outloc,
                                        const std::string& UNUSED(method)) {
  TRACE("DifferentialOperators::Grad_par( Field2D )");

  return DDY(var, dy) * invSg(covariantMetricTensor);
}

Field3D DifferentialOperators::Grad_par(const Field3D& var, const Field3D& dy,
                                        const MetricTensor& covariantMetricTensor,
                                        CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Grad_par( Field3D )");

  return DDY(var, dy, outloc, method) * invSg(covariantMetricTensor);
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

Field2D DifferentialOperators::Div_par(const Field2D& f, const Field2D& Bxy,
                                       const Field2D& dy,
                                       const MetricTensor& covariantMetricTensor,
                                       CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Div_par( Field2D )");

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy();

  return Bxy * Grad_par(f / Bxy_floc, dy, covariantMetricTensor, outloc, method);
}

Field3D DifferentialOperators::Div_par(const Field3D& f, const Field3D& Bxy,
                                       const Field3D& dy,
                                       const MetricTensor& covariantMetricTensor,
                                       CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Div_par( Field3D )");

  // Need Bxy at location of f, which might be different from location of this
  // Coordinates object
  auto Bxy_floc = f.getCoordinates()->Bxy();

  if (!f.hasParallelSlices()) {
    // No yup/ydown fields. The Grad_par operator will
    // shift to field aligned coordinates
    return Bxy * Grad_par(f / Bxy_floc, dy, covariantMetricTensor, outloc, method);
  }

  // Need to modify yup and ydown fields
  Field3D f_B = f / Bxy_floc;
  f_B.splitParallelSlices();
  for (int i = 0; i < f.getMesh()->ystart; ++i) {
    f_B.yup(i) = f.yup(i) / Bxy_floc.yup(i);
    f_B.ydown(i) = f.ydown(i) / Bxy_floc.ydown(i);
  }
  return Bxy * Grad_par(f_B, dy, covariantMetricTensor, outloc, method);
}

/////////////////////////////////////////////////////////
// second parallel derivative (b dot Grad)(b dot Grad)
// Note: For parallel Laplacian use Laplace_par

Field2D DifferentialOperators::Grad2_par2(const Field2D& f, const Field2D& dy,
                                          MetricTensor& covariantMetricTensor,
                                          CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Grad2_par2( Field2D )");

  return Grad2_par2_DDY_invSg(covariantMetricTensor, dy, outloc, method)
             * DDY(f, dy, outloc, method)
         + D2DY2(f, outloc, method) / covariantMetricTensor.Getg22();
}

Field3D DifferentialOperators::Grad2_par2(const Field3D& f, const FieldMetric& dy,
                                          MetricTensor& covariantMetricTensor,
                                          CELL_LOC outloc, const std::string& method) {
  TRACE("DifferentialOperators::Grad2_par2( Field3D )");
  if (outloc == CELL_DEFAULT) {
    outloc = f.getLocation();
  }

  Field3D result = ::DDY(f, outloc, method);

  Field3D r2 = D2DY2(f, outloc, method) / covariantMetricTensor.Getg22();

  result = Grad2_par2_DDY_invSg(covariantMetricTensor, dy, outloc, method) * result + r2;

  ASSERT2(result.getLocation() == outloc)

  return result;
}

FieldMetric DifferentialOperators::Laplace_par(const Field2D& f, const Field2D& g22,
                                               const Field2D& J, const Field2D& dy,
                                               CELL_LOC outloc) {
  return D2DY2(f, outloc) / g22 + DDY(J / g22, dy, outloc) * DDY(f, dy, outloc) / J;
}

Field3D DifferentialOperators::Laplace_par(const Field3D& f, const Field3D& g22,
                                           const Field3D& J, const Field2D& dy,
                                           CELL_LOC outloc) {
  return D2DY2(f, outloc) / g22 + DDY(J / g22, dy, outloc) * ::DDY(f, outloc) / J;
}

// Full Laplacian operator on scalar field

FieldMetric DifferentialOperators::Laplace(
    const Field2D& f, MetricTensor& covariantMetricTensor, const Field2D& dy,
    const Field2D& G1, const Field2D& G2, const Field2D& G3, CELL_LOC outloc,
    const std::string& dfdy_boundary_conditions, const std::string& dfdy_dy_region) {
  TRACE("DifferentialOperators::Laplace( Field2D )");

  return G1 * DDX(f, dy, outloc) + G2 * DDY(f, dy, outloc)
         + covariantMetricTensor.Getg11() * D2DX2(f, outloc)
         + covariantMetricTensor.Getg22() * D2DY2(f, outloc)
         + 2.0 * covariantMetricTensor.Getg12()
               * D2DXDY(f, outloc, "DEFAULT", "RGN_NOBNDRY", dfdy_boundary_conditions,
                        dfdy_dy_region);
}

Field3D DifferentialOperators::Laplace(const Field3D& f,
                                       MetricTensor& covariantMetricTensor,
                                       const Field3D& G1, const Field3D& G2,
                                       const Field3D& G3, CELL_LOC outloc,
                                       const std::string& dfdy_boundary_conditions,
                                       const std::string& dfdy_dy_region) {
  TRACE("DifferentialOperators::Laplace( Field3D )");

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
                                              const Field2D& J, const Field2D& dx,
                                              const Field2D& dy) {
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
                                            const FieldMetric& dy, CELL_LOC outloc,
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
  *ptr = DDY(*invSgCache, dy, outloc, method) * invSg(covariantMetricTensor);
  Grad2_par2_DDY_invSgCache[method] = std::move(ptr);
  return *Grad2_par2_DDY_invSgCache[method];
}

void DifferentialOperators::invalidateAndRecalculateCachedVariables() {
  Grad2_par2_DDY_invSgCache.clear();
  invSgCache.reset();
}
