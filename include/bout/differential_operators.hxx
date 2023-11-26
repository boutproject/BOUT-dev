
#ifndef BOUT_DIFFERENTIALOPERATORS_HXX
#define BOUT_DIFFERENTIALOPERATORS_HXX

#include "bout/metricTensor.hxx"
//#include "bout/index_derivs_interface.hxx"

class DifferentialOperators {

  using FieldMetric = MetricTensor::FieldMetric;

public:
  DifferentialOperators(Mesh* mesh, FieldMetric intShiftTorsion, CELL_LOC location,
                        FieldMetric& dx, FieldMetric& dy, FieldMetric& dz);

  //  DifferentialOperators(DifferentialOperators operators,
  //                        DifferentialOperators::FieldMetric& dx);

  FieldMetric DDX(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY");

  FieldMetric DDY(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY") const;

  FieldMetric DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& method = "DEFAULT",
                  const std::string& region = "RGN_NOBNDRY");

  Field3D DDX(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  Field3D DDY(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDZ(const Field3D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY");

  /// Gradient along magnetic field  b.Grad(f)
  FieldMetric Grad_par(const Field2D& var, const MetricTensor& covariantMetricTensor,
                       CELL_LOC outloc = CELL_DEFAULT,
                       const std::string& method = "DEFAULT");

  Field3D Grad_par(const Field3D& var, const MetricTensor& covariantMetricTensor,
                   CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT");

  /// Advection along magnetic field V*b.Grad(f)
  Field2D Vpar_Grad_par(const Field2D& v, const Field2D& f,
                        const MetricTensor& covariantMetricTensor,
                        CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                        const MetricTensor& covariantMetricTensor,
                        CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  FieldMetric Div_par(const Field2D& f, const Field2D& Bxy,
                      const MetricTensor& covariantMetricTensor,
                      CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& method = "DEFAULT");

  Field3D Div_par(const Field3D& f, const Field2D& Bxy,
                  const MetricTensor& covariantMetricTensor,
                  CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT");

  // Second derivative along magnetic field
  FieldMetric Grad2_par2(const Field2D& f, MetricTensor& covariantMetricTensor,
                         CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT");

  Field3D Grad2_par2(const Field3D& f, MetricTensor& covariantMetricTensor,
                     CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");

  //  // Perpendicular Laplacian operator, using only X-Z derivatives
  //  // NOTE: This might be better bundled with the Laplacian inversion code
  //  // since it makes use of the same coefficients and FFT routines
  //  FieldMetric Delp2(const Field2D& f, const Field2D& g11, const Field2D& G1,
  //                    CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  //
  //  Field3D Delp2(const Field3D& f, MetricTensor& covariantMetricTensor, const Field3D& G1,
  //                const Field3D& G3, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);
  //
  //  FieldPerp Delp2(const FieldPerp& f, CELL_LOC outloc = CELL_DEFAULT, bool useFFT = true);

  // Full parallel Laplacian operator on scalar field
  // Laplace_par(f) = Div( b (b dot Grad(f)) )
  FieldMetric Laplace_par(const Field2D& f, const Field2D& g22, const Field2D& J,
                          CELL_LOC outloc = CELL_DEFAULT);

  Field3D Laplace_par(const Field3D& f, const Field3D& g22, const Field3D& J,
                      CELL_LOC outloc = CELL_DEFAULT);

  // Full Laplacian operator on scalar field
  FieldMetric Laplace(const Field2D& f, MetricTensor& covariantMetricTensor,
                      const Field2D& G1, const Field2D& G2, const Field2D& G3,
                      CELL_LOC outloc = CELL_DEFAULT,
                      const std::string& dfdy_boundary_conditions = "free_o3",
                      const std::string& dfdy_dy_region = "");

  Field3D Laplace(const Field3D& f, MetricTensor& covariantMetricTensor,
                  const Field3D& G1, const Field3D& G2, const Field3D& G3,
                  CELL_LOC outloc = CELL_DEFAULT,
                  const std::string& dfdy_boundary_conditions = "free_o3",
                  const std::string& dfdy_dy_region = "");

  // Full perpendicular Laplacian, in form of inverse of Laplacian operator in LaplaceXY
  // solver
  Field2D Laplace_perpXY(const Field2D& A, const Field2D& f,
                         MetricTensor& covariantMetricTensor, const Field2D& J);

  void invalidateAndRecalculateCachedVariables();

  const FieldMetric& Grad2_par2_DDY_invSg(const MetricTensor& covariantMetricTensor,
                                          CELL_LOC outloc,
                                          const std::string& method) const;

private:
  Mesh* mesh;
  FieldMetric intShiftTorsion;
  CELL_LOC location;
  FieldMetric& dx;
  FieldMetric& dy;
  FieldMetric& dz;

  /// Cache variable for Grad2_par2
  mutable std::map<std::string, std::unique_ptr<FieldMetric>> Grad2_par2_DDY_invSgCache;
  mutable std::unique_ptr<FieldMetric> invSgCache{nullptr};

  FieldMetric& invSg(const MetricTensor& covariantMetricTensor) const;
};

#endif //BOUT_DIFFERENTIALOPERATORS_HXX
