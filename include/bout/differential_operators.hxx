
#ifndef BOUT_DIFFERENTIALOPERATORS_HXX
#define BOUT_DIFFERENTIALOPERATORS_HXX

#include "bout/metricTensor.hxx"

class DifferentialOperators {

  using FieldMetric = MetricTensor::FieldMetric;

public:
  DifferentialOperators(Mesh* mesh);

  Field2D DDX(const Field2D& f, const Field2D& dx, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field2D DDY(const Field2D& f, const Field2D& dy, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field2D DDZ(const Field2D& f, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDX(const Field3D& f, const Field3D& dx, const Field3D& dz,
              const FieldMetric& intShiftTorsion, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDY(const Field3D& f, const Field3D& dy, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  Field3D DDZ(const Field3D& f, const Field3D& dz, CELL_LOC outloc = CELL_DEFAULT,
              const std::string& method = "DEFAULT",
              const std::string& region = "RGN_NOBNDRY") const;

  /// Gradient along magnetic field  b.Grad(f)
  Field2D Grad_par(const Field2D& var, const Field2D& dy, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT");

  Field3D Grad_par(const Field3D& var, const Field3D& dy, CELL_LOC outloc = CELL_DEFAULT,
                   const std::string& method = "DEFAULT");

  /// Advection along magnetic field V*b.Grad(f)
  Field2D Vpar_Grad_par(const Field2D& v, const Field2D& f,
                        const Coordinates* coordinates, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  Field3D Vpar_Grad_par(const Field3D& v, const Field3D& f,
                        const Coordinates* coordinates, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT");

  /// Divergence along magnetic field  Div(b*f) = B.Grad(f/B)
  Field2D Div_par(const Field2D& f, const Field2D& Bxy, const Field2D& dy,
                  CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT");

  Field3D Div_par(const Field3D& f, const Field3D& Bxy, const Field3D& dy,
                  CELL_LOC outloc = CELL_DEFAULT, const std::string& method = "DEFAULT");

  // Second derivative along magnetic field
  Field2D Grad2_par2(const Field2D& f, const Field2D& dy,
                     MetricTensor& covariantMetricTensor, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");

  Field3D Grad2_par2(const Field3D& f, const FieldMetric& dy,
                     MetricTensor& covariantMetricTensor, CELL_LOC outloc = CELL_DEFAULT,
                     const std::string& method = "DEFAULT");

  Field2D Laplace_par(const Field2D& f, const Field2D& g22, const Field2D& J,
                      const Field2D& dy, CELL_LOC outloc = CELL_DEFAULT);

  Field3D Laplace_par(const Field3D& f, const Field3D& g22, const Field3D& J,
                      const Field2D& dy, CELL_LOC outloc = CELL_DEFAULT);

  // Full Laplacian operator on scalar field
  Field2D Laplace(const Field2D& f, MetricTensor& covariantMetricTensor,
                  const Field2D& dy, const Field2D& G1, const Field2D& G2,
                  const Field2D& G3, CELL_LOC outloc = CELL_DEFAULT,
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
                         MetricTensor& covariantMetricTensor, const Field2D& J,
                         const Field2D& dx, const Field2D& dy);

private:
  Mesh* mesh;
};

#endif //BOUT_DIFFERENTIALOPERATORS_HXX
