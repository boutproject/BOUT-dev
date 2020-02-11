#include <bout/snb.hxx>
#include <bout.hxx>
#include <boutexception.hxx>
#include <field_factory.hxx>
#include <output.hxx>

// Convert __LINE__ to string S__LINE__
#define S(x) #x
#define S_(x) S(x)
#define S__LINE__ S_(__LINE__)

#define EXPECT_TRUE(expr)                                                       \
  if (!expr) {                                                                  \
    throw BoutException("Line " S__LINE__ " Expected true, got false: " #expr); \
  }

#define EXPECT_FALSE(expr)                                                      \
  if (expr) {                                                                   \
    throw BoutException("Line " S__LINE__ " Expected false, got true: " #expr); \
  }

/// Is \p field equal to \p reference, with a tolerance of \p tolerance?
template <class T, class U>
bool IsFieldEqual(const T& field, const U& reference,
                  const std::string& region = "RGN_ALL", BoutReal tolerance = 1e-10) {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference[i]) > tolerance) {
      output.write("Field: {:e}, reference: {:e}, tolerance: {:e}\n", field[i], reference[i],
                   tolerance);
      return false;
    }
  }
  return true;
}

/// Is \p field equal to \p reference, with a tolerance of \p tolerance?
/// Overload for BoutReals
template <class T>
bool IsFieldEqual(const T& field, BoutReal reference,
                  const std::string& region = "RGN_ALL", BoutReal tolerance = 1e-10) {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference) > tolerance) {
      output.write("Field: {:e}, reference: {:e}, tolerance: {:e}\n", field[i], reference,
                   tolerance);
      return false;
    }
  }
  return true;
}

/// Is \p field close to \p reference, with a relative tolerance of \p tolerance?
template <class T, class U>
bool IsFieldClose(const T& field, const U& reference,
                  const std::string& region = "RGN_ALL", BoutReal tolerance = 1e-4) {
  for (auto i : field.getRegion(region)) {
    if (fabs(field[i] - reference[i])
        > tolerance * (fabs(reference[i]) + fabs(field[i]))) {
      output.write("Field: {:e}, reference: {:e}, tolerance: {:e}\n", field[i], reference[i],
                   tolerance * (fabs(reference[i]) + fabs(field[i])));
      return false;
    }
  }
  return true;
}

int main(int argc, char** argv) {
  using bout::HeatFluxSNB;

  BoutInitialise(argc, argv);

  ///////////////////////////////////////////////////////////
  // When there is a temperature gradient the flux is nonzero

  {
    FieldFactory factory;
    auto Te = factory.create3D("5 + cos(y)");
    auto Ne = factory.create3D("1e18 * (1 + 0.5*sin(y))");

    mesh->communicate(Te, Ne);

    HeatFluxSNB snb;

    Field3D Div_q_SH;
    Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

    // Check that flux is not zero
    EXPECT_FALSE(IsFieldEqual(Div_q_SH, 0.0, "RGN_NOBNDRY"));
    EXPECT_FALSE(IsFieldEqual(Div_q, 0.0, "RGN_NOBNDRY"));
  }

  ///////////////////////////////////////////////////////////
  // When the temperature is constant there is no flux

  {
    FieldFactory factory;
    auto Te = factory.create3D("1.5");
    auto Ne = factory.create3D("1e18 * (1 + 0.5*sin(y))");

    mesh->communicate(Te, Ne);

    HeatFluxSNB snb;

    Field3D Div_q_SH;
    Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

    // Check that flux is zero
    EXPECT_TRUE(IsFieldEqual(Div_q_SH, 0.0, "RGN_NOBNDRY"));
    EXPECT_TRUE(IsFieldEqual(Div_q, 0.0, "RGN_NOBNDRY"));
  }

  ///////////////////////////////////////////////////////////
  // In the collisional limit the SH and SNB fluxes are close

  {
    FieldFactory factory;
    auto Te = factory.create3D("1 + 0.01*sin(y)");
    auto Ne = factory.create3D("1e20 * (1 + 0.5*sin(y))");
    mesh->communicate(Te, Ne);

    HeatFluxSNB snb;

    Field3D Div_q_SH;
    Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

    // Check that flux is zero
    EXPECT_TRUE(IsFieldClose(Div_q, Div_q_SH, "RGN_NOBNDRY"));
  }

  ///////////////////////////////////////////////////////////
  // In the collisionless limit the SH and SNB fluxes are different

  {
    FieldFactory factory;
    auto Te = factory.create3D("1e3 + 0.01*sin(y)");
    auto Ne = factory.create3D("1e19 * (1 + 0.5*sin(y))");
    mesh->communicate(Te, Ne);

    HeatFluxSNB snb;

    Field3D Div_q_SH;
    Field3D Div_q = snb.divHeatFlux(Te, Ne, &Div_q_SH);

    // Check that fluxes are not equal
    EXPECT_FALSE(IsFieldClose(Div_q, Div_q_SH, "RGN_NOBNDRY"));
  }

  ///////////////////////////////////////////////////////////
  // Reversing the temperature gradient reverses the flux

  {
    Field3D Ne = 1e19;

    FieldFactory factory;
    auto Te = factory.create3D("10 + 0.01*sin(y)");
    mesh->communicate(Te);

    HeatFluxSNB snb;

    Field3D Div_q_SH_1;
    Field3D Div_q_1 = snb.divHeatFlux(Te, Ne, &Div_q_SH_1);

    auto Te2 = factory.create3D("10 - 0.01*sin(y)");
    mesh->communicate(Te2);

    Field3D Div_q_SH_2;
    Field3D Div_q_2 = snb.divHeatFlux(Te2, Ne, &Div_q_SH_2);

    // Check that fluxes are reversed in y
    for (int y = mesh->ystart; y <= mesh->yend; y++) {
      if (fabs(Div_q_SH_2(0, y, 0) - Div_q_SH_1(0, mesh->yend - y + mesh->ystart, 0))
          > 1e-6 * (fabs(Div_q_SH_2(0, y, 0)) + fabs(Div_q_SH_1(0, y, 0)))) {
        throw BoutException("SH: y = {:d}: {:e} != {:e}", y, Div_q_SH_2(0, y, 0),
                            Div_q_SH_1(0, mesh->yend - y + mesh->ystart, 0));
      }

      if (fabs(Div_q_2(0, y, 0) - Div_q_1(0, mesh->yend - y + mesh->ystart, 0))
          > 1e-6 * (fabs(Div_q_2(0, y, 0)) + fabs(Div_q_1(0, y, 0)))) {
        throw BoutException("SNB: y = {:d}: {:e} != {:e}", y, Div_q_2(0, y, 0),
                            Div_q_1(0, mesh->yend - y + mesh->ystart, 0));
      }
    }
  }

  BoutFinalise();

  output << "All tests passed\n";

  return 0;
}
