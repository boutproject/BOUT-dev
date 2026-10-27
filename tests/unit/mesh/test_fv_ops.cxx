#include "gtest/gtest.h"
#include <cmath>

#include "bout/build_defines.hxx"
#include "bout/constants.hxx"
#include "bout/coordinates.hxx"
#include "bout/fv_ops.hxx"
#include "bout/fv_ops_impl.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx"

#include "fake_mesh.hxx"
#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx"

using bout::globals::mesh;

using FVOpsTest = FakeMeshFixture;

TEST_F(FVOpsTest, Div_a_Grad_perp_FCI_no_coeffs) {
  EXPECT_THROW(FV::Div_a_Grad_perp_FCI::create(mesh, 1), BoutException);
}

TEST_F(FVOpsTest, Div_a_Grad_perp_FCI_coeffs) {
  static_cast<FakeMesh*>(bout::globals::mesh)
      ->setGridDataSource(new FakeGridDataSource({{"dagp_fv_XX", 1.0},
                                                  {"dagp_fv_XZ", 1.0},
                                                  {"dagp_fv_ZX", 1.0},
                                                  {"dagp_fv_ZZ", 1.0},
                                                  {"dagp_fv_volume", 1.0}}));

  auto op = FV::Div_a_Grad_perp_FCI::create(mesh, 1);
  Field3D a{1.0}, f{1.0};
  (*op)(a, f, true);
  (*op)(a, f, false);
}

namespace {

// Two Y guards let the parallel fluxes reconstruct their neighboring guard cell.
// Keep the existing one-guard fixture as the default for other unit tests.
template <typename Limiter>
class FVVHeatingTest : public FakeMeshFixture_tmpl<5, 20, 3, false, 2> {
public:
  struct Fields {
    Field3D density{1.0}, velocity{0.0}, wave{1.0};
  };

  Fields fields(BoutReal flow_sign, bool uniform_velocity = false) {
    Fields result;
    const int ny_physical = mesh->yend - mesh->ystart + 1;
    Coordinates::FieldMetric jacobian{1.0}, dy{1.0};
    for (int x = 0; x < mesh->LocalNx; ++x) {
      for (int y = 0; y < mesh->LocalNy; ++y) {
        const BoutReal theta = TWOPI * (y - mesh->ystart) / ny_physical;
        for (int z = 0; z < mesh->LocalNz; ++z) {
          const BoutReal phase = theta + 0.2 * x + 0.3 * z;
          result.density(x, y, z) = 1.0 + 0.3 * cos(phase);
          result.velocity(x, y, z) =
              uniform_velocity ? 0.6 : flow_sign * (0.4 + 0.6 * sin(2.0 * phase));
          result.wave(x, y, z) = 0.8 + 0.1 * cos(phase + 0.2);
          jacobian(x, y, z) = 1.0 + 0.2 * cos(theta);
          dy(x, y, z) = 0.7 + 0.1 * sin(theta);
        }
      }
    }
    mesh->getCoordinates()->setJ(jacobian, false);
    mesh->getCoordinates()->setDy(dy, false);
    return result;
  }

  void checkEnergyBalance(bool periodic, bool fix_momentum_flux, BoutReal flow_sign) {
    static_cast<FakeMesh*>(mesh)->ix_separatrix = periodic ? mesh->LocalNx : 0;
    const auto f = fields(flow_sign);
    const auto* coord = mesh->getCoordinates();
    Field3D particle_flow, kinetic_flow;
    // The density operator fixes its physical boundary flux to midpoint n*v.
    const auto div_particle =
        FV::Div_par_mod<Limiter>(f.density, f.velocity, f.wave, particle_flow, true);
    const auto div_momentum =
        FV::Div_par_fvv<Limiter>(f.density, f.velocity, f.wave, fix_momentum_flux);
    const auto heating = FV::Div_par_fvv_heating<Limiter>(
        f.density, f.velocity, f.wave, kinetic_flow, fix_momentum_flux);

    BoutReal integrated_change = 0.0, outward = 0.0, heating_magnitude = 0.0;
    for (int x = mesh->xstart; x <= mesh->xend; ++x) {
      for (int z = 0; z < mesh->LocalNz; ++z) {
        outward += kinetic_flow(x, mesh->yend + 1, z) - kinetic_flow(x, mesh->ystart, z);
        if (!periodic) {
          // Both momentum boundary modes should give the prescribed midpoint
          // kinetic energy flux, with heating compensating their different work.
          for (int face : {mesh->ystart, mesh->yend + 1}) {
            const BoutReal n_mid =
                0.5 * (f.density(x, face - 1, z) + f.density(x, face, z));
            const BoutReal v_mid =
                0.5 * (f.velocity(x, face - 1, z) + f.velocity(x, face, z));
            const BoutReal area = face == mesh->ystart
                                      ? coord->cell_area_ylow()(x, face, z)
                                      : coord->cell_area_yhigh()(x, face - 1, z);
            EXPECT_NEAR(kinetic_flow(x, face, z),
                        0.5 * n_mid * v_mid * v_mid * v_mid * area, 1e-13);
          }
        }
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          SCOPED_TRACE(::testing::Message()
                       << "x=" << x << " y=" << y << " z=" << z << " sign=" << flow_sign);
          const BoutReal v = f.velocity(x, y, z);
          // d(n*v^2/2)/dt = v*d(n*v)/dt - v^2*dn/dt/2, for unit mass.
          const BoutReal kinetic_rhs =
              -v * div_momentum(x, y, z) + 0.5 * v * v * div_particle(x, y, z);
          const BoutReal volume = coord->cell_volume()(x, y, z);
          // kinetic_flow is already area integrated, and is indexed at y-low.
          const BoutReal divergence =
              (kinetic_flow(x, y + 1, z) - kinetic_flow(x, y, z)) / volume;
          EXPECT_NEAR(kinetic_rhs + heating(x, y, z) + divergence, 0.0, 1e-13);
          integrated_change += (kinetic_rhs + heating(x, y, z)) * volume;
          heating_magnitude += std::abs(heating(x, y, z)) * volume;
        }
      }
    }
    EXPECT_GT(heating_magnitude, 1e-4);
    EXPECT_NEAR(integrated_change + outward, 0.0, 2e-12);
    if (periodic) {
      EXPECT_NEAR(outward, 0.0, 1e-13);
    }
  }
};

using HeatingLimiters = ::testing::Types<FV::MC, FV::MinMod, FV::Upwind>;
TYPED_TEST_SUITE(FVVHeatingTest, HeatingLimiters);

TYPED_TEST(FVVHeatingTest, PeriodicKineticAndThermalEnergyBalance) {
  this->checkEnergyBalance(true, false, 1.0);
}

TYPED_TEST(FVVHeatingTest, OpenMidpointMomentumEnergyBalance) {
  for (BoutReal sign : {1.0, -1.0}) {
    this->checkEnergyBalance(false, true, sign);
  }
}

TYPED_TEST(FVVHeatingTest, OpenReconstructedMomentumEnergyBalance) {
  for (BoutReal sign : {1.0, -1.0}) {
    this->checkEnergyBalance(false, false, sign);
  }
}

TYPED_TEST(FVVHeatingTest, UniformVelocityHasNoPeriodicHeating) {
  const auto f = this->fields(1.0, true);
  Field3D flow;
  const auto heating =
      FV::Div_par_fvv_heating<TypeParam>(f.density, f.velocity, f.wave, flow);
  BOUT_FOR_SERIAL(i, heating.getRegion("RGN_NOBNDRY")) {
    EXPECT_DOUBLE_EQ(heating[i], 0.0);
  }
}

} // namespace
