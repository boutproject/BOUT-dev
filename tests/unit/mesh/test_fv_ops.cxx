#include "gtest/gtest.h"
#include <cmath>

#include "bout/build_defines.hxx"
#include "bout/constants.hxx"
#include "bout/coordinates.hxx"
#include "bout/fv_ops.hxx"
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
