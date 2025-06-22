/**************************************************************
 * Gyro-averaging operators
 *
 *
 * 2010-09-03 Ben Dudson <bd512@york.ac.uk>
 *    * Initial version, simple averaging operator
 *
 **************************************************************************
 * Copyright 2010 - 2025 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************/

#include <bout/difops.hxx>
#include <bout/globals.hxx>
#include <bout/gyro_average.hxx>
#include <bout/invert_laplace.hxx>
#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>

Field3D gyroTaylor0(const Field3D& f, const Field3D& rho) {
  return f + SQ(rho) * Delp2(f);
}

Field3D gyroPade0(const Field3D& f, BoutReal rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  const Field2D a{1.0, f.getMesh()};
  const Field2D d{-rho * rho, f.getMesh()};

  // Invert, leaving boundaries unchanged

  Timer timer("invert");

  auto lap = Laplacian::create(f.getMesh());

  lap->setCoefA(a);
  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(inner_boundary_flags);
  lap->setOuterBoundaryFlags(outer_boundary_flags);

  return lap->solve(f).setLocation(f.getLocation());
}

Field3D gyroPade0(const Field3D& f, const Field2D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  const Field2D a{1.0, f.getMesh()};
  const Field2D d{-rho * rho};

  // Invert, leaving boundaries unchanged
  Timer timer("invert");

  auto lap = Laplacian::create(f.getMesh());

  lap->setCoefA(a);
  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(inner_boundary_flags);
  lap->setOuterBoundaryFlags(outer_boundary_flags);

  return lap->solve(f).setLocation(f.getLocation());
}

Field3D gyroPade0(const Field3D& f, const Field3D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  // Have to use Z average of rho for efficient inversion
  return gyroPade0(f, DC(rho), inner_boundary_flags, outer_boundary_flags);
}

Field3D gyroPade1(const Field3D& f, BoutReal rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  const Field2D a{1.0, f.getMesh()};
  const Field2D d{-0.5 * rho * rho, f.getMesh()};

  // Invert, leaving boundaries unchanged
  Timer timer("invert");

  auto lap = Laplacian::create(f.getMesh());

  lap->setCoefA(a);
  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(inner_boundary_flags);
  lap->setOuterBoundaryFlags(outer_boundary_flags);

  return lap->solve(f).setLocation(f.getLocation());
}

Field3D gyroPade1(const Field3D& f, const Field2D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  const Field2D a{1.0, f.getMesh()};
  const Field2D d{-0.5 * rho * rho};

  // Invert, leaving boundaries unchanged
  Timer timer("invert");

  auto lap = Laplacian::create(f.getMesh());

  lap->setCoefA(a);
  lap->setCoefC(1.0);
  lap->setCoefD(d);
  lap->setInnerBoundaryFlags(inner_boundary_flags);
  lap->setOuterBoundaryFlags(outer_boundary_flags);

  return lap->solve(f).setLocation(f.getLocation());
}

Field3D gyroPade1(const Field3D& f, const Field3D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  return gyroPade1(f, DC(rho), inner_boundary_flags, outer_boundary_flags);
}

Field2D gyroPade1(const Field2D& f, const Field2D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  // Very inefficient implementation
  Field3D tmp = f;
  tmp = gyroPade1(tmp, rho, inner_boundary_flags, outer_boundary_flags);
  return DC(tmp);
}

Field3D gyroPade2(const Field3D& f, BoutReal rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  Field3D result =
      gyroPade1(gyroPade1(f, rho, inner_boundary_flags, outer_boundary_flags), rho,
                inner_boundary_flags, outer_boundary_flags);

  result.getMesh()->communicate(result);
  result = 0.5 * rho * rho * Delp2(result);
  result.applyBoundary("dirichlet");
  return result;
}

Field3D gyroPade2(const Field3D& f, const Field2D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  Field3D result =
      gyroPade1(gyroPade1(f, rho, inner_boundary_flags, outer_boundary_flags), rho,
                inner_boundary_flags, outer_boundary_flags);
  result.getMesh()->communicate(result);
  result = 0.5 * rho * rho * Delp2(result);
  result.applyBoundary("dirichlet");
  return result;
}

Field3D gyroPade2(const Field3D& f, const Field3D& rho, int inner_boundary_flags,
                  int outer_boundary_flags) {
  // Have to use Z average of rho for efficient inversion
  return gyroPade2(f, DC(rho), inner_boundary_flags, outer_boundary_flags);
}
