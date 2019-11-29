/**************************************************************************
 * Implements the shifted metric method for parallel derivatives
 *
 * By default fields are stored so that X-Z are orthogonal, and so not aligned
 * in Y. This implementation uses Interpolation objects for interpolation
 * rather than FFTs.
 *
 **************************************************************************
 * Copyright 2018 B.D.Dudson, P. Hill, J. Omotani, J. Parker
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
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
 **************************************************************************/

#include "shiftedmetricinterp.hxx"
#include "mask.hxx"
#include "bout/constants.hxx"
#include <interpolation_factory.hxx>

ShiftedMetricInterp::ShiftedMetricInterp(Mesh& mesh, CELL_LOC location_in,
                                         Field2D zShift_in)
    : ParallelTransform(mesh), location(location_in), zShift(std::move(zShift_in)) {
  // check the coordinate system used for the grid data source
  ShiftedMetricInterp::checkInputGrid();

  // Create the Interpolation objects and set whether they go up or down the
  // magnetic field
  interp_yup = InterpolationFactory::getInstance().create(&mesh);
  interp_yup->setYOffset(1);

  interp_ydown = InterpolationFactory::getInstance().create(&mesh);
  interp_ydown->setYOffset(-1);

  // Find the index positions where the magnetic field line intersects the next
  // x-z plane
  Field3D xt_prime(&mesh), zt_prime_up(&mesh), zt_prime_down(&mesh);
  xt_prime.allocate();
  zt_prime_up.allocate();
  zt_prime_down.allocate();

  for (const auto& i : xt_prime) {
    // no interpolation in x, all field lines stay at constant x
    xt_prime[i] = i.x();
  }

  for (const auto& i : zt_prime_up.getRegion(RGN_NOY)) {
    // Field line moves in z by an angle zShift(i,j+1)-zShift(i,j) when going
    // from j to j+1, but we want the shift in index-space
    zt_prime_up[i] =
        static_cast<BoutReal>(i.z())
        + (zShift[i.yp()] - zShift[i]) * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
  }

  // Make a mask, to skip interpolating to points inside y-boundaries that need to be set
  // by boundary conditions.
  auto mask_up = BoutMask(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);
  for (auto it = mesh.iterateBndryUpperY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      mask_up(it.ind, mesh.yend, z) = true;
    }
  }

  interp_yup->calcWeights(xt_prime, zt_prime_up, mask_up);

  for (const auto& i : xt_prime) {
    // no interpolation in x, all field lines stay at constant x
    xt_prime[i] = i.x();
  }

  for (const auto& i : zt_prime_down.getRegion(RGN_NOY)) {
    // Field line moves in z by an angle -(zShift(i,j)-zShift(i,j-1)) when going
    // from j to j-1, but we want the shift in index-space
    zt_prime_down[i] =
        static_cast<BoutReal>(i.z())
        - (zShift[i] - zShift[i.ym()]) * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
  }

  auto mask_down = BoutMask(mesh.LocalNx, mesh.LocalNy, mesh.LocalNz);
  for (auto it = mesh.iterateBndryLowerY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      mask_down(it.ind, mesh.ystart, z) = true;
    }
  }

  interp_ydown->calcWeights(xt_prime, zt_prime_down, mask_down);

  // Set up interpolation to/from field-aligned coordinates
  interp_to_aligned = InterpolationFactory::getInstance().create(&mesh);
  interp_from_aligned = InterpolationFactory::getInstance().create(&mesh);

  Field3D zt_prime_to(&mesh), zt_prime_from(&mesh);
  zt_prime_to.allocate();
  zt_prime_from.allocate();

  for (const auto& i : zt_prime_to) {
    // Field line moves in z by an angle zShift(i,j) when going
    // from y0 to y(j), but we want the shift in index-space
    zt_prime_to[i] = static_cast<BoutReal>(i.z())
                     + zShift[i] * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
  }

  interp_to_aligned->calcWeights(xt_prime, zt_prime_to);

  for (const auto& i : zt_prime_from) {
    // Field line moves in z by an angle zShift(i,j) when going
    // from y0 to y(j), but we want the shift in index-space.
    // Here we reverse the shift, so subtract zShift
    zt_prime_from[i] = static_cast<BoutReal>(i.z())
                       - zShift[i] * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
  }

  interp_from_aligned->calcWeights(xt_prime, zt_prime_from);

  // Create regions for parallel boundary conditions
  Field2D dy;
  mesh.get(dy, "dy", 1.);
  auto forward_boundary =
      new BoundaryRegionPar("parallel_forward", BNDRY_PAR_FWD, +1, &mesh);
  for (auto it = mesh.iterateBndryUpperY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      forward_boundary->add_point(
          it.ind, mesh.yend, z,
          mesh.GlobalX(it.ind),                           // x
          2. * PI * mesh.GlobalY(mesh.yend + 0.5),        // y
          2. * PI * BoutReal(z) / BoutReal(mesh.GlobalNz) // z
              + 0.5 * (zShift(it.ind, mesh.yend + 1) - zShift(it.ind, mesh.yend)),
          0.25
              * (dy(it.ind, mesh.yend) // dy/2
                 + dy(it.ind, mesh.yend + 1)),
          0. // angle?
      );
    }
  }
  auto backward_boundary =
      new BoundaryRegionPar("parallel_backward", BNDRY_PAR_BKWD, -1, &mesh);
  for (auto it = mesh.iterateBndryLowerY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      backward_boundary->add_point(
          it.ind, mesh.ystart, z,
          mesh.GlobalX(it.ind),                           // x
          2. * PI * mesh.GlobalY(mesh.ystart - 0.5),      // y
          2. * PI * BoutReal(z) / BoutReal(mesh.GlobalNz) // z
              + 0.5 * (zShift(it.ind, mesh.ystart) - zShift(it.ind, mesh.ystart - 1)),
          0.25
              * (dy(it.ind, mesh.ystart - 1) // dy/2
                 + dy(it.ind, mesh.ystart)),
          0. // angle?
      );
    }
  }

  // Add the boundary region to the mesh's vector of parallel boundaries
  mesh.addBoundaryPar(forward_boundary);
  mesh.addBoundaryPar(backward_boundary);
}

void ShiftedMetricInterp::checkInputGrid() {
  std::string coordinates_type = "";
  if (!mesh.get(coordinates_type, "coordinates_type")) {
    if (coordinates_type != "orthogonal") {
      throw BoutException("Incorrect coordinate system type '" + coordinates_type
                          + "' used to generate metric components for ShiftedMetric. "
                            "Should be 'orthogonal'.");
    }
  } // else: coordinate_system variable not found in grid input, indicates older input
    //       file so must rely on the user having ensured the type is correct
}

/*!
 * Calculate the Y up and down fields
 */
void ShiftedMetricInterp::calcParallelSlices(Field3D& f) {
  AUTO_TRACE();

  // Ensure that yup and ydown are different fields
  f.splitParallelSlices();

  // Interpolate f onto yup and ydown fields
  f.yup() = interp_yup->interpolate(f);
  f.ydown() = interp_ydown->interpolate(f);
}

/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetricInterp::toFieldAligned(const Field3D& f,
                                                  const std::string& region) {
  return interp_to_aligned->interpolate(f, region);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetricInterp::fromFieldAligned(const Field3D& f,
                                                    const std::string& region) {
  return interp_from_aligned->interpolate(f, region);
}
