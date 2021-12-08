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

ShiftedMetricInterp::ShiftedMetricInterp(Mesh& mesh, CELL_LOC location_in,
                                         Field2D zShift_in, Options* opt)
    : ParallelTransform(mesh, opt), location(location_in), zShift(std::move(zShift_in)),
      ydown_index(mesh.ystart) {
  // check the coordinate system used for the grid data source
  ShiftedMetricInterp::checkInputGrid();

  // Allocate space for interpolator cache: y-guard cells in each direction
  parallel_slice_interpolators.resize(mesh.ystart * 2);

  // Create the Interpolation objects and set whether they go up or down the
  // magnetic field
  auto& interp_options = options["zinterpolation"];
  // Careful with the indices/offsets! Offsets are 1-indexed (as 0 would be the original
  // slice), and Mesh::ystart is the number of guard cells. The parallel slice vector
  // stores the offsets as
  //    {+1, ..., +n, -1, ..., -n}
  // Once parallel_slice_interpolators is initialised though, each interpolator stores its
  // offset, so we don't need to faff about after this
  for (int y_offset = 0; y_offset < mesh.ystart; ++y_offset) {
    parallel_slice_interpolators[yup_index + y_offset] =
      ZInterpolationFactory::getInstance().create(&interp_options, y_offset + 1, &mesh);
    parallel_slice_interpolators[ydown_index + y_offset] =
      ZInterpolationFactory::getInstance().create(&interp_options, -y_offset - 1, &mesh);

    // Find the index positions where the magnetic field line intersects the x-z plane
    // y_offset points up
    Field3D zt_prime_up(&mesh), zt_prime_down(&mesh);
    zt_prime_up.allocate();
    zt_prime_down.allocate();

    for (const auto& i : zt_prime_up.getRegion(RGN_NOY)) {
      // Field line moves in z by an angle zShift(i,j+1)-zShift(i,j) when going
      // from j to j+1, but we want the shift in index-space
      zt_prime_up[i] =
          static_cast<BoutReal>(i.z())
          + (zShift[i.yp(y_offset + 1)] - zShift[i])
            * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
    }

    parallel_slice_interpolators[yup_index + y_offset]->calcWeights(zt_prime_up);

    for (const auto& i : zt_prime_down.getRegion(RGN_NOY)) {
      // Field line moves in z by an angle -(zShift(i,j)-zShift(i,j-1)) when going
      // from j to j-1, but we want the shift in index-space
      zt_prime_down[i] =
          static_cast<BoutReal>(i.z())
          - (zShift[i] - zShift[i.ym(y_offset + 1)])
            * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
    }

    parallel_slice_interpolators[ydown_index + y_offset]->calcWeights(zt_prime_down);
  }

  // Set up interpolation to/from field-aligned coordinates
  interp_to_aligned = ZInterpolationFactory::getInstance().create(&interp_options, 0, &mesh);
  interp_from_aligned =
      ZInterpolationFactory::getInstance().create(&interp_options, 0, &mesh);

  Field3D zt_prime_to(&mesh), zt_prime_from(&mesh);
  zt_prime_to.allocate();
  zt_prime_from.allocate();

  for (const auto& i : zt_prime_to) {
    // Field line moves in z by an angle zShift(i,j) when going
    // from y0 to y(j), but we want the shift in index-space
    zt_prime_to[i] = static_cast<BoutReal>(i.z())
                     + zShift[i] * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
  }

  interp_to_aligned->calcWeights(zt_prime_to);

  for (const auto& i : zt_prime_from) {
    // Field line moves in z by an angle zShift(i,j) when going
    // from y0 to y(j), but we want the shift in index-space.
    // Here we reverse the shift, so subtract zShift
    zt_prime_from[i] = static_cast<BoutReal>(i.z())
                       - zShift[i] * static_cast<BoutReal>(mesh.GlobalNz) / TWOPI;
  }

  interp_from_aligned->calcWeights(zt_prime_from);

  // Create regions for parallel boundary conditions
  Field2D dy;
  mesh.get(dy, "dy", 1.);
  auto forward_boundary_xin =
      new BoundaryRegionPar("parallel_forward_xin", BNDRY_PAR_FWD_XIN, +1, &mesh);
  for (auto it = mesh.iterateBndryUpperY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      forward_boundary_xin->add_point(
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
  auto backward_boundary_xin =
      new BoundaryRegionPar("parallel_backward_xin", BNDRY_PAR_BKWD_XIN, -1, &mesh);
  for (auto it = mesh.iterateBndryLowerY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      backward_boundary_xin->add_point(
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
  // Create regions for parallel boundary conditions
  auto forward_boundary_xout =
      new BoundaryRegionPar("parallel_forward_xout", BNDRY_PAR_FWD_XOUT, +1, &mesh);
  for (auto it = mesh.iterateBndryUpperY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      forward_boundary_xout->add_point(
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
  auto backward_boundary_xout =
      new BoundaryRegionPar("parallel_backward_xout", BNDRY_PAR_BKWD_XOUT, -1, &mesh);
  for (auto it = mesh.iterateBndryLowerY(); not it.isDone(); it.next()) {
    for (int z = mesh.zstart; z <= mesh.zend; z++) {
      backward_boundary_xout->add_point(
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
  mesh.addBoundaryPar(forward_boundary_xin);
  mesh.addBoundaryPar(backward_boundary_xin);
  mesh.addBoundaryPar(forward_boundary_xout);
  mesh.addBoundaryPar(backward_boundary_xout);
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
  for (const auto& interp : parallel_slice_interpolators) {
    f.ynext(interp->y_offset) = interp->interpolate(f);
  }
}

/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
Field3D ShiftedMetricInterp::toFieldAligned(const Field3D& f, const std::string& region) {
  ASSERT2(f.getDirectionY() == YDirectionType::Standard);
  return interp_to_aligned->interpolate(f, region).setDirectionY(YDirectionType::Aligned);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
Field3D ShiftedMetricInterp::fromFieldAligned(const Field3D& f,
                                              const std::string& region) {
  ASSERT2(f.getDirectionY() == YDirectionType::Aligned);
  return interp_from_aligned->interpolate(f, region).setDirectionY(
      YDirectionType::Standard);
}
