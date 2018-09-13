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
#include "bout/constants.hxx"
#include <interpolation_factory.hxx>

ShiftedMetricInterp::ShiftedMetricInterp(Mesh& mesh) : localmesh(mesh) {

  // Read the zShift angle from the mesh
  if (localmesh.get(zShift, "zShift")) {
    // No zShift variable. Try qinty in BOUT grid files
    localmesh.get(zShift, "qinty");
  }

  // Create the Interpolation objects and set whether they go up or down the
  // magnetic field
  interp_yup = InterpolationFactory::getInstance()->create(&localmesh);
  interp_yup->setYOffset(1);

  interp_ydown = InterpolationFactory::getInstance()->create(&localmesh);
  interp_ydown->setYOffset(-1);

  // Find the index positions where the magnetic field line intersects the next
  // x-z plane
  Field3D xt_prime(&localmesh), zt_prime_up(&localmesh), zt_prime_down(&localmesh);
  xt_prime.allocate();
  zt_prime_up.allocate();
  zt_prime_down.allocate();

  for (const auto& i : xt_prime) {
    // no interpolation in x, all field lines stay at constant x
    xt_prime[i] = i.x;
  }

  for (const auto& i : zt_prime_up.region(RGN_NOY)) {
    // Field line moves in z by an angle zShift(i,j+1)-zShift(i,j) when going
    // from j to j+1, but we want the shift in index-space
    zt_prime_up[i] = static_cast<BoutReal>(i.z)
                     + (zShift[i.yp()] - zShift[i])
                           * static_cast<BoutReal>(localmesh.GlobalNz) / TWOPI;
  }

  interp_yup->calcWeights(xt_prime, zt_prime_up);

  for (const auto& i : xt_prime) {
    // no interpolation in x, all field lines stay at constant x
    xt_prime[i] = i.x;
  }

  for (const auto& i : zt_prime_down.region(RGN_NOY)) {
    // Field line moves in z by an angle -(zShift(i,j)-zShift(i,j-1)) when going
    // from j to j-1, but we want the shift in index-space
    zt_prime_down[i] = static_cast<BoutReal>(i.z)
                       - (zShift[i] - zShift[i.ym()])
                             * static_cast<BoutReal>(localmesh.GlobalNz) / TWOPI;
  }

  interp_ydown->calcWeights(xt_prime, zt_prime_down);

  // Set up interpolation to/from field-aligned coordinates
  interp_to_aligned = InterpolationFactory::getInstance()->create(&localmesh);
  interp_from_aligned = InterpolationFactory::getInstance()->create(&localmesh);

  Field3D zt_prime_to(&localmesh), zt_prime_from(&localmesh);
  zt_prime_to.allocate();
  zt_prime_from.allocate();

  for (const auto& i : zt_prime_to) {
    // Field line moves in z by an angle zShift(i,j) when going
    // from y0 to y(j), but we want the shift in index-space
    zt_prime_to[i] = static_cast<BoutReal>(i.z)
                     + zShift[i] * static_cast<BoutReal>(localmesh.GlobalNz) / TWOPI;
  }

  interp_to_aligned->calcWeights(xt_prime, zt_prime_to);

  for (const auto& i : zt_prime_from) {
    // Field line moves in z by an angle zShift(i,j) when going
    // from y0 to y(j), but we want the shift in index-space.
    // Here we reverse the shift, so subtract zShift
    zt_prime_from[i] = static_cast<BoutReal>(i.z)
                       - zShift[i] * static_cast<BoutReal>(localmesh.GlobalNz) / TWOPI;
  }

  interp_from_aligned->calcWeights(xt_prime, zt_prime_from);
}

/*!
 * Calculate the Y up and down fields
 */
void ShiftedMetricInterp::calcYUpDown(Field3D& f) {
  TRACE("ShiftedMetricInterp::calcYUpDown");

  // Ensure that yup and ydown are different fields
  f.splitYupYdown();

  // Interpolate f onto yup and ydown fields
  f.yup() = interp_yup->interpolate(f);
  f.ydown() = interp_ydown->interpolate(f);
}

/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
const Field3D ShiftedMetricInterp::toFieldAligned(const Field3D& f) {
  return interp_to_aligned->interpolate(f);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
const Field3D ShiftedMetricInterp::fromFieldAligned(const Field3D& f) {
  return interp_from_aligned->interpolate(f);
}
