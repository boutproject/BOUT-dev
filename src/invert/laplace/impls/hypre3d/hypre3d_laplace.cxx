/**************************************************************************
 * 3D Laplacian Solver
 *                           Using Hypre Solvers
 *
 **************************************************************************
 * Copyright 2021 J.Omotani, C.MacMackin
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
#include "bout/build_config.hxx"

#if BOUT_HAS_HYPRE

#include "hypre3d_laplace.hxx"

#include <bout/mesh.hxx>
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <bout/assert.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <bout/hypre_interface.hxx>
#include <bout/operatorstencil.hxx>

LaplaceHypre3d::LaplaceHypre3d(Options *opt, const CELL_LOC loc, Mesh *mesh_in) :
  Laplacian(opt, loc, mesh_in),
  A(0.0), C1(1.0), C2(1.0), D(1.0), Ex(0.0), Ez(0.0),
  opts(opt == nullptr ? Options::getRoot()->getSection("laplace") : opt),
  lowerY(localmesh->iterateBndryLowerY()), upperY(localmesh->iterateBndryUpperY()),
  indexer(std::make_shared<GlobalIndexer<Field3D>>(localmesh,
						   getStencil(localmesh, lowerY, upperY))),
  operator3D(indexer), solution(indexer), rhs(indexer), linearSystem(*localmesh, *opts)
{
  // Provide basic initialisation of field coefficients, etc.
  // Get relevent options from user input
  A.setLocation(location);
  C1.setLocation(location);
  C2.setLocation(location);
  D.setLocation(location);
  Ex.setLocation(location);
  Ez.setLocation(location);

  // Initialise Hypre objects
  linearSystem.setMatrix(&operator3D);
  linearSystem.setRHSVector(&rhs);
  linearSystem.setSolutionVector(&solution);

  // Get y boundary flags
  lower_boundary_flags = (*opts)["lower_boundary_flags"].withDefault(0);
  upper_boundary_flags = (*opts)["upper_boundary_flags"].withDefault(0);

  // Checking flags are set to something which is not implemented
  // This is done binary (which is possible as each flag is a power of 2)
  if ( global_flags & ~implemented_flags ) {
    throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in LaplaceHypre3d");
  }
  if ( inner_boundary_flags & ~implemented_boundary_flags ) {
    throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in LaplaceHypre3d");
  }
  if ( outer_boundary_flags & ~implemented_boundary_flags ) {
    throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in LaplaceHypre3d");
  }
  if ( lower_boundary_flags & ~implemented_boundary_flags ) {
    throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in LaplaceHypre3d");
  }
  if ( upper_boundary_flags & ~implemented_boundary_flags ) {
    throw BoutException("Attempted to set Laplacian inversion boundary flag that is not implemented in LaplaceHypre3d");
  }    
  if(localmesh->periodicX) {
    throw BoutException("LaplaceHypre3d does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
  }

  // Set up boundary conditions in operator
  BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
    if(inner_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on inner X boundary
      operator3D(i, i) = -1./coords->dx[i]/sqrt(coords->g_11[i]);
      operator3D(i, i.xp()) = 1./coords->dx[i]/sqrt(coords->g_11[i]);
    } else {
      // Dirichlet on inner X boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.xp()) = 0.5;
    }
  }

  BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
    if(outer_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on outer X boundary
      operator3D(i, i) = 1./coords->dx[i]/sqrt(coords->g_11[i]);
      operator3D(i, i.xm()) = -1./coords->dx[i]/sqrt(coords->g_11[i]);
    } else {
      // Dirichlet on outer X boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.xm()) = 0.5;
    }
  }

  BOUT_FOR_SERIAL(i, indexer->getRegionLowerY()) {
    if(lower_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on lower Y boundary
      operator3D(i, i) = -1./coords->dy[i]/sqrt(coords->g_22[i]);
      operator3D(i, i.yp()) = 1./coords->dy[i]/sqrt(coords->g_22[i]);
    } else {
      // Dirichlet on lower Y boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.yp()) = 0.5;
    }
  }

  BOUT_FOR_SERIAL(i, indexer->getRegionUpperY()) {
    if(upper_boundary_flags & INVERT_AC_GRAD) {
      // Neumann on upper Y boundary
      operator3D(i, i) = 1./coords->dy[i]/sqrt(coords->g_22[i]);
      operator3D(i, i.ym()) = -1./coords->dy[i]/sqrt(coords->g_22[i]);
    } else {
      // Dirichlet on upper Y boundary
      operator3D(i, i) = 0.5;
      operator3D(i, i.ym()) = 0.5;
    }
  }
}


LaplaceHypre3d::~LaplaceHypre3d() {
}


Field3D LaplaceHypre3d::solve(const Field3D &b_in, const Field3D &x0) {
  // If necessary, update the values in the matrix operator
  if (updateRequired) {
    updateMatrix3D();
  }

  auto b = b_in;
  // Make sure b has a unique copy of the data
  b.allocate();

  // Adjust vectors to represent boundary conditions and check that
  // boundary cells are finite
  BOUT_FOR_SERIAL(i, indexer->getRegionInnerX()) {
    const BoutReal val = (inner_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    ASSERT1(finite(val));
    if (!(inner_boundary_flags & INVERT_RHS)) {
      b[i] = val;
    }
    else {
      ASSERT1(finite(b[i]));
    }
  }

  BOUT_FOR_SERIAL(i, indexer->getRegionOuterX()) {
    const BoutReal val = (outer_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    ASSERT1(finite(val));
    if (!(outer_boundary_flags & INVERT_RHS)) {
      b[i] = val;
    }
    else {
      ASSERT1(finite(b[i]));
    }
  }

  BOUT_FOR_SERIAL(i, indexer->getRegionLowerY()) {
    const BoutReal val = (lower_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    ASSERT1(finite(val));
    if (!(lower_boundary_flags & INVERT_RHS)) {
      b[i] = val;
    }
    else {
      ASSERT1(finite(b[i]));
    }
  }

  BOUT_FOR_SERIAL(i, indexer->getRegionUpperY()) {
    const BoutReal val = (upper_boundary_flags & INVERT_SET) ? x0[i] : 0.;
    ASSERT1(finite(val));
    if (!(upper_boundary_flags & INVERT_RHS)) {
      b[i] = val;
    }
    else {
      ASSERT1(finite(b[i]));
    }
  }

  rhs.importValuesFromField(b);
  solution.importValuesFromField(x0);
  rhs.assemble();
  solution.assemble();

  // Invoke solver
  { Timer timer("hypresolve");
    linearSystem.solve();
  }

  // Create field from solution
  Field3D result = solution.toField();
  localmesh->communicate(result);
  if (result.hasParallelSlices()) {
    BOUT_FOR(i, indexer->getRegionLowerY()) {
      result.ydown()[i] = result[i];
    }
    BOUT_FOR(i, indexer->getRegionUpperY()) {
      result.yup()[i] = result[i];
    }
    for (int b = 1; b < localmesh->ystart; b++) {
      BOUT_FOR(i, indexer->getRegionLowerY()) {
        result.ydown(b)[i.ym(b)] = result[i];
      }
      BOUT_FOR(i, indexer->getRegionUpperY()) {
        result.yup(b)[i.yp(b)] = result[i];
      }
    }
  }
  for (int b = 1; b < localmesh->xstart; b++) {
    BOUT_FOR(i, indexer->getRegionInnerX()) {
      result[i.xm(b)] = result[i];
    }
    BOUT_FOR(i, indexer->getRegionOuterX()) {
      result[i.xp(b)] = result[i];
    }
  }

  return result;
}

Field2D LaplaceHypre3d::solve(const Field2D &b) {
  return Laplacian::solve(b);
}

bout::HypreMatrix<Field3D>& LaplaceHypre3d::getMatrix3D() {
  if (updateRequired) {
    updateMatrix3D();
  }
  return operator3D;
}

void LaplaceHypre3d::updateMatrix3D() {
  const Field3D dc_dx = issetC ? DDX(C2) : Field3D();
  const Field3D dc_dy = issetC ? DDY(C2) : Field3D();
  const Field3D dc_dz = issetC ? DDZ(C2) : Field3D();
  const Field2D dJ_dy = DDY(coords->J/coords->g_22);
  
  // Set up the matrix for the internal points on the grid.
  // Boundary conditions were set in the constructor.
  BOUT_FOR_SERIAL(l, indexer->getRegionNobndry()) {
    // Index is called l for "location". It is not called i so as to
    // avoid confusing it with the x-index.

    // Calculate coefficients for the terms in the differential operator
    BoutReal C_df_dx = coords->G1[l], C_df_dz = coords->G3[l];
    if (issetD) {
      C_df_dx *= D[l];
      C_df_dz *= D[l];
    }
    if (issetC) {
      C_df_dx += (coords->g11[l]*dc_dx[l] + coords->g12[l]*dc_dy[l] +
		  coords->g13[l]*dc_dz[l])/C1[l];
      C_df_dz += (coords->g13[l]*dc_dx[l] + coords->g23[l]*dc_dy[l] +
		  coords->g33[l]*dc_dz[l])/C1[l];
    }
    if (issetE) {
      C_df_dx += Ex[l];
      C_df_dz += Ez[l];
    }

    BoutReal C_d2f_dx2 = coords->g11[l],
      C_d2f_dy2 = (coords->g22[l] - 1.0/coords->g_22[l]),
      C_d2f_dz2 = coords->g33[l];
    if (issetD) {
      C_d2f_dx2 *= D[l];
      C_d2f_dy2 *= D[l];
      C_d2f_dz2 *= D[l];
    }
  
    BoutReal C_d2f_dxdz = 2*coords->g13[l];
    if (issetD) {
      C_d2f_dxdz *= D[l];
    }
  
    // Adjust the coefficients to include finite-difference factors
    if (nonuniform) {
      C_df_dx += C_d2f_dx2 * coords->d1_dx[l];
    }
    C_df_dx /= 2 * coords->dx[l];
    C_df_dz /= 2 * coords->dz;

    C_d2f_dx2 /= SQ(coords->dx[l]);
    C_d2f_dy2 /= SQ(coords->dy[l]);
    C_d2f_dz2 /= SQ(coords->dz);

    C_d2f_dxdz /= 4 * coords->dx[l] * coords->dz;

    operator3D(l, l) = -2 * (C_d2f_dx2 + C_d2f_dy2 + C_d2f_dz2) + A[l];
    operator3D(l, l.xp()) = C_df_dx + C_d2f_dx2;
    operator3D(l, l.xm()) = -C_df_dx + C_d2f_dx2;
    operator3D(l, l.zp()) = C_df_dz + C_d2f_dz2;
    operator3D(l, l.zm()) = -C_df_dz + C_d2f_dz2;
    operator3D(l, l.xp().zp()) = C_d2f_dxdz;
    operator3D(l, l.xp().zm()) = -C_d2f_dxdz;
    operator3D(l, l.xm().zp()) = -C_d2f_dxdz;
    operator3D(l, l.xm().zm()) = C_d2f_dxdz;
    // The values stored in the y-boundary are already interpolated
    // up/down, so we don't want the matrix to do any such
    // interpolation there.
    const int yup = (l.y() == localmesh->yend && upperY.intersects(l.x())) ? -1 : 0,
              ydown = (l.y() == localmesh->ystart && lowerY.intersects(l.x())) ? -1 : 0;
    operator3D.yup(yup)(l, l.yp()) = 0.0;
    operator3D.ydown(ydown)(l, l.ym()) = 0.0;
    operator3D.yup(yup)(l, l.xp().yp()) = 0.0;
    operator3D.ydown(ydown)(l, l.xp().ym()) = 0.0;
    operator3D.yup(yup)(l, l.xm().yp()) = 0.0;
    operator3D.ydown(ydown)(l, l.xm().ym()) = 0.0;
    operator3D.yup(yup)(l, l.yp().zp()) = 0.0;
    operator3D.yup(yup)(l, l.yp().zm()) = 0.0;
    operator3D.ydown(ydown)(l, l.ym().zp()) = 0.0;
    operator3D.ydown(ydown)(l, l.ym().zm()) = 0.0;
  }

  // Must add these (rather than assign) so that elements used in
  // interpolation don't overwrite each other.
  BOUT_FOR_SERIAL(l, indexer->getRegionNobndry()) {
    BoutReal C_df_dy = (coords->G2[l] - dJ_dy[l]/coords->J[l]);
    if (issetD) {
      C_df_dy *= D[l];
    }
    if (issetC) {
      C_df_dy += (coords->g12[l]*dc_dx[l] + (coords->g22[l] - 1./coords->g_22[l])*dc_dy[l] +
		  coords->g23[l]*dc_dz[l])/C1[l];
    }

    BoutReal C_d2f_dy2 = (coords->g22[l] - 1.0/coords->g_22[l]);
    if (issetD) {
      C_d2f_dy2 *= D[l];
    }

    BoutReal C_d2f_dxdy = 2*coords->g12[l],
      C_d2f_dydz = 2*coords->g23[l];
    if (issetD) {
      C_d2f_dxdy *= D[l];
      C_d2f_dydz *= D[l];
    }
  
    // Adjust the coefficients to include finite-difference factors
    if (nonuniform) {
      C_df_dy += C_d2f_dy2*coords->d1_dy[l];
    }
    C_df_dy /= 2*coords->dy[l];
    C_d2f_dy2 /= SQ(coords->dy[l]);
    C_d2f_dxdy /= 4*coords->dx[l]; // NOTE: This value is not completed here. It needs to
                                   // be divide by dx(i +/- 1, j, k) when using to set a
                                   // matrix element
    C_d2f_dydz /= 4*coords->dy[l]*coords->dz;

    // The values stored in the y-boundary are already interpolated
    // up/down, so we don't want the matrix to do any such
    // interpolation there.
    const int yup = (l.y() == localmesh->yend && upperY.intersects(l.x())) ? -1 : 0,
      ydown = (l.y() == localmesh->ystart && lowerY.intersects(l.x())) ? -1 : 0;
    
    operator3D.yup(yup)(l, l.yp()) += C_df_dy + C_d2f_dy2;
    operator3D.ydown(ydown)(l, l.ym()) += -C_df_dy + C_d2f_dy2;
    operator3D.yup(yup)(l, l.xp().yp()) += C_d2f_dxdy/coords->dy[l.xp()];
    operator3D.ydown(ydown)(l, l.xp().ym()) += -C_d2f_dxdy/coords->dy[l.xp()];
    operator3D.yup(yup)(l, l.xm().yp()) += -C_d2f_dxdy/coords->dy[l.xm()];
    operator3D.ydown(ydown)(l, l.xm().ym()) += C_d2f_dxdy/coords->dy[l.xm()];
    operator3D.yup(yup)(l, l.yp().zp()) += C_d2f_dydz;
    operator3D.yup(yup)(l, l.yp().zm()) += -C_d2f_dydz;
    operator3D.ydown(ydown)(l, l.ym().zp()) += -C_d2f_dydz;
    operator3D.ydown(ydown)(l, l.ym().zm()) += C_d2f_dydz;    
  }
  operator3D.assemble();
  linearSystem.setupAMG(&operator3D);

  updateRequired = false;
}

OperatorStencil<Ind3D> LaplaceHypre3d::getStencil(Mesh* localmesh,
						  const RangeIterator &lowerYBound,
						  const RangeIterator &upperYBound) {
  OperatorStencil<Ind3D> stencil;

  // Get the pattern used for interpolation. This is assumed to be the
  // same across the whole grid.
  const auto pw = localmesh->getCoordinates()->getParallelTransform().getWeightsForYDownApproximation(localmesh->xstart, localmesh->ystart + 1, localmesh->zstart);
  std::vector<OffsetInd3D> interpPattern;
  std::transform(pw.begin(), pw.end(), std::back_inserter(interpPattern),
          [localmesh](ParallelTransform::PositionsAndWeights p) -> OffsetInd3D {
            return {localmesh->xstart - p.i, localmesh->ystart - p.j,
		    localmesh->LocalNz - p.k < p.k ? p.k - localmesh->LocalNz : p.k};
          });

  OffsetInd3D zero;

  // Add interior cells
  const std::vector<OffsetInd3D> interpolatedUpElements = {zero.yp(), zero.xp().yp(), zero.xm().yp(),
						     zero.yp().zp(), zero.yp().zm()},
    interpolatedDownElements = {zero.ym(), zero.xp().ym(), zero.xm().ym(), zero.ym().zp(),
				zero.ym().zm()};
  std::set<OffsetInd3D> interiorStencil = {zero, zero.xp(), zero.xm(),
					   zero.zp(), zero.zm(),
					   zero.xp().zp(), zero.xp().zm(),
					   zero.xm().zp(), zero.xm().zm()},
    lowerEdgeStencil = interiorStencil, upperEdgeStencil = interiorStencil;

  for (const auto& i : interpolatedDownElements) {
    for (auto& j : interpPattern) {
      interiorStencil.insert(i + j);
      upperEdgeStencil.insert(i + j);
    }
    lowerEdgeStencil.insert(i);
  }
  for (const auto& i : interpolatedUpElements) {
    for (auto& j : interpPattern) {
      interiorStencil.insert(i + j);
      lowerEdgeStencil.insert(i + j);
    }
    upperEdgeStencil.insert(i);
  }
  const std::vector<OffsetInd3D> interiorStencilVector(interiorStencil.begin(), interiorStencil.end()),
    lowerEdgeStencilVector(lowerEdgeStencil.begin(), lowerEdgeStencil.end()),
    upperEdgeStencilVector(upperEdgeStencil.begin(), upperEdgeStencil.end());

  // If there is a lower y-boundary then create a part of the stencil
  // for cells immediately adjacent to it.
  if (lowerYBound.max() - lowerYBound.min() > 0) {
    stencil.add([index = localmesh->ystart, lowerYBound](Ind3D ind) -> bool {
		  return index == ind.y() && lowerYBound.intersects(ind.x()); },
      lowerEdgeStencilVector);
  }

  // If there is an upper y-boundary then create a part of the stencil
  // for cells immediately adjacent to it.
  if (upperYBound.max() - upperYBound.min() > 0) {
    stencil.add([index = localmesh->yend, upperYBound](Ind3D ind) -> bool {
		  return index == ind.y() && upperYBound.intersects(ind.x()); },
      upperEdgeStencilVector);
  }

  // Create a part of the stencil for the interior cells. Although the
  // test here would also pass for the edge-cells immediately adjacent
  // to upper/lower y-boundaries, because those tests are run first
  // the cells will be assigned to those regions.
  stencil.add([localmesh](Ind3D ind) -> bool {
		return (localmesh->xstart <= ind.x() &&
			ind.x() <= localmesh->xend &&
			localmesh->ystart <= ind.y() &&
			ind.y() <= localmesh->yend &&
			localmesh->zstart <= ind.z() &&
			ind.z() <= localmesh->zend); },
    interiorStencilVector);

  // Add Y boundaries before X boundaries so corners are assigned to
  // the former.

  // Note that the tests for whether a cell is in the boundary
  // actually are not exclusive enough. They really only correspond to
  // whether the cell is in the guard regions. However, the tests are
  // much simpler to implement this way and other checks will prevent
  // guard-cells which are not boundaries from having memory
  // pre-allocated.
  
  // Add lower Y boundary.
  stencil.add([index = localmesh->ystart - 1](Ind3D ind) -> bool {
		return ind.y() == index;
	      }, {zero, zero.yp()});
  // Add upper Y boundary
  stencil.add([index = localmesh->yend + 1](Ind3D ind) -> bool {
		return ind.y() == index;
	      }, {zero, zero.ym()});
  // Add inner X boundary
  if (localmesh->firstX()) {
    stencil.add(
        [index = localmesh->xstart - 1](Ind3D ind) -> bool { return ind.x() == index; },
        {zero, zero.xp()});
  }
  // Add outer X boundary
  if (localmesh->lastX()) {
    stencil.add([index = localmesh->xend + 1](Ind3D ind) -> bool { return ind.x() == index; },
                {zero, zero.xm()});
  }

  return stencil;
}

#endif // BOUT_HAS_HYPRE
