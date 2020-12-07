/**************************************************************************
 * Perpendicular Laplacian inversion. Serial code using FFT
 * and band-diagonal solver
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#include <globals.hxx>
#include "serial_band.hxx"

#include <bout/mesh.hxx>
#include <fft.hxx>
#include <utils.hxx>
#include <boutexception.hxx>
#include <lapack_routines.hxx>
#include <bout/constants.hxx>
#include <bout/openmpwrap.hxx>

#include <output.hxx>

//#define SECONDORDER // Define to use 2nd order differencing

LaplaceSerialBand::LaplaceSerialBand(Options *opt, const CELL_LOC loc, Mesh *mesh_in)
    : Laplacian(opt, loc, mesh_in), Acoef(0.0), Ccoef(1.0), Dcoef(1.0) {
#if BOUT_USE_METRIC_3D
  throw BoutException("Laplacian serial band solver does not support 3D metric yet.");
#endif
  Acoef.setLocation(location);
  Ccoef.setLocation(location);
  Dcoef.setLocation(location);

  if(!localmesh->firstX() || !localmesh->lastX())
    throw BoutException("LaplaceSerialBand only works for localmesh->NXPE = 1");
  if(localmesh->periodicX) {
      throw BoutException("LaplaceSerialBand does not work with periodicity in the x direction (localmesh->PeriodicX == true). Change boundary conditions or use serial-tri or cyclic solver instead");
    }
  // Allocate memory

  int ncz = localmesh->LocalNz;
  bk.reallocate(localmesh->LocalNx, ncz / 2 + 1);
  bk1d.reallocate(localmesh->LocalNx);

  //Initialise bk to 0 as we only visit 0<= kz <= maxmode in solve
  for(int kz=maxmode+1; kz < ncz/2 + 1; kz++){
    for (int ix=0; ix<localmesh->LocalNx; ix++){
      bk(ix, kz) = 0.0;
    }
  }

  xk.reallocate(localmesh->LocalNx, ncz / 2 + 1);
  xk1d.reallocate(localmesh->LocalNx);

  // Initialise xk to 0 as we only visit 0<= kz <= maxmode in solve
  for (int kz = maxmode + 1; kz < ncz / 2 + 1; kz++) {
    for (int ix = 0; ix < localmesh->LocalNx; ix++) {
      xk(ix, kz) = 0.0;
    }
  }

  A.reallocate(localmesh->LocalNx, 5);
}

FieldPerp LaplaceSerialBand::solve(const FieldPerp& b) { return solve(b, b); }

FieldPerp LaplaceSerialBand::solve(const FieldPerp& b, const FieldPerp& x0) {
#if not(BOUT_USE_METRIC_3D)
  ASSERT1(localmesh == b.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  FieldPerp x{emptyFrom(b)};

  int jy = b.getIndex();

  int ncz = localmesh->LocalNz;
  int ncx = localmesh->LocalNx-1;

  int xbndry = localmesh->xstart; // Width of the x boundary
  // If the flags to assign that only one guard cell should be used is set
  if((global_flags & INVERT_BOTH_BNDRY_ONE) || (localmesh->xstart < 2))
    xbndry = 1;

  BOUT_OMP(parallel for)
  for(int ix=0;ix<localmesh->LocalNx;ix++) {
    // for fixed ix,jy set a complex vector rho(z)
    
    if(((ix < xbndry) && (inner_boundary_flags & INVERT_SET)) ||
       ((ncx-ix < xbndry) && (outer_boundary_flags & INVERT_SET))) {
      // Use the values in x0 in the boundary
      rfft(x0[ix], ncz, &bk(ix, 0));
    }else
      rfft(b[ix], ncz, &bk(ix, 0));
  }
  
  int xstart, xend;
  // Get range for 4th order: Need at least 2 each side
  if(xbndry > 1) {
    xstart = xbndry;
    xend = ncx-xbndry;
  }else {
    xstart = 2;
    xend = localmesh->LocalNx-2;
  }

  const auto kwave_fac = TWOPI / coords->zlength(); // wave number is 1/[rad]
  for(int iz=0;iz<=maxmode;iz++) {
    // solve differential equation in x

    BoutReal coef1 = 0.0, coef2 = 0.0, coef3 = 0.0, coef4 = 0.0, coef5 = 0.0, coef6 = 0.0;
    ///////// PERFORM INVERSION /////////
      
    // shift freqs according to FFT convention
    const Field2D kwave_ = iz * kwave_fac;

    // set bk1d
    for(int ix=0;ix<localmesh->LocalNx;ix++)
      bk1d[ix] = bk(ix, iz);

    // Fill in interior points

    for(int ix=xstart;ix<=xend;ix++) {
      BoutReal kwave = kwave_(ix, jy);
#ifdef SECONDORDER 
      // Use second-order differencing. Useful for testing the tridiagonal solver
      // with different boundary conditions
      dcomplex a,b,c;
      tridagCoefs(ix, jy, iz, a, b, c, &Ccoef, &Dcoef);

      A(ix, 0) = 0.;
      A(ix, 1) = a;
      A(ix, 2) = b + Acoef(ix, jy);
      A(ix, 3) = c;
      A(ix, 4) = 0.;
#else
      // Set coefficients
      coef1 = coords->g11(ix,jy);  // X 2nd derivative
      coef2 = coords->g33(ix,jy);  // Z 2nd derivative
      coef3 = coords->g13(ix,jy);  // X-Z mixed derivatives
      coef4 = 0.0;          // X 1st derivative
      coef5 = 0.0;          // Z 1st derivative
      coef6 = Acoef(ix,jy); // Constant

      // Multiply Delp2 component by a factor
      coef1 *= Dcoef(ix,jy);
      coef2 *= Dcoef(ix,jy);
      coef3 *= Dcoef(ix,jy);

      if(all_terms) {
        coef4 = coords->G1(ix,jy);
        coef5 = coords->G3(ix,jy);
      }

      if(nonuniform) {
        // non-uniform localmesh correction
        if((ix != 0) && (ix != ncx))
          coef4 += coords->g11(ix,jy)*( (1.0/coords->dx(ix+1,jy)) - (1.0/coords->dx(ix-1,jy)) )/(2.0*coords->dx(ix,jy));
      }

      // A first order derivative term (1/c)\nabla_perp c\cdot\nabla_\perp x
    
      if((ix > 1) && (ix < (localmesh->LocalNx-2)))
        coef4 += coords->g11(ix,jy) * (Ccoef(ix-2,jy) - 8.*Ccoef(ix-1,jy) + 8.*Ccoef(ix+1,jy) - Ccoef(ix+2,jy)) / (12.*coords->dx(ix,jy)*(Ccoef(ix,jy)));

      // Put into matrix
      coef1 /= 12.* SQ(coords->dx(ix,jy));
      coef2 *= SQ(kwave);
      coef3 *= kwave / (12. * coords->dx(ix,jy));
      coef4 /= 12. * coords->dx(ix,jy);
      coef5 *= kwave;

      A(ix, 0) = dcomplex(-coef1 + coef4, coef3);
      A(ix, 1) = dcomplex(16. * coef1 - 8 * coef4, -8. * coef3);
      A(ix, 2) = dcomplex(-30. * coef1 - coef2 + coef6, coef5);
      A(ix, 3) = dcomplex(16. * coef1 + 8 * coef4, 8. * coef3);
      A(ix, 4) = dcomplex(-coef1 - coef4, -coef3);
#endif
    }

    if(xbndry < 2) {
      // Use 2nd order near edges

      int ix = 1;

      auto kwave = kwave_(ix, jy);
      coef1=coords->g11(ix,jy)/(SQ(coords->dx(ix,jy)));
      coef2=coords->g33(ix,jy);
      coef3= kwave * coords->g13(ix,jy)/(2. * coords->dx(ix,jy));
        
      // Multiply Delp2 component by a factor
      coef1 *= Dcoef(ix,jy);
      coef2 *= Dcoef(ix,jy);
      coef3 *= Dcoef(ix,jy);

      A(ix, 0) = 0.0; // Should never be used
      A(ix, 1) = dcomplex(coef1, -coef3);
      A(ix, 2) = dcomplex(-2.0 * coef1 - SQ(kwave) * coef2 + coef4, 0.0);
      A(ix, 3) = dcomplex(coef1, coef3);
      A(ix, 4) = 0.0;

      ix = ncx-1;

      coef1=coords->g11(ix,jy)/(SQ(coords->dx(ix,jy)));
      coef2=coords->g33(ix,jy);
      coef3= kwave * coords->g13(ix,jy)/(2. * coords->dx(ix,jy));

      A(ix, 0) = 0.0;
      A(ix, 1) = dcomplex(coef1, -coef3);
      A(ix, 2) = dcomplex(-2.0 * coef1 - SQ(kwave) * coef2 + coef4, 0.0);
      A(ix, 3) = dcomplex(coef1, coef3);
      A(ix, 4) = 0.0; // Should never be used
    }

    // Boundary conditions

    for(int ix=0;ix<xbndry;ix++) {
      // Set zero-value. Change to zero-gradient if needed

      if(!(inner_boundary_flags & (INVERT_RHS|INVERT_SET)))
        bk1d[ix] = 0.0;
      if(!(outer_boundary_flags & (INVERT_RHS|INVERT_SET)))
        bk1d[ncx-ix] = 0.0;

      A(ix, 0) = A(ix, 1) = A(ix, 3) = A(ix, 4) = 0.0;
      A(ix, 2) = 1.0;

      A(ncx - ix, 0) = A(ncx - ix, 1) = A(ncx - ix, 3) = A(ncx - ix, 4) = 0.0;
      A(ncx - ix, 2) = 1.0;
    }

    if(iz == 0) {
      // DC
	
      // Inner boundary
      if(inner_boundary_flags & (INVERT_DC_GRAD+INVERT_SET) || inner_boundary_flags & (INVERT_DC_GRAD+INVERT_RHS)) {
        // Zero gradient at inner boundary. 2nd-order accurate
        // Boundary at midpoint
        for (int ix=0;ix<xbndry;ix++) {
          A(ix, 0) = 0.;
          A(ix, 1) = 0.;
          A(ix, 2) = -.5 / sqrt(coords->g_11(ix, jy)) / coords->dx(ix, jy);
          A(ix, 3) = .5 / sqrt(coords->g_11(ix, jy)) / coords->dx(ix, jy);
          A(ix, 4) = 0.;
        }
        
      }
      else if(inner_boundary_flags & INVERT_DC_GRAD) {
        // Zero gradient at inner boundary. 2nd-order accurate
        // Boundary at midpoint
        for (int ix=0;ix<xbndry;ix++) {
          A(ix, 0) = 0.;
          A(ix, 1) = 0.;
          A(ix, 2) = -.5;
          A(ix, 3) = .5;
          A(ix, 4) = 0.;
        }
        
      }
      else if(inner_boundary_flags & INVERT_DC_GRADPAR) {
        for (int ix=0;ix<xbndry;ix++) {
          A(ix, 0) = 0.;
          A(ix, 1) = 0.;
          A(ix, 2) = -3. / sqrt(coords->g_22(ix, jy));
          A(ix, 3) = 4. / sqrt(coords->g_22(ix + 1, jy));
          A(ix, 4) = -1. / sqrt(coords->g_22(ix + 2, jy));
        }
      }
      else if(inner_boundary_flags & INVERT_DC_GRADPARINV) {
        for (int ix=0;ix<xbndry;ix++) {
          A(ix, 0) = 0.;
          A(ix, 1) = 0.;
          A(ix, 2) = -3. * sqrt(coords->g_22(ix, jy));
          A(ix, 3) = 4. * sqrt(coords->g_22(ix + 1, jy));
          A(ix, 4) = -sqrt(coords->g_22(ix + 2, jy));
        }
      }
      else if (inner_boundary_flags & INVERT_DC_LAP) {
        for (int ix=0;ix<xbndry;ix++) {
          A(ix, 0) = 0.;
          A(ix, 1) = 0.;
          A(ix, 2) = 1.;
          A(ix, 3) = -2;
          A(ix, 4) = 1.;
        }
      }
      
      // Outer boundary
      if(outer_boundary_flags & INVERT_DC_GRAD) {
        // Zero gradient at outer boundary
        for (int ix=0;ix<xbndry;ix++)
          A(ncx - ix, 1) = -1.0;
      }
	
    }else {
      // AC
	
      // Inner boundarySQ(kwave)*coef2
      if(inner_boundary_flags & INVERT_AC_GRAD) {
        // Zero gradient at inner boundary
        for (int ix=0;ix<xbndry;ix++)
          A(ix, 3) = -1.0;
      }else if(inner_boundary_flags & INVERT_AC_LAP) {
        // Enforce zero laplacian for 2nd and 4th-order
	  
        int ix = 1;

        auto kwave = kwave_(ix, jy);
        coef1=coords->g11(ix,jy)/(12.* SQ(coords->dx(ix,jy)));
	
        coef2=coords->g33(ix,jy);
	
        coef3= kwave * coords->g13(ix,jy)/(2. * coords->dx(ix,jy));
        
        coef4 = Acoef(ix,jy);
	  
        // Combine 4th order at 1 with 2nd order at 0
        A(1, 0) = 0.0; // Not used
        A(1, 1) = dcomplex(
            (14. - SQ(coords->dx(0, jy) * kwave) * coords->g33(0, jy) / coords->g11(0, jy)) *
                coef1,
            -coef3);
        A(1, 2) = dcomplex(-29. * coef1 - SQ(kwave) * coef2 + coef4, 0.0);
        A(1, 3) = dcomplex(16. * coef1, coef3);
        A(1, 4) = dcomplex(-coef1, 0.0);

        coef1=coords->g11(ix,jy)/(SQ(coords->dx(ix,jy)));
        coef2=coords->g33(ix,jy);
        coef3= kwave * coords->g13(ix,jy)/(2. * coords->dx(ix,jy));

        // Use 2nd order at 1
        A(0, 0) = 0.0; // Should never be used
        A(0, 1) = 0.0;
        A(0, 2) = dcomplex(coef1, -coef3);
        A(0, 3) = dcomplex(-2.0 * coef1 - SQ(kwave) * coef2 + coef4, 0.0);
        A(0, 4) = dcomplex(coef1, coef3);
      }
	
      // Outer boundary
      if(outer_boundary_flags & INVERT_AC_GRAD) {
        // Zero gradient at outer boundary
        for (int ix=0;ix<xbndry;ix++)
          A(ncx - ix, 1) = -1.0;
      }else if(outer_boundary_flags & INVERT_AC_LAP) {
        // Enforce zero laplacian for 2nd and 4th-order
        // NOTE: Currently ignoring XZ term and coef4 assumed zero on boundary
        // FIX THIS IF IT WORKS

        int ix = ncx-1;
	  
        coef1=coords->g11(ix,jy)/(12.* SQ(coords->dx(ix,jy)));
	
        coef2=coords->g33(ix,jy);

        auto kwave = kwave_(ix, jy);
        coef3= kwave * coords->g13(ix,jy)/(2. * coords->dx(ix,jy));
        
        coef4 = Acoef(ix,jy);
	  
        // Combine 4th order at ncx-1 with 2nd order at ncx
        A(ix, 0) = dcomplex(-coef1, 0.0);
        A(ix, 1) = dcomplex(16. * coef1, -coef3);
        A(ix, 2) = dcomplex(-29. * coef1 - SQ(kwave) * coef2 + coef4, 0.0);
        A(ix, 3) = dcomplex(
            (14. -
             SQ(coords->dx(ncx, jy) * kwave) * coords->g33(ncx, jy) / coords->g11(ncx, jy)) *
                coef1,
            coef3);
        A(ix, 4) = 0.0; // Not used

        coef1=coords->g11(ix,jy)/(SQ(coords->dx(ix,jy)));
        coef2=coords->g33(ix,jy);
        coef3= kwave * coords->g13(ix,jy)/(2. * coords->dx(ix,jy));

        // Use 2nd order at ncx - 1
        A(ncx, 0) = dcomplex(coef1, -coef3);
        A(ncx, 1) = dcomplex(-2.0 * coef1 - SQ(kwave) * coef2 + coef4, 0.0);
        A(ncx, 2) = dcomplex(coef1, coef3);
        A(ncx, 3) = 0.0; // Should never be used
        A(ncx, 4) = 0.0;
      }
    }
    
    // Perform inversion
    cband_solve(A, localmesh->LocalNx, 2, 2, bk1d);

    if((global_flags & INVERT_KX_ZERO) && (iz == 0)) {
      // Set the Kx = 0, n = 0 component to zero. For now just subtract
      // Should do in the inversion e.g. Sherman-Morrison formula
        
      dcomplex offset(0.0);
      for(int ix=0;ix<=ncx;ix++)
        offset += bk1d[ix];
      offset /= static_cast<BoutReal>(ncx + 1);
      for(int ix=0;ix<=ncx;ix++)
        bk1d[ix] -= offset;
    }
      
    // Fill xk
    for (int ix=0; ix<=ncx; ix++)
      xk(ix, iz) = bk1d[ix];
  }
  
  // Done inversion, transform back

  for(int ix=0; ix<=ncx; ix++){
    if(global_flags & INVERT_ZERO_DC)
      xk(ix, 0) = 0.0;

    irfft(&xk(ix, 0), ncz, x[ix]);
  }

  checkData(x);
  return x;
#else
  return FieldPerp{};
#endif
}
