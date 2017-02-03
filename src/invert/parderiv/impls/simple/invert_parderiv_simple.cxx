/************************************************************************
 * Inversion of parallel derivatives
 * Intended for use in preconditioner for reduced MHD
 * 
 * Inverts a matrix of the form 
 *
 * (A + B * Grad2_par2) x = r
 * 
 * Stages:
 * - Problem trivially parallel in X, so gather all data for fixed X onto
 *   a single processor. Split MXSUB locations between NYPE processors
 * - Use nearest neighbour for twist-shift location.
 *   This splits the problem into one or more (cyclic) tridiagonal problems
 *   (number depends on the q value)
 * - Solve each of these tridiagonal systems O(Nz*Ny)
 * - Scatter data back 
 *
 * Author: Ben Dudson, University of York, June 2009
 * 
 * Known issues:
 * ------------
 *
 * - Assumes all points are in the core: no handling of SOL / PF regions
 * - This algorithm can only use MXSUB processors, so if NYPE > MXSUB
 *   (i.e. Np > Nx) then efficiency will fall.
 * - When higher-order interpolation is used at the twist-shift location,
 *   errors can become large
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
 ************************************************************************/

#include "mpi.h"

#include <invert_parderiv_simple.hxx>

#include <globals.hxx>
#include <utils.hxx>
#include <comm_group.hxx> // Gather/scatter operations

#include <lapack_routines.hxx> // For tridiagonal inversions
 
/// Parallel inversion routine
const Field3D InvertParSimple::solve(const Field3D &rc) {
  static BoutReal *senddata;
  static BoutReal *recvdata;
  static int max_size = 0;
  static BoutReal *resultdata;

#ifdef CHECK
  msg_stack.push("invert_parderiv");
#endif

  // Copy (to get rid of const)
  Field3D r = rc;

  // Decide on x range. Solve in boundary conditions
  int xs = (mesh->firstX()) ? 0 : 2;
  int xe = (mesh->lastX()) ? mesh->LocalNx-1 : (mesh->LocalNx-3);

  int nxsolve = xe - xs + 1; // Number of X points to solve    
  int nylocal = mesh->yend - mesh->ystart + 1;

  if(max_size == 0) {
    // allocate working memory
    senddata = new BoutReal[nxsolve * nylocal * (2 + mesh->LocalNz) ]; // Problem data sent out
  }
    
  // coefficients for derivative term
  Field2D coeff1;
  coeff1.allocate();
  Field2D sg = sqrt(mesh->g_22); // Needed for first Y derivative
  for(int i=xs;i<=xe;i++)
    for(int j=mesh->ystart;j<=mesh->yend;j++) {
      // See Grad2_par2 in difops.cpp for these coefficients
      coeff1[i][j] = (1./sg[i][j+1] - 1./sg[i][j-1])/(4.*SQ(mesh->dy[i][j])) / sg[i][j];
    }
    
  Field2D coeff2 = 1. / (mesh->g_22 * SQ(mesh->dy)); // Second derivative
    
  // coefficients in tridiagonal matrix
  Field2D acoeff = B * (coeff2 - coeff1); // a coefficient (y-1)
  Field2D bcoeff = A + -2.*B*coeff2; // b coefficient (diagonal)
  Field2D ccoeff = B * (coeff2 + coeff1); // c coefficient (y+1)

  // Create a group of fields to send together
  FieldGroup sendfields;
  sendfields.add(acoeff);
  sendfields.add(bcoeff);
  sendfields.add(ccoeff);
  sendfields.add(r);
    
  // Create a field for the result
  Field3D result;

  // Create an iterator over surfaces
  DistribSurfaceIter* surf = mesh->iterateSurfacesDistrib();

  // Iterate over the surfaces
  for(surf->first(); !surf->isDone(); surf->next()) {
    int ysize = surf->ySize(); // Size of this surface
      
    // Make sure our arrays are big enough
    if(ysize > max_size) {
      // Need to allocate more memory
      if(max_size > 0) {
        delete[] recvdata;
        delete[] resultdata;
      }
      recvdata = new BoutReal[ysize * (3+mesh->LocalNz)];  // Problem data received (to be solved)
      resultdata = new BoutReal[ysize*mesh->LocalNz];  // Inverted result
      max_size = ysize;
    }
      
    // Gather data
    surf->gather(sendfields, recvdata);
      
    // Check that this processor still has something to do
    if(ysize > 0) {
      // Perform inversion
      BoutReal ts; // Twist-shift matching (if closed)
      if(surf->closed(ts)) {
        cyclicSolve(ysize,  // Number of y locations
                    surf->xpos,  // The x index being solved
                    recvdata,    // Interleaved coefficients and data
                    resultdata);
      }else {
        // Open field-lines
        throw BoutException("Sorry; invertParSimple can't cope with open field lines yet\n");
      }
    }
      
    // Scatter the result back
    surf->scatter(resultdata, result);
  }
    
#ifdef CHECK
  msg_stack.pop();
#endif

  // done
  return result;
}

/// Solve the cyclic problem for core points
/*!
 * Assumes that the twist-shift location is between the last and first point
 *
 * @param[in]    ysize    Number of y points
 * @param[in]    coeffs   Matrix coefficients
 * @param[in]    xpos     X location. Needed for shift
 * @param[in]    data     Data to be inverted. Input is interleaved
 * @param[out]   result
 */
void InvertParSimple::cyclicSolve(int ysize, int xpos, BoutReal *data, BoutReal *result, bool coef3d) {
#ifdef CHECK
  msg_stack.push("cyclic_solve(%d, %d)", ysize, xpos);
#endif

  BoutReal ts; // Twist-shift angle
  if(!mesh->periodicY(xpos,ts))
    ts = 0.; // Should be an error, but just separate field-lines
    
  int tshift = ROUND(ts / mesh->dz); // Nearest neighbour

  static int ylen = 0;
  static bool *done; // Record whether a particular location has been inverted
  static BoutReal *avec, *bvec, *cvec; // Matrix coefficients
  static BoutReal *rvec, *xvec;
    
  int ncz = mesh->LocalNz;
    
  if(ysize > ylen) {
    // Initialise
    if(ylen > 0) {
      delete[] done;
      delete[] avec;
      delete[] bvec;
      delete[] cvec;
      delete[] rvec;
      delete[] xvec;
    }
    done = new bool[ncz];
    //.allocate largest possible array
    avec = new BoutReal[ysize * ncz];
    bvec = new BoutReal[ysize * ncz];
    cvec = new BoutReal[ysize * ncz];
    rvec = new BoutReal[ysize * ncz];
    xvec = new BoutReal[ysize * ncz];

    ylen = ysize;
  }

  for(int i=0;i<ncz;i++)
    done[i] = false;
    
  int z0 = 0; // Starting z location
  do {
    /////////////////////////////////////////////
    // Get the matrix coefficients. Could invert
    // everything in-place, but simpler to use new arrays
      
    int zpos = z0;
    int ind = 0; // Index in matrix
    do {
      if(coef3d) {
        // Matrix coefficients depend on Z
        for(int y=0;y<ysize;y++) {
          avec[ind] = data[(4*ncz)*y         + zpos]; 
          bvec[ind] = data[(4*ncz)*y + ncz   + zpos];
          cvec[ind] = data[(4*ncz)*y + 2*ncz + zpos];
	    
          rvec[ind] = data[(4*ncz)*y + 3*ncz + zpos];
	    
          ind++;
        }
      }else {
        // Matrix coefficients 2D. Each y location has [a,b,c,<RHS>]
        for(int y=0;y<ysize;y++) {
          avec[ind] = data[(3+ncz)*y]; 
          bvec[ind] = data[(3+ncz)*y+1];
          cvec[ind] = data[(3+ncz)*y+2];
	    
          rvec[ind] = data[(3+ncz)*y + 3 + zpos];
          ind++;
        }
      }

      done[zpos] = true; // Mark as solved

      // Get the next z index
      zpos = (((zpos + tshift) % ncz) + ncz) % ncz; // Next location

      if((zpos != z0) && done[zpos]) {
        // Somehow hit a different fieldline. Should never happen
        throw BoutException("ERROR: Crossed streams in invpar::cyclic_solve!\n");
      }
	
    }while(zpos != z0);
      
    /////////////////////////////////////////////
    // Solve cyclic tridiagonal system
      
    cyclic_tridag(avec, bvec, cvec, rvec, xvec, ind);
      
    /////////////////////////////////////////////
    // Copy result back
      
    zpos = z0;
    ind = 0;
    do {
      for(int y=0;y<ysize;y++) {
        result[y*ncz + zpos] = xvec[ind];
        ind++;
      }
      zpos = (((zpos + tshift) % ncz) + ncz) % ncz; // Next location
    }while(zpos != z0);

    /////////////////////////////////////////////

    // Find the next place to start
    z0 = -1;
    for(int i=0;i<ncz;i++) {
      if(!done[i])
        z0 = i;
    }
  }while(z0 != -1); // Keep going until everything's been inverted
    
#ifdef CHECK
  msg_stack.pop();
#endif
}

