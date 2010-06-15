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

#include "invert_parderiv.h"

#include "globals.h"
#include "utils.h"
#include "mesh_topology.h"
#include "comm_group.h" // Gather/scatter operations

#include "lapack_routines.h" // For tridiagonal inversions

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

namespace invpar { 

  /***********************************************************************
   *                          COMMUNICATIONS
   * 
   * Need to perform a load of gather and scatter operations simultaneously
   * For now, use MPI_Gather etc. to do each operation in turn.
   *
   ***********************************************************************/
  
  static bool comms_init = false; ///< Signals whether the comms are set up
  MPI_Comm comm_in; ///< Communicator

  /// Initialise communicators
  bool init_comms()
  {
    if(comms_init)
      return true; // Already initialised
    
    int ype0 = YPROC(jyseps1_1+1); // The processor at beginning of twist shift
    
    if(NXPE > 1) {
      /// Create a group and communicator for this X location
      MPI_Group group_world;
      
      MPI_Comm_group(MPI_COMM_WORLD, &group_world); // Get the entire group
      
      bout_error("SORRY: invert_parderivs can't cope with NXPE > 1 yet\n");
    }else {
      // All processors. May need to re-number 
      
      if(ype0 == 0) { 
	comm_in = MPI_COMM_WORLD;
      }else {
	// Rotate around so rank0 is at beginning of twist shift
	
	bout_error("SORRY: invert_parderivs can't cope with this topology yet\n");
      }
    }

    comms_init = true;

    return true;
  }
  
  /***********************************************************************
   *                          INVERSION ROUTINES
   * 
   * 
   ***********************************************************************/

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
  void cyclic_solve(int ysize, int xpos, real *data, real *result, bool coef3d = false)
  {
#ifdef CHECK
    msg_stack.push("cyclic_solve(%d, %d)", ysize, xpos);
#endif

    int tshift = ROUND(ShiftAngle[xpos] / mesh->dz); // Nearest neighbour

    static int ylen = 0;
    static bool *done; // Record whether a particular location has been inverted
    static real *avec, *bvec, *cvec; // Matrix coefficients
    static real *rvec, *xvec;
    
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
      // Allocate largest possible array
      avec = new real[ysize * ncz];
      bvec = new real[ysize * ncz];
      cvec = new real[ysize * ncz];
      rvec = new real[ysize * ncz];
      xvec = new real[ysize * ncz];

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
	  bout_error("ERROR: Crossed streams in invpar::cyclic_solve!\n");
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

  /***********************************************************************
   *                         EXTERNAL INTERFACE
   * 
   * 
   ***********************************************************************/
  
  /// Parallel inversion routine
  const Field3D invert_parderiv(const Field2D &A, const Field2D &B, const Field3D &r)
  {
    static real *senddata;
    static real *recvdata;
    static real *resultdata;
    static Comm_handle_t *handles; // Gather and scatter handles
    static bool initialised = false;

#ifdef CHECK
    msg_stack.push("invert_parderiv");
#endif

    // Decide on x range. Solve in boundary conditions
    int xs = (IDATA_DEST < 0) ? 0 : MXG;
    int xe = (ODATA_DEST < 0) ? mesh->ngx-1 : (mesh->ngx-1 - MXG);

    int nxsolve = xe - xs + 1; // Number of X points to solve
    
    if(!initialised) {
      // Allocate working memory
      senddata = new real[nxsolve * MYSUB * (3 + ncz) ]; // Problem data sent out
      recvdata = new real[MY * (3+ncz)];  // Problem data received (to be solved)
      resultdata = new real[MY*ncz];  // Inverted result
      handles = new Comm_handle_t[NYPE];
      initialised = true;

      // Initialise communications
      init_comms();
    }
    
    Field2D sg = sqrt(mesh->g_22); // Needed for first Y derivative

    int x0 = xs; ///< Starting x of this round of inversions
    int offset  = 0; // Location in data array
    for(int xpos=xs; xpos <= xe; xpos++) {
      
      /// Calculate matrix coefficients
      real *dptr = senddata + offset; // Start of this chunk of data
      
      int off0 = offset; // This for debugging only

      // First all the matrix coefficients (2D only in this case)
      for(int j=jstart;j<=jend;j++) {
	// See Grad2_par2 in difops.cpp for these coefficients
	real coeff1 = (1./sg[xpos][j+1] - 1./sg[xpos][j-1])/(4.*SQ(mesh->dy[xpos][j])) / sg[xpos][j];
	real coeff2 = 1. / (mesh->g_22[xpos][j] * SQ(mesh->dy[xpos][j])); // Second derivative
	
	senddata[offset] = B[xpos][j] * (coeff2 - coeff1); // a coefficient (y-1)
	offset++;
	
	senddata[offset] = A[xpos][j] + -2.*B[xpos][j]*coeff2;  // b coefficient (diagonal)
	offset++;
	
	senddata[offset] = B[xpos][j] * (coeff2 + coeff1); // c coefficient (y+1);
	offset++;

	// Then the vector to be inverted (3D)
	for(int k=0;k<ncz;k++) {
	  senddata[offset] = r[xpos][j][k];
	  offset++;
	}
      }

      // Check the amount of data is correct
      if(offset - off0 != MYSUB * (3 + ncz)) {
	output.write("In invert_parderiv: wrong amount of data: %d, %d\n", offset - off0, MYSUB * (3 + ncz));
	bout_error("aborting\n");
      }
      
      /// Gather data onto processors
      
      int yproc = (xpos-x0) % NYPE; // the destination processor

      // Start a gather operation (blocking or nonblocking)
      Comm_gather_start(dptr, MYSUB * (3 + ncz), PVEC_REAL_MPI_TYPE,
			recvdata,
			yproc, comm_in,
			handles+yproc);

      if((yproc == (NYPE-1)) || (xpos == xe)) {
	// Either each processor has a chunk of data, or run out of data
	
	int nsolve = yproc + 1; // Number of x slice being solved

	/// Perform inversion if data is available
	if(nsolve >= PE_YIND) {
	  // Wait for the gather to finish
	  if(!Comm_wait(handles+PE_YIND)) {
	    bout_error("Gather failed\n");
	  }
	  
	  cyclic_solve(NYPE * MYSUB,  // Number of y locations
		       x0 + PE_YIND,  // The x index being solved
		       recvdata,      // Interleaved coefficients and data
		       resultdata);
	}

	// Need to wait for all the gathers to finish (frees memory)
	Comm_wait_all(nsolve, handles);

	// Scatter result back
	for(int xrec = x0; xrec <= xpos; xrec++) { // Loop over the current range
	  yproc = (xrec-x0) % NYPE;
	  Comm_scatter_start(resultdata, MYSUB*ncz, PVEC_REAL_MPI_TYPE,
			     senddata + (xrec-xs)*MYSUB*ncz, // Put result back into senddata
			     yproc, comm_in,
			     handles+yproc);
	}

	// Wait for scatters to finish
	Comm_wait_all(nsolve, handles);

	x0 += NYPE; // Shift starting place for next time
      }
    }
    
    Field3D result;
    result.Allocate();

    // Result is now in senddata. Copy across
    int ind = 0;
    for(int xpos=xs; xpos <= xe; xpos++) {
      for(int j=jstart;j<=jend;j++) {
	for(int k=0;k<ncz;k++) {
	  result[xpos][j][k] = senddata[ind];
	  ind++;
	}
      }
    }
    
#ifdef CHECK
    msg_stack.pop();
#endif

    // done
    return result;
  }

  const Field3D invert_parderiv(real val, const Field2D &B, const Field3D &r)
  {
    Field2D A;
    A = val;
    return invert_parderiv(A, B, r);
  }
  
  const Field3D invert_parderiv(const Field2D &A, real val, const Field3D &r)
  {
    Field2D B;
    B = val;
    return invert_parderiv(A, B, r);
  }
  
  const Field3D invert_parderiv(real val, real val2, const Field3D &r)
  {
    Field2D A, B;
    A = val;
    B = val2;
    return invert_parderiv(A, B, r);
  }
  
} // End of namespace invpar
