/**************************************************************************
 * Differential geometry 
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
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

#include "communicator.h"
#include "geometry.h"
#include "globals.h"
#include "utils.h"
#include "derivs.h"
#include "grid.h"

#include <math.h>
#include <stdlib.h>

int gaussj(real **a, int n);

int geometry()
{
  return geometry(0);
}

int geometry(int readJB)
{
  int jx, jy;
  real **a;

  output.write("Calculating differential geometry terms\n");

  // Check input metrics
  if((!finite(g11)) || (!finite(g22)) || (!finite(g33))) {
    output.write("\tERROR: Diagonal metrics are not finite!\n");
    exit(1);
  }
  if((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    output.write("\tERROR: Diagonal metrics are negative!\n");
  }
  if((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    output.write("\tERROR: Off-diagonal metrics are not finite!\n");
    exit(1);
  }

  // Make sure elements are allocated
  g_11.Allocate();
  g_22.Allocate();
  g_33.Allocate();
  g_12.Allocate();
  g_13.Allocate();
  g_23.Allocate();
  
  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  a = rmatrix(3, 3);

  for(jx=0;jx<ngx;jx++) {
    for(jy=0;jy<ngy;jy++) {
      // set elements of g
      a[0][0] = g11[jx][jy];
      a[1][1] = g22[jx][jy];
      a[2][2] = g33[jx][jy];
      
      a[0][1] = a[1][0] = g12[jx][jy];
      a[1][2] = a[2][1] = g23[jx][jy];
      a[0][2] = a[2][0] = g13[jx][jy];

      
      /*
      output.write("Inverting %e, %e, %e, %e, %e, %e\n", a[0][0], a[1][1], a[2][2],
		   a[0][1], a[1][2], a[0][2]);
      */
      
      
      // invert
      if(gaussj(a, 3)) {
	output.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
	return 1;
      }

      /*
      output.write("       => %e, %e, %e, %e, %e, %e\n", a[0][0], a[1][1], a[2][2],
		   a[0][1], a[1][2], a[0][2]);
      
      if( (a[0][0] <= 0.0) || (a[1][1] <= 0.0) || (a[2][2] <= 0.0) )
	exit(1);
      */

      // put elements into g_{ij}
      g_11[jx][jy] = a[0][0];
      g_22[jx][jy] = a[1][1];
      g_33[jx][jy] = a[2][2];
      
      g_12[jx][jy] = a[0][1];
      g_13[jx][jy] = a[0][2];
      g_23[jx][jy] = a[1][2];
    }
  }

  free_rmatrix(a);

  // Check result metrics
  if((!finite(g_11)) || (!finite(g_22)) || (!finite(g_33))) {
    output.write("\tERROR: Diagonal g_ij metrics are not finite!\n");
    exit(1);
  }
  if((min(g_11) <= 0.0) || (min(g_22) <= 0.0) || (min(g_33) <= 0.0)) {
    output.write("\tERROR: Diagonal g_ij metrics are negative!\n");
    exit(1);
  }
  if((!finite(g_12)) || (!finite(g_13)) || (!finite(g_23))) {
    output.write("\tERROR: Off-diagonal g_ij metrics are not finite!\n");
    exit(1);
  }
  
  // Test the inversion: multiply g_ai*g^ib = delta_a^b
  real err, maxerr;
  
  maxerr = max(abs( (g_11*g11 + g_12*g12 + g_13*g13) - 1 ));
  if((err = max(abs( (g_12*g12 + g_22*g22 + g_23*g23) - 1 ))) > maxerr)
    maxerr = err;
  if((err = max(abs( (g_13*g13 + g_23*g23 + g_33*g33) - 1 ))) > maxerr)
    maxerr = err; 
  output.write("\tMaximum error in diagonal inversion is %e\n", maxerr);
  
  maxerr = max(abs(g_11*g12 + g_12*g22 + g_13*g23));
  if((err = max(abs(g_11*g13 + g_12*g23 + g_13*g33))) > maxerr)
    maxerr = err;
  if((err = max(abs(g_12*g13 + g_22*g23 + g_23*g33))) > maxerr)
    maxerr = err;
  output.write("\tMaximum error in off-diagonal inversion is %e\n", maxerr);

  // calculate using g^-1 = det[g^ij], J = sqrt(g)
  J = (g11*g22*g33 + 2.0*g12*g13*g23 - g11*g23*g23 - g22*g13*g13 - g33*g12*g12)^(-0.5);
  
  if(readJB) {
    // Attempt to read J from the grid file
    
    Field2D Jcalc = J;
    if(grid_load2d(J, "J")) {
      output.write("\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
      J = Jcalc;
    }else {
      // Compare calculated and loaded values
      
      output.write("\tMaximum difference in J is %e\n", max(abs(J - Jcalc)));
    }   
  }
  
  // Check jacobian
  if(!finite(J)) {
    output.write("\tERROR: Jacobian not finite everywhere!\n");
    exit(1);
  }
  if(min(abs(J)) < 1.0e-10) {
    output.write("\tERROR: Jacobian becomes very small\n");
    exit(0);
  }

  Bxy = sqrt(g_22)/J;

  if(readJB) {
    // Attempt to read Bxy from the grid file

    Field2D Bcalc = Bxy;

    if(grid_load2d(Bxy, "Bxy")) {
      output.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from metric tensor\n");
      Bxy = Bcalc;
    }else {
      output.write("\tMaximum difference in Bxy is %e\n", max(abs(Bxy - Bcalc)));
    }
  }

  // Check Bxy
  if(!finite(Bxy)) {
    output.write("\tERROR: Bxy not finite everywhere!\n");
    exit(1);
  }

  geometry_derivs();
  
  return 0;
}

void geometry_derivs()
{
  // Calculate Christoffel symbol terms (15 independent values)
  // Note: This calculation is completely general: metric 
  // tensor can be 2D or 3D. For 2D, all DDZ terms are zero

  G1_11 = 0.5*g11*DDX(g_11) 
    + g12*(DDX(g_12) - 0.5*DDY(g_11))
    + g13*(DDX(g_13) - 0.5*DDZ(g_11));
  G1_22 = g11*(DDY(g_12) - 0.5*DDX(g_22))
    + 0.5*g12*DDY(g_22)
    + g13*(DDY(g_23) - 0.5*DDZ(g_22));
  G1_33 = g11*(DDZ(g_13) - 0.5*DDX(g_33))
    + g12*(DDZ(g_23) - 0.5*DDY(g_33))
    + 0.5*g13*DDZ(g_33);
  G1_12 = 0.5*g11*DDY(g_11)
    + 0.5*g12*DDX(g_22)
    + 0.5*g13*(DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G1_13 = 0.5*g11*DDZ(g_11)
    + 0.5*g12*(DDZ(g_12) + DDX(g_23) - DDY(g_13))
    + 0.5*g13*DDX(g_33);

  G2_11 = 0.5*g12*DDX(g_11)
    + g22*(DDX(g_12) - 0.5*DDY(g_11))
    + g23*(DDX(g_13) - 0.5*DDZ(g_11));
  G2_22 = g12*(DDY(g_12) - 0.5*DDX(g_22))
    + 0.5*g22*DDY(g_22)
    + g23*(DDY(g23) - 0.5*DDZ(g_22));
  G2_33 = g12*(DDZ(g_13) - 0.5*DDX(g_33))
    + g22*(DDZ(g_23) - 0.5*DDY(g_33))
    + 0.5*g23*DDZ(g_33);
  G2_12 = 0.5*g12*DDY(g_11)
    + 0.5*g22*DDX(g_22)
    + 0.5*g23*(DDY(g_13) + DDX(g_23) - DDZ(g_12));
  G2_23 = 0.5*g12*(DDZ(g_12) + DDY(g_13) - DDX(g_23))
    + 0.5*g22*DDZ(g_22)
    + 0.5*g23*DDY(g_33);
  
  G3_11 = 0.5*g13*DDX(g_11)
    + g23*(DDX(g_12) - 0.5*DDY(g_11))
    + g33*(DDX(g_13) - 0.5*DDZ(g_11));
  G3_22 = g13*(DDY(g_12) - 0.5*DDX(g_22))
    + 0.5*g23*DDY(g_22)
    + g33*(DDY(g_23) - 0.5*DDZ(g_22));
  G3_33 = g13*(DDZ(g_13) - 0.5*DDX(g_33))
    + g23*(DDZ(g_23) - 0.5*DDY(g_33))
    + 0.5*g33*DDZ(g_33);
  G3_13 = 0.5*g13*DDZ(g_11)
    + 0.5*g23*(DDZ(g_12) + DDX(g_23) - DDY(g_13))
    + 0.5*g33*DDX(g_33);
  G3_23 = 0.5*g13*(DDZ(g_12) + DDY(g_13) - DDX(g_23))
    + 0.5*g23*DDZ(g_22)
    + 0.5*g33*DDY(g_33);

  G1 = (DDX(J*g11) + DDY(J*g12) + DDZ(J*g13))/J;
  G2 = (DDX(J*g12) + DDY(J*g22) + DDZ(J*g23))/J;
  G3 = (DDX(J*g13) + DDY(J*g23) + DDZ(J*g33))/J;

  // Communicate christoffel symbol terms
  output.write("\tCommunicating connection terms\n");

  Communicator com;

  com.add(G1_11);
  com.add(G1_22);
  com.add(G1_33);
  com.add(G1_12);
  com.add(G1_13);
  
  com.add(G2_11);
  com.add(G2_22);
  com.add(G2_33);
  com.add(G2_12);
  com.add(G2_23);
  
  com.add(G3_11);
  com.add(G3_22);
  com.add(G3_33);
  com.add(G3_13);
  com.add(G3_23);

  com.add(G1);
  com.add(G2);
  com.add(G3);

  com.run();
}

/*******************************************************************************
 * Gauss-Jordan matrix inversion
 * used to invert metric tensor
 *******************************************************************************/

#define SWAP(a,b) {temp=(a); (a) = (b); (b) = temp;}

// Invert an nxn matrix using Gauss-Jordan elimination with full pivoting
int gaussj(real **a, int n)
{
  static int *indxc, *indxr, *ipiv, len = 0;
  int i, icol, irow, j, k, l, ll;
  float big, dum, pivinv, temp;

  // Make sure enough temporary memory is allocated
  if(n > len) {
    if(len == 0) {
      indxc = ivector(n);
      indxr = ivector(n);
      ipiv = ivector(n);
    }else {
      indxc = ivresize(indxc, n);
      indxr = ivresize(indxr, n);
      ipiv = ivresize(ipiv, n);
    }
    len = n;
  }

  for(i=0;i<n;i++)
    ipiv[i] = 0;

  for(i=0;i<n;i++) { // Main loop over columns
    big = 0.0;
    irow = icol = -1;
    for(j=0;j<n;j++) { // search for pivot element
      if(ipiv[j] != 1) {
	for(k=0;k<n;k++) {
	  if(ipiv[k] == 0) {
	    if(fabs(a[j][k]) >= big) {
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }else if(ipiv[k] > 1) {
	    output.write("Error in GaussJ: Singular matrix-1\n");
	    return 1;
	  }
	}
      }
    }
    
    if(irow == -1) {
      // All elements zero!!
      output.write("Error in GaussJ: Singular matrix-3\n");
      return 3;
    }

    ++(ipiv[icol]);
    // Now have pivot element, so interchange rows to put pivot
    // on the diagonal
    if(irow != icol) {
      for(l=0;l<n;l++)
	SWAP(a[irow][l],a[icol][l]);
    }
    indxr[i] = irow;
    indxc[i] = icol;

    if(a[icol][icol] == 0.0) {
      output.write("Error in GaussJ: Singular matrix-2\n");
      return 2;
    }
    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<n;l++)
      a[icol][l] *= pivinv;

    for(ll=0;ll<n;ll++) { // reduce rows
      if(ll != icol) {    // except for the pivot one
	dum = a[ll][icol];
	a[ll][icol] = 0.0;
	for(l=0;l<n;l++)
	  a[ll][l] -= a[icol][l]*dum;
	
      }
    }
  }
  // end of main loop. Unscramble solution due to column interchanges
  for(l=n-1;l>=0;l--) {
    if(indxr[l] != indxc[l])
      for(k=0;k<n;k++)
	SWAP(a[k][indxr[l]], a[k][indxc[l]]);
  }
  // done.

  return 0;
}
