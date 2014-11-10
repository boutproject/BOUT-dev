/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

#include <bout/coordinates.hxx>
#include <utils.hxx>
#include <msg_stack.hxx>

Coordinates::Coordinates(Mesh *mesh) {
  dx = 1.0; dy = 1.0; dz = 1.0;
  
  J = 1.0;
  Bxy = 1.0;
    
  // Identity metric tensor
  
  g11 = 1.0; g22 = 1.0; g33 = 1.0;
  g12 = 0.0; g13 = 0.0; g23 = 0.0;
  
  g_11 = 1.0; g_22 = 1.0; g_33 = 1.0;
  g_12 = 0.0; g_13 = 0.0; g_23 = 0.0;
  
  if(mesh->get(dx, "dx")) {
    output.write("\tWARNING: differencing quantity 'dx' not found. Set to 1.0\n");
    dx = 1.0;
  }

  if(mesh->periodicX())
    mesh->communicate(dx);

  if(mesh->get(dy, "dy")) {
    output.write("\tWARNING: differencing quantity 'dy' not found. Set to 1.0\n");
    dy = 1.0;
  }
  
  int zperiod;
  BoutReal ZMIN, ZMAX;
  Options* options = Options::getRoot();
  if(options->isSet("zperiod")) {
    OPTION(options, zperiod, 1);
    ZMIN = 0.0;
    ZMAX = 1.0 / (double) zperiod;
  }else {
    OPTION(options, ZMIN, 0.0);
    OPTION(options, ZMAX, 1.0);
    
    zperiod = ROUND(1.0 / (ZMAX - ZMIN));
  }

  zlength = (ZMAX-ZMIN)*TWOPI;
  dz = zlength/(ngz-1);
  
  // Diagonal components of metric tensor g^{ij} (default to 1)
  mesh->get(g11, "g11", 1.0);
  mesh->get(g22, "g22", 1.0);
  mesh->get(g33, "g33", 1.0);
  
  // Off-diagonal elements. Default to 0
  mesh->get(g12, "g12", 0.0);
  mesh->get(g13, "g13", 0.0);
  mesh->get(g23, "g23", 0.0);
  
  // Check input metrics
  if((!finite(g11)) || (!finite(g22)) || (!finite(g33))) {
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }
  
  /// Calculate contravariant metric components
  if(calcCovariant())
    throw BoutException("Error in calcCovariant call");

  /// Calculate Jacobian and Bxy
  if(jacobian())
    throw BoutException("Error in jacobian call");
  
  // Attempt to read J from the grid file
  Field2D Jcalc = J;
  if(mesh->get(J, "J")) {
    output.write("\tWARNING: Jacobian 'J' not found. Calculating from metric tensor\n");
    J = Jcalc;
  }else {
    // Compare calculated and loaded values  
    output.write("\tMaximum difference in J is %e\n", max(abs(J - Jcalc)));
    
    // Re-evaluate Bxy using new J
    Bxy = sqrt(g_22)/J;
  }

  // Attempt to read Bxy from the grid file
  Field2D Bcalc = Bxy;
  if(mesh->get(Bxy, "Bxy")) {
    output.write("\tWARNING: Magnitude of B field 'Bxy' not found. Calculating from metric tensor\n");
    Bxy = Bcalc;
  }else {
    output.write("\tMaximum difference in Bxy is %e\n", max(abs(Bxy - Bcalc)));
    // Check Bxy
    if(!finite(Bxy))
      throw BoutException("\tERROR: Bxy not finite everywhere!\n");
  }
  
  //////////////////////////////////////////////////////
  /// Calculate Christoffel symbols. Needs communication
  if(geometry()) {
    throw BoutException("Differential geometry failed\n");
  }
  
  //////////////////////////////////////////////////////
  /// Non-uniform meshes. Need to use DDX, DDY
  
  Field2D d2x, d2y; // d^2 x / d i^2
  // Read correction for non-uniform meshes
  if(mesh->get(d2x, "d2x")) {
    output.write("\tWARNING: differencing quantity 'd2x' not found. Calculating from dx\n");
    d1_dx = DDX(1./dx)*dx; // d/di(1/dx)
  }else
    d1_dx = -d2x / (dx*dx);
  
  if(mesh->get(d2y, "d2y")) {
    output.write("\tWARNING: differencing quantity 'd2y' not found. Calculating from dy\n");
    d1_dy = DDY(1./dy)*dy; // d/di(1/dy)
  }else
    d1_dy = -d2y / (dy*dy);
}

Coordinates::~Coordinates() {
  // Gaussj working arrays
  if(ilen > 0) {
    ivfree(indxc);
    ivfree(indxr);
    ivfree(ipiv);
  }
}

void Coordinates::outputVars(Datafile &file) {
  file.add(dx,    "dx",    0);
  file.add(dy,    "dy",    0);
  file.add(dz,    "dz",    0);
  
  file.add(g11,   "g11",   0);
  file.add(g22,   "g22",   0);
  file.add(g33,   "g33",   0);
  file.add(g12,   "g12",   0);
  file.add(g13,   "g13",   0);
  file.add(g23,   "g23",   0);
  
  file.add(g_11,  "g_11",  0);
  file.add(g_22,  "g_22",  0);
  file.add(g_33,  "g_33",  0);
  file.add(g_12,  "g_12",  0);
  file.add(g_13,  "g_13",  0);
  file.add(g_23,  "g_23",  0);
  
  file.add(J,     "J",     0);
}


int Coordinates::geometry() {
#ifdef CHECK
  msg_stack.push("Coordinates::geometry");
#endif
  
  output.write("Calculating differential geometry terms\n");

  if(min(abs(dx)) < 1e-8)
    throw BoutException("dx magnitude less than 1e-8");

  if(min(abs(dy)) < 1e-8)
    throw BoutException("dy magnitude less than 1e-8");

  if(fabs(dz) < 1e-8)
    throw BoutException("dz magnitude less than 1e-8");

  // Check input metrics
  if((!finite(g11)) || (!finite(g22)) || (!finite(g33))) {
    throw BoutException("\tERROR: Diagonal metrics are not finite!\n");
  }
  if((min(g11) <= 0.0) || (min(g22) <= 0.0) || (min(g33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal metrics are negative!\n");
  }
  if((!finite(g12)) || (!finite(g13)) || (!finite(g23))) {
    throw BoutException("\tERROR: Off-diagonal metrics are not finite!\n");
  }
  
  if((!finite(g_11)) || (!finite(g_22)) || (!finite(g_33))) {
    throw BoutException("\tERROR: Diagonal g_ij metrics are not finite!\n");
  }
  if((min(g_11) <= 0.0) || (min(g_22) <= 0.0) || (min(g_33) <= 0.0)) {
    throw BoutException("\tERROR: Diagonal g_ij metrics are negative!\n");
  }
  if((!finite(g_12)) || (!finite(g_13)) || (!finite(g_23))) {
    throw BoutException("\tERROR: Off-diagonal g_ij metrics are not finite!\n");
  }
  
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

  FieldGroup com;

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

  communicate(com);
  
#ifdef CHECK
  msg_stack.pop();
#endif
  
  return 0;
}

int Coordinates::calcCovariant() {
#ifdef CHECK
  msg_stack.push("Coordinates::calcCovariant");
#endif
  
  // Make sure metric elements are allocated
  g_11.allocate();
  g_22.allocate();
  g_33.allocate();
  g_12.allocate();
  g_13.allocate();
  g_23.allocate();
  
  // Perform inversion of g^{ij} to get g_{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects

  BoutReal** a = rmatrix(3, 3);
  
  for(int jx=0;jx<ngx;jx++) {
    for(int jy=0;jy<ngy;jy++) {
      // set elements of g
      a[0][0] = g11(jx, jy);
      a[1][1] = g22(jx, jy);
      a[2][2] = g33(jx, jy);
      
      a[0][1] = a[1][0] = g12(jx, jy);
      a[1][2] = a[2][1] = g23(jx, jy);
      a[0][2] = a[2][0] = g13(jx, jy);
      
      // invert
      if(gaussj(a, 3)) {
	output.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
	return 1;
      }
      
      // put elements into g_{ij}
      g_11(jx, jy) = a[0][0];
      g_22(jx, jy) = a[1][1];
      g_33(jx, jy) = a[2][2];
      
      g_12(jx, jy) = a[0][1];
      g_13(jx, jy) = a[0][2];
      g_23(jx, jy) = a[1][2];
    }
  }

  free_rmatrix(a);
  
  BoutReal maxerr, err;
  maxerr = max(abs( (g_11*g11 +
		     g_12*g12 + 
		     g_13*g13)- 1 ));
  if((err = max(abs( (g_12*g12 +
		      g_22*g22 +
		      g_23*g23) - 1 ))) > maxerr)
    maxerr = err;
  
  if((err = max(abs( (g_13*g13 + 
		      g_23*g23 + 
		      g_33*g33) - 1 ))) > maxerr)
    maxerr = err; 
  output.write("\tMaximum error in diagonal inversion is %e\n", maxerr);
  
  
  maxerr = max(abs(g_11*g12 + 
		   g_12*g22 + 
		   g_13*g23));
  
  if((err = max(abs(g_11*g13 + 
		    g_12*g23 + 
		    g_13*g33))) > maxerr)
    maxerr = err;
  
  if((err = max(abs(g_12*g13 + 
		    g_22*g23 + 
		    g_23*g33))) > maxerr)
    maxerr = err;
  
  output.write("\tMaximum error in off-diagonal inversion is %e\n", maxerr);
  
#ifdef CHECK
  msg_stack.pop();
#endif

  return 0;
}

int Coordinates::calcContravariant() {
  // Make sure metric elements are allocated
  g11.allocate();
  g22.allocate();
  g33.allocate();
  g12.allocate();
  g13.allocate();
  g23.allocate();
  
  // Perform inversion of g_{ij} to get g^{ij}
  // NOTE: Currently this bit assumes that metric terms are Field2D objects
  
  BoutReal** a = rmatrix(3, 3);
  
  for(int jx=0;jx<ngx;jx++) {
    for(int jy=0;jy<ngy;jy++) {
      // set elements of g
      a[0][0] = g_11(jx, jy);
      a[1][1] = g_22(jx, jy);
      a[2][2] = g_33(jx, jy);
      
      a[0][1] = a[1][0] = g_12(jx, jy);
      a[1][2] = a[2][1] = g_23(jx, jy);
      a[0][2] = a[2][0] = g_13(jx, jy);
      
      // invert
      if(gaussj(a, 3)) {
	output.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
	return 1;
      }
      
      // put elements into g_{ij}
      g11(jx, jy) = a[0][0];
      g22(jx, jy) = a[1][1];
      g33(jx, jy) = a[2][2];
      
      g12(jx, jy) = a[0][1];
      g13(jx, jy) = a[0][2];
      g23(jx, jy) = a[1][2];
    }
  }

  free_rmatrix(a);

  BoutReal maxerr, err;
  maxerr = max(abs( (g_11*g11 +
		     g_12*g12 + 
		     g_13*g13)- 1 ));
  if((err = max(abs( (g_12*g12 +
		      g_22*g22 +
		      g_23*g23) - 1 ))) > maxerr)
    maxerr = err;
  
  if((err = max(abs( (g_13*g13 + 
		      g_23*g23 + 
		      g_33*g33) - 1 ))) > maxerr)
    maxerr = err; 
  output.write("\tMaximum error in diagonal inversion is %e\n", maxerr);
  
  
  maxerr = max(abs(g_11*g12 + 
		   g_12*g22 + 
		   g_13*g23));
  
  if((err = max(abs(g_11*g13 + 
		    g_12*g23 + 
		    g_13*g33))) > maxerr)
    maxerr = err;
  
  if((err = max(abs(g_12*g13 + 
		    g_22*g23 + 
		    g_23*g33))) > maxerr)
    maxerr = err;
  
  output.write("\tMaximum error in off-diagonal inversion is %e\n", maxerr);
  return 0;
}

int Coordinates::jacobian() {
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  J = 1. / sqrt(g11*g22*g33 + 
                2.0*g12*g13*g23 - 
                g11*g23*g23 - 
                g22*g13*g13 - 
                g33*g12*g12);
  
  // Check jacobian
  if(!finite(mesh->J)) {
    output.write("\tERROR: Jacobian not finite everywhere!\n");
    return 1;
  }
  if(min(abs(mesh->J)) < 1.0e-10) {
    output.write("\tERROR: Jacobian becomes very small\n");
    return 1;
  }
  
  Bxy = sqrt(mesh->g_22)/mesh->J;
  
  return 0;
}



/*******************************************************************************
 * Gauss-Jordan matrix inversion
 * used to invert metric tensor
 *******************************************************************************/

// Invert an nxn matrix using Gauss-Jordan elimination with full pivoting
int Coordinates::gaussj(BoutReal **a, int n) {
  int i, icol, irow, j, k, l, ll;
  float big, dum, pivinv;

  // Make sure enough temporary memory is allocated
  if(n > ilen) {
    if(ilen == 0) {
      indxc = ivector(n);
      indxr = ivector(n);
      ipiv = ivector(n);
    }else {
      indxc = ivresize(indxc, n);
      indxr = ivresize(indxr, n);
      ipiv = ivresize(ipiv, n);
    }
    ilen = n;
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
	swap(a[irow][l],a[icol][l]);
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
	swap(a[k][indxr[l]], a[k][indxc[l]]);
  }
  // done.

  return 0;
}
