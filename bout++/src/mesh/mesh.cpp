
#include "mesh.h"
#include "globals.h"
#include "utils.h"
#include "derivs.h"

#include <cmath>

BoutReal Mesh::wtime_comms = 0.0;

/**************************************************************************
 * Data sources
 **************************************************************************/

int Mesh::addSource(GridDataSource &source)
{
  return addSource(&source);
}

int Mesh::addSource(GridDataSource *source)
{
  source_list.push_front(source);
  return 0;
}

int Mesh::load(GridDataSource &source)
{
  std::list<GridDataSource*> old_list;
  
  old_list = source_list;
  source_list.clear();
  
  if(addSource(source))
    return 1;
  
  return load();
}

GridDataSource* Mesh::findSource(const char *name)
{
#ifdef CHECK
  int msg_point = msg_stack.push("Finding source");
#endif
  for(std::list<GridDataSource*>::iterator it = source_list.begin(); 
      it != source_list.end(); it++) {
    
    // Query this source
    if((*it) != NULL)
      if((*it)->hasVar(name)) {
#ifdef CHECK
	msg_stack.pop(msg_point);
#endif
	return *it;
      }
  }
#ifdef CHECK
  msg_stack.pop(msg_point);
#endif
  return NULL;
}

/**************************************************************************
 * Data get routines
 **************************************************************************/

int Mesh::get(Vector2D &var, const char *name)
{
  return get(var, string(name));
}

int Mesh::get(Vector3D &var, const char *name)
{
  return get(var, string(name));
}

int Mesh::get(Vector2D &var, const string &name)
{
#ifdef CHECK
  msg_stack.push("Loading 2D vector: Mesh::get(Vector2D, %s)", name.c_str());
#endif

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;
    
    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");
    
  }else {
    output << "\tReading contravariant vector " << name << endl;
    
    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif  

  return 0;
}

int Mesh::get(Vector3D &var, const string &name)
{
#ifdef CHECK
  msg_stack.push("Loading 3D vector: Mesh::get(Vector3D, %s)", name.c_str());
#endif

  if(var.covariant) {
    output << "\tReading covariant vector " << name << endl;
    
    get(var.x, name+"_x");
    get(var.y, name+"_y");
    get(var.z, name+"_z");
    
  }else {
    output << "\tReading contravariant vector " << name << endl;
    
    get(var.x, name+"x");
    get(var.y, name+"y");
    get(var.z, name+"z");
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif  

  return 0;
}

/**************************************************************************
 * Communications
 **************************************************************************/

int Mesh::communicate(FieldData &f)
{
  FieldGroup group;
  group.add(f);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2)
{
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3)
{
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4)
{
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  group.add(f4);
  return communicate(group);
}

comm_handle Mesh::send(FieldData &f)
{
  FieldGroup group;
  group.add(f);
  return send(group);
}

/// This is a bit of a hack for now to get FieldPerp communications
/// The FieldData class needs to be changed to accomodate FieldPerp objects
int Mesh::communicate(FieldPerp &f)
{
  comm_handle recv[2];
  
  BoutReal **fd = f.getData();
  
  int nin = xstart; // Number of x points in inner guard cell
  int nout = ngx-xend-1; // Number of x points in outer guard cell
  
  // Post receives for guard cell regions
  recv[0] = irecvXIn(fd[0],       nin*ngz, 0);
  recv[1] = irecvXOut(fd[xend+1], nout*ngz, 1);
  
  // Send data
  sendXIn(fd[xstart], nin*ngz, 1);
  sendXOut(fd[xend-nout+1], nout*ngz, 0);
  
  return 0;
}

/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

int Mesh::geometry()
{
#ifdef CHECK
  msg_stack.push("Mesh::geometry");
#endif
  
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

int Mesh::calcCovariant()
{
#ifdef CHECK
  msg_stack.push("Mesh::calcCovariant");
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
      a[0][0] = g11[jx][jy];
      a[1][1] = g22[jx][jy];
      a[2][2] = g33[jx][jy];
      
      a[0][1] = a[1][0] = g12[jx][jy];
      a[1][2] = a[2][1] = g23[jx][jy];
      a[0][2] = a[2][0] = g13[jx][jy];
      
      // invert
      if(gaussj(a, 3)) {
	output.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
	return 1;
      }
      
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

int Mesh::calcContravariant()
{
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
      a[0][0] = g_11[jx][jy];
      a[1][1] = g_22[jx][jy];
      a[2][2] = g_33[jx][jy];
      
      a[0][1] = a[1][0] = g_12[jx][jy];
      a[1][2] = a[2][1] = g_23[jx][jy];
      a[0][2] = a[2][0] = g_13[jx][jy];
      
      // invert
      if(gaussj(a, 3)) {
	output.write("\tERROR: metric tensor is singular at (%d, %d)\n", jx, jy);
	return 1;
      }
      
      // put elements into g_{ij}
      g11[jx][jy] = a[0][0];
      g22[jx][jy] = a[1][1];
      g33[jx][jy] = a[2][2];
      
      g12[jx][jy] = a[0][1];
      g13[jx][jy] = a[0][2];
      g23[jx][jy] = a[1][2];
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

int Mesh::jacobian()
{
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)
  J = sqrt(g11*g22*g33 + 
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
int Mesh::gaussj(BoutReal **a, int n)
{
  static int *indxc, *indxr, *ipiv, len = 0;
  int i, icol, irow, j, k, l, ll;
  float big, dum, pivinv;

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
