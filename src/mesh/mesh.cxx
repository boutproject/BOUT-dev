
#include <globals.hxx>
#include <bout/mesh.hxx>
#include <utils.hxx>
#include <derivs.hxx>
#include <msg_stack.hxx>

#include <cmath>

#include "meshfactory.hxx"

#include <output.hxx>

Mesh* Mesh::create(GridDataSource *s, Options *opt) {
  return MeshFactory::getInstance()->createMesh(s, opt);
}

Mesh* Mesh::create(Options *opt) {
  return create(NULL, opt);
}

Mesh::Mesh(GridDataSource *s) : source(s) {
  if(s == NULL)
    throw BoutException("GridDataSource passed to Mesh::Mesh() is NULL");

  ilen = 0; // For gaussj routine

  // Will be set to true if any variable has a free boundary condition applied to the corresponding boundary
  freeboundary_xin = false;
  freeboundary_xout = false;
  freeboundary_ydown = false;
  freeboundary_yup = false;
}

Mesh::~Mesh() {
  delete source;

  // Gaussj working arrays
  if(ilen > 0) {
    ivfree(indxc);
    ivfree(indxr);
    ivfree(ipiv);
  }
}

/**************************************************************************
 * Functions for reading data from external sources
 *
 * These functions are delegated to a GridDataSource object,
 * which may then read from a file, options, or other sources.
 **************************************************************************/

/// Get an integer
int Mesh::get(int &ival, const string &name) {
  MsgStackItem msg("Mesh::get(ival)");

  if(!source->get(this, ival, name))
    return 1;

  return 0;
}

/// A BoutReal number
int Mesh::get(BoutReal &rval, const string &name) {
  MsgStackItem msg("Mesh::get(rval)");

  if(!source->get(this, rval, name))
    return 1;

  return 0;
}

int Mesh::get(Field2D &var, const string &name, BoutReal def) {
  MsgStackItem msg("Loading 2D field: Mesh::get(Field2D)");

  // Ensure data allocated
  var.allocate();

  if(!source->get(this, var, name, def))
    return 1;

  // Communicate to get guard cell data
  Mesh::communicate(var);

  // Check that the data is valid
  var.checkData(true);

  return 0;
}

int Mesh::get(Field3D &var, const string &name, BoutReal def) {
  MsgStackItem msg("Loading 3D field: Mesh::get(Field3D)");

  // Ensure data allocated
  var.allocate();

  if(!source->get(this, var, name, def))
    return 1;

  // Communicate to get guard cell data
  Mesh::communicate(var);

  // Check that the data is valid
  var.checkData(true);

  return 0;
}

/**************************************************************************
 * Data get routines
 **************************************************************************/

int Mesh::get(Vector2D &var, const string &name) {
  msg_stack.push("Loading 2D vector: Mesh::get(Vector2D, %s)", name.c_str());

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

  msg_stack.pop();

  return 0;
}

int Mesh::get(Vector3D &var, const string &name) {
  msg_stack.push("Loading 3D vector: Mesh::get(Vector3D, %s)", name.c_str());

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

  msg_stack.pop();

  return 0;
}

/**************************************************************************
 * Communications
 **************************************************************************/

int Mesh::communicate(FieldData &f) {
  FieldGroup group;
  group.add(f);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2) {
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3) {
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  return communicate(group);
}

int Mesh::communicate(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4) {
  FieldGroup group;
  group.add(f1);
  group.add(f2);
  group.add(f3);
  group.add(f4);
  return communicate(group);
}

comm_handle Mesh::send(FieldData &f) {
  FieldGroup group;
  group.add(f);
  return send(group);
}

/// This is a bit of a hack for now to get FieldPerp communications
/// The FieldData class needs to be changed to accomodate FieldPerp objects
int Mesh::communicate(FieldPerp &f) {
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

  // Wait for receive
  wait(recv[0]);
  wait(recv[1]);

  return 0;
}

int Mesh::msg_len(const vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt) {
  int len = 0;

  /// Loop over variables
  for(std::vector<FieldData*>::const_iterator it = var_list.begin(); it != var_list.end(); it++) {
    if((*it)->is3D()) {
      len += (xlt - xge) * (ylt - yge) * (ngz-1) * (*it)->BoutRealSize();
    }else
      len += (xlt - xge) * (ylt - yge) * (*it)->BoutRealSize();
  }

  return len;
}

bool Mesh::periodicY(int jx) const {
  BoutReal ts; return periodicY(jx, ts);
}

int Mesh::ySize(int jx) const {
  // Get the size of a surface in Y using MPI communicator
  MPI_Comm comm = getYcomm(jx);

  int local = yend - ystart + 1;
  int all;
  MPI_Allreduce(&local, &all, 1, MPI_INT, MPI_SUM, comm);
  return all;
}

bool Mesh::hasBndryLowerY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = (int) !(iterateBndryLowerY().isDone());
  int allbndry;
  MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm());
  answer = (bool) allbndry;
  calc = true;
  return answer;
}

bool Mesh::hasBndryUpperY() {
  static bool calc = false, answer;
  if(calc) return answer; // Already calculated

  int mybndry = (int) !(iterateBndryUpperY().isDone());
  int allbndry;
  MPI_Allreduce(&mybndry, &allbndry, 1, MPI_INT, MPI_BOR, getXcomm());
  answer = (bool) allbndry;
  calc = true;
  return answer;
}

/**************************************************************************
 * Differential geometry
 * Calculates the covariant metric tensor, and christoffel symbol terms
 * given the contravariant metric tensor terms
 **************************************************************************/

int Mesh::geometry() {
  TRACE("Mesh::geometry");

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
    + g23*(DDY(g_23) - 0.5*DDZ(g_22));
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

  // Note that G1, G2 and G3 has nothing to do with the Cristoffel symbols, but
  // serve as variables which only has to be calculated once
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

  return 0;
}

int Mesh::calcCovariant() {
  TRACE("Mesh::calcCovariant");

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

  return 0;
}

int Mesh::calcContravariant() {
  TRACE("Mesh::calcContravariant");

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

int Mesh::jacobian() {
  TRACE("Mesh::jacobian");
  // calculate Jacobian using g^-1 = det[g^ij], J = sqrt(g)

  Field2D g = g11*g22*g33 +
    2.0*g12*g13*g23 -
    g11*g23*g23 -
    g22*g13*g13 -
    g33*g12*g12;

  // Check that g is positive
  if(min(g) < 0.0) {
    throw BoutException("The determinant of g^ij is somewhere less than 0.0");
  }
  J = 1. / sqrt(g);

  // Check jacobian
  if(!finite(J)) {
    throw BoutException("\tERROR: Jacobian not finite everywhere!\n");
  }
  if(min(abs(J)) < 1.0e-10) {
    throw BoutException("\tERROR: Jacobian becomes very small\n");
  }

  if(min(g_22) < 0.0)
    throw BoutException("g_22 is somewhere less than 0.0");
  Bxy = sqrt(g_22)/J;

  return 0;
}

const vector<int> Mesh::readInts(const string &name, int n) {
  vector<int> result;

  if(source->hasVar(name)) {
    if(!source->get(this, result, name, n, 0)) {
      // Error reading
      throw BoutException("Could not read integer array '%s'\n", name.c_str());
    }
  }else {
    // Not found
    throw BoutException("Missing integer array %s\n", name.c_str());
  }

  return result;
}

/*******************************************************************************
 * Gauss-Jordan matrix inversion
 * used to invert metric tensor
 *******************************************************************************/

// Invert an nxn matrix using Gauss-Jordan elimination with full pivoting
int Mesh::gaussj(BoutReal **a, int n) {
  TRACE("Mesh::gaussj");
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
            throw BoutException("Error in GaussJ: Singular matrix-1\n");
          }
        }
      }
    }

    if(irow == -1) {
      // All elements zero!!
      throw BoutException("Error in GaussJ: Singular matrix-3\n");
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
      throw BoutException("Error in GaussJ: Singular matrix-2\n");
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

/*******************************************************************************
 * AverageY
 *******************************************************************************/

/// Not very efficient version, as a fallback
const Field3D Mesh::averageY(const Field3D &f) {
  TRACE("Mesh::averageY");

  Field3D result;
  result.allocate();

  Field2D xy;
  xy.allocate();

  for(int jz=0; jz<ngz;jz++) {

    for(int jx=0; jx<ngx;jx++) {
      for(int jy=0; jy<ngy;jy++) {
        xy(jx, jy) = f(jx, jy, jz);
      }
    }

    xy = averageY(xy);

    for(int jx=0; jx<ngx;jx++) {
      for(int jy=0; jy<ngy;jy++) {
        result(jx, jy, jz) = xy(jx, jy);
      }
    }
  }

  return result;
}

/*******************************************************************************
 * AverageX
 *******************************************************************************/

/// Not very efficient version, as a fallback
const Field3D Mesh::averageX(const Field3D &f) {
  TRACE("Mesh::averageX");

  Field3D result;
  result.allocate();

  Field2D xy;
  xy.allocate();

  for(int jz=0; jz<ngz;jz++) {

    for(int jx=0; jx<ngx;jx++) {
      for(int jy=0; jy<ngy;jy++) {
        xy(jx, jy) = f(jx, jy, jz);
      }
    }

    xy = averageX(xy);

    for(int jx=0; jx<ngx;jx++) {
      for(int jy=0; jy<ngy;jy++) {
        result(jx, jy, jz) = xy(jx, jy);
      }
    }
  }

  return result;
}
