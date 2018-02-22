/*!
 * \file lapack_routines.cxx
 *
 * Serial code to invert a complex tridiagonal system
 *
 * Complex banded matrix solver
 *
 * Solves a banded matrix given the matrix in compact form
 * a[0...(n-1)][0...(m1+m2)]
 * and the rhs vector
 * b[0...(n-1)]
 * 
 * a is overwritten, and b is replaced by the solution
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
 */


#include <bout/globals.hxx>
#include <bout/dcomplex.hxx>
#include <bout/boutexception.hxx>

#ifdef LAPACK

// LAPACK prototypes
extern "C" {
  /// Complex tridiagonal inversion
  void zgtsv_(int *n, int *nrhs, fcmplx *dl, fcmplx *d, fcmplx *du, fcmplx * b, int *ldb, int *info);
  /// BoutReal (double) tridiagonal inversion
  void dgtsv_(int *n, int *nrhs, BoutReal *dl, BoutReal *d, BoutReal *du, BoutReal *b, int *ldb, int *info); 
  /// Complex band solver
  void zgbsv_(int *n, int *kl, int *ku, int *nrhs, fcmplx *ab, int *ldab, int *ipiv, fcmplx *b, int *ldb, int *info);
}

/// Use LAPACK routine ZGTSV. About 25% slower than the simple NR routine for 260 points
int tridag(const dcomplex *a, const dcomplex *b, const dcomplex *c, const dcomplex *r, dcomplex *u, int n) {
  
  // Lapack routines overwrite their inputs, so need to copy

  // allocate memory
  Array<fcmplx> dl(n), d(n), du(n), x(n);

  for (int i = 0; i < n; i++) {
    // Diagonal
    d[i].r = b[i].real();
    d[i].i = b[i].imag();

    // Off-diagonal terms
    if (i != (n - 1)) {
      dl[i].r = a[i + 1].real();
      dl[i].i = a[i + 1].imag();
      
      du[i].r = c[i].real();
      du[i].i = c[i].imag();
    }

    x[i].r = r[i].real();
    x[i].i = r[i].imag();
  }

  /* LAPACK ZGTSV routine.

     n - Size of the array
     nrhs - Number of RHS vectors to solve
     dl - lower band
     d - diagonal values
     du - upper band
     x  - input and output values
     info - output status
  */
  int nrhs = 1;
  int info;
  zgtsv_(&n, &nrhs, dl.begin(), d.begin(), du.begin(), x.begin(), &n, &info);
  
  if (info != 0) {
    // Some sort of problem
    throw BoutException("Problem in LAPACK ZGTSV routine\n");
  }

  // Copy result back
  for (int i = 0; i < n; i++) {
    u[i] = dcomplex(x[i].r, x[i].i);
  }

  return 0;
}

/* Real tridiagonal solver
 * 
 * Returns true on success
 */
bool tridag(const BoutReal *a, const BoutReal *b, const BoutReal *c, const BoutReal *r, BoutReal *u, int n) {

  // Lapack routines overwrite their inputs, so need to copy
  Array<BoutReal> dl(n), d(n), du(n), x(n);

  for (int i = 0; i < n; i++) {
    // Diagonal
    d[i] = b[i];

    // Off-diagonal terms
    if (i != (n - 1)) {
      dl[i] = a[i + 1];

      du[i] = c[i];
    }

    x[i] = r[i];
  }

  /* LAPACK DGTSV routine.

     n - Size of the array
     nrhs - Number of RHS vectors to solve
     dl - lower band
     d - diagonal values
     du - upper band
     x  - input and output values
     info - output status
  */
  int nrhs = 1;
  int info;
  dgtsv_(&n, &nrhs, dl.begin(), d.begin(), du.begin(), x.begin(), &n, &info);
  
  if (info != 0) {
    // Some sort of problem
    throw BoutException("Problem in LAPACK DGTSV routine\n");
  }

  // Copy result back
  for (int i = 0; i < n; i++) {
    u[i] = x[i];
  }

  return true;
}

/* Real cyclic tridiagonal solver
 * 
 * Uses Sherman-Morrison formula
 */
void cyclic_tridag(BoutReal *a, BoutReal *b, BoutReal *c, BoutReal *r, BoutReal *x, int n) {
  if (n <= 2) {
    throw BoutException("n too small in cyclic_tridag");
  }

  Array<BoutReal> u(n), z(n);

  BoutReal gamma = -b[0];

  // Save original values of b (restore after)
  BoutReal b0 = b[0];
  BoutReal bn = b[n - 1];

  // Modify b
  b[0] = b[0] - gamma;
  b[n - 1] = b[n - 1] - c[n - 1] * a[0] / gamma;

  // Solve tridiagonal system Ax=r
  if (!tridag(a, b, c, r, x, n))
    throw BoutException("ERROR: first tridag call failed in cyclic_tridag\n");

  u[0] = gamma;
  u[n - 1] = c[n - 1];
  for (int i = 1; i < (n - 1); i++) {
    u[i] = 0.;
  }

  // Solve Az = u
  if (!tridag(a, b, c, u.begin(), z.begin(), n)) {
    throw BoutException("ERROR: second tridag call failed in cyclic_tridag\n");
  }
    
  BoutReal fact =
      (x[0] + a[0] * x[n - 1] / gamma) / (1.0 + z[0] + a[0] * z[n - 1] / gamma);

  for (int i = 0; i < n; i++) {
    x[i] -= fact * z[i];
  }

  // Restore coefficients
  b[0] = b0;
  b[n - 1] = bn;
}

/*! Complex band solver using ZGBSV
 * 
 * n      Size of the matrix (number of equations)
 * kl     Number of subdiagonals
 * ku     Number of superdiagonals
 * nrhs   Number of RHS vectors to solve
 * ab     Array values (2D array)
 * ldab   Leading dim. size of ab = 2*KL+KU+1
 * ipiv   output integer array containing pivot permutation
 * b      RHS vectors, and solution
 * ldb    length of b (= n)
 * info   output status
 *
 * This code is about 25% faster than NR
 * Timing for 260 points (s): NR: 5.698204e-05, LAPACK: 4.410744e-05
 */
void cband_solve(dcomplex **a, int n, int m1, int m2, dcomplex *b) {
  int nrhs = 1;
  int kl = m1;
  int ku = m2;
  int ldab = 2 * kl + ku + 1;
  int ldb = n;
  
  Array<fcmplx> AB(ldab * n);
  Array<int> ipiv(n);
  Array<fcmplx> x(n);

  // Copy RHS data
  for (int i = 0; i < n; i++) {
    x[i].r = b[i].real();
    x[i].i = b[i].imag();
  }
  
  // Put matrix elements into AB(ldab, N) (FORTRAN -> loops over ldab fastest)
  // A is organised into rows, but AB is in columns. First kl not set
  for (int j = 0; j < n; j++) {
    for (int i = 0; i <= (ku + kl); i++) {
      // AB(kl + i, j) = A[j - ku + i][kl+ku - i]

      if (((j - ku + i) >= 0) && ((j - ku + i) < n)) {
        AB[j * ldab + kl + i].r = a[j - ku + i][kl + ku - i].real();
        AB[j * ldab + kl + i].i = a[j - ku + i][kl + ku - i].imag();
      }
    }
  }

  int info;
  zgbsv_(&n, &kl, &ku, &nrhs, AB.begin(), &ldab, ipiv.begin(), x.begin(), &ldb, &info);

  // Copy result back
  for (int i = 0; i < n; i++) {
    b[i] = dcomplex(x[i].r, x[i].i);
  }
}

#else
///////////////////////////////////////////////////////////////////////
// No LAPACK used. Instead fall back on Numerical Recipes routine.
//
// NOTE: THESE MUST BE REMOVED FOR PUBLIC RELEASE
///////////////////////////////////////////////////////////////////////

#include <bout/utils.hxx>

/// Tri-diagonal complex matrix inversion (from Numerical Recipes)
int tridag(const dcomplex *a, const dcomplex *b, const dcomplex *c, const dcomplex *r,
           dcomplex *u, int n) {
  dcomplex bet;
  Array<dcomplex> gam(n);

  if (b[0] == 0.0) {
    throw BoutException("Tridag: Rewrite equations\n");
  }

  bet = b[0];
  u[0] = r[0] / bet;

  for (int j = 1; j < n; j++) {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0) {
      throw BoutException("Tridag: Zero pivot\n");
    }
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }

  for (int j = n - 2; j >= 0; j--) {
    u[j] = u[j] - gam[j + 1] * u[j + 1];
  }

  return 0;
}

/// Tri-diagonal matrix inversion (BoutReal)
bool tridag(const BoutReal *a, const BoutReal *b, const BoutReal *c, const BoutReal *r, BoutReal *x, int n) {
  int j;  
  
  BoutReal bet;
  Array<BoutReal> gam(n);
  
  if(b[0] == 0.0) {
    throw BoutException("Tridag: Rewrite equations\n");
  }
  
  bet = b[0];
  x[0] = r[0] / bet;
  
  for(j=1;j<n;j++) {
    gam[j] = c[j-1]/bet;
    bet = b[j]-a[j]*gam[j];
    if(bet == 0.0) {
      throw BoutException("Tridag: Zero pivot\n");
    }
    x[j] = (r[j]-a[j]*x[j-1])/bet;
  }
  
  for(j=n-2;j>=0;j--) {
    x[j] = x[j]-gam[j+1]*x[j+1];
  }
  
  return true;
}

/// Solve a cyclic tridiagonal matrix
void cyclic_tridag(BoutReal *a, BoutReal *b, BoutReal *c, BoutReal *r, BoutReal *x, int n) {
  if (n <= 2)
    throw BoutException("ERROR: n too small in cyclic_tridag(BoutReal)");
  
  Array<BoutReal> u(n), z(n);
  
  BoutReal gamma = -b[0];
  
  // Save original values of b (restore after)
  BoutReal b0 = b[0];
  BoutReal bn = b[n-1];
  
  // Modify b
  b[0] = b[0] - gamma;
  b[n-1] = b[n-1] - c[n-1]*a[0]/gamma;
  
  // Solve tridiagonal system Ax=r
  if(!tridag(a, b, c, r, x, n))
    throw BoutException("First tridag call failed in cyclic_tridag(BoutReal)");
  
  u[0] = gamma;
  u[n-1] = c[n-1];
  for(int i=1;i<(n-1);i++)
    u[i] = 0.;
  
  // Solve Az = u
  if(!tridag(a, b, c, u.begin(), z.begin(), n))
    throw BoutException("ERROR: second tridag call failed in cyclic_tridag(BoutReal)\n");
  
  BoutReal fact = (x[0] + a[0]*x[n-1]/gamma) / // v.x / (1 + v.z)
    (1.0 + z[0] + a[0]*z[n-1]/gamma); 
  
  for(int i=0;i<n;i++)
    x[i] -= fact*z[i];
  
  // Restore coefficients
  b[0] = b0;
  b[n-1] = bn;
}


const BoutReal TINY = 1.0e-20;

void cbandec(dcomplex **a, unsigned long n, unsigned int m1, unsigned int m2,
	     dcomplex **al, unsigned long indx[], dcomplex *d) {
  unsigned long i,j,k,l;
  unsigned int mm;
  dcomplex dum;

  mm=m1+m2+1;
  l=m1;

  for (i=0;i<m1;i++) {
    for (j=m1-i;j<mm;j++) a[i][j-l]=a[i][j];
    l--;
    for (j=mm-l-1;j<mm;j++) a[i][j]=0.0;
  }
  *d = 1.0;
  l=m1;
  for (k=1;k<=n;k++) {
    dum=a[k-1][0];
    i=k;
    if (l < n) l++;
    for (j=k+1;j<=l;j++) {
      if (abs(a[j-1][0]) > abs(dum)) {
	dum=a[j-1][0];
	i=j;
      }
    }
    indx[k-1]=i;
    if (dum == 0.0) a[k-1][0]=TINY;
    if (i != k) {
      *d = -(*d);
      for (j=0;j<mm;j++) swap(a[k-1][j],a[i-1][j]);
    }
    for (i=k;i<l;i++) {
      dum=a[i][0]/a[k-1][0];
      al[k-1][i-k]=dum;
      for (j=1;j<mm;j++) a[i][j-1]=a[i][j]-dum*a[k-1][j];
      a[i][mm-1]=0.0;
    }
  }
}

void cbanbks(dcomplex **a, unsigned long n, unsigned int m1, unsigned int m2,
	     dcomplex **al, unsigned long indx[], dcomplex b[]) {
  unsigned long i,k,l;
  unsigned int mm;
  dcomplex dum;
  
  mm=m1+m2+1;
  l=m1;
  for (k=1;k<=n;k++) {
    i=indx[k-1];
    if (i != k) swap(b[k-1],b[i-1]);
    if (l < n) l++;
    for (i=k+1;i<=l;i++) b[i-1] -= al[k-1][i-k-1]*b[k-1];
  }
  l=1;
  for (i=n;i>=1;i--) {
    dum=b[i-1];
    for (k=2;k<=l;k++) dum -= a[i-1][k-1]*b[k+i-2];
    b[i-1]=dum/a[i-1][0];
    if (l < mm) l++;
  }
}

void cband_solve(dcomplex **a, int n, int m1, int m2, dcomplex *b) {
  static dcomplex **al;
  static int an = 0, am1 = 0; // Allocated sizes
  dcomplex d;
  
  if(an < n) {
    if(an != 0) {
      free_matrix(al);
    }
    al = matrix<dcomplex>(n, m1); //Never freed
    an = n;
    am1 = m1;
  }
  
  if(am1 < m1) {
    if(am1 != 0)
      free_matrix(al);
    al =  matrix<dcomplex>(an, m1);
    am1 = m1;
  }

  Array<unsigned long> indx(n);
  
  // LU decompose matrix
  cbandec(a, n, m1, m2, al, indx.begin(), &d);

  // Solve
  cbanbks(a, n, m1, m2, al, indx.begin(), b);
}

#endif // LAPACK

// Common functions

/// Solve a cyclic tridiagonal matrix
void cyclic_tridag(dcomplex *a, dcomplex *b, dcomplex *c, dcomplex *r, dcomplex *x, int n) {
  if (n <= 2)
    throw BoutException("n too small in cyclic_tridag(dcomplex)");
  
  Array<dcomplex> u(n), z(n);
  
  dcomplex gamma = -b[0];
  
  // Save original values of b (restore after)
  dcomplex b0 = b[0];
  dcomplex bn = b[n-1];
  
  // Modify b
  b[0] = b[0] - gamma;
  b[n-1] = b[n-1] - c[n-1]*a[0]/gamma;
  
  // Solve tridiagonal system Ax=r
  if(tridag(a, b, c, r, x, n))
    throw BoutException("First tridag call failed in cyclic_tridag(dcomplex)");
  
  u[0] = gamma;
  u[n-1] = c[n-1];
  for(int i=1;i<(n-1);i++)
    u[i] = 0.;
  
  // Solve Az = u
  if(tridag(a, b, c, u.begin(), z.begin(), n))
    throw BoutException("Second tridag call failed in cyclic_tridag(dcomplex)\n");
  
  dcomplex fact = (x[0] + a[0]*x[n-1]/gamma) / // v.x / (1 + v.z)
    (1.0 + z[0] + a[0]*z[n-1]/gamma); 
  
  for(int i=0;i<n;i++)
    x[i] -= fact*z[i];
  
  // Restore coefficients
  b[0] = b0;
  b[n-1] = bn;
}
