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


#include <globals.hxx>
#include <dcomplex.hxx>
#include <boutexception.hxx>
#include <utils.hxx>

// for dev
#include <output.hxx>
#include "boutcomm.hxx"

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

/// Use LAPACK routine ZGTSV
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

//  if( BoutComm::rank() == 0 ) {
//    std::cout << "proc " << BoutComm::rank() <<  " vn " << d[n-2].r << endl;
//  }
///
///  BoutReal S = a[2].real()/d[2].r;
///  BoutReal fac = S;
///  for (int i = 2; i < n-1; i++) {
///    fac = fac*c[i-1].real()*a[i].real()/(d[i-1].r*d[i].r);
///    S = S + fac;
///  }
///
//  if( BoutComm::rank() == 1 ) {
//    std::cout << "proc " << BoutComm::rank() <<  " S " << S << endl;
//  }

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
 */
void cband_solve(Matrix<dcomplex> &a, int n, int m1, int m2, Array<dcomplex> &b) {
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
        AB[j * ldab + kl + i].r = a(j - ku + i, kl + ku - i).real();
        AB[j * ldab + kl + i].i = a(j - ku + i, kl + ku - i).imag();
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
// No LAPACK available. Routines throw exceptions

/// Tri-diagonal complex matrix inversion
int tridag(const dcomplex*, const dcomplex*, const dcomplex*, const dcomplex*, dcomplex*, int) {
  throw BoutException("complex tridag function not available. Compile BOUT++ with Lapack support.");
}

/// Tri-diagonal matrix inversion (BoutReal)
bool tridag(const BoutReal*, const BoutReal*, const BoutReal*, const BoutReal*, BoutReal*, int) {
  throw BoutException("tridag function not available. Compile BOUT++ with Lapack support.");
}

/// Solve a cyclic tridiagonal matrix
void cyclic_tridag(BoutReal*, BoutReal*, BoutReal*, BoutReal*, BoutReal*, int) {
  throw BoutException("cyclic_tridag function not available. Compile BOUT++ with Lapack support.");
}

void cband_solve(Matrix<dcomplex>&, int, int, int, Array<dcomplex>&) {
  throw BoutException("cband_solve function not available. Compile BOUT++ with Lapack support.");
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
