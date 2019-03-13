/************************************************************************
 * Cyclic reduction for direct solution of a complex tridiagonal system
 *
 *
 * (b0  c0      a0)
 * (a1  b1  c1    )
 * (    a2  b2  c2)
 * (c3      a3  b3)

 ic    = (a1  b1  c1)
 alpha = a2 / b1
 ic    = (a1  (b2 - (a2/b1)*c1)  c2)
 alpha = a3 / b1
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

#ifndef __CYCLIC_REDUCE_H__
#define __CYCLIC_REDUCE_H__

#ifdef DIAGNOSE
#undef DIAGNOSE
#endif
//#define DIAGNOSE 1

#include "mpi.h"
#include "utils.hxx"
#include "msg_stack.hxx"
#include <lapack_routines.hxx>

#include "bout/assert.hxx"
#include "boutexception.hxx"

#include "output.hxx"

#include "bout/openmpwrap.hxx"

template <class T> class CyclicReduce {
public:
  CyclicReduce() = default;

  CyclicReduce(MPI_Comm c, int size) : comm(c), N(size), Nsys(0) {
    MPI_Comm_size(c, &nprocs);
    MPI_Comm_rank(c, &myproc);
  }

  /// Set parameters
  /// @param[in] c  The communicator of all processors involved in the solve
  /// @param[in] size  The number of rows on this processor
  void setup(MPI_Comm c, int size) {
    comm = c;

    int np, myp;
    MPI_Comm_size(c, &np);
    MPI_Comm_rank(c, &myp);
    if ((size != N) || (np != nprocs) || (myp != myproc))
      Nsys = 0; // Need to re-size
    N = size;
    periodic = false;
    nprocs = np;
    myproc = myp;
  }

  ~CyclicReduce() = default;

  /// Specify that the tridiagonal system is periodic
  /// By default not periodic
  void setPeriodic(bool p = true) { periodic = p; }

  void setCoefs(const Array<T> &a, const Array<T> &b, const Array<T> &c) {
    ASSERT2(a.size() == b.size());
    ASSERT2(a.size() == c.size());
    ASSERT2(a.size() == N);

    Matrix<T> aMatrix(1, N);
    Matrix<T> bMatrix(1, N);
    Matrix<T> cMatrix(1, N);

    BOUT_OMP(parallel for)    
    for (int i = 0; i < N; ++i) {
      aMatrix(0, i) = a[i];
      bMatrix(0, i) = b[i];
      cMatrix(0, i) = c[i];
    }

    setCoefs(aMatrix, bMatrix, cMatrix);
  }

  /// Set the entries in the matrix to be inverted
  ///
  /// @param[in] a   Left diagonal. Should have size [nsys][N]
  ///                where N is set in the constructor or setup
  /// @param[in] b   Diagonal values. Should have size [nsys][N]
  /// @param[in] c   Right diagonal. Should have size [nsys][N]
  void setCoefs(const Matrix<T>& a, const Matrix<T>& b, const Matrix<T>& c) {
    TRACE("CyclicReduce::setCoefs");

    int nsys = std::get<0>(a.shape());

    // Make sure correct memory arrays allocated
    allocMemory(nprocs, nsys, N);

    // Fill coefficient array
    BOUT_OMP(parallel for)
    for (int j = 0; j < Nsys; j++) {
      for (int i = 0; i < N; i++) {
        coefs(j, 4 * i) = a(j, i);
        coefs(j, 4 * i + 1) = b(j, i);
        coefs(j, 4 * i + 2) = c(j, i);
        // 4*i + 3 will contain RHS
      }
    }
  }

  /// Solve a set of tridiagonal systems
  ///
  /// @param[in] rhs Array storing Values of the rhs for a single system
  /// @param[out] x  Array storing the result for a single system
  void solve(const Array<T> &rhs, Array<T> &x) {
    ASSERT2(rhs.size() == x.size());
    ASSERT2(rhs.size() == N);

    Matrix<T> rhsMatrix(1, N);
    Matrix<T> xMatrix(1, N);

    // Copy input data into matrix
    BOUT_OMP(parallel for)
    for (int i = 0; i < N; ++i) {
      rhsMatrix(0, i) = rhs[i];
    }

    // Solve
    solve(rhsMatrix, xMatrix);

    // Copy result back into argument
    BOUT_OMP(parallel for)
    for (int i = 0; i < N; ++i) {
      x[i] = xMatrix(0, i);
    }
  };

  /// Solve a set of tridiagonal systems
  ///
  /// @param[in] rhs Matrix storing Values of the rhs for each system
  /// @param[out] x  Matrix storing the result for each system
  void solve(const Matrix<T> &rhs, Matrix<T> &x) {
    TRACE("CyclicReduce::solve");
    ASSERT2(static_cast<int>(std::get<0>(rhs.shape())) == Nsys);
    ASSERT2(static_cast<int>(std::get<0>(x.shape())) == Nsys);
    ASSERT2(static_cast<int>(std::get<1>(rhs.shape())) == N);
    ASSERT2(static_cast<int>(std::get<1>(x.shape())) == N);

    // Multiple RHS
    int nrhs = std::get<0>(rhs.shape());

    if (nrhs != Nsys)
      throw BoutException("Sorry, can't yet handle nrhs != nsys");

    // Insert RHS into coefs array. Ordered to allow efficient partitioning
    // for MPI send/receives
    BOUT_OMP(parallel for)
    for (int j = 0; j < Nsys; j++) {
      for (int i = 0; i < N; i++) {
        coefs(j, 4 * i + 3) = rhs(j, i);
      }
    }

    ///////////////////////////////////////
    // Reduce local part of the matrix to interface equations
    reduce(Nsys, N, coefs, myif);

    ///////////////////////////////////////
    // Gather all interface equations onto single processor
    // NOTE: Need to replace with divide-and-conquer at some point
    //
    // There are Nsys sets of equations to gather, and nprocs processors
    // which can be used. Each processor therefore sends interface equations
    // to several different processors for gathering
    //
    // e.g. 3 processors (nproc=3), 4 sets of equations (Nsys=4):
    //
    // PE 0: [1a 2a 3a 4a]  PE 1: [1b 2b 3b 4b]  PE 2: [1c 2c 3c 4c]
    //
    // Gathered into:
    //
    // PE 0: [1a 1b 1c]  PE 1: [2a 2b 2c]   PE 2: [3a 3b 3c]
    //       [4a 4b 4c]
    //
    // Here PE 0 would have myns=2, PE 1 and 2 would have myns=1

    int ns = Nsys / nprocs;      // Number of systems to assign to all processors
    int nsextra = Nsys % nprocs; // Number of processors with 1 extra

    MPI_Request *req = new MPI_Request[nprocs];

    if (myns > 0) {
      // Post receives from all other processors
      req[myproc] = MPI_REQUEST_NULL;
      for (int p = 0; p < nprocs; p++) { // Loop over processor
        // 2 interface equations per processor
        // myns systems to solve
        // 3 coefficients + 1 RHS value
        int len = 2 * myns * 4 * sizeof(T); // Length of data in bytes

        if (p == myproc) {
          // Just copy the data
	  BOUT_OMP(parallel for)
          for (int i = 0; i < myns; i++)
            for (int j = 0; j < 8; j++)
              ifcs(i, 8 * p + j) = myif(sys0 + i, j);
        } else {
#ifdef DIAGNOSE
          output << "Expecting to receive " << len << " from " << p << endl;
#endif
          MPI_Irecv(&recvbuffer(p, 0), len,
                    MPI_BYTE, // Just sending raw data, unknown type
                    p,        // Destination processor
                    p,        // Identifier
                    comm,     // Communicator
                    &req[p]); // Request
        }
      }
    }

    // Send data
    int s0 = 0;
    for (int p = 0; p < nprocs; p++) { // Loop over processor
      int nsp = ns;
      if (p < nsextra)
        nsp++;
      if ((p != myproc) && (nsp > 0)) {
#ifdef DIAGNOSE
        output << "Sending to " << p << endl;
        for (int i = 0; i < 8; i++)
          output << "value " << i << " : " << myif(s0, i) << endl;
#endif
        MPI_Send(&myif(s0, 0),        // Data pointer
                 8 * nsp * sizeof(T), // Number
                 MPI_BYTE,            // Type
                 p,                   // Destination
                 myproc,              // Message identifier
                 comm);               // Communicator
      }
      s0 += nsp;
    }

    if (myns > 0) {
      // Wait for data
      int p;
      do {
        MPI_Status stat;
        MPI_Waitany(nprocs, req, &p, &stat);
        if (p != MPI_UNDEFINED) {
// p is the processor number. Copy data
#ifdef DIAGNOSE
          output << "Copying received data from " << p << endl;
#endif
	  BOUT_OMP(parallel for)
          for (int i = 0; i < myns; i++)
            for (int j = 0; j < 8; j++) {
#ifdef DIAGNOSE
              output << "Value " << j << " : " << recvbuffer(p, 8 * i + j) << endl;
#endif
              ifcs(i, 8 * p + j) = recvbuffer(p, 8 * i + j);
            }
          req[p] = MPI_REQUEST_NULL;
        }
      } while (p != MPI_UNDEFINED);

      ///////////////////////////////////////
      if (nprocs > 1) {
// Reduce the interface equations to a pair of equations
#ifdef DIAGNOSE
        output << "Reducing again\n";
#endif
        reduce(myns, 2 * nprocs, ifcs, if2x2);
      } else {
        // Already just a pair of equations
        if2x2 = ifcs;
      }
      ///////////////////////////////////////
      // Solve the 2x2 system directly

      // For OpenMP, ensure that memory won't be modified inside parallel loop
      if2x2.ensureUnique();
      x1.ensureUnique();
      xn.ensureUnique();
      
      BOUT_OMP(parallel for)
      for (int i = 0; i < myns; ++i) {
        //  (a  b) (x1) = (b1)
        //  (c  d) (xn)   (bn)

        T a, b, c, d;
        a = if2x2(i, 1);
        b = if2x2(i, 2);
        c = if2x2(i, 4);
        d = if2x2(i, 5);
        if (periodic) {
          b += if2x2(i, 0);
          c += if2x2(i, 6);
        }
        T b1 = if2x2(i, 3);
        T bn = if2x2(i, 7);

        // Solve
        T det = a * d - b * c; // Determinant
        x1[i] = (d * b1 - b * bn) / det;
        xn[i] = (-c * b1 + a * bn) / det;

#ifdef DIAGNOSE
        output << "system " << i << endl;
        output << "(" << a << ", " << b << ") (" << x1[i] << ") = (" << b1 << ")\n";
        output << "(" << c << ", " << d << ") (" << xn[i] << ")   (" << bn << ")\n\n";
#endif
      }

      // Solve the interface equations
      back_solve(myns, 2 * nprocs, ifcs, x1, xn, ifx);
    }

    if (nprocs > 1) {
      ///////////////////////////////////////
      // Scatter back solution

      // Post receives
      for (int p = 0; p < nprocs; p++) { // Loop over processor
        int nsp = ns;
        if (p < nsextra)
          nsp++;
        int len = 2 * nsp * sizeof(T); // 2 values per system

        if (p == myproc) {
          // Just copy the data
	  BOUT_OMP(parallel for)	  
          for (int i = 0; i < myns; i++) {
            x1[sys0 + i] = ifx(i, 2 * p);
            xn[sys0 + i] = ifx(i, 2 * p + 1);
          }
          req[p] = MPI_REQUEST_NULL;
        } else if (nsp > 0) {
#ifdef DIAGNOSE
          output << "Expecting receive from " << p << " of size " << len << endl;
#endif
          MPI_Irecv(&recvbuffer(p, 0), len,
                    MPI_BYTE, // Just sending raw data, unknown type
                    p,        // Destination processor
                    p,        // Identifier
                    comm,     // Communicator
                    &req[p]); // Request
        } else
          req[p] = MPI_REQUEST_NULL;
      }

      if (myns > 0) {
        // Send data
        for (int p = 0; p < nprocs; p++) { // Loop over processor
          if (p != myproc) {
	    BOUT_OMP(parallel for)	    
            for (int i = 0; i < myns; i++) {
              ifp[2 * i] = ifx(i, 2 * p);
              ifp[2 * i + 1] = ifx(i, 2 * p + 1);
#ifdef DIAGNOSE
              output << "Returning: " << ifp[2 * i] << ", " << ifp[2 * i + 1] << " to "
                     << p << endl;
#endif
            }
            MPI_Send(std::begin(ifp), 2 * myns * sizeof(T), MPI_BYTE, p,
                     myproc, // Message identifier
                     comm);
          }
        }
      }

      // Wait for data
      int fromproc;
      int nsp;
      do {
        MPI_Status stat;
        MPI_Waitany(nprocs, req, &fromproc, &stat);

        if (fromproc != MPI_UNDEFINED) {
          // fromproc is the processor number. Copy data

          int s0 = fromproc * ns;
          if (fromproc > nsextra) {
            s0 += nsextra;
          } else
            s0 += fromproc;

          nsp = ns;
          if (fromproc < nsextra)
            nsp++;
	  
	  BOUT_OMP(parallel for)
          for (int i = 0; i < nsp; i++) {
            x1[s0 + i] = recvbuffer(fromproc, 2 * i);
            xn[s0 + i] = recvbuffer(fromproc, 2 * i + 1);
#ifdef DIAGNOSE
            output << "Received x1,xn[" << s0 + i << "] = " << x1[s0 + i] << ", "
                   << xn[s0 + i] << " from " << fromproc << endl;
#endif
          }
          req[fromproc] = MPI_REQUEST_NULL;
        }
      } while (fromproc != MPI_UNDEFINED);
    }

    ///////////////////////////////////////
    // Solve local equations
    back_solve(Nsys, N, coefs, x1, xn, x);
    delete[] req;
  }

private:
  MPI_Comm comm;             ///< Communicator
  int nprocs{0}, myproc{-1}; ///< Number of processors and ID of my processor

  int N{0};    ///< Total size of the problem
  int Nsys{0}; ///< Number of independent systems to solve
  int myns;    ///< Number of systems for interface solve on this processor
  int sys0;    ///< Starting system index for interface solve

  bool periodic{false}; ///< Is the domain periodic?

  Matrix<T> coefs; ///< Starting coefficients, rhs [Nsys, {3*coef,rhs}*N]
  Matrix<T> myif;  ///< Interface equations for this processor

  Matrix<T> recvbuffer; ///< Buffer for receiving from other processors
  Matrix<T> ifcs;       ///< Coefficients for interface solve
  Matrix<T> if2x2;      ///< 2x2 interface equations on this processor
  Matrix<T> ifx;        ///< Solution of interface equations
  Array<T> ifp;         ///< Interface equations returned to processor p
  Array<T> x1, xn;      ///< Interface solutions for back-solving

  /// Allocate memory arrays
  /// @param[in] np   Number of processors
  /// @param[in] nsys  Number of independent systems to solve
  /// @param[in] n     Size of each system of equations
  void allocMemory(int np, int nsys, int n) {
    if ((nsys == Nsys) && (n == N) && (np == nprocs))
      return; // No need to allocate memory

    nprocs = np;
    Nsys = nsys;
    N = n;

    // Work out how many systems are going to be solved on this processor
    int ns = nsys / nprocs;      // Number of systems to assign to all processors
    int nsextra = nsys % nprocs; // Number of processors with 1 extra

    myns = ns;          // Number of systems to gather onto this processor
    sys0 = ns * myproc; // Starting system number
    if (myproc < nsextra) {
      myns++;
      sys0 += myproc;
    } else {
      sys0 += nsextra;
    }

    coefs.reallocate(Nsys, 4 * N);
    myif.reallocate(Nsys, 8);

    // Note: The recvbuffer is used to receive data in both stages of the solve:
    //  1. In the gather step, this processor will receive myns interface equations
    //     from each processor.
    //  2. In the scatter step, this processor receives the solved interface values
    //     from each processor. The number of systems of equations received will
    //     vary from myns to myns+1 (if myproc >= nsextra).
    // The size of the array reserved is therefore (myns+1)
    recvbuffer.reallocate(nprocs, (myns + 1) * 8);

    // Some interface systems to be solved on this processor
    // Note that the interface equations are organised by system (myns as first argument)
    // but communication buffers are organised by processor (nprocs first).
    ifcs.reallocate(myns, 2 * 4 * nprocs); // Coefficients for interface solve
    if (nprocs > 1) {
      if2x2.reallocate(myns, 2 * 4); // 2x2 interface equations on this processor
    }
    ifx.reallocate(myns, 2 * nprocs); // Solution of interface equations
    ifp.reallocate(myns * 2);         // Solution to be sent to processor p
    // Each system to be solved on this processor has two interface equations from each
    // processor

    x1.reallocate(Nsys);
    xn.reallocate(Nsys);
  }

  /// Calculate interface equations
  ///
  /// This reduces ns separate systems of equations, each consisting
  /// of nloc rows on this processor, to two interface rows for each system,
  /// which are stored in ifc.
  ///
  /// (a1 b1 c1                  )
  /// (   a2 b2 c2               )       (A1 B1 C1   )
  /// (      a3 b3 c3            )   =>  (   A2 B2 C2)
  /// (              ...         )
  /// (                  an bn cn)
  void reduce(int ns, int nloc, Matrix<T> &co, Matrix<T> &ifc) {
#ifdef DIAGNOSE
    if (nloc < 2)
      throw BoutException("CyclicReduce::reduce nloc < 2");
#endif
    
    BOUT_OMP(parallel for)
    for (int j = 0; j < ns; j++) {
      // Calculate upper interface equation

      // v_l <- v_(k+N-2)
      // b_u <- b_{k+N-2}
      for (int i = 0; i < 4; i++) {
        ifc(j, i) = co(j, 4 * (nloc - 2) + i);
      }

      for (int i = nloc - 3; i >= 0; i--) {
        // Check for zero pivot
        if (std::abs(ifc(j, 1)) < 1e-10)
          throw BoutException("Zero pivot in CyclicReduce::reduce");

        // beta <- v_{i,i+1} / v_u,i
        T beta = co(j, 4 * i + 2) / ifc(j, 1);

        // v_u <- v_i - beta * v_u
        ifc(j, 1) = co(j, 4 * i + 1) - beta * ifc(j, 0);
        ifc(j, 0) = co(j, 4 * i);
        ifc(j, 2) *= -beta;
        // ic columns  {i-1, i, N-1}

        // b_u <- b_i - beta*b_u
        ifc(j, 3) = co(j, 4 * i + 3) - beta * ifc(j, 3);
      }

      // Calculate lower interface equation
      // Uses next 4 ifc values for this system so have +4 in indexinf

      // v_l <- v_(k+1)
      // b_l <- b_{k+1}
      for (int i = 0; i < 4; i++)
        ifc(j, 4 + i) = co(j, 4 + i);

      for (int i = 2; i < nloc; i++) {

        if (std::abs(ifc(j, 4 + 1)) < 1e-10)
          throw BoutException("Zero pivot in CyclicReduce::reduce");

        // alpha <- v_{i,i-1} / v_l,i-1
        T alpha = co(j, 4 * i) / ifc(j, 4 + 1);

        // v_l <- v_i - alpha*v_l
        ifc(j, 4 + 0) *= -alpha;
        ifc(j, 4 + 1) = co(j, 4 * i + 1) - alpha * ifc(j, 4 + 2);
        ifc(j, 4 + 2) = co(j, 4 * i + 2);
        // columns of ic are {0, i, i + 1}

        // b_l <- b_{k + i} - alpha*b_l
        ifc(j, 4 + 3) = co(j, 4 * i + 3) - alpha * ifc(j, 4 + 3);
      }

#ifdef DIAGNOSE
      output << "Lower: " << ifc(j, 4 + 0) << ", " << ifc(j, 4 + 1) << ", "
             << ifc(j, 4 + 2) << " : " << ifc(j, 4 + 3) << endl;
      output << "Upper: " << ifc(j, 0) << ", " << ifc(j, 1) << ", " << ifc(j, 2) << " : "
             << ifc(j, 3) << endl;
#endif
    }

    // Lower system couples {0, N-1, N}
    // Upper system couples {-1. 0, N-1}
  }

  /// Back-solve from x at ends (x1, xn) to obtain remaining values
  /// Coefficients ordered [ns, nloc*(a,b,c,r)]
  void back_solve(int ns, int nloc, const Matrix<T>& co, const Array<T>& x1,
                  const Array<T>& xn, Matrix<T>& xa) {

    xa.ensureUnique(); // Going to be modified, so call this outside parallel region
    
    // Tridiagonal system, solve using serial Thomas algorithm
    // xa -- Result for each system
    // co -- Coefficients & rhs for each system
    BOUT_OMP(parallel for)
    for (int i = 0; i < ns; i++) { // Loop over systems
      Array<T> gam(nloc); // Thread-local array
      T bet = 1.0;
      xa(i, 0) = x1[i]; // Already know the first
      gam[1] = 0.;
      for (int j = 1; j < nloc - 1; j++) {
        bet = co(i, 4 * j + 1) - co(i, 4 * j) * gam[j]; // bet = b[j]-a[j]*gam[j]
        xa(i, j) = (co(i, 4 * j + 3) - co(i, 4 * j) * xa(i, j - 1)) /
                   bet;                      // x[j] = (r[j]-a[j]*x[j-1])/bet;
        gam[j + 1] = co(i, 4 * j + 2) / bet; // gam[j+1] = c[j]/bet
      }
      xa(i, nloc - 1) = xn[i]; // Know the last value

      for (int j = nloc - 2; j > 0; j--) {
        xa(i, j) = xa(i, j) - gam[j + 1] * xa(i, j + 1);
      }
    }
  }
};

#endif // __CYCLIC_REDUCE_H__
