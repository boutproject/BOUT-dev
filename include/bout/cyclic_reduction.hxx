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
#include "bout/utils.hxx"
#include <bout/lapack_routines.hxx>
#include "bout/boutexception.hxx"

template <class T>
class CyclicReduce {
public:
  
  CyclicReduce() {
    nprocs = 0;
    myproc = -1;
    N = 0;
    Nsys = 0;
  }
  
  CyclicReduce(MPI_Comm c, int size) : comm(c), N(size), Nsys(0), periodic(false) {
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
    if((size != N) || (np != nprocs) || (myp != myproc))
      freeMemory(); // Need to re-size
    N = size;
    periodic = false;
    nprocs = np;
    myproc = myp;
  }

  ~CyclicReduce() {
    freeMemory();
  }

  /// Specify that the tridiagonal system is periodic
  /// By default not periodic
  void setPeriodic(bool p=true) {periodic=p;}
  
  /// Set up a single equation to invert
  void setCoefs(T a[], T b[], T c[]) {
    setCoefs(1, &a, &b, &c);
  }

  /// Set the entries in the matrix to be inverted
  ///
  /// @param[in] nsys   The number of independent matrices to be solved
  /// @param[in] a   Left diagonal. Should have size [nsys][N]
  ///                where N is set in the constructor or setup
  /// @param[in] b   Diagonal values. Should have size [nsys][N]
  /// @param[in] c   Right diagonal. Should have size [nsys][N]
  void setCoefs(int nsys, T **a, T **b, T **c) {
    // Make sure correct memory arrays allocated
    allocMemory(nprocs, nsys, N);

    // Fill coefficient array
    for(int j=0;j<Nsys;j++)
      for(int i=0;i<N;i++) {
        coefs[j][4*i] = a[j][i];
        coefs[j][4*i + 1] = b[j][i];
        coefs[j][4*i + 2] = c[j][i];
        // 4*i + 3 will contain RHS
      }
  }

  /// Solve a single triadiagonal system
  /// 
  void solve(T rhs[], T x[]) {
    // Solving single system
    solve(1, &rhs, &x);
  }

  /// Solve a set of tridiagonal systems
  /// 
  void solve(int nrhs, T **rhs, T **x) {
    // Multiple RHS
    
    if(nrhs != Nsys)
      throw BoutException("Sorry, can't yet handle nrhs != nsys");
    
    // Insert RHS into coefs array. Ordered to allow efficient partitioning
    // for MPI send/receives
    for(int j=0;j<Nsys;j++)
      for(int i=0;i<N;i++) {
        coefs[j][4*i + 3] = rhs[j][i];
      }

    ///////////////////////////////////////
    // Reduce local part of the matrix to interface equations
    reduce(Nsys, N, coefs, myif);
    
    ///////////////////////////////////////
    // Gather all interface equations onto single processor
    // NOTE: Need to replace with divide-and-conquer at some point
    
    int ns = Nsys / nprocs; // Number of systems to assign to all processors
    int nsextra = Nsys % nprocs;  // Number of processors with 1 extra 
    
    MPI_Request* req = new MPI_Request[nprocs];

    if(myns > 0) {
      // Post receives from all other processors
      req[myproc] = MPI_REQUEST_NULL;
      for(int p=0;p<nprocs;p++) { // Loop over processor
        // 2 interface equations per processor
        // myns systems to solve
        // 3 coefficients + 1 RHS value
        int len = 2 * myns * 4 * sizeof(T); // Length of data in bytes
        
        if(p == myproc) {
          // Just copy the data
          for(int i=0;i<myns; i++)
            for(int j=0;j<8;j++)
              ifcs[i][8*p + j] = myif[sys0+i][j];
        }else {
#ifdef DIAGNOSE
          output << "Expecting to receive " << len << " from " << p << endl;
#endif
          MPI_Irecv(recvbuffer[p], 
                    len, 
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
    for(int p=0;p<nprocs;p++) { // Loop over processor
      int nsp = ns;
      if(p < nsextra)
        nsp++;
      if((p != myproc) && (nsp > 0)) {
#ifdef DIAGNOSE
        output << "Sending to " << p << endl;
        for(int i=0;i<8;i++)
          output << "value " << i << " : " << myif[s0][i] << endl;
#endif
        MPI_Send(myif[s0],        // Data pointer
                 8*nsp*sizeof(T), // Number
                 MPI_BYTE,        // Type
		 p,               // Destination
                 myproc,          // Message identifier
                 comm);           // Communicator
      }
      s0 += nsp;
    }
    
    if(myns > 0) {
      // Wait for data
      int p;
      do {
        MPI_Status stat;
        MPI_Waitany(nprocs, req, &p, &stat);
        if(p != MPI_UNDEFINED) {
          // p is the processor number. Copy data
#ifdef DIAGNOSE
          output << "Copying received data from " << p << endl;
#endif
          for(int i=0;i<myns; i++)
            for(int j=0;j<8;j++) {
#ifdef DIAGNOSE
              output << "Value " << j << " : " << recvbuffer[p][8*i + j] << endl;
#endif
              ifcs[i][8*p + j] = recvbuffer[p][8*i + j];
            }
          req[p] = MPI_REQUEST_NULL;
        }
      }while(p != MPI_UNDEFINED);

      ///////////////////////////////////////
      if(nprocs > 1) {
	// Reduce the interface equations to a pair of equations
#ifdef DIAGNOSE
	output << "Reducing again\n";
#endif
	reduce(myns, 2*nprocs, ifcs, if2x2);
      }else {
	// Already just a pair of equations
	if2x2 = ifcs;
      }
      ///////////////////////////////////////
      // Solve the 2x2 system directly
        
      for(int i=0;i<myns;i++) {
	//  (a  b) (x1) = (b1)
	//  (c  d) (xn)   (bn)
          
	T a, b, c, d;
	a = if2x2[i][1];
	b = if2x2[i][2];
	c = if2x2[i][4];
	d = if2x2[i][5];
	if(periodic) {
	  b += if2x2[i][0];
	  c += if2x2[i][6];
	}
	T b1 = if2x2[i][3];
	T bn = if2x2[i][7];
          
	// Solve
	T det = a*d - b*c; // Determinant
	x1[i] = (d*b1 - b*bn) / det;
	xn[i] = (-c*b1 + a*bn) / det;
          
#ifdef DIAGNOSE    
	output << "system " << i << endl;
	output << "(" << a << ", " << b << ") ("<<x1[i]<<") = (" << b1 << ")\n";
	output << "(" << c << ", " << d << ") ("<<xn[i]<<")   (" << bn << ")\n\n";
#endif
      }
      
      // Solve the interface equations
      back_solve(myns, 2*nprocs, ifcs, x1, xn, ifx);
    }
    
    if(nprocs > 1) { 
      ///////////////////////////////////////
      // Scatter back solution
      
      // Post receives
      for(int p=0;p<nprocs;p++) { // Loop over processor
        int nsp = ns;
	if(p < nsextra)
	  nsp++;
	int len = 2 * nsp * sizeof(T); // 2 values per system
        
	if(p == myproc) {
	  // Just copy the data
	  for(int i=0;i<myns; i++) {
	    x1[sys0+i] = ifx[i][2*p];
	    xn[sys0+i] = ifx[i][2*p+1];
	  }
          req[p] = MPI_REQUEST_NULL;
	}else if(nsp > 0) {
#ifdef DIAGNOSE
          output << "Expecting receive from " << p << " of size " << len << endl;
#endif
	  MPI_Irecv(recvbuffer[p],
		    len,
		    MPI_BYTE, // Just sending raw data, unknown type
		    p,        // Destination processor
		    p,        // Identifier
		    comm,     // Communicator
		    &req[p]); // Request
	}else
          req[p] = MPI_REQUEST_NULL;
      }
      
      if(myns > 0) {
        // Send data
        for(int p=0;p<nprocs;p++) { // Loop over processor
          if(p != myproc) {
            for(int i=0;i<myns;i++) {
              ifp[2*i]   = ifx[i][2*p];
              ifp[2*i+1] = ifx[i][2*p+1];
#ifdef DIAGNOSE
              output << "Returning: " << ifp[2*i] 
                     << ", " << ifp[2*i+1] << " to " << p << endl;
#endif
            }
            MPI_Send(ifp,
                     2*myns*sizeof(T),
                     MPI_BYTE,
                     p,
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
	if(fromproc != MPI_UNDEFINED) {
	  // fromproc is the processor number. Copy data
	  
	  int s0 = fromproc*ns;
	  if(fromproc > nsextra) {
	    s0 += nsextra;
	  }else
	    s0 += fromproc;
	    
	  nsp = ns;
	  if(fromproc < nsextra)
	    nsp++;
	  
	  for(int i=0;i<nsp; i++) {
	    x1[s0+i] = recvbuffer[fromproc][2*i];
	    xn[s0+i] = recvbuffer[fromproc][2*i+1];
#ifdef DIAGNOSE
            output << "Received x1,xn[" << s0+i << "] = " << x1[s0+i] << ", "
		   << xn[s0+i] << " from " << fromproc << endl;
#endif
	  }
	  req[fromproc] = MPI_REQUEST_NULL;
	}
      }while(fromproc != MPI_UNDEFINED);
    }
    
    ///////////////////////////////////////
    // Solve local equations
    back_solve(Nsys, N, coefs, x1, xn, x);
    delete[] req;
  }
  
private:
  MPI_Comm comm;      ///< Communicator
  int nprocs, myproc; ///< Number of processors and ID of my processor
  
  int N;         ///< Total size of the problem
  int Nsys;      ///< Number of independent systems to solve
  int myns;      ///< Number of systems for interface solve on this processor
  int sys0;      ///< Starting system index for interface solve
  
  bool periodic; ///< Is the domain periodic?

  T **coefs;  ///< Starting coefficients, rhs [Nsys, {3*coef,rhs}*N]
  T **myif;   ///< Interface equations for this processor
  
  T **recvbuffer; ///< Buffer for receiving from other processors
  T **ifcs;   ///< Coefficients for interface solve
  T **if2x2;  ///< 2x2 interface equations on this processor
  T **ifx;    ///< Solution of interface equations
  T *ifp;     ///< Interface equations returned to processor p
  T *x1, *xn; ///< Interface solutions for back-solving

  /// Allocate memory arrays
  /// @param[in[ np   Number of processors
  /// @param[in] nsys  Number of independent systems to solve
  /// @param[in] n     Size of each system of equations
  void allocMemory(int np, int nsys, int n) {
    if( (nsys == Nsys) && (n == N) && (np == nprocs))
      return; // No need to allocate memory
    
    if(Nsys > 0)
      freeMemory(); // Free existing memory
    
    nprocs = np;
    Nsys   = nsys;
    N      = n;

    // Work out how many systems are going to be solved on this processor
    int ns = nsys / nprocs; // Number of systems to assign to all processors
    int nsextra = nsys % nprocs;  // Number of processors with 1 extra 
    
    myns = ns; // Number of systems to gather onto this processor
    sys0 = ns * myproc; // Starting system number
    if(myproc < nsextra) {
      myns++;
      sys0 += myproc;
    }else
      sys0 += nsextra;
    
    int my = myns;
    if(my == 0)
      my = 0;

    coefs = matrix<T>(Nsys, 4*N);
      
    myif = matrix<T>(Nsys, 8);
    
    recvbuffer = matrix<T>(nprocs, my*8); // Buffer for receiving from other processors
    ifcs = matrix<T>(my, 2*4*nprocs);     // Coefficients for interface solve
    if(nprocs > 1)
      if2x2 = matrix<T>(my, 2*4);         // 2x2 interface equations on this processor
    ifx  = matrix<T>(my, 2*nprocs);       // Solution of interface equations
    ifp = new T[my*2];     // Solution to be sent to processor p
    x1 = new T[Nsys];
    xn = new T[Nsys];
    
  }

  /// Free all memory arrays allocated by allocMemory()
  void freeMemory() {
    if(Nsys == 0)
      return;
    
    // Free all working memory
    free_matrix(coefs);
    free_matrix(myif);
    free_matrix(recvbuffer);
    free_matrix(ifcs);
    if(nprocs > 1)
      free_matrix(if2x2);
    free_matrix(ifx);
    delete[] ifp;
    delete[] x1;
    delete[] xn;
    
    N = Nsys = 0;
  }

  /// Calculate interface equations
  void reduce(int ns, int nloc, T **co, T **ifc) {
#ifdef DIAGNOSE
    if(nloc < 2)
      throw BoutException("CyclicReduce::reduce nloc < 2");
#endif
    for(int j=0;j<ns;j++) {
      T *c = co[j];    // Coefficients for system j
      T *ic = ifc[j];  // Interface equations for system j
        
      // Calculate upper interface equation
      
      // v_l <- v_(k+N-2)
      // b_u <- b_{k+N-2}
      for(int i=0;i<4;i++) {
        ic[i] = c[4*(nloc-2)+i];
      }
      
      for(int i=nloc-3;i>=0;i--) {
	// Check for zero pivot
	if(abs(ic[1]) < 1e-10)
	  throw BoutException("Zero pivot in CyclicReduce::reduce");
	
        // beta <- v_{i,i+1} / v_u,i
        T beta = c[4*i+2] / ic[1];
          
        // v_u <- v_i - beta * v_u
        ic[1] = c[4*i + 1] - beta * ic[0];
        ic[0] = c[4*i];
	ic[2] *= -beta;
        // ic columns  {i-1, i, N-1}
        
        // b_u <- b_i - beta*b_u
        ic[3] = c[4*i + 3] - beta*ic[3];
      }
      
      ic += 4;
      
      // Calculate lower interface equation
      
      // v_l <- v_(k+1)
      // b_l <- b_{k+1}
      for(int i=0;i<4;i++)
        ic[i] = c[4+i];
        
      for(int i=2;i<nloc;i++) {
	
	if(abs(ic[1]) < 1e-10)
	  throw BoutException("Zero pivot in CyclicReduce::reduce");
	
        // alpha <- v_{i,i-1} / v_l,i-1
        T alpha = c[4*i] / ic[1];
          
        // v_l <- v_i - alpha*v_l
	ic[0] *= -alpha;
        ic[1] = c[4*i + 1] - alpha*ic[2];
        ic[2] = c[4*i + 2];
        // columns of ic are {0, i, i+1}
          
        // b_l <- b_{k+i} - alpha*b_l
        ic[3] = c[4*i + 3] - alpha * ic[3];   
      }
      
#ifdef DIAGNOSE
      output << "Lower: " << ic[0] << ", " << ic[1] << ", " << ic[2] << " : " << ic[3] << endl;
      output << "Upper: " << ic[0] << ", " << ic[1] << ", " << ic[2] << " : " << ic[3] << endl;
#endif
    }

    // Lower system couples {0, N-1, N}
    // Upper system couples {-1. 0, N-1}
  }
  
  /// Back-solve from x at ends (x1, xn) to obtain remaining values
  /// Coefficients ordered [ns, nloc*(a,b,c,r)]
  void back_solve(int ns, int nloc, T **co, T *x1, T *xn, T **xa) {
    // Tridiagonal system, solve using serial Thomas algorithm
    
    T *gam = new T[nloc];
    for(int i=0;i<ns;i++) { // Loop over systems
      T *c = co[i]; // Coefficients & rhs for this system
      T *x = xa[i]; // Result for this system
      T bet = 1.0;
      x[0] = x1[i]; // Already know the first 
      gam[1] = 0.;
      for(int j=1;j<nloc-1;j++) {
        bet = c[4*j+1] - c[4*j]*gam[j]; // bet = b[j]-a[j]*gam[j]
        x[j] = (c[4*j+3] - c[4*j]*x[j-1])/bet;  // x[j] = (r[j]-a[j]*x[j-1])/bet;
        gam[j+1] = c[4*j+2] / bet;    // gam[j+1] = c[j]/bet
      }
      x[nloc-1] = xn[i]; // Know the last value
      
      for(int j=nloc-2;j>0;j--) {
        x[j] = x[j]-gam[j+1]*x[j+1];
      }
    }
    delete[] gam;
  }
};

#endif // __CYCLIC_REDUCE_H__
