/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Geometrical Multigrid Solver
 *
 * Equation solved is:
 *  d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2016 K.S. Kang
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

#include "bout/build_config.hxx"
#include "multigrid_laplace.hxx"

#if not BOUT_USE_METRIC_3D

#include <bout/mesh.hxx>
#include <msg_stack.hxx>
#include <bout/openmpwrap.hxx>

#ifdef _OPENMP
#include <omp.h>
#endif

BoutReal soltime=0.0,settime=0.0;

LaplaceMultigrid::LaplaceMultigrid(Options *opt, const CELL_LOC loc, Mesh *mesh_in) :
  Laplacian(opt, loc, mesh_in),
  A(0.0), C1(1.0), C2(1.0), D(1.0) {

  TRACE("LaplaceMultigrid::LaplaceMultigrid(Options *opt)");
  
  // periodic x-direction not handled: see MultigridAlg::communications
  ASSERT1(!localmesh->periodicX);

  A.setLocation(location);
  C1.setLocation(location);
  C2.setLocation(location);
  D.setLocation(location);

  // Get Options in Laplace Section
  if (!opt) opts = Options::getRoot()->getSection("laplace");
  else opts=opt;
  opts->get("multigridlevel",mglevel,100,true);
  opts->get("rtol",rtol,pow(10.0,-8),true);
  opts->get("atol",atol,pow(10.0,-20),true);
  opts->get("dtol",dtol,pow(10.0,5),true);
  opts->get("smtype",mgsm,1,true);
#ifdef _OPENMP
  if (mgsm != 0 && omp_get_max_threads()>1) {
    output_warn << "WARNING: in multigrid Laplace solver, for smtype!=0 the smoothing cannot be parallelised with OpenMP threads."<<endl
                << "         Consider using smtype=0 instead when using OpenMP threads."<<endl;
  }
#endif
  opts->get("jacomega",omega,0.8,true);
  opts->get("solvertype",mgplag,1,true);
  opts->get("cftype",cftype,0,true);
  opts->get("mergempi",mgmpi,63,true);
  opts->get("checking",pcheck,0,true);
  mgcount = 0;

  // Initialize, allocate memory, etc.
  comms_tagbase = 385; // Some random number
  
  int implemented_global_flags = INVERT_START_NEW;
  if ( global_flags & ~implemented_global_flags ) {
    throw BoutException("Attempted to set Laplacian inversion flag that is not implemented in LaplaceMultigrid.");
  }
  int implemented_boundary_flags = INVERT_AC_GRAD + INVERT_SET + INVERT_DC_GRAD; // INVERT_DC_GRAD does not actually do anything, but harmless to set while comparing to Fourier solver with Neumann boundary conditions
  if ( inner_boundary_flags & ~implemented_boundary_flags ) {
    throw BoutException("Attempted to set Laplacian inner boundary inversion flag that is not implemented in LaplaceMultigrid.");
  }
  if ( outer_boundary_flags & ~implemented_boundary_flags ) {
    throw BoutException("Attempted to set Laplacian outer boundary inversion flag that is not implemented in LaplaceMultigrid.");
  }
  
  commX = localmesh->getXcomm();
  
  Nx_local = localmesh->xend - localmesh->xstart + 1; // excluding guard cells
  Nx_global = localmesh->GlobalNx - 2*localmesh->xstart; // excluding guard cells
  
  if (mgcount == 0) {
    output <<"Nx="<<Nx_global<<"("<<Nx_local<<")"<<endl;
  }
  Nz_global = localmesh->GlobalNz;
  Nz_local = Nz_global; // No parallelization in z-direction (for now)
  // 
  //else {
  //  Nz_local = localmesh->zend - localmesh->zstart + 1; // excluding guard cells
  //  Nz_global = localmesh->GlobalNz - 2*localmesh->zstart; // excluding guard cells
  // }
  if (mgcount==0) {
    output <<"Nz="<<Nz_global<<"("<<Nz_local<<")"<<endl;
  }

  // Compute available levels along x-direction
  if (mglevel >1) {
    int nn = Nx_global;
    for (int n = mglevel;n > 1; n--) {
      if ( nn%2 != 0 )  {
	output<<"Size of global x-domain is not a multiple of 2^"<<mglevel-1<<" mglevel is changed to "<<mglevel-n+1<<endl;
        mglevel = mglevel - n + 1;
        break;
      }
      nn = nn/2;
    }
    // ... and check the same for z-direction
    nn = Nz_global;
    for (int n = mglevel;n > 1; n--) {
      if ( nn%2 != 0 )  {
	output<<"Size of global z-domain is not a multiple of 2^ "<<mglevel-1<<" mglevel is changed to "<<mglevel-n+1<<endl;
        mglevel = mglevel - n + 1;
        break;
      }
      nn = nn/2;
    }
  }
  else mglevel = 1;

  // Compute available levels on each processor along x-direction
  // aclevel is the number of levels that can be used in parallel, i.e. set by
  // the grid size on a single processor
  // If the number of levels is higher than aclevel, then the grid is collected
  // to a single processor, and a new multigrid solver (called sMG) is created
  // to run in serial to compute the coarsest (mglevel-aclevel) levels
  int aclevel,adlevel;
  if (mglevel >1) {
    int nn = Nx_local;
    aclevel = mglevel;
    for (int n = aclevel;n > 1; n--) {
      if ( nn%2 != 0 )  {
	output<<"Size of local x-domain is not a multiple of 2^"<<aclevel<<" aclevel is changed to "<<aclevel-n+1<<endl;
        aclevel = aclevel - n + 1;
        break;
      }
      nn = nn/2;
    }
    // ... and check the same for z-direction
    nn = Nz_local;
    for (int n = aclevel;n > 1; n--) {
      if ( nn%2 != 0 )  {
	output<<"Size of local z-domain is not a multiple of 2^ "<<aclevel<<" aclevel is changed to "<<aclevel-n<<endl;
        aclevel = aclevel - n + 1;
        break;
      }
      nn = nn/2;
    }
  }
  else aclevel = 1;
  adlevel = mglevel - aclevel;

  kMG = bout::utils::make_unique<Multigrid1DP>(aclevel, Nx_local, Nz_local, Nx_global,
                                               adlevel, mgmpi, commX, pcheck);
  kMG->mgplag = mgplag;
  kMG->mgsm = mgsm; 
  kMG->cftype = cftype;
  kMG->rtol = rtol;
  kMG->atol = atol;
  kMG->dtol = dtol;
  kMG->omega = omega;
  kMG->setValueS();

  // Set up Multigrid Cycle

  x.reallocate((Nx_local + 2) * (Nz_local + 2));
  b.reallocate((Nx_local + 2) * (Nz_local + 2));

  if (mgcount == 0) {  
    output<<" Smoothing type is ";
    if (mgsm == 0) {
      output<<"Jacobi smoother";
      output<<"with omega = "<<omega<<endl;
    }
    else if(mgsm ==1) output<<" Gauss-Seidel smoother"<<endl;
    else throw BoutException("Undefined smoother");
    output<<"Solver type is ";
    if (mglevel == 1) output<<"PGMRES with simple Preconditioner"<<endl;
    else if(mgplag == 1) output<<"PGMRES with multigrid Preconditioner"<<endl;
    else output<<"Multigrid solver with merging "<<mgmpi<<endl;
#ifdef OPENMP
BOUT_OMP(parallel)
BOUT_OMP(master)
    {
      output<<"Num threads = "<<omp_get_num_threads()<<endl;
    } 
#endif
  }  
}

FieldPerp LaplaceMultigrid::solve(const FieldPerp& b_in, const FieldPerp& x0) {

  TRACE("LaplaceMultigrid::solve(const FieldPerp, const FieldPerp)");

  ASSERT1(localmesh == b_in.getMesh() && localmesh == x0.getMesh());
  ASSERT1(b_in.getLocation() == location);
  ASSERT1(x0.getLocation() == location);

  checkData(b_in);
  checkData(x0);
  ASSERT3(b_in.getIndex() == x0.getIndex());

  BoutReal t0,t1;
  
  yindex = b_in.getIndex();
  int level = kMG->mglevel-1;
  int lzz = kMG->lnz[level];
  int lz2 = lzz+2;
  int lxx = kMG->lnx[level];

  if ( global_flags & INVERT_START_NEW ) {
    // set initial guess to zero
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for collapse(2))
    for (int i=1; i<lxx+1; i++) {
      for (int k=1; k<lzz+1; k++) {
        x[i*lz2+k] = 0.;
      }
    }
  } else {
    // Read initial guess into local array, ignoring guard cells
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for collapse(2))
    for (int i=1; i<lxx+1; i++) {
      for (int k=1; k<lzz+1; k++) {
        int i2 = i-1+localmesh->xstart;
        int k2 = k-1;
        x[i*lz2+k] = x0[i2][k2];
      }
    }
  }
  
  // Read RHS into local array
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for collapse(2))
  for (int i=1; i<lxx+1; i++) {
    for (int k=1; k<lzz+1; k++) {
      int i2 = i-1+localmesh->xstart;
      int k2 = k-1;
      b[i*lz2+k] = b_in(i2, k2);
    }
  }
  
  if (localmesh->firstX()) {
    if ( inner_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify gradient to set at inner boundary
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
	  x[k] = -x0(localmesh->xstart-1, k2)*sqrt(coords->g_11(localmesh->xstart, yindex))*coords->dx(localmesh->xstart, yindex); 
        }
      } else {
        // zero gradient inner boundary condition
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          // set inner guard cells
          x[k] = 0.0;
        }
      }
    } else {
      // Dirichlet boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify value to set at inner boundary
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          x[k] = 2.*x0(localmesh->xstart-1, k2); 
        // this is the value to set at the inner boundary
        }
      }
      else {
        // zero value inner boundary condition
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          // set inner guard cells
          x[k] = 0.;
        }
      }
    }
  }
  if (localmesh->lastX()) {
    if ( outer_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify gradient to set at outer boundary
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
        x[(lxx+1)*lz2+k] = x0(localmesh->xend+1, k2)*sqrt(coords->g_11(localmesh->xend, yindex))*coords->dx(localmesh->xend, yindex); 
        // this is the value to set the gradient to at the outer boundary
        }
      }
      else {
        // zero gradient outer boundary condition
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          // set outer guard cells
          x[(lxx+1)*lz2+k] = 0.;
        }
      }
    }
    else {
      // Dirichlet boundary condition
      if ( outer_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify value to set at outer boundary
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          x[(lxx+1)*lz2+k]=2.*x0(localmesh->xend+1, k2); 
          // this is the value to set at the outer boundary
        }
      }
      else {
        // zero value inner boundary condition
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          // set outer guard cells
          x[(lxx+1)*lz2+k] = 0.;
        }
      }
    }
  }

  // Exchange ghost cells of initial guess
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
  for(int i=0;i<lxx+2;i++) {
    x[i*lz2] = x[(i+1)*lz2-2];
    x[(i+1)*lz2-1] = x[i*lz2+1];
  }
   
  if(mgcount == 0) {
    soltime = 0.0;
    settime = 0.0;
  }
  else {
    if(kMG->pcheck > 0) {
      kMG->setPcheck(0);
    }
  }

  t0 = bout::globals::mpi->MPI_Wtime();
  generateMatrixF(level);

  if (kMG->xNP > 1)
    bout::globals::mpi->MPI_Barrier(commX);

  if ((pcheck == 3) && (mgcount == 0)) {
    FILE *outf;
    char outfile[256];
    sprintf(outfile,"test_matF_%d.mat",kMG->rProcI);
    output<<"Out file= "<<outfile<<endl;
    outf = fopen(outfile,"w");
    int dim =  (lxx+2)*(lzz+2);
    fprintf(outf,"dim = %d (%d, %d)\n",dim,lxx,lzz);

    for(int i = 0;i<dim;i++) {
      fprintf(outf,"%d ==",i);
      for(int j=0;j<9;j++) fprintf(outf,"%12.6f,",kMG->matmg[level][i*9+j]);
      fprintf(outf,"\n");
    }  
    fclose(outf);
  }

  if (level > 0) kMG->setMultigridC(0);

  if((pcheck == 3) && (mgcount == 0)) {
    for(int i = level; i> 0;i--) {
      output<<i<<"dimension= "<<kMG->lnx[i-1]<<"("<<kMG->gnx[i-1]<<"),"<<kMG->lnz[i-1]<<endl;
      
      FILE *outf;
      char outfile[256];
      sprintf(outfile,"test_matC%1d_%d.mat",i,kMG->rProcI);
      output<<"Out file= "<<outfile<<endl;
      outf = fopen(outfile,"w");
      int dim =  (kMG->lnx[i-1]+2)*(kMG->lnz[i-1]+2);
      fprintf(outf,"dim = %d (%d,%d)\n",dim,kMG->lnx[i-1],kMG->lnz[i-1]);
  
      for(int ii = 0;ii<dim;ii++) {
        fprintf(outf,"%d ==",ii);
        for(int j=0;j<9;j++) fprintf(outf,"%12.6f,",kMG->matmg[i-1][ii*9+j]);
        fprintf(outf,"\n");
      }  
      fclose(outf);
    }
  }

  t1 = bout::globals::mpi->MPI_Wtime();
  settime += t1-t0;

  // Compute solution.

  mgcount++;
  if (pcheck > 0)
    t0 = bout::globals::mpi->MPI_Wtime();

  kMG->getSolution(std::begin(x), std::begin(b), 0);

  if (pcheck > 0) {
    t1 = bout::globals::mpi->MPI_Wtime();
    soltime += t1-t0;
    if(mgcount%300 == 0) {
      output<<"Accumulated execution time at "<<mgcount<<" Sol "<<soltime<<" ( "<<settime<<" )"<<endl;
      settime = 0.;
      soltime = 0.;
    }
  }

  FieldPerp result{emptyFrom(b_in)};

#if CHECK > 2
  // Make any unused elements NaN so that user does not try to do calculations with them
  BOUT_FOR(i, result.getRegion("RGN_ALL")) { result[i] = BoutNaN; }
#endif

  // Copy solution into a FieldPerp to return
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for collapse(2))
  for (int i=1; i<lxx+1; i++) {
    for (int k=1; k<lzz+1; k++) {
      int i2 = i-1+localmesh->xstart;
      int k2 = k-1;
      result(i2, k2) = x[i*lz2+k];
    }
  }
  if (localmesh->firstX()) {
    if ( inner_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify gradient to set at inner boundary
        int i2 = -1+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = x[lz2+k] - x0(localmesh->xstart-1, k2)*sqrt(coords->g_11(localmesh->xstart, yindex))*coords->dx(localmesh->xstart, yindex);
        }
      }
      else {
        // zero gradient inner boundary condition
        int i2 = -1+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = x[lz2+k];
        }
      }
    }
    else {
      // Dirichlet boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify value to set at inner boundary
        int i2 = -1+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = 2.*x0(localmesh->xstart-1,k2) - x[lz2+k];
        }
      }
      else {
        // zero value inner boundary condition
        int i2 = -1+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = -x[lz2+k];
        }
      }
    }
  }
  if (localmesh->lastX()) {
    if ( outer_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify gradient to set at outer boundary
        int i2 = lxx+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = x[lxx*lz2+k] + x0(localmesh->xend+1, k2)*sqrt(coords->g_11(localmesh->xend, yindex))*coords->dx(localmesh->xend, yindex);
        }
      }
      else {
        // zero gradient outer boundary condition
        int i2 = lxx+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = x[lxx*lz2+k];
        }
      }
    }
    else {
      // Dirichlet boundary condition
      if ( outer_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify value to set at outer boundary
        int i2 = lxx+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = 2.*x0(localmesh->xend+1,k2) - x[lxx*lz2+k];
        }
      }
      else {
        // zero value inner boundary condition
        int i2 = lxx+localmesh->xstart;
BOUT_OMP(parallel default(shared) )
BOUT_OMP(for)
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          result(i2, k2) = -x[lxx*lz2+k];
        }
      }
    }
  }

  checkData(result);

  return result;
}

void LaplaceMultigrid::generateMatrixF(int level) {
  TRACE("LaplaceMultigrid::generateMatrixF(int)");
  
  // Set (fine-level) matrix entries

  BoutReal *mat;
  mat = kMG->matmg[level];
  int llx = kMG->lnx[level];
  int llz = kMG->lnz[level];

BOUT_OMP(parallel default(shared))
BOUT_OMP(for collapse(2))
  for (int i=1; i<llx+1; i++) {
    for (int k=1; k<llz+1; k++) {
      int i2 = i-1+localmesh->xstart;
      int k2 = k-1;
      int k2p  = (k2+1)%Nz_global;
      int k2m  = (k2+Nz_global-1)%Nz_global;

      BoutReal dz = coords->dz(i2, yindex);
      BoutReal ddx_C = (C2(i2+1, yindex, k2) - C2(i2-1, yindex, k2))/2./coords->dx(i2, yindex)/C1(i2, yindex, k2);
      BoutReal ddz_C =
          (C2(i2, yindex, k2p) - C2(i2, yindex, k2m)) / 2. / dz / C1(i2, yindex, k2);

      BoutReal ddx = D(i2, yindex, k2)*coords->g11(i2, yindex)/coords->dx(i2, yindex)/coords->dx(i2, yindex); 
               // coefficient of 2nd derivative stencil (x-direction)

      BoutReal ddz = D(i2, yindex, k2) * coords->g33(i2, yindex) / SQ(dz);
      // coefficient of 2nd derivative stencil (z-direction)

      BoutReal dxdz =
          D(i2, yindex, k2) * 2. * coords->g13(i2, yindex) / coords->dx(i2, yindex) / dz;
      // coefficient of mixed derivative stencil (could assume zero, at least initially,
      // if easier; then check this is true in constructor)

      BoutReal dxd = (D(i2, yindex, k2)*coords->G1(i2, yindex)
        + coords->g11(i2, yindex)*ddx_C
        + coords->g13(i2, yindex)*ddz_C // (could assume zero, at least initially, if easier; then check this is true in constructor)
      )/coords->dx(i2, yindex); // coefficient of 1st derivative stencil (x-direction)
      if (nonuniform) {
        // add correction for non-uniform dx
        dxd += D(i2, yindex, k2)*coords->d1_dx(i2, yindex);
      }

      BoutReal dzd =
          (D(i2, yindex, k2) * coords->G3(i2, yindex) + coords->g33(i2, yindex) * ddz_C
           + coords->g13(i2, yindex)
                 * ddx_C // (could assume zero, at least initially, if easier; then check
                         // this is true in constructor)
           )
          / dz; // coefficient of 1st derivative stencil (z-direction)

      int ic = i*(llz+2)+k;
      mat[ic*9] = dxdz/4.;
      mat[ic*9+1] = ddx - dxd/2.;
      mat[ic*9+2] = -dxdz/4.;
      mat[ic*9+3] = ddz - dzd/2.;
      mat[ic*9+4] = A(i2, yindex, k2) - 2.*(ddx+ddz); // coefficient of no-derivative component
      mat[ic*9+5] = ddz + dzd/2.;
      mat[ic*9+6] = -dxdz/4.;
      mat[ic*9+7] = ddx+dxd/2.;
      mat[ic*9+8] = dxdz/4.;
    }
  }

  // Here put boundary conditions

  if (kMG->rProcI == 0) {
    if ( inner_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int k = 1;k<llz+1; k++) {
        int ic = llz+2 +k;
        mat[ic*9+3] += mat[ic*9];
        mat[ic*9+4] += mat[ic*9+1];
        mat[ic*9+5] += mat[ic*9+2];
        b[ic] -= mat[ic*9]*x[k-1];
        b[ic] -= mat[ic*9+1]*x[k];
        b[ic] -= mat[ic*9+2]*x[k+1];
        mat[ic*9] = 0.;
        mat[ic*9+1] = 0.;
        mat[ic*9+2] = 0.;
      }
    }
    else {
      // Dirichlet boundary condition
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int k = 1;k<llz+1; k++) {
        int ic = llz+2 +k;
        mat[ic*9+3] -= mat[ic*9];
        mat[ic*9+4] -= mat[ic*9+1];
        mat[ic*9+5] -= mat[ic*9+2];
        b[ic] -= mat[ic*9]*x[k-1];
        b[ic] -= mat[ic*9+1]*x[k];
        b[ic] -= mat[ic*9+2]*x[k+1];
        mat[ic*9] = 0.;
        mat[ic*9+1] = 0.;
        mat[ic*9+2] = 0.;
      }
    }
  }
  if (kMG->rProcI == kMG->xNP-1) {
    if ( outer_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int k = 1;k<llz+1; k++) {
        int ic = llx*(llz+2)+k;
        mat[ic*9+3] += mat[ic*9+6];
        mat[ic*9+4] += mat[ic*9+7];
        mat[ic*9+5] += mat[ic*9+8];
        b[ic] -= mat[ic*9+6]*x[(llx+1)*(llz+2)+k-1];
        b[ic] -= mat[ic*9+7]*x[(llx+1)*(llz+2)+k];
        b[ic] -= mat[ic*9+8]*x[(llx+1)*(llz+2)+k+1];
        mat[ic*9+6] = 0.;
        mat[ic*9+7] = 0.;
        mat[ic*9+8] = 0.;
      }
    }
    else {
      // Dirichlet boundary condition
BOUT_OMP(parallel default(shared))
BOUT_OMP(for)
      for(int k = 1;k<llz+1; k++) {
        int ic = llx*(llz+2)+k;
        mat[ic*9+3] -= mat[ic*9+6];
        mat[ic*9+4] -= mat[ic*9+7];
        mat[ic*9+5] -= mat[ic*9+8];
        b[ic] -= mat[ic*9+6]*x[(llx+1)*(llz+2)+k-1];
        b[ic] -= mat[ic*9+7]*x[(llx+1)*(llz+2)+k];
        b[ic] -= mat[ic*9+8]*x[(llx+1)*(llz+2)+k+1];
        mat[ic*9+6] = 0.;
        mat[ic*9+7] = 0.;
        mat[ic*9+8] = 0.;
      }
    }
  }
}

#endif // BOUT_USE_METRIC_3D
