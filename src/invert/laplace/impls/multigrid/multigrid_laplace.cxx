/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Geometrical Multigrid Solver
 *
 * Equation solved is:
 *  d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2015 K.S. Kang
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

#include "multigrid_laplace.hxx"
int mgcount = 0;
BoutReal soltime=0.0,settime=0.0;


LaplaceMultigrid::LaplaceMultigrid(Options *opt) : 
  Laplacian(opt),
  A(0.0), C1(1.0), C2(1.0), D(1.0)
{

  // Get Options in Laplace Section
  if (!opt) opts = Options::getRoot()->getSection("laplace");
  else opts=opt;
  opts->get("multigridlevel",mglevel,7,true);
  opts->get("rtol",rtol,pow(10.0,-8),true);
  opts->get("atol",atol,pow(10.0,-20),true);
  opts->get("dtol",dtol,pow(10.0,5),true);
  opts->get("smtype",mgsm,1,true);
  opts->get("solvertype",mgplag,1,true);
  opts->get("cftype",cftype,0,true);
  opts->get("checking",pcheck,0,true);

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
  if (nonuniform) {
    throw BoutException("nonuniform option is not implemented in LaplaceMultigrid.");
  }
  
  commX = mesh->getXcomm();
  MPI_Comm_size(commX,&xNP);
  if(xNP > 1) {
    xProcI = mesh->getXProcIndex();
    if(xProcI == 0) xProcM = xNP-1;
    else xProcM = xProcI-1;
    if(xProcI == xNP-1) xProcP = 0;
    else xProcP = xProcI+1;
  }
  else {
    xProcI = 0;
    xProcM = 0;
    xProcP = 0;
  }
  
  Nx_local = mesh->xend - mesh->xstart + 1; // excluding guard cells
  Nx_global = mesh->GlobalNx - 2*mesh->xstart; // excluding guard cells
  
  if(mgcount == 0) {
    output <<"xNP="<<xNP<<"("<<xProcI<<","<<xProcM<<","<<xProcP<<"), Nx="<<Nx_global<<"("<<Nx_local<<")"<<endl;
  }
  zNP = 1;
  commXZ = commX;
  commZ = MPI_COMM_WORLD;
  if(zNP == 1) {
    zProcI = 0;
    zProcM = 0;
    zProcP = 0;
    Nz_global = mesh->GlobalNz - 1;
           // z-direction has one extra, unused point for historical reasons
    Nz_local = Nz_global; // No parallelization in z-direction (for now)
  }
  // 
  //else {
  //  Nz_local = mesh->zend - mesh->zstart + 1; // excluding guard cells
  //  Nz_global = mesh->GlobalNz - 2*mesh->zstart; // excluding guard cells
  // }
  if(mgcount==0) {
    output <<"zNP="<<zNP<<"("<<zProcI<<"), Nz="<<Nz_global<<"("<<Nz_local<<")"<<endl;
  }
  // Compute available levels along x-direction
  if(mglevel >1) {
    int nn = Nx_local;
    for(int n = mglevel;n > 1; n--) {
      if ( nn%2 != 0 )  {
	output<<"Size of local x-domain is not a power of 2^"<<mglevel<<" mglevel is changed to"<<mglevel-n+1<<endl;
        mglevel = mglevel - n + 1;
        n = 1;
      }
      nn = nn/2;
    }
  }
  if(mglevel >1) {
    int nn = Nz_local;
    for(int n = mglevel;n > 1; n--) {
      if ( nn%2 != 0 )  {
	output<<"Size of local z-domain is not a power of 2"<<mglevel <<" mglevel is changed to "<<mglevel - n + 1<<endl;
        mglevel = mglevel - n + 1;
        n = 1;
      }
      nn = nn/2;
    }
  }
  
  /* Momory allocate for Multigrid */
  gnx = new int[mglevel];
  gnz = new int[mglevel];
  lnx = new int[mglevel];
  lnz = new int[mglevel];
  gnx[mglevel-1] = Nx_global;
  gnz[mglevel-1] = Nz_global;
  lnx[mglevel-1] = Nx_local;
  lnz[mglevel-1] = Nz_local;
  if(mglevel > 1) {
    for(int i=mglevel-1;i>0;i--) {
      gnx[i-1] = gnx[i]/2;
      gnz[i-1] = gnz[i]/2;
      lnx[i-1] = lnx[i]/2;
      lnz[i-1] = lnz[i]/2;
    }
  }
  if(mgcount == 0) { 
    output<<"Before alloc mat "<<mglevel<<"("<<gnx[mglevel-1]<<","<<gnz[mglevel-1]<<"="<<lnx[mglevel-1]<<","<<lnz[mglevel-1]<<")"<<endl;
  }
  matmg = new BoutReal *[mglevel];

  for(int i = 0;i<mglevel;i++) {
    if(mgcount == 0) {
      output<<"In alloc mat "<<i<<"("<<gnx[i]<<","<<gnz[i]<<"="<<lnx[i]<<","<<lnz[i]<<")"<<(lnx[i]+2)*(lnz[i]+2)*9<<endl;
    }
    matmg[i] = new BoutReal[(lnx[i]+2)*(lnz[i]+2)*9];
  }
  
  // Set up Multigrid Cycle

  x = new BoutReal[(lnx[mglevel-1]+2)*(lnz[mglevel-1]+2)];
  b = new BoutReal[(lnx[mglevel-1]+2)*(lnz[mglevel-1]+2)];

  if(mgcount == 0) {  
    output<<" Smoothing type is ";
    if(mgsm ==0) output<<"Jacobi smoother"<<endl;
    else if(mgsm ==1) output<<" Gauss-Seidel smoother"<<endl;
    else throw BoutException("Undefined smoother");
    output<<"Solver type is ";
    if(mglevel == 1) output<<"PGMRES with simple Preconditioner"<<endl;
    else if(mgplag == 1) output<<"PGMRES with multigrid Preconditioner"<<endl;
    else output<<"Multigrid solver"<<endl;
  }
  
}

LaplaceMultigrid::~LaplaceMultigrid() {
  // Finalize, deallocate memory, etc.
  delete [] x;
  delete [] b;
  for(int i = 0;i<mglevel;i++) delete [] matmg[i];
  delete [] matmg;
  delete [] lnz;
  delete [] lnx;
  delete [] gnz;
  delete [] gnx;
}

const FieldPerp LaplaceMultigrid::solve(const FieldPerp &b_in, const FieldPerp &x0) {

  BoutReal t0,t1;
  
  yindex = b_in.getIndex();
  
  int lzz = lnz[mglevel-1];
  int lz2 = lzz+2;
  int lxx = lnx[mglevel-1];
  if(pcheck == 1) {
    output<<mgcount<<"Start of solve at "<<yindex<<endl;
  }


  if ( global_flags & INVERT_START_NEW ) {
    // set initial guess to zero
    for (int i=1; i<lxx+1; i++) {
      int i2 = i-1+mesh->xstart;
      for (int k=1; k<lzz+1; k++) {
        int k2 = k;
        x[i*lz2+k] = 0.;
      }
    }
  }
  else {
    // Read initial guess into local array, ignoring guard cells
    for (int i=1; i<lxx+1; i++) {
      int i2 = i-1+mesh->xstart;
      for (int k=1; k<lzz+1; k++) {
        int k2 = k-1;
        x[i*lz2+k] = x0[i2][k2];
      }
    }
  }
  if(pcheck == 1) {
    output<<"End of Initial condition at "<<xProcI<<endl;
  }
  
  // Read RHS into local array
  for (int i=1; i<lxx+1; i++) {
    int i2 = i-1+mesh->xstart;
    for (int k=1; k<lzz+1; k++) {
      int k2 = k-1;
      b[i*lz2+k] = b_in[i2][k2];
    }
  }
  
  if (mesh->firstX()) {
    if ( inner_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify gradient to set at inner boundary
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
	  x[k] = -x0[mesh->xstart-1][k2]*sqrt(mesh->g_11[mesh->xstart][yindex])*mesh->dx[mesh->xstart][yindex]; 
        }
      }
      else {
        // zero gradient inner boundary condition
        for (int k=1; k<lzz+1; k++) {
          // set inner guard cells
          x[k] = 0.0;
        }
      }
    }
    else {
      // Dirichlet boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify value to set at inner boundary
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          x[k] = 2.*x0[mesh->xstart-1][k2]; 
        // this is the value to set at the inner boundary
        }
      }
      else {
        // zero value inner boundary condition
        for (int k=1; k<lzz+1; k++) {
          // set inner guard cells
          x[k] = 0.;
        }
      }
    }
  }
  if (mesh->lastX()) {
    if ( outer_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      if ( inner_boundary_flags & INVERT_SET ) {
        // guard cells of x0 specify gradient to set at outer boundary
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
        x[(lxx+1)*lz2+k] = x0[mesh->xend+1][k2]*sqrt(mesh->g_11[mesh->xend][yindex])*mesh->dx[mesh->xend][yindex]; 
        // this is the value to set the gradient to at the outer boundary
        }
      }
      else {
        // zero gradient outer boundary condition
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
        for (int k=1; k<lzz+1; k++) {
          int k2 = k-1;
          x[(lxx+1)*lz2+k]=2.*x0[mesh->xend+1][k2]; 
          // this is the value to set at the outer boundary
        }
      }
      else {
        // zero value inner boundary condition
        for (int k=1; k<lzz+1; k++) {
          // set outer guard cells
          x[(lxx+1)*lz2+k] = 0.;
        }
      }
    }
  }

  // Exchange ghost cells of initial guess
  if(zNP > 1) {
    MPI_Status  status[4];
    MPI_Request requests[4];
    MPI_Datatype xvector;
    int stag,rtag,ierr;
    ierr = MPI_Type_vector(lxx+2, 1, lz2, MPI_DOUBLE, &xvector);
    ierr = MPI_Type_commit(&xvector);
    // Send to z+ and recieve from z-
    ierr = MPI_Sendrecv(&x[lz2-2],1,xvector,zProcP,stag,
                        &x[0],1,xvector,zProcM,rtag,commZ,status);
    // Send to z- and recieve from z+
    ierr = MPI_Sendrecv(&x[1],1,xvector,zProcM,stag,&x[lz2-1],
                        1,xvector,zProcP,rtag,commZ,status);
    ierr = MPI_Type_free(&xvector);
  }
  else {
    for(int i=0;i<lxx+2;i++) {
      x[i*lz2] = x[(i+1)*lz2-2];
      x[(i+1)*lz2-1] = x[i*lz2+1];
    }
  }
   
  if(mgcount == 0) {
    soltime = 0.0;
    settime = 0.0;
  }

  t0 = MPI_Wtime();
  generateMatrixF(mglevel-1);  
  if(pcheck == 1) {
    output<<"End of Matrix F generation at "<<xProcI<<endl;
  }

  if(xNP >1) MPI_Barrier(commX);

  /*
  FILE *outf;
  outf = fopen("test_matF.mat","w");
  int dim =  (lnx[mglevel-1]+2)*(lnz[mglevel-1]+2);
  fprintf(outf,"dim = %d\n",dim);

  for(int i = 0;i<dim;i++) {
    fprintf(outf,"%d ==",i);
    for(int j=0;j<9;j++) fprintf(outf,"%12.6f,",matmg[mglevel-1][i*9+j]);
    fprintf(outf,"\n");
  }  
  fclose(outf);
  */
  for(int i = mglevel-1; i> 0;i--) {
    setMatrixC(i);

    /*
    char outfile[256];
    sprintf(outfile,"test_matC%1d.mat",i);
    outf = fopen(outfile,"w");
    int dim =  (lnx[i-1]+2)*(lnz[i-1]+2);
    fprintf(outf,"dim = %d\n",dim);
  
    for(int ii = 0;ii<dim;ii++) {
      fprintf(outf,"%d ==",ii);
      for(int j=0;j<9;j++) fprintf(outf,"%12.6f,",matmg[i-1][ii*9+j]);
      fprintf(outf,"\n");
    }  
    fclose(outf);
    */
  }
  if(pcheck == 1) {
    output<<"End of all matrix at "<<xProcI<<endl;
  }

  t1 = MPI_Wtime();
  settime += t1-t0;
  // Compute solution.. 
  communications(x,mglevel-1);

  mgcount++;
  t0 = MPI_Wtime();
  if(mglevel == 1) pGMRES(x,b,mglevel-1,1);
  else if(mgplag == 1) pGMRES(x,b,mglevel-1,1);
  else solveMG(x,b,mglevel-1);
  t1 = MPI_Wtime();
  soltime += t1-t0;
  //if(mgcount%100 == 0) {
  //  output<<"Accumulated execution time at "<<mgcount<<" Sol "<<soltime<<" ( "<<settime<<" )"<<endl;
  //}
  
  FieldPerp result;
  result.allocate();
  #if CHECK>2
    // Make any unused elements NaN so that user does not try to do calculations with them
    result = 1./0.;
  #endif
  // Copy solution into a FieldPerp to return
  for (int i=0; i<lxx+2; i++) {
    int i2 = i-1+mesh->xstart;
    for (int k=1; k<lzz+1; k++) {
      int k2 = k-1;
      result[i2][k2] = x[i*lz2+k];
    }
  }
  result.setIndex(yindex); // Set the index of the FieldPerp to be returned
  
  return result;
  
}

void LaplaceMultigrid::communications(BoutReal* x, int level) {
  // nx does not include guard cells
  
  MPI_Status  status[4];
  MPI_Request requests[4];
  int stag,rtag,ierr;

  if(zNP > 1) {
    MPI_Datatype xvector;
    ierr = MPI_Type_vector(lnx[level], 1, lnz[level]+2, MPI_DOUBLE, &xvector);
    ierr = MPI_Type_commit(&xvector);
    // Send to z+ and recieve from z-
    stag = zProcI;
    rtag = zProcM;
    ierr = MPI_Sendrecv(&x[2*(lnz[level]+2)-2],1,xvector,zProcP,stag,
                        &x[lnz[level]+2],1,xvector,zProcM,rtag,commZ,status);
    // Send to z- and recieve from z+
    stag = zProcI+zNP;
    rtag = zProcM+zNP;
    ierr = MPI_Sendrecv(&x[1],1,xvector,zProcM,stag,&x[2*(lnz[level]+2)-1],
                        1,xvector,zProcP,rtag,commZ,status);
    ierr = MPI_Type_free(&xvector);
  }
  else {
    for(int i=1;i<lnx[level]+1;i++) {
      x[i*(lnz[level]+2)] = x[(i+1)*(lnz[level]+2)-2];
      x[(i+1)*(lnz[level]+2)-1] = x[i*(lnz[level]+2)+1];
    }
  }
  if(xNP > 1) {
    // Send to x+ and recieve from x-
    stag = xProcI; 
    rtag = xProcM;
    ierr = MPI_Sendrecv(&x[lnx[level]*(lnz[level]+2)],lnz[level]+2,
                MPI_DOUBLE,xProcP,stag,&x[0],lnz[level]+2,MPI_DOUBLE,
                xProcM,rtag,commX,status);
    // Send to x- and recieve from x+
    stag = xProcI+xNP;
    rtag = xProcP+xNP;;
    ierr = MPI_Sendrecv(&x[lnz[level]+2],lnz[level]+2,MPI_DOUBLE,xProcM,stag,
                        &x[(lnx[level]+1)*(lnz[level]+2)],lnz[level]+2,
                        MPI_DOUBLE,xProcP,rtag,commX,status);

  }

}

void LaplaceMultigrid::generateMatrixF(int level) {

  // Set (fine-level) matrix entries

  int i2,k2;


  for (int i=1; i<lnx[level]+1; i++) {
    i2 = i-1+mesh->xstart;
    for (int k=1; k<lnz[level]+1; k++) {
      k2 = k-1;
      int k2p  = (k2+1)%Nz_global;
      int k2m  = (k2+Nz_global-1)%Nz_global;
      
      BoutReal ddx_C = (C2[i2+1][yindex][k2] - C2[i2-1][yindex][k2])/2./mesh->dx[i2][yindex]/C1[i2][yindex][k2];
      BoutReal ddz_C = (C2[i2][yindex][k2p] - C2[i2][yindex][k2m]) /2./mesh->dz/C1[i2][yindex][k2];
      
      BoutReal ddx = D[i2][yindex][k2]*mesh->g11[i2][yindex]/mesh->dx[i2][yindex]/mesh->dx[i2][yindex]; 
               // coefficient of 2nd derivative stencil (x-direction)
      
      BoutReal ddz = D[i2][yindex][k2]*mesh->g33[i2][yindex]/mesh->dz/mesh->dz; 
              // coefficient of 2nd derivative stencil (z-direction)
      
      BoutReal dxdz = D[i2][yindex][k2]*mesh->g13[i2][yindex]/mesh->dx[i2][yindex]/mesh->dz/2.; 
              // coefficient of mixed derivative stencil (could assume zero, at least initially, 
              // if easier; then check this is true in constructor)
      
      BoutReal dxd = (D[i2][yindex][k2]*2.*mesh->G1[i2][yindex]
        + mesh->g11[i2][yindex]*ddx_C
        + mesh->g13[i2][yindex]*ddz_C // (could assume zero, at least initially, if easier; then check this is true in constructor)
      )/mesh->dx[i2][yindex]; // coefficient of 1st derivative stencil (x-direction)
      
      BoutReal dzd = (D[i2][yindex][k2]*2.*mesh->G3[i2][yindex]
        + mesh->g33[i2][yindex]*ddz_C
        + mesh->g13[i2][yindex]*ddx_C // (could assume zero, at least initially, if easier; then check this is true in constructor)
      )/mesh->dz; // coefficient of 1st derivative stencil (z-direction)
      
      int ic = i*(lnz[level]+2)+k;
      matmg[level][ic*9] = dxdz/4.;
      matmg[level][ic*9+1] = ddx - dxd/2.;
      matmg[level][ic*9+2] = -dxdz/4.;
      matmg[level][ic*9+3] = ddz - dzd/2.;
      matmg[level][ic*9+4] = A[i2][yindex][k2] - 2.*(ddx+ddz); // coefficient of no-derivative component
      matmg[level][ic*9+5] = ddz + dzd/2.;
      matmg[level][ic*9+6] = -dxdz/4.;
      matmg[level][ic*9+7] = ddx+dxd/2.;
      matmg[level][ic*9+8] = dxdz/4.;
    }
  }

  // Here put boundary conditions

  if (xProcI == 0) {
    if ( inner_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      for(int k = 1;k<lnz[level]+1; k++) {
        int ic = lnz[level]+2 +k;
        matmg[level][ic*9+3] += matmg[level][ic*9];
        matmg[level][ic*9+4] += matmg[level][ic*9+1];
        matmg[level][ic*9+5] += matmg[level][ic*9+2];
        b[ic] -= matmg[level][ic*9]*x[k-1];
        b[ic] -= matmg[level][ic*9+1]*x[k];
        b[ic] -= matmg[level][ic*9+2]*x[k+1];
        matmg[level][ic*9] = 0.;
        matmg[level][ic*9+1] = 0.;
        matmg[level][ic*9+2] = 0.;
      }
    }
    else {
      // Dirichlet boundary condition
      for(int k = 1;k<lnz[level]+1; k++) {
        int ic = lnz[level]+2 +k;
        matmg[level][ic*9+3] -= matmg[level][ic*9];
        matmg[level][ic*9+4] -= matmg[level][ic*9+1];
        matmg[level][ic*9+5] -= matmg[level][ic*9+2];
        b[ic] -= matmg[level][ic*9]*x[k-1];
        b[ic] -= matmg[level][ic*9+1]*x[k];
        b[ic] -= matmg[level][ic*9+2]*x[k+1];
        matmg[level][ic*9] = 0.;
        matmg[level][ic*9+1] = 0.;
        matmg[level][ic*9+2] = 0.;
      }
    }
  }
  if (xProcI == xNP-1) {
    if ( outer_boundary_flags & INVERT_AC_GRAD ) {
      // Neumann boundary condition
      for(int k = 1;k<lnz[level]+1; k++) {
        int ic = lnx[level]*(lnz[level]+2)+k;
        matmg[level][ic*9+3] += matmg[level][ic*9+6];
        matmg[level][ic*9+4] += matmg[level][ic*9+7];
        matmg[level][ic*9+5] += matmg[level][ic*9+8];
        b[ic] -= matmg[level][ic*9+6]*x[(lnx[level]+1)*(lnz[level]+2)+k-1];
        b[ic] -= matmg[level][ic*9+7]*x[(lnx[level]+1)*(lnz[level]+2)+k];
        b[ic] -= matmg[level][ic*9+8]*x[(lnx[level]+1)*(lnz[level]+2)+k+1];
        matmg[level][ic*9+6] = 0.;
        matmg[level][ic*9+7] = 0.;
        matmg[level][ic*9+8] = 0.;
      }
    }
    else {
      // Dirichlet boundary condition
      for(int k = 1;k<lnz[level]+1; k++) {
        int ic = lnx[level]*(lnz[level]+2)+k;
        matmg[level][ic*9+3] -= matmg[level][ic*9+6];
        matmg[level][ic*9+4] -= matmg[level][ic*9+7];
        matmg[level][ic*9+5] -= matmg[level][ic*9+8];
        b[ic] -= matmg[level][ic*9+6]*x[(lnx[level]+1)*(lnz[level]+2)+k-1];
        b[ic] -= matmg[level][ic*9+7]*x[(lnx[level]+1)*(lnz[level]+2)+k];
        b[ic] -= matmg[level][ic*9+8]*x[(lnx[level]+1)*(lnz[level]+2)+k+1];
        matmg[level][ic*9+6] = 0.;
        matmg[level][ic*9+7] = 0.;
        matmg[level][ic*9+8] = 0.;
      }
    }
  }
}

void LaplaceMultigrid::setMatrixC(int level) {

  BoutReal ratio=8.0,val; 

  for(int i = 1;i<lnx[level-1]+1;i++) {
    int i2 = 2*i-1;
    for(int k = 1;k<lnz[level-1]+1;k++) {
      int k2 = 2*k-1;
      int mm = i*(lnz[level-1]+2)+k;
      int m0 = i2*(lnz[level]+2)+k2;
      int m1 = i2*(lnz[level]+2)+k2+1;
      int m2 = (i2+1)*(lnz[level]+2)+k2;
      int m3 = (i2+1)*(lnz[level]+2)+k2+1;
      val = matmg[level][m0*9+4]+matmg[level][m1*9+4];
      val += matmg[level][m2*9+4] + matmg[level][m3*9+4];
      val += matmg[level][m0*9+5] + matmg[level][m1*9+3];
      val += matmg[level][m2*9+5] + matmg[level][m3*9+3];
      val += matmg[level][m0*9+7] + matmg[level][m2*9+1];
      val += matmg[level][m1*9+7] + matmg[level][m3*9+1];
      val += matmg[level][m0*9+8] + matmg[level][m3*9];
      val += matmg[level][m1*9+6] + matmg[level][m2*9+2];
      matmg[level-1][mm*9+4] = val/ratio;
      val = matmg[level][m0*9+1]+matmg[level][m1*9+1];
      val += matmg[level][m0*9+2]+matmg[level][m1*9];
      matmg[level-1][mm*9+1] = val/ratio;
      val = matmg[level][m0*9+3]+matmg[level][m2*9+3];
      val += matmg[level][m0*9+6]+matmg[level][m2*9];
      matmg[level-1][mm*9+3] = val/ratio;
      val = matmg[level][m1*9+5]+matmg[level][m3*9+5];
      val += matmg[level][m1*9+8]+matmg[level][m3*9+2];
      matmg[level-1][mm*9+5] = val/ratio;
      val = matmg[level][m2*9+7]+matmg[level][m3*9+7];
      val += matmg[level][m2*9+8]+matmg[level][m3*9+6];
      matmg[level-1][mm*9+7] = val/ratio;
      matmg[level-1][mm*9] = matmg[level][m0*9]/ratio;
      matmg[level-1][mm*9+2] = matmg[level][m1*9+2]/ratio;
      matmg[level-1][mm*9+6] = matmg[level][m2*9+6]/ratio;
      matmg[level-1][mm*9+8] = matmg[level][m3*9+8]/ratio;      
    }
  }

}
