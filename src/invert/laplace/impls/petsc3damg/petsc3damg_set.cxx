 /**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Algebraic multigrid Solver
 *                              with PETSc library
 *
 * Equation solved is:
 *  d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 *
 **************************************************************************
 * Copyright 2018 K.S. Kang kskang@ipp.mpg.de
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

#include "bout/mesh.hxx"
#include "boutcomm.hxx"

#include "petsc3damg.hxx"

void LaplacePetsc3DAmg::settingSolver(int kflag){
  //  Timer timer("invert");
  TRACE("LaplacePetsc3DAmg::settingSolver(int)");
  
  if(!opts) {
    // If no options supplied, use default
    opts = Options::getRoot()->getSection("petsc3damg");
  }

  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate(BoutComm::get(), &ksp ); 
  KSPGetPC(ksp, &pc);
  
  // Configure Linear Solver
  if(kflag == 0) KSPSetOperators( ksp, MatA, MatA );
  else KSPSetOperators( ksp, MatA, MatP );
   // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
    // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
  std::string solt;
  solt.assign(soltype,0,2);
  if(solt == "di") {
    if(xNP*zNP == 1) {
      KSPSetType(ksp,KSPPREONLY);
      PCSetType(pc,PCLU);
    }
    else {
      KSPSetType( ksp, KSPGMRES );
      KSPSetInitialGuessNonzero(ksp, (PetscBool) true );
      PCSetType(pc,PCBJACOBI);
      if(soltype == "direct1") {
        int *blks = new int[Nx_global*Ny_global];
	for(int k=0;k<Ny_global;k++) {
          for(int i=0;i<Nx_global;i++) blks[k*Nx_global+i] = Nz_local;
	}
        PCBJacobiSetTotalBlocks(pc,Nx_global*Ny_global,blks);
        delete [] blks;
      }
    }
  }
  else {
    KSPSetType( ksp, KSPGMRES );
    KSPSetInitialGuessNonzero(ksp, (PetscBool) true );
    if(soltype == "gmres") {
      PCSetType(pc,PCBJACOBI);
      int *blks = new int[Ny_global*Nx_global*2];
      for(int k = 0;k<Ny_global;k++) {
        for(int i=0;i<Nx_global;i++) {
  	  blks[2*(k*Nx_global+i)] = Nz_local/2;
	  blks[2*(k*Nx_global+i)+1] = Nz_local - Nz_local/2;
        }
      }
      PCBJacobiSetTotalBlocks(pc,Ny_global*Nx_global*2,blks);
      delete [] blks;
    }
    else {
      if(solt == "ml") {
	PCSetType(pc,PCML);
        PCMGSetLevels(pc,mglevel,NULL);
	PCMGSetCycleType(pc,PC_MG_CYCLE_V);
	PCMGSetNumberSmooth(pc,2);
      }
      else if(solt == "hy") {
	PCSetType(pc,PCHYPRE);
	PCHYPRESetType(pc,"boomeramg");
        if(soltype == "hypre0") {
  	  char mclev[3];
	  if(mglevel > 9) sprintf(mclev,"%2d",mglevel);
	  else sprintf(mclev,"0%1d",mglevel);
	  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_max_levels",mclev);
	  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_down","2");
	  PetscOptionsSetValue(NULL,"-pc_hypre_boomeramg_grid_sweeps_up","2");
        }
      }
      else { // For gamg 
        PCSetType(pc,PCGAMG);
        PCGAMGSetType(pc,PCGAMGAGG);
        if(soltype == "gamgag") {
	  PCGAMGSetNSmooths(pc,2);
        }
        else if(soltype == "gamgag0") {
  	  PCGAMGSetNSmooths(pc,3);
        }
        else if(soltype == "gamgag1") {
  	  PCGAMGSetNSmooths(pc,4);
        }
      //      PCGAMGSetCycleType(pc,PC_MG_CYCLE_V);
        PCMGSetLevels(pc,mglevel,NULL);
	PCMGSetCycleType(pc,PC_MG_CYCLE_V);
	PCMGSetNumberSmooth(pc,2);
      }
    }
    if(rightpre) KSPSetPCSide(ksp, PC_RIGHT); // Right preconditioning
    else         KSPSetPCSide(ksp, PC_LEFT);  // Left preconditioning
  }
  KSPSetTolerances(ksp,rtol,atol,dtol,maxits);  
  KSPSetFromOptions(ksp);
}

const Field3D LaplacePetsc3DAmg::solve(const Field3D &rhs, const Field3D &x0) {
  // Timer timer("invert");
  TRACE("LaplacePetsc3DAmg::solve(const Field3D, const Field3D)");
  
  // Load initial guess x0 into xs and rhs into bs
  Mesh *mesh = rhs.getMesh();  // Where to get initializing LaplacePetscAmg
  Coordinates *coords = mesh->getCoordinates();
  BoutReal tmss, tms, tmf, tmG, tmS, tmR, tsol;
  if(fcheck) tmss = MPI_Wtime();
  generateMatrixA(elemf);
  if(diffpre > 0) generateMatrixP(elemf);
  if(fcheck) {
    tms = MPI_Wtime();
    tmG = tms - tmss;
  }

  settingSolver(diffpre);
  if(fcheck) {
    tmf = MPI_Wtime();
    tmS = tmf - tms;
  }
  // MPI_Barrier(MPI_COMM_WORLD);
  
  int ind,i2,i,j,j2,j2p,j2m,k,k2,nn,nxzt;
  PetscScalar val,volm;
  nxzt = nzt*nxt;
  if(fcheck) tms = MPI_Wtime();
  if ( global_flags & INVERT_START_NEW ) {
    // set initial guess to zero
    for (k=0; k < Ny_local; k++) {
      k2 = k + mystart;
      for (i=0; i < Nx_local; i++) {
        i2 = i + mxstart;
        volm = coords->dy(i2,k2)*coords->dx(i2, k2)*coords->dz;
        for (j=0; j < Nz_local; j++) {
          ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
          val = 0.;
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
          val = rhs(i2, k2, j+mzstart)*volm;
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
        }
      }
    }
  }
  else {
    // Read initial guess into local array, ignoring guard cells
    for (k=0; k < Ny_local; k++) {
      k2 = k + mystart;
      for (i=0; i < Nx_local; i++) {
        i2 = i + mxstart;
        for (j=0; j < Nz_local; j++) {
          ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
          val = x0(i2, k2, j+mzstart);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
          val = rhs(i2, k2, j+mzstart)*volm;
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
        }
      }
    }
  }
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);
  // For the boundary conditions
  int tmax;
  tmax = nyt;
  if(nxt > tmax) tmax = nxt;
  BoutReal tval[tmax*nzt];
  if(yProcI == 0) {
    BoutReal dhx,dhy,dhz,volm,gt12,gt22,gt23,ddJ,ddx_C,ddy_C,ddz_C;
    BoutReal ddy,dxdy,dzdy,dyd,dydz;
    if(ybdcon != 0) { // Boundary value along with y-direction, i.e. x and z index
      k2 = mystart;
      if(ybdcon%3 == 1) { // Neumann
        for(i = 0; i< nxt; i++) {
	  i2 = i + mxstart - lxs;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = i*nzt + j;
            tval[nn] = -x0(i2, k2-1, j2)*sqrt(coords->g_22(i2, k2))*coords->dy(i2, k2); 
	  }
	}
      }
      else { // Dirichlet
        for(i = 0; i< nxt; i++) {
	  i2 = i + mxstart - lxs;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = i*nzt + j;
            tval[nn] = 2.*x0(i2, k2-1, j2); 
	  }
	}
      }
      dhz = coords->dz;
      for(i = 0;i < Nx_local;i++) {
	i2 = i + mxstart;
	dhx = coords->dx(i2,k2);
	dhy = coords->dy(i2,k2);
        for(j = 0;j < Nz_local;j++) {
          j2 = j + mzstart;
	  j2m = (j2-1+Nz_global)%Nz_global;
	  j2p = (j2+1)%Nz_global;
  	  gt12 = coords->g12(i2,k2);
	  gt22 = coords->g22(i2,k2) - 1.0/coords->g22(i2,k2);
	  gt23 = coords->g23(i2,k2);

	  ddJ = (coords->J(i2,k2+1)/coords->g22(i2,k2+1) - coords->J(i2,k2-1)/coords->g22(i2,k2-1))
	    /2./dhy/coords->J(i2,k2);
          ddx_C = (C2(i2+1,k2,j2)-C2(i2-1,k2,j2))/2./dhx/C1(i2,k2,j2);
          ddy_C = (C2(i2,k2+1,j2)-C2(i2,k2-1,j2))/2./dhy/C1(i2,k2,j2);
          ddz_C = (C2(i2,k2,j2p)-C2(i2,k2,j2m))/2./dhz/C1(i2,k2,j2);
      
        // ddx = D(i2,k2,j2)*gt11/dhx/dhx; 
          ddy = D(i2,k2,j2)*gt22/dhy/dhy;    
        // ddz = D(i2,k2,j2)*gt33/dhz/dhz; 
      
          dxdy = 2.0*D(i2,k2,j2)*gt12/dhx/dhy; 
        // dxdz = 2.0*D(i2,k2,j2)*gt13/dhx/dhz; 
          dydz = 2.0*D(i2,k2,j2)*gt23/dhy/dhz; 
       
        // dxd = (D(i2,k2,j2)*coords->G1(i2,k2) + gt11*ddx_C + gt12*ddy_C + gt13*ddz_C)/dhx;
          dyd = (D(i2,k2,j2)*(coords->G2(i2,k2) - ddJ)+gt12*ddx_C + gt22*ddy_C + gt23*ddz_C)/dhy;
        // dzd = (D(i2,k2,j2)*coords->G3(i2,k2) + gt33*ddz_C +gt13*ddx_C + gt23*ddy_C)/dhz;
          volm = dhx*dhy*dhz;

          val = -(ddy - dyd/2.0)*volm*tval[i*nzt + j]; // Mat_2
	  if((i2 == 0) && (lxs == 0)) { // i == 0
	    if(xbdcon == 0)  val -= dydz/4.0*volm*tval[(nxt-1)*nzt+j]; // Periodic
	    else  val -= dydz/8.0*volm*tval[i*nzt+j]; // Dirichlet and Neumann
	  }
	  else val -= dydz/4.0*volm*tval[(i-1)*nzt+j]; // mat_0
	  if((i2+1 == Nx_global) && (Nx_local + lxs == nxt)) {
	    if(xbdcon == 0) val += dydz/4.0*volm*tval[j];
	    else val += dydz/4.0*volm*tval[i*nzt+j]; // Dirichlet and Neumann 
	  }
          else val += dydz/4.0*volm*tval[(i+1)*nzt+j]; // mat_4
	  
	  if((j2 == 0) && (lzs == 0)) { // j == 0 always Periodic
	      val -= dxdy/4.0*volm*tval[i*nzt+nzt-1];
	  }
	  else val -= dxdy/4.0*volm*tval[i*nzt+j-1];  // mat_1
	  if((j2+1 == Nz_global) && (Nz_local + lzs == nzt)) {
            val += dxdy/4.0*volm*tval[i*nzt];
	  }
	  else val += dxdy/4.0*volm*tval[i*nzt+j+1];   //mat_2
          ind = gindices[lys*nxzt + (i+lxs)*nzt+j+lzs];
          VecSetValues( bs, 1, &ind, &val, ADD_VALUES );
        }
      }
    }
  }

  if(yProcI == yNP - 1) {
    BoutReal dhx,dhy,dhz,volm,gt12,gt22,gt23,ddJ,ddx_C,ddy_C,ddz_C;
    BoutReal ddy,dxdy,dzdy,dyd,dydz;
    if(ybdcon != 0) { // Boundary value along with y-direction, i.e. x and z index
      k2 = Ny_global - 1;
      if(ybdcon%5 == 1) { // Neumann
        for(i = 0; i< nxt; i++) {
	  i2 = i + mxstart - lxs;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = i*nzt + j;
            tval[nn] = x0(i2, k2+1, j2)*sqrt(coords->g_22(i2, k2))*coords->dy(i2, k2); 
	  }
	}
      }
      else { // Dirichlet
        for(i = 0; i< nxt; i++) {
	  i2 = i + mxstart - lxs;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = i*nzt + j;
            tval[nn] = 2.*x0(i2, k2+1, j2); 
	  }
	}
      }
      dhz = coords->dz;
      for(i = 0;i < Nx_local;i++) {
	i2 = i + mxstart;
	dhx = coords->dx(i2,k2);
	dhy = coords->dy(i2,k2);
        for(j = 0;j < Nz_local;j++) {
          j2 = j + mzstart;
	  j2m = (j2-1+Nz_global)%Nz_global;
	  j2p = (j2+1)%Nz_global;
  	  gt12 = coords->g12(i2,k2);
	  gt22 = coords->g22(i2,k2) - 1.0/coords->g22(i2,k2);
	  gt23 = coords->g23(i2,k2);

	  ddJ = (coords->J(i2,k2+1)/coords->g22(i2,k2+1) - coords->J(i2,k2-1)/coords->g22(i2,k2-1))
	    /2./dhy/coords->J(i2,k2);
          ddx_C = (C2(i2+1,k2,j2)-C2(i2-1,k2,j2))/2./dhx/C1(i2,k2,j2);
          ddy_C = (C2(i2,k2+1,j2)-C2(i2,k2-1,j2))/2./dhy/C1(i2,k2,j2);
          ddz_C = (C2(i2,k2,j2p)-C2(i2,k2,j2m))/2./dhz/C1(i2,k2,j2); 
      
        // ddx = D(i2,k2,j2)*gt11/dhx/dhx; 
          ddy = D(i2,k2,j2)*gt22/dhy/dhy;    
        // ddz = D(i2,k2,j2)*gt33/dhz/dhz; 
      
          dxdy = 2.0*D(i2,k2,j2)*gt12/dhx/dhy; 
        // dxdz = 2.0*D(i2,k2,j2)*gt13/dhx/dhz; 
          dydz = 2.0*D(i2,k2,j2)*gt23/dhy/dhz; 
       
        // dxd = (D(i2,k2,j2)*coords->G1(i2,k2) + gt11*ddx_C + gt12*ddy_C + gt13*ddz_C)/dhx;
          dyd = (D(i2,k2,j2)*(coords->G2(i2,k2) - ddJ)+gt12*ddx_C + gt22*ddy_C + gt23*ddz_C)/dhy;
        // dzd = (D(i2,k2,j2)*coords->G3(i2,k2) + gt33*ddz_C +gt13*ddx_C + gt23*ddy_C)/dhz;
          volm = dhx*dhy*dhz;
	  
          val = -(ddy + dyd/2.0)*volm*tval[i*nzt + j];  // mat_16
	  if((i2 == 0) && (lxs == 0)) { // i == 0
	    if(xbdcon == 0)  val += dydz/4.0*volm*tval[(nxt-1)*nzt+j]; // Periodic
	    else  val += dydz/8.0*volm*tval[i*nzt+j]; // Dirichlet and Neumann
	  }
	  else val += dydz/4.0*volm*tval[(i-1)*nzt+j]; // mat_14
	  if((i2+1 == Nx_global) && (Nx_local + lxs == nxt)) {
	    if(xbdcon == 0) val -= dydz/4.0*volm*tval[j];
	    else val -= dydz/4.0*volm*tval[i*nzt+j]; // Dirichlet and Neumann 
	  }
          else val -= dydz/4.0*volm*tval[(i+1)*nzt+j];// mat_18
	  
	  if((j2 == 0) && (lzs == 0)) { // j == 0 always Periodic
	      val += dxdy/4.0*volm*tval[i*nzt+nzt-1];
	  }
	  else val += dxdy/4.0*volm*tval[i*nzt+j-1]; // mat_15
	  if((j2+1 == Nz_global) && (Nz_local + lzs == nzt)) {
            val -= dxdy/4.0*volm*tval[i*nzt];
	  }
	  else val -= dxdy/4.0*volm*tval[i*nzt+j+1];    // mat_17
          ind = gindices[(nyt-lys-1)*nxzt + (i+lxs)*nzt+j+lzs];
          VecSetValues( bs, 1, &ind, &val, ADD_VALUES );
        }
      }
    }
  }
  
  if(xProcI == 0) { 
    BoutReal dhx,dhy,dhz,volm,gt11,gt12,gt13,gt22,gt23,ddJ,ddx_C,ddy_C,ddz_C;
    BoutReal ddx,dxdz,dxdy,dxd,dydz;
    if(xbdcon != 0) { // Boundary value along with x-direction, i.e. y and z index
      i2 = mxstart;
      if(xbdcon%3 == 1) { // Neumann
        for(k = 0; k< nyt; k++) {
	  k2 = k + mystart - lys;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = k*nzt + j;
            tval[nn] = -x0(i2-1, k2, j2)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
	  }
	}
      }
      else { // Dirichlet
        for(k = 0; k< nyt; k++) {
	  k2 = i + mystart - lys;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = k*nzt + j;
            tval[nn] = 2.*x0(i2-1, k2, j2); 
	  }
	}
      }
      dhz = coords->dz;
      for(k = 0;k < Ny_local;k++) {
	k2 = k + mystart;
	dhx = coords->dx(i2,k2);
	dhy = coords->dy(i2,k2);
        for(j = 0;j < Nz_local;j++) {
          j2 = j + mzstart;
	  j2m = (j2-1+Nz_global)%Nz_global;
	  j2p = (j2+1)%Nz_global;
  	  gt11 = coords->g11(i2,k2);
  	  gt12 = coords->g12(i2,k2);
  	  gt13 = coords->g13(i2,k2);
	  gt22 = coords->g22(i2,k2) - 1.0/coords->g22(i2,k2);
	  gt23 = coords->g23(i2,k2);

	  ddJ = (coords->J(i2+1,k2)/coords->g22(i2+1,k2) - coords->J(i2-1,k2)/coords->g22(i2-1,k2))
	    /2./dhy/coords->J(i2,k2);
          ddx_C = (C2(i2+1,k2,j2)-C2(i2-1,k2,j2))/2./dhx/C1(i2,k2,j2);
          ddy_C = (C2(i2,k2+1,j2)-C2(i2,k2-1,j2))/2./dhy/C1(i2,k2,j2);
          ddz_C = (C2(i2,k2,j2p)-C2(i2,k2,j2m))/2./dhz/C1(i2,k2,j2);
      
          ddx = D(i2,k2,j2)*gt11/dhx/dhx; 
	//  ddy = D(i2,k2,j2)*gt22/dhy/dhy;    
        //  ddz = D(i2,k2,j2)*gt33/dhz/dhz; 
      
          dxdy = 2.0*D(i2,k2,j2)*gt12/dhx/dhy; 
          dxdz = 2.0*D(i2,k2,j2)*gt13/dhx/dhz; 
	  dydz = 2.0*D(i2,k2,j2)*gt23/dhy/dhz; 
       
          dxd = (D(i2,k2,j2)*coords->G1(i2,k2) + gt11*ddx_C + gt12*ddy_C + gt13*ddz_C)/dhx;
	//  dyd = (D(i2,k2,j2)*(coords->G2(i2,k2) - ddJ)+gt12*ddx_C + gt22*ddy_C + gt23*ddz_C)/dhy;
        //  dzd = (D(i2,k2,j2)*coords->G3(i2,k2) + gt33*ddz_C +gt13*ddx_C + gt23*ddy_C)/dhz;
          volm = dhx*dhy*dhz;  

          val = -(ddx - dxd/2.0)*volm*tval[k*nzt + j]; // mat_6
	  if((k2 == 0) && (lys == 0)) { // k == 0
	    if(ybdcon == 0)  val -= dydz/4.0*volm*tval[(nyt-1)*nzt+j]; // Periodic
	    else  val -= dydz/8.0*volm*tval[k*nzt+j]; // Dirichlet and Neumann
	  }
	  else val -= dydz/4.0*volm*tval[(k-1)*nzt+j]; // mat_0
	  if((k2+1 == Ny_global) && (Ny_local + lys == nyt)) {
	    if(ybdcon == 0) val += dydz/4.0*volm*tval[j];
	    else val += dydz/4.0*volm*tval[k*nzt+j]; // Dirichlet and Neumann 
	  }
          else val += dydz/4.0*volm*tval[(k+1)*nzt+j]; //mat_14
	  
	  if((j2 == 0) && (lzs == 0)) { // j == 0 always Periodic
	      val -= dxdz/4.0*volm*tval[k*nzt+nzt-1];
	  }
	  else val -= dxdz/4.0*volm*tval[k*nzt+j-1]; //mat_5
	  if((j2+1 == Nz_global) && (Nz_local + lzs == nzt)) {
            val += dxdz/4.0*volm*tval[k*nzt];
	  }
	  else val += dxdz/4.0*volm*tval[k*nzt+j+1]; // mat_7
          ind = gindices[(k+lys)*nxzt + lxs*nzt+j+lzs];
          VecSetValues( bs, 1, &ind, &val, ADD_VALUES );
        }
      }
    }
  }

  if(xProcI == xNP - 1) { 
    BoutReal dhx,dhy,dhz,volm,gt11,gt12,gt13,gt22,gt23,ddJ,ddx_C,ddy_C,ddz_C;
    BoutReal ddx,dxdz,dxdy,dxd,dydz;
    if(xbdcon != 0) { // Boundary value along with x-direction, i.e. y and z index
      i2 = Nx_global - 1;
      if(xbdcon%5 == 1) { // Neumann
        for(k = 0; k< nyt; k++) {
	  k2 = k + mystart - lys;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = k*nzt + j;
            tval[nn] = x0(i2+1, k2, j2)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
	  }
	}
      }
      else { // Dirichlet
        for(k = 0; k< nyt; k++) {
	  k2 = k + mystart - lys;
	  for(j = 0;j < nzt; j++) {
	    j2 = j + mzstart - lzs;
	    nn = k*nzt + j;
            tval[nn] = 2.*x0(i2+1, k2, j2); 
	  }
	}
      }
      dhz = coords->dz;
      for(k = 0;k < Ny_local;k++) {
	k2 = k + mystart;
	dhx = coords->dx(i2,k2);
	dhy = coords->dy(i2,k2);
        for(j = 0;j < Nz_local;j++) {
          j2 = j + mzstart;
	  j2m = (j2-1+Nz_global)%Nz_global;
	  j2p = (j2+1)%Nz_global;
  	  gt11 = coords->g11(i2,k2);
  	  gt12 = coords->g12(i2,k2);
  	  gt13 = coords->g13(i2,k2);
	  gt22 = coords->g22(i2,k2) - 1.0/coords->g22(i2,k2);
	  gt23 = coords->g23(i2,k2);

	  ddJ = (coords->J(i2,k2+1)/coords->g22(i2,k2+1) - coords->J(i2,k2-1)/coords->g22(i2,k2-1))
	    /2./dhy/coords->J(i2,k2);
          ddx_C = (C2(i2+1,k2,j2)-C2(i2-1,k2,j2))/2./dhx/C1(i2,k2,j2);
          ddy_C = (C2(i2,k2+1,j2)-C2(i2,k2-1,j2))/2./dhy/C1(i2,k2,j2);
          ddz_C = (C2(i2,k2,j2p)-C2(i2,k2,j2m))/2./dhz/C1(i2,k2,j2);
      
          ddx = D(i2,k2,j2)*gt11/dhx/dhx; 
	//   ddy = D(i2,k2,j2)*gt22/dhy/dhy;    
        // ddz = D(i2,k2,j2)*gt33/dhz/dhz; 
      
          dxdy = 2.0*D(i2,k2,j2)*gt12/dhx/dhy; 
          dxdz = 2.0*D(i2,k2,j2)*gt13/dhx/dhz; 
	  dydz = 2.0*D(i2,k2,j2)*gt23/dhy/dhz; 
       
          dxd = (D(i2,k2,j2)*coords->G1(i2,k2) + gt11*ddx_C + gt12*ddy_C + gt13*ddz_C)/dhx;
	//  dyd = (D(i2,k2,j2)*(coords->G2(i2,k2) - ddJ)+gt12*ddx_C + gt22*ddy_C + gt23*ddz_C)/dhy;
        //  dzd = (D(i2,k2,j2)*coords->G3(i2,k2) + gt33*ddz_C +gt13*ddx_C + gt23ddy_C)/dhz;
          volm = dhx*dhy*dhz;
	  	  
          val = -(ddx + dxd/2.0)*volm*tval[k*nzt + j];  // mat_12
	  if((k2 == 0) && (lys == 0)) { // k == 0
	    if(ybdcon == 0)  val += dydz/4.0*volm*tval[(nyt-1)*nzt+j]; // Periodic
	    else  val += dydz/8.0*volm*tval[k*nzt+j]; // Dirichlet and Neumann
	  }
	  else val += dydz/4.0*volm*tval[(k-1)*nzt+j]; // mat_4
	  if((k2+1 == Ny_global) && (Ny_local + lys == nyt)) {
	    if(ybdcon == 0) val -= dydz/4.0*volm*tval[j];
	    else val -= dydz/4.0*volm*tval[k*nzt+j]; // Dirichlet and Neumann 
	  }
          else val -= dydz/4.0*volm*tval[(k+1)*nzt+j];// mat_18
	  
	  if((j2 == 0) && (lzs == 0)) { // j == 0 always Periodic
	      val += dxdz/4.0*volm*tval[k*nzt+nzt-1];
	  }
	  else val += dxdz/4.0*volm*tval[k*nzt+j-1]; // mat_11
	  if((j2+1 == Nz_global) && (Nz_local + lzs == nzt)) {
            val -= dxdz/4.0*volm*tval[k*nzt];
	  }
	  else val -= dxdz/4.0*volm*tval[k*nzt+j+1];    // mat_13
          ind = gindices[(k+lys)*nxzt + (nxt-lxs-1)*nzt+j+lzs];
          VecSetValues( bs, 1, &ind, &val, ADD_VALUES );
        }
      }
    }
  }
  
  // Assemble RHS Vector
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);
  if(fcheck) {
    tmf = MPI_Wtime();
    tmR = tmf - tms;
  }
  //  VecView(bs,PETSC_VIEWER_STDOUT_WORLD);

  
  // Solve the system
  PCType typepc;
  PCGetType(pc,&typepc);
  if(fcheck) {
    Vec rs;
    PetscScalar norm;
    VecDuplicate(xs,&rs);
    MatResidual(MatA,bs,xs,rs);
    VecNorm(rs,NORM_2,&norm);
    output<<"Norm of rhs "<< (double)norm<< ":: "<<yNP<<":"<<xNP<<endl;
    VecDestroy(&rs);
  }
  if(fcheck) tms = MPI_Wtime();
  KSPSetUp(ksp);
  KSPSolve(ksp,bs,xs);
  if(fcheck) {
    tmf = MPI_Wtime();
    tsol = tmf - tms;
  }
  //  MPI_Barrier(MPI_COMM_WORLD);
  // output <<"After solvs"<<yindex<<"(After set)"<<endl;

  if(fcheck) {
    int its;
    Vec rs;
    PetscScalar norm;
    PCType typepc;
    PCGetType(pc,&typepc);
    VecDuplicate(xs,&rs);
    //    VecGhostUpdateBegin(xs,INSERT_VALUES,SCATTER_FORWARD);
    // VecGhostUpdateEnd(xs,INSERT_VALUES,SCATTER_FORWARD);
    MatResidual(MatA,bs,xs,rs);
    //    VecView(rs,PETSC_VIEWER_STDOUT_WORLD);
    VecNorm(rs,NORM_2,&norm);
    KSPGetIterationNumber(ksp,&its);
    output<<"Timing:"<< tmf - tmss<<" Sol: "<<tsol<<":"<<tmG<<","<<tmS<<","<<tmR<<endl;
    output<<"Norm of error "<< (double)norm<< " iterations "<<its<<":"<<typepc<<endl;
    VecDestroy(&rs);
  }
  KSPConvergedReason reason;
  KSPGetConvergedReason( ksp, &reason );
  
  if(reason <= 0) {
    throw BoutException("LaplacePetscAmg failed to converge. Reason %d", reason);
  }
  
  //////////////////////////
  // Copy data into result
  
  //  MPI_Barrier(MPI_COMM_WORLD);
  Field3D result(mesh);
  result.allocate();
  
  for (k=0; k<Ny_local; k++) {
    k2 = k + mystart;
    for(i = 0;i < Nx_local;i++) {
      i2 = i+mxstart;
      for(j= 0;j < Nz_local;j++) {
        ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
        VecGetValues(xs, 1, &ind, &val );
        result(i2, k2, j+mzstart) = val;
      }
    }
  }
  
  // X and Y boundary approximations on guard cells
  if(yProcI == 0) {
    if(ybdcon != 0) { // Boundary value along with y-direction, i.e. x and z index
      k2 = mystart;
      if(ybdcon%3 == 1) {
        for(i = 0;i < Nx_local;i++) {
          i2 = i+mxstart;
          for(j= 0;j < Nz_local;j++) {
            val = -x0(i2, k2-1, j+mzstart)*sqrt(coords->g_22(i2, k2))*coords->dy(i2, k2); 
            result(i2, k2-1, j+mzstart) = val + result(i2,k2,j+mzstart);
          }
	}
      }
      else {
        for(i = 0;i < Nx_local;i++) {
          i2 = i+mxstart;
          for(j= 0;j < Nz_local;j++) {
            result(i2, k2, j+mzstart) = 2.*x0(i2, k2, j+mzstart) - result(i2, k2+1, j+mzstart); 
          }
	}
      }
    }
  }
  if(yProcI == yNP-1) {
    if(ybdcon != 0) { // Boundary value along with y-direction, i.e. x and z index
      k2 = Ny_global-1;
      if(ybdcon%5 == 1) {
        for(i = 0;i < Nx_local;i++) {
          i2 = i+mxstart;
          for(j= 0;j < Nz_local;j++) {
            val = x0(i2, k2+1, j+mzstart)*sqrt(coords->g_22(i2, k2))*coords->dy(i2, k2); 
            result(i2, k2+1, j+mzstart) = val + result(i2,k2,j+mzstart);
          }
	}
      }
      else {
        for(i = 0;i < Nx_local;i++) {
          i2 = i+mxstart;
          for(j= 0;j < Nz_local;j++) {
            result(i2, k2+1, j+mzstart) = 2.*x0(i2, k2+1, j+mzstart) - result(i2, k2, j+mzstart); 
          }
	}
      }
    }
  }
    
  if(xProcI == 0) {
    if(xbdcon != 0) {  // Boundary value along with x-direction, i.e. y and z index
      i2 = mxstart;
      if(xbdcon%3 == 1) {
	for(k = 0;k < Ny_local;k++) {
	  k2 = k + mystart;
          for (j = 0; j < Nz_local; j++) {
            val = -x0(i2-1, k2, j+mzstart)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
            result(i2-1, k2, j+mzstart) = val + result(i2, k2, j+mzstart);
          }
	}
      }
      else {      // Dirichlet boundary condition
	for(k = 0;k < Ny_local;k++) {
	  k2 = k + mystart;
          for (j = 0; j < Nz_local; j++) {
            result(i2-1, k2, j+mzstart) = 2.*x0(i2-1, k2, j+mzstart) - result(i2, k2, j+mzstart); 
          }
        }
      }
    }
  }

  // Outer X boundary
  if(xProcI == xNP -1) {
    if(xbdcon != 0) {  // Boundary value along with x-direction, i.e. y and z index
      i2 = Nx_global - 1;
      if(xbdcon%5 == 1) {
	for(k = 0;k < Ny_local;k++) {
	  k2 = k + mystart;
          for (j = 0; j < Nz_local; j++) {
            val = x0(i2+1, k2, j+mzstart)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
            result(i2+1, k2, j+mzstart) = val + result(i2, k2, j+mzstart);
          }
	}
      }
      else {      // Dirichlet boundary condition
	for(k = 0;k < Ny_local;k++) {
	  k2 = k + mystart;
          for (j = 0; j < Nz_local; j++) {
            result(i2+1, k2, j+mzstart) = 2.*x0(i2+1, k2, j+mzstart) - result(i2, k2, j+mzstart); 
          }
        }
      }
    }
  }

  MatDestroy( &MatA );
  if(diffpre > 0) MatDestroy( &MatP );
  KSPDestroy( &ksp );
  // Set the index of the FieldPerp to be returned
  // MPI_Barrier(MPI_COMM_WORLD);
  return result;
}

