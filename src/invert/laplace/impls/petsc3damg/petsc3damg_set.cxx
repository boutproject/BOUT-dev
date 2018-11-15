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

#include "petsc3damg.hxx"
#include <boutcomm.hxx>

void LaplacePetsc3DAmg::settingSolver(int kflag){
  //  Timer timer("invert");
  TRACE("LaplacePetsc3DAmg::settingSolver(int)");
  
  if(!opts) {
    // If no options supplied, use default
    opts = Options::getRoot()->getSection("petscamg");
  }

  //////////////////////////////////////////////////
  // Set up KSP
  
  // Declare KSP Context 
  KSPCreate(BoutComm::get(), &ksp ); 
  KSPGetPC(ksp, &pc);
  
  // Configure Linear Solver
  if(kflag == 0) KSPSetOperators( ksp, MatA, MatA );
  else KSPSetOperators( ksp, MatA, MatP );
   // Convergence Parameters. Solution is considered converged if |r_k| < max( rtol * |b| , atol )
    // where r_k = b - Ax_k. The solution is considered diverged if |r_k| > dtol * |b|.
  string solt;
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
        int *blks = new int[Nx_global];
        for(int i=0;i<Nx_global;i++) blks[i] = Nz_local;
        PCBJacobiSetTotalBlocks(pc,Nx_global,blks);
        delete [] blks;
      }
    }
  }
  else {
    KSPSetType( ksp, KSPGMRES );
    KSPSetInitialGuessNonzero(ksp, (PetscBool) true );
    if(soltype == "gmres") {
      PCSetType(pc,PCBJACOBI);
      int *blks = new int[Nx_global*2];
      for(int i=0;i<Nx_global;i++) {
	blks[2*i] = Nz_local/2;
	blks[2*i+1] = Nz_local - Nz_local/2;
      }
      PCBJacobiSetTotalBlocks(pc,Nx_global*2,blks);
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
      /*
        else if(soltype == "gamggeo" ) {
	  PCSetType(pc,PCGAMG);
          PCGAMGSetType(pc,PCGAMGGEO);
        }
      */
        if(soltype == "gamgag") {
	  PCGAMGSetNSmooths(pc,2);
        }
        else if(soltype == "gamgag0") {
  	  PCGAMGSetNSmooths(pc,3);
        }
        else if(soltype == "gamgag1") {
  	  PCGAMGSetNSmooths(pc,4);
        }
      /*
        else if(soltype == "gamgcla") {
          PCSetType(pc,PCGAMG);
          PCGAMGSetType(pc,PCGAMGCLASSICAL);
	//	PCGAMGSetNSmooths(pc,1);
      }
      */
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
  Coordinates *coords = mesh->coordinates();
  BoutReal tmss, tms, tmf, tmG, tmS, tmR, tsol, tmB;
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
  
  int ind,i2,i,j,j2,k,k2,nxzt;
  PetscScalar val,area;
  nxzt = nzt*nxt;
  if(fcheck) tms = MPI_Wtime();
  if ( global_flags & INVERT_START_NEW ) {
    // set initial guess to zero
    for (k=0; k < Ny_local; k++) {
      k2 = k + mystart;
      for (i=0; i < Nx_local; i++) {
        i2 = i + mxstart;
        area = coords->dx(i2, k2)*coords->dz;
        for (j=0; j < Nz_local; j++) {
          ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
          val = 0.;
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
          val = rhs(i2, k2, j+mzstart)*area;
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
        area = coords->dx(i2, k2)*coords->dz;
        for (j=0; j < Nz_local; j++) {
          ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
          val = x0(i2, k2, j+mzstart);
          VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
        
          val = rhs(i2, k2, j+mzstart)*area;
          VecSetValues( bs, 1, &ind, &val, INSERT_VALUES );
        }
      }
    }
  }
  VecAssemblyBegin(bs);
  VecAssemblyEnd(bs);
  //  VecView(bs,PETSC_VIEWER_STDOUT_WORLD);
  // Assemble Trial Solution Vector
  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);
  // For the boundary conditions 
  // MPI_Barrier(MPI_COMM_WORLD);
  BoutReal tval[nzt],dval,ddx_C,ddz_C;
  if (mesh->firstX()) {
    for (k=0; k<Ny_local; k++) {
      k2 = k + mystart;

      i2 = mxstart;
      if ( inner_boundary_flags & INVERT_AC_GRAD ) {
        // Neumann boundary condition
        //      if ( inner_boundary_flags & INVERT_SET ) {
          // Set the value of guard cells specify gradient to set at inner boundary
          // tval = df/dn = (v_ghost - v_in)/distance 
          for (j = 0; j < nzt; j++) {
            tval[j] = -x0(i2-1, k2, j+mzstart-lzs)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
          }
          //   }
          //else {
          // zero gradient inner boundary condition
          //for (int j = 0; j<nzt; j++) {
            // set inner guard cells
          //  tval[j] = 0.0;
          // }
          // }
      }
      else {      // Dirichlet boundary condition
        // if ( inner_boundary_flags & INVERT_SET ) {
          // guard cells of x0 specify value to set at inner boundary
          // tval = f = (v_ghost + v_in)/2.0
          for (int j = 0; j < nzt; j++) {
            tval[j] = 2.*x0(i2-1, k2, j+mzstart-lzs); 
          // this is the value to set at the inner boundary
          }
          // }
          //else {
          // zero value inner boundary condition
          // for (int j=0; j < nzt; j++) {
            // set inner guard cells
          //  tval[j] = 0.;
          // }
          // }
      }
      for(j = 0;j < Nz_local;j++) {
        j2 = j + mzstart;
        area = coords->dx(i2, k2)*coords->dz;
        ddx_C = (C2(i2+1, k2, j2) - C2(i2-1, k2, j2))/2./coords->dx(i2, k2)/C1(i2, k2, j2);
        ddz_C = (C2(i2, k2, (j2+1)%nzt) - C2(i2, k2, (j2-1+nzt)%nzt)) /2./coords->dz/C1(i2, k2, j2);
        dval = D(i2, k2, j2)*coords->g11(i2, k2)/coords->dx(i2, k2)/coords->dx(i2, k2);
        dval -= (D(i2, k2, j2)*2.*coords->G1(i2, k2) + coords->g11(i2, k2)*ddx_C
                 + coords->g13(i2, k2)*ddz_C)/coords->dx(i2, k2)/2.0;
        dval *= area;
        val = -tval[j+lzs]*dval;
        dval = D(i2, k2, j2)*coords->g13(i2, k2)/coords->dx(i2, k2)/coords->dz/8.;
        dval *= area;
        if(lzs == 0 && j == 0) val -= dval*tval[nzt-1];
        else val -= dval*tval[j-1];
        if(lzs == 0 && j == nzt-1) val += dval*tval[0];
        else val += dval*tval[j+1];
        ind = gindices[(k+lys)*nxzt + (0+lxs)*nzt+j+lzs];
        VecSetValues( bs, 1, &ind, &val, ADD_VALUES );
      }  
    }
  }
  if (mesh->lastX()) {
    for (k=0; k<Ny_local; k++) {
      k2 = k + mystart;

      i2 = mesh->xend;
      if ( outer_boundary_flags & INVERT_AC_GRAD ) {
        // Neumann boundary condition
        //      if ( inner_boundary_flags & INVERT_SET ) {
          // guard cells of x0 specify gradient to set at outer boundary
          // tval = df/dn = (v_ghost - v_in)/distance 
          for (j= 0; j < nzt; j++) {
            tval[j] = x0(i2+1, k2, j+mzstart-lzs)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
          // this is the value to set the gradient to at the outer boundary
          }
          //}
          //else {
          // zero gradient outer boundary condition
          //for (j=0; j<nzt; j++) {
            // set outer guard cells
            //tval[j] = 0.;
          //}
          //}
      }
      else {
        // Dirichlet boundary condition
        // if ( outer_boundary_flags & INVERT_SET ) {
          // guard cells of x0 specify value to set at outer boundary
          for (j=0; j< nzt; j++) {
            tval[j]=2.*x0(i2+1, k2, j+mzstart-lzs); 
            // this is the value to set at the outer boundary
          }
          //}
          //else {
          // zero value inner boundary condition
          //for (int j=0; j<nzt; j++) {
            // set outer guard cells
            //tval[j] = 0.;
          //}
          // }
      }
      for(j = 0;j < Nz_local;j++) {
        j2 = j+mzstart;
        area = coords->dx(i2, k2)*coords->dz;
        ddx_C = (C2(i2+1, k2, j2) - C2(i2-1, k2, j2))/2./coords->dx(i2, k2)/C1(i2, k2, j2);
        ddz_C = (C2(i2, k2, (j2+1)%nzt) - C2(i2, k2, (j2-1+nzt)%nzt)) /2./coords->dz/C1(i2, k2, j2);
        dval = D(i2, k2, j2)*coords->g11(i2, k2)/coords->dx(i2, k2)/coords->dx(i2, k2);
        dval += (D(i2, k2, j2)*2.*coords->G1(i2, k2) + coords->g11(i2, k2)*ddx_C
                     + coords->g13(i2, k2)*ddz_C)/coords->dx(i2, k2)/2.0;
        dval *= area;
        val = -tval[j+lzs]*dval;
        dval = D(i2, k2, j2)*coords->g13(i2, k2)/coords->dx(i2, k2)/coords->dz/8.;
        dval *= area;
        if(lzs == 0 && j == 0) val += dval*tval[nzt-1];
        else val += dval*tval[j-1];
        if(lzs == 0 && j == nzt-1) val -= dval*tval[0];
        else val -= dval*tval[j+1];
        ind = gindices[(k+lys)*nxzt + (nxt-1)*nzt+j+lzs];
        VecSetValues( bs, 1, &ind, &val, ADD_VALUES );
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
    int its;
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
  
  // Inner X boundary approximations on guard cells
  // Need to modify
  
  if(mesh->firstX()) {
    for (k=0; k<Ny_local; k++) {
      k2 = k + mystart;

      i2 = mxstart;
      if ( inner_boundary_flags & INVERT_AC_GRAD ) {
        // Set THE VALUE guard cells specify gradient to set at inner boundary
        // tval = df/dn = (v_ghost - v_in)/distance 
        for (j = 0; j < Nz_local; j++) {
          val = -x0(i2-1, k2, j+mzstart)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
          result(i2-1, k2, j+mzstart) = val + result(i2, k2, j+mzstart);
        }
      }
      else {      // Dirichlet boundary condition
          // guard cells of x0 specify value to set at inner boundary
          // tval = f = (v_ghost + v_in)/2.0
        for (int j = 0; j < Nz_local; j++) {
            result(i2-1, k2, j+mzstart) = 2.*x0(i2-1, k2, j+mzstart) - result(i2, k2, j+mzstart); 
          // this is the value to set at the inner boundary
        }
      }
      output<<"Put results first "<<reason<<endl;
    }
  }

  // Outer X boundary
  if (mesh->lastX()) {
    for (k=0; k<Ny_local; k++) {
      k2 = k + mystart;

      i2 = mesh->xend;
      if ( outer_boundary_flags & INVERT_AC_GRAD ) {
        // Neumann boundary condition
          // guard cells of x0 specify gradient to set at outer boundary
          // tval = df/dn = (v_ghost - v_in)/distance 
        for (j = 0; j < Nz_local; j++) {
          val = x0(i2+1, k2, j+mzstart)*sqrt(coords->g_11(i2, k2))*coords->dx(i2, k2); 
          result(i2+1, k2, j+mzstart) = val + result(i2, k2, j+mzstart);
          // this is the value to set the gradient to at the outer boundary
        }
      }
      else {
        // Dirichlet boundary condition
        // guard cells of x0 specify value to set at outer boundary
        for (j = 0; j < Nz_local; j++) {
          result(i2+1, k2, j+mzstart) = 2.*x0(i2+1, k2, j+mzstart) - result(i2, k2, j+mzstart); 
          // this is the value to set at the outer boundary
        }
      }
    }
  }
  // MPI_Barrier(MPI_COMM_WORLD);

  MatDestroy( &MatA );
  if(diffpre > 0) MatDestroy( &MatP );
  KSPDestroy( &ksp );
  // Set the index of the FieldPerp to be returned
  // MPI_Barrier(MPI_COMM_WORLD);
  return result;
}



