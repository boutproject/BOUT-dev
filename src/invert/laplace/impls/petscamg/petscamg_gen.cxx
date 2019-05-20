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

#include "petscamg.hxx"

void LaplacePetscAmg::generateMatrixA(int kflag) {

  TRACE("LaplacePetscAmg::generateMatrixA(int)");
  
  // Set (fine-level) matrix entries
  Coordinates *coords = localmesh->getCoordinates();
  int i,k,i2,k2,k2p,k2m,icc,irow,icol,nn,dz,*dzz,oz,*ozz;
  BoutReal ddx_C,ddz_C,ddx,ddz,dxdz,dxd,dzd,area;
  PetscScalar lval[9],val;
  nn = Nx_local*Nz_local;
  dz = 9;
  oz = 3; // For x-directional decomposition only
  // oz = 5 //For 2-d decompositions
  if(kflag >1) {
    dzz = new int[nn];
    ozz = new int[nn];
    for(i=0;i<Nx_local;i++) {
      for(k = 0;k<Nz_local;k++) {
	icc = i*Nz_local + k;
	if(i == 0) {
	  dzz[icc] = 6;
          if(xProcI == 0) ozz[icc] = 0;
	  else ozz[icc] = 3;
	}
	else if(i == nxt-1) {
	  dzz[icc] = 6;
	  if(xProcI == xNP-1) ozz[icc] = 0;
	  else ozz[icc] = 3;
	}
	else {
	  dzz[icc] = 9;
          ozz[i] = 0;
	}
	if(zNP > 1) {
	  if((k == 0) || (k == Nz_local-1)) {
	    k2 = dzz[icc]/3;
	    dzz[icc] -= k2;
	    ozz[icc] += k2;
	  }
	}
      }
    }
    MatCreateAIJ(commX,nn,nn,PETSC_DETERMINE,PETSC_DETERMINE,dz,dzz,oz,ozz,&MatA );
    delete [] dzz;
    delete [] ozz;
  }
  else { 
    MatCreateAIJ(commX,nn,nn,PETSC_DETERMINE,PETSC_DETERMINE,dz,NULL,oz,NULL,&MatA );
  }
  MatSetFromOptions(MatA);

  
  for (i = 0; i < Nx_local; i++) {
    i2 = i+localmesh->xstart;
BOUT_OMP(parallel default(shared) private(k2))
BOUT_OMP(for)
    for (k = 0; k < Nz_local; k++) {
      k2 = k+mzstart;
      k2p  = (k2+1)%nzt;
      k2m  = (k2+Nz_global-1)%nzt;
      
      ddx_C = (C2(i2+1, yindex, k2) - C2(i2-1, yindex, k2))/2./coords->dx(i2, yindex)/C1(i2, yindex, k2);
      ddz_C = (C2(i2, yindex, k2p) - C2(i2, yindex, k2m))/2./coords->dz/C1(i2, yindex, k2);
      
      ddx = D(i2, yindex, k2)*coords->g11(i2, yindex)/coords->dx(i2, yindex)/coords->dx(i2, yindex); 
               // coefficient of 2nd derivative stencil (x-direction)
      
      ddz = D(i2, yindex, k2)*coords->g33(i2, yindex)/coords->dz/coords->dz; 
              // coefficient of 2nd derivative stencil (z-direction)
      
      dxdz = 2.0*D(i2, yindex, k2)*coords->g13(i2, yindex)/coords->dx(i2, yindex)/coords->dz; 
              // coefficient of mixed derivative stencil (could assume zero, at least initially, 
              // if easier; then check this is true in constructor)
      
      dxd = (D(i2, yindex, k2)*coords->G1(i2, yindex) + coords->g11(i2, yindex)*ddx_C
        + coords->g13(i2, yindex)*ddz_C)/coords->dx(i2, yindex);
             // (could assume zero, at least initially, if easier; then check this is true in constructor)
             // coefficient of 1st derivative stencil (x-direction)
      
      dzd = (D(i2, yindex, k2)*coords->G3(i2, yindex) + coords->g33(i2, yindex)*ddz_C
        + coords->g13(i2, yindex)*ddx_C)/coords->dz;
             // (could assume zero, at least initially, if easier; then check this is true in constructor)
             // coefficient of 1st derivative stencil (z-direction)
      area = coords->dx(i2, yindex)*coords->dz;

      // Put Matrix element with global numbering
      lval[0] = area*dxdz/4.;
      lval[1] = (ddx - dxd/2.)*area;
      lval[2] = -area*dxdz/4.;
      lval[3] = (ddz - dzd/2.)*area;
      lval[4] = (A(i2, yindex, k2) - 2.*(ddx+ddz))*area;
      lval[5] = (ddz + dzd/2.)*area;
      lval[6] = -area*dxdz/4.;
      lval[7] = (ddx+dxd/2.)*area;
      lval[8] = area*dxdz/4.;
      
      icc = (i+lxs)*nzt+k+lzs;
      irow = gindices[icc];
      
      if((xProcI == 0) && (i == 0))  {
        if( inner_boundary_flags & INVERT_AC_GRAD ) {
            // Neumann boundary condition
	  //          output <<"NS"<<irow<<":"<<lval[0]<<","<<lval[1]<<","<<lval[2]<<endl;
          lval[3] += lval[0];
          lval[4] += lval[1];
          lval[5] += lval[2];
          lval[0] = 0.;
          lval[1] = 0.;
          lval[2] = 0.;
        }
        else {
           // Dirichlet boundary condition
	  //          output <<"DS"<<irow<<":"<<lval[0]<<","<<lval[1]<<","<<lval[2]<<endl;
          lval[3] -= lval[0];
          lval[4] -= lval[1];
          lval[5] -= lval[2];
          lval[0] = 0.;
          lval[1] = 0.;
          lval[2] = 0.;
        }
      }
      else {
        icc = (i+lxs-1)*nzt+k+lzs;
        icol = gindices[icc];
        val = lval[1];
	// output<<"1V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
        if((k == 0) && (lzs == 0))  icc = (i+lxs)*nzt-1;
        else icc = (i-1+lxs)*nzt+k+lzs-1;
        icol = gindices[icc];
        val = lval[0];
	// output <<"0V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);

        if((k == Nz_local-1) && (lzs == 0))  icc = (i-1+lxs)*nzt;
        else icc = (i-1+lxs)*nzt+k+lzs+1;
        icol = gindices[icc];
        val = lval[2];
	// output <<"2V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);     
      }
      if((xProcI == xNP-1) && (i == Nx_local-1)) {
        if ( outer_boundary_flags & INVERT_AC_GRAD ) {
             // Neumann boundary condition
	  //          output <<"NF"<<irow<<":"<<lval[6]<<","<<lval[7]<<","<<lval[8]<<endl;
          lval[3] += lval[6];
          lval[4] += lval[7];
          lval[5] += lval[8];
          lval[6] = 0.;
          lval[7] = 0.;
          lval[8] = 0.;
        }
        else {
          // Dirichlet boundary condition
	  //          output <<"DF"<<irow<<":"<<lval[6]<<","<<lval[7]<<","<<lval[8]<<endl;
          lval[3] -= lval[6];
          lval[4] -= lval[7];
          lval[5] -= lval[8];
          lval[6] = 0.;
          lval[7] = 0.;
          lval[8] = 0.;
        }
      }
      else {
        icc = (i+lxs+1)*nzt+k+lzs;
        icol = gindices[icc];
        val = lval[7];
	// output <<"7V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
        if((k == 0) && (lzs == 0))  icc = (i+lxs+2)*nzt-1;
        else icc = (i+1+lxs)*nzt+k+lzs-1;
        icol = gindices[icc];
        val = lval[6];
	// output <<"6V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);

        if((k == Nz_local-1) && (lzs == 0))  icc = (i+1+lxs)*nzt;
        else icc = (i+1+lxs)*nzt+k+lzs+1;
        icol = gindices[icc];
        val = lval[8];
	// output <<"8V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);     
      }
      val = lval[4];
      // output <<"="<<i<<"N="<<k<<"("<<icc<<","<<irow<<")"<<val<<endl;
      MatSetValues(MatA,1,&irow,1,&irow,&val,INSERT_VALUES);
      
      if((k == 0) && (lzs == 0))  icc = (i+lxs+1)*nzt-1;
      else icc = (i+lxs)*nzt+k+lzs-1;
      icol = gindices[icc];
      val = lval[3];
      // output <<"3V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
      MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);

      if((k == Nz_local-1) && (lzs == 0))  icc = (i+lxs)*nzt;
      else icc = (i+lxs)*nzt+k+lzs+1;
      icol = gindices[icc];
      val = lval[5];
      // output <<"5V"<<i<<"N="<<k<<"("<<icc<<","<<icol<<")"<<val<<endl;
      MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);     
    }
  }
  // Assemble Matrix
  MatAssemblyBegin( MatA, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatA, MAT_FINAL_ASSEMBLY );
}

void LaplacePetscAmg::generateMatrixP(int kflag) {

  TRACE("LaplacePetscAmg::generateMatrixP(int)");
  
  // Set (fine-level) matrix entries

  Coordinates *coords = localmesh->getCoordinates();
  int i,k,i2,k2,k2p,k2m,icc,irow,icol,nn,dz,*dzz,oz,*ozz;
  BoutReal ddx_C,ddz_C,ddx,ddz,dxdz,dxd,dzd,area;
  PetscScalar lval[5],val;
  nn = Nx_local*Nz_local;
  dz = 5;
  oz = 1; // For x-directional decomposition only
  // oz = 2 //For 2-d decompositions
  if(kflag >1) {
    dzz = new int[nn];
    ozz = new int[nn];
    for(i=0;i<Nx_local;i++) {
      for(k = 0;k<Nz_local;k++) {
	icc = i*Nz_local + k;
	if(i == 0) {
	  dzz[icc] = 4;
          if(xProcI == 0) ozz[icc] = 0;
	  else ozz[icc] = 1;
	}
	else if(i == nxt-1) {
	  dzz[icc] = 4;
	  if(xProcI == xNP-1) ozz[icc] = 0;
	  else ozz[icc] = 1;
	}
	else {
	  dzz[icc] = 5;
          ozz[i] = 0;
	}
	if(zNP > 1) {
	  if((k == 0) || (k == Nz_local-1)) {
	    dzz[icc] -= 1;
	    ozz[icc] += 1;
	  }
	}
      }
    }
  }
  
  MatCreateAIJ( commX,nn,nn,Nglobal,Nglobal,dz,dzz,oz,ozz, &MatP );                                
  MatSetFromOptions(MatP);

  
  for (i = 0; i < Nx_local; i++) {
    i2 = i+mxstart;
BOUT_OMP(parallel default(shared) private(k2))
BOUT_OMP(for)
    for (k = 0; k < Nz_local; k++) {
      k2 = k+mzstart;
      k2p  = (k2+1)%nzt;
      k2m  = (k2+Nz_global-1)%nzt;
      
      ddx_C = (C2(i2+1, yindex, k2) - C2(i2-1, yindex, k2))/2./coords->dx(i2, yindex)/C1(i2, yindex, k2);
      ddz_C = (C2(i2, yindex, k2p) - C2(i2, yindex, k2m)) /2./coords->dz/C1(i2, yindex, k2);
      
      ddx = D(i2, yindex, k2)*coords->g11(i2, yindex)/coords->dx(i2, yindex)/coords->dx(i2, yindex); 
               // coefficient of 2nd derivative stencil (x-direction)
      
      ddz = D(i2, yindex, k2)*coords->g33(i2, yindex)/coords->dz/coords->dz; 
              // coefficient of 2nd derivative stencil (z-direction)
      
      dxdz = D(i2, yindex, k2)*coords->g13(i2, yindex)/coords->dx(i2, yindex)/coords->dz/2.; 
              // coefficient of mixed derivative stencil (could assume zero, at least initially, 
              // if easier; then check this is true in constructor)
      
      dxd = (D(i2, yindex, k2)*2.*coords->G1(i2, yindex) + coords->g11(i2, yindex)*ddx_C
        + coords->g13(i2, yindex)*ddz_C)/coords->dx(i2, yindex);
             // (could assume zero, at least initially, if easier; then check this is true in constructor)
             // coefficient of 1st derivative stencil (x-direction)
      
      dzd = (D(i2, yindex, k2)*2.*coords->G3(i2, yindex) + coords->g33(i2, yindex)*ddz_C
        + coords->g13(i2, yindex)*ddx_C)/coords->dz;
             // (could assume zero, at least initially, if easier; then check this is true in constructor)
             // coefficient of 1st derivative stencil (z-direction)
      area = coords->dx(i2, yindex)*coords->dz;


      // Put Matrix element with global numbering
      lval[0] = (ddx - dxd/2.)*area;
      lval[1] = (ddz - dzd/2.)*area;
      lval[2] = (A(i2, yindex, k2) - 2.*(ddx+ddz))*area;
      lval[3] = (ddz + dzd/2.)*area;
      lval[4] = (ddx+dxd/2.)*area;
      
      icc = (i+lxs)*nzt+k+lzs;
      irow = gindices[icc];
      if((xProcI == 0) && (i == 0))  {
        if( inner_boundary_flags & INVERT_AC_GRAD ) {
            // Neumann boundary condition
          lval[2] += lval[0];
          lval[0] = 0.;
        }
        else {
           // Dirichlet boundary condition
          lval[2] -= lval[0];
          lval[0] = 0.;
        }
      }
      else {
        icc = (i+lxs-1)*nzt+k+lzs;
        icol = gindices[icc];
        val = lval[0];
        MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
      }
      if((xProcI == xNP-1) && (i == nxt-1)) {
        if ( outer_boundary_flags & INVERT_AC_GRAD ) {
             // Neumann boundary condition
          lval[2] += lval[4];
          lval[4] = 0.;
        }
        else {
          // Dirichlet boundary condition
          lval[2] -= lval[4];
          lval[4] = 0.;
        }
      }
      else {
        icc = (i+lxs+1)*nzt+k+lzs;
        icol = gindices[icc];
        val = lval[4];
        MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
      }
      val = lval[2];
      MatSetValues(MatP,1,&irow,1,&irow,&val,INSERT_VALUES);
      
      if((k == 0) && (lzs == 0))  icc = (i+lxs+1)*nzt-1;
      else icc = (i+lxs)*nzt+k+lzs-1;
      icol = gindices[icc];
      val = lval[1];
      MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);

      if((k == nzt-1) && (lzs == 0))  icc = (i+lxs)*nzt;
      else icc = (i+lxs)*nzt+k+lzs+1;
      icol = gindices[icc];
      val = lval[3];
      MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);     
    }
  }
  // Assemble Matrix
  MatAssemblyBegin( MatP, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatP, MAT_FINAL_ASSEMBLY );
}

