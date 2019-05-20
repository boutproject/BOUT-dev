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

void LaplacePetsc3DAmg::generateMatrixA(int kflag) {

  TRACE("LaplacePetsc3DAmg::generateMatrixA(int)");
  
  // Set (fine-level) matrix entries
  Coordinates *coords = mesh->coordinates();
  int i,j,k,i2,k2,j2,i2m,i2p,j2m,j2p,k2p,k2m,icc,irow,icol,nn,dz,*dzz,oz,*ozz;
  BoutReal ddx_C,ddz_C,ddy_C,ddx,ddz,ddy,dxdz,dxdy,dydz,dxd,dzd,dyd,ddJ,area,volm;
  PetscScalar lval[19],val,dhz,dhx,dhy,gt11,gt12,gt13,gt22,gt23,gt33;
  nn = Nx_local*Nz_local*Ny_local;
  int nxz = Nx_local*Nz_local;
  dz = 19;
  oz = 12; // For x- and y- directional decomposition only
  // oz = 5 //For 3-d decompositions
  if(kflag >1) {
    int bcase,tdz,toz;
    dzz = new int[nn];
    ozz = new int[nn];
    for(k=0;k<Ny_local;k++) {
      for(i=0;i<Nx_local;i++) {
        for(j = 0;j<Nz_local;j++) {
	  icc = k*nxz + i*Nz_local + j;
	  bcase = 0;
	  if(k == 0) {
	    if(lys == 1) bcase += 1;
	    else if((yProcI == 0) && (ybdcon != 0)) bcase += 10;
	  }
	  else if(k == Ny_local -1) {
	    if(nyt == Ny_local+lys+1) bcase += 1;
	    else if((yProcI == yNP-1) && (ybdcon != 0)) bcase += 10;
	  }  
	  if(i == 0) {
	    if(lxs == 1) bcase += 1;
	    else if((xProcI == 0) && (xbdcon != 0)) bcase += 10;
	  }
	  else if(i == Nx_local -1) {
	    if(nxt == Nx_local + lxs + 1) bcase += 1;
	    else if((xProcI == xNP-1) && (xbdcon != 0)) bcase += 10;
	  }
	  if(j == 0) {
            if(lzs == 1) bcase += 1;
	    else if((zProcI == 0) && (zbdcon != 0)) bcase += 10;
	  }
	  else if(j == Nz_local - 1) {
	    if(lzs == Nz_local + lzs + 1) bcase += 1;
	    else if((zProcI == zNP-1) && (zbdcon != 0)) bcase += 10;
	  }
          int odz, tdz;
	  switch(bcase) {
	    case 0: tdz = 19;
	      odz = 0;
	      break;
	    case 1: tdz = 14;
	      odz = 5;
	      break; 
	    case 2: tdz = 10;
	      odz = 9;
	      break; 
	    case 3: tdz = 7;
	      odz = 12;
	      break; 
	    case 10: tdz = 14;
	      odz = 0;
	      break; 
	    case 11: tdz = 10;
	      odz = 4;
	      break; 
	    case 12: tdz = 7;
	      odz = 7;
	      break; 
	    case 20: tdz = 10;
	      odz = 0;
	      break; 
	    case 21: tdz = 7;
	      odz = 3;
	      break; 
	    case 30: tdz = 7;
	      odz = 0;
	      break; 
	    default: output<<"Error in BD at "<<icc<<":"<<bcase<<endl;
	      break;
	  }
	  dzz[icc] = tdz;
	  ozz[icc] = odz;
	}
      }
    }
    MatCreateAIJ(BoutComm::get(),nn,nn,PETSC_DETERMINE,PETSC_DETERMINE,dz,dzz,oz,ozz,&MatA );
    delete [] dzz;
    delete [] ozz;
  }
  else { 
    MatCreateAIJ(BoutComm::get(),nn,nn,PETSC_DETERMINE,PETSC_DETERMINE,dz,NULL,oz,NULL,&MatA );
  }
  MatSetFromOptions(MatA);
  
  dhz = coords->dz;
  int ybds,ybde,xbds,xbde;
  for (k = 0; k < Ny_local; k++) {
    k2 = k + mesh->ystart;
    k2p  = k2+1;            // How to handel when k2p = Ny_global and k2m = -1
    k2m  = k2-1; 
    for (i = 0; i < Nx_local; i++) {
      i2 = i+mesh->xstart;
      i2p  = i2+1;            // How to handel when i2p = Ny_global and i2m = -1
      i2m  = i2-1;
      dhx = coords->dx(i2,k2);
      dhy = coords->dy(i2,k2);
BOUT_OMP(parallel default(shared) private(j,j2,j2p,j2m))
BOUT_OMP(for)
      for (j = 0; j < Nz_local; j++) {
        j2 = j+mzstart;
        j2p  = (j2+1)%nzt;
        j2m  = (j2+Nz_global-1)%nzt;
        gt11 = coords->g11(i2,k2);
	gt12 = coords->g12(i2,k2);
	gt13 = coords->g13(i2,k2);
	gt22 = coords->g22(i2,k2) - 1.0/coords->g22(i2,k2);
	gt23 = coords->g23(i2,k2);
	gt33 = coords->g33(i2,k2);

	ddJ = (coords->J(i2,k2p)/coords->g22(i2,k2p) - coords->J(i2,k2m)/coords->g22(i2,k2m))
	    /2./dhy/coords->J(i2,k2);
        ddx_C = (C2(i2p,k2,j2)-C2(i2m,k2,j2))/2./dhx/C1(i2,k2,j2);
        ddy_C = (C2(i2,k2p,j2)-C2(i2,k2m,j2))/2./dhy/C1(i2,k2,j2);
        ddz_C = (C2(i2,k2,j2p)-C2(i2,k2,j2m))/2./dhz/C1(i2,k2,j2);
      
        ddx = D(i2,k2,j2)*gt11/dhx/dhx; 
        ddy = D(i2,k2,j2)*gt22/dhy/dhy;    
        ddz = D(i2,k2,j2)*gt33/dhz/dhz; 
      
        dxdy = 2.0*D(i2,k2,j2)*gt12/dhx/dhy; 
        dxdz = 2.0*D(i2,k2,j2)*gt13/dhx/dhz; 
        dydz = 2.0*D(i2,k2,j2)*gt23/dhy/dhz; 
       
        dxd = (D(i2,k2,j2)*coords->G1(i2,k2) + gt11*ddx_C + gt12*ddy_C + gt13*ddz_C)/dhx;
        dyd = (D(i2,k2,j2)*(coords->G2(i2,k2) - ddJ)+gt12*ddx_C + gt22*ddy_C + gt23*ddz_C)/dhy;
        dzd = (D(i2,k2,j2)*coords->G3(i2,k2) + gt33*ddz_C +gt13*ddx_C + gt23*ddy_C)/dhz;
        volm = dhx*dhy*dhz;
       
      // Put Matrix element with global numbering
        lval[0] = dydz*volm/4.0;
        lval[1] = volm*dxdy/4.0;
        lval[2] = (ddy - dyd/2.0)*volm;
        lval[3] = -dxdy*volm/4.0;
        lval[4] = -dydz/4.0*volm;
        lval[5] = volm*dxdz/4.0;
        lval[6] = (ddx - dxd/2.0)*volm;
        lval[7] = -volm*dxdz/4.0;
        lval[8] = (ddz - dzd/2.0)*volm;
        lval[9] = (A(i2,k2,j2) - 2.0*(ddx+ddy+ddz))*volm;
        lval[10] = (ddz + dzd/2.0)*volm;
        lval[11] = -volm*dxdz/4.0;
        lval[12] = (ddx+dxd/2.0)*volm;
        lval[13] = volm*dxdz/4.0;
        lval[14] = -dydz*volm/4.0;
        lval[15] = -dxdy/4.0*volm;
        lval[16] = (ddy+dyd/2.0)*volm;
        lval[17] = dxdy/4.0*volm;
        lval[18] = dydz/4.0*volm;
        ybds = 0;
	ybde = 0;
	xbds = 0;
	xbde = 0;
	if((yProcI == 0) && (k == 0)) {
	  if(ybdcon > 0) {
	    if(ybdcon%3 == 1) { // Neumann BD
	      lval[6] += lval[0];
	      lval[8] += lval[1];
	      lval[9] += lval[2];
	      lval[10] += lval[3];
	      lval[12] += lval[4];
	    }
	    else { // Dirichlet BD
	      lval[6] -= lval[0];
	      lval[8] -= lval[1];
	      lval[9] -= lval[2];
	      lval[10] -= lval[3];
	      lval[12] -= lval[4];
	    }
	    lval[0] = 0.0;
	    lval[1] = 0.0;
	    lval[2] = 0.0;
	    lval[3] = 0.0;
	    lval[4] = 0.0;
	    ybds = 1;
	  }
	}
        if((yProcI == yNP-1) && (k == Ny_local-1)) {
	  if(ybdcon > 0) {
	    if(ybdcon%5 == 1) { // Neumann BD
	      lval[6] += lval[14];
	      lval[8] += lval[15];
	      lval[9] += lval[16];
	      lval[10] += lval[17];
	      lval[12] += lval[18];
	    }
	    else { // Dirichlet BD
	      lval[6] -= lval[14];
	      lval[8] -= lval[15];
	      lval[9] -= lval[16];
	      lval[10] -= lval[17];
	      lval[12] -= lval[18];
	    }
	    lval[14] = 0.0;
	    lval[15] = 0.0;
	    lval[16] = 0.0;
	    lval[17] = 0.0;
	    lval[18] = 0.0;
	    ybde = 1;
	  }
	}
        if((xProcI == 0) && (i == 0)) {
	  if(xbdcon > 0) {
	    if(xbdcon%3 == 1) {
	      lval[2] += lval[0];
	      lval[8] += lval[5];
	      lval[9] += lval[6];
	      lval[10] += lval[7];
	      lval[16] += lval[14];
	    }
	    else {
	      lval[2] -= lval[0];
	      lval[8] -= lval[5];
	      lval[9] -= lval[6];
	      lval[10] -= lval[7];
	      lval[16] -= lval[14];
	    }
	    lval[0] = 0.0;
	    lval[5] = 0.0;
	    lval[6] = 0.0;
	    lval[7] = 0.0;
	    lval[14] = 0.0;
	    xbds = 1;
	  }
	}
        if((xProcI == xNP-1) && (i == Nx_local-1)) {
	  if(xbdcon > 0) {
	    if(xbdcon%5 == 1) {
	      lval[2] += lval[4];
	      lval[8] += lval[11];
	      lval[9] += lval[12];
	      lval[10] += lval[13];
	      lval[16] += lval[18];
	    }
	    else {
	      lval[2] -= lval[4];
	      lval[8] -= lval[11];
	      lval[9] -= lval[12];
	      lval[10] -= lval[13];
	      lval[16] -= lval[18];
	    }
	    lval[4] = 0.0;
	    lval[11] = 0.0;
	    lval[12] = 0.0;
	    lval[13] = 0.0;
	    lval[18] = 0.0;
	    xbde = 1;
	  }
	}
	  // Finish set BD condition
        icc = (k+lys)*nxzt + (i+lxs)*nzt+j+lzs;
        irow = gindices[icc];
        icol = irow;  
        val = lval[9];
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
        if((j == 0) && (lzs == 0)) icc = (k+lys)*nxzt + (i+lxs+1)*nzt-1;
	else icc = (k+lys)*nxzt + (i+lxs)*nzt+j+lzs -1;
        icol = gindices[icc];
        val = lval[8];
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
        if((k == Nz_local-1) && (lzs == 0))  icc = (k+lys)*nxzt + (i+lxs)*nzt;
        else icc = (k+lys)*nxzt + (i+lxs)*nzt+j+lzs+1;
        icol = gindices[icc];
        val = lval[10];
        MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
       
        if(xbds == 0) {
	  icc = (k+lys)*nxzt + (i+lxs-1)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[6];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((j == 0) && (lzs == 0)) icc = (k+lys)*nxzt + (i+lxs)*nzt-1;
	  else icc = (k+lys)*nxzt + (i+lxs-1)*nzt+j+lzs -1;
          icol = gindices[icc];
          val = lval[5];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((k == Nz_local-1) && (lzs == 0))  icc = (k+lys)*nxzt + (i+lxs-1)*nzt;
          else icc = (k+lys)*nxzt + (i+lxs-1)*nzt+j+lzs+1;
          icol = gindices[icc];
          val = lval[7];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);     
	}
	if(xbde == 0) {
          icc = (k+lys)*nxzt + (i+lxs+1)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[12];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((j == 0) && (lzs == 0)) icc = (k+lys)*nxzt + (i+lxs+2)*nzt-1;
	  else icc = (k+lys)*nxzt + (i+lxs+1)*nzt+j+lzs -1;
          icol = gindices[icc];
          val = lval[11];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((k == Nz_local-1) && (lzs == 0))  icc = (k+lys)*nxzt + (i+lxs+1)*nzt;
          else icc = (k+lys)*nxzt + (i+lxs+1)*nzt+j+lzs+1;
          icol = gindices[icc];
          val = lval[13];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);     
        }
        if(ybds == 0) {
	  icc = (k+lys-1)*nxzt + (i+lxs)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[2];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((j == 0) && (lzs == 0)) icc = (k+lys-1)*nxzt + (i+lxs+1)*nzt-1;
	  else icc = (k+lys-1)*nxzt + (i+lxs)*nzt+j+lzs -1;
          icol = gindices[icc];
          val = lval[1];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((k == Nz_local-1) && (lzs == 0))  icc = (k+lys-1)*nxzt + (i+lxs)*nzt;
          else icc = (k+lys-1)*nxzt + (i+lxs)*nzt+j+lzs+1;
          icol = gindices[icc];
          val = lval[3];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
	  if(xbds == 0) {
	    icc = (k+lys-1)*nxzt + (i+lxs-1)*nzt+j+lzs;
            icol = gindices[icc];
            val = lval[0];
            MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
	  }
	  if(xbde == 0) {
	    icc = (k+lys-1)*nxzt + (i+lxs+1)*nzt+j+lzs;
            icol = gindices[icc];
            val = lval[4];
            MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
	  }
	}
        if(ybde == 0) {
	  icc = (k+lys+1)*nxzt + (i+lxs)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[16];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((j == 0) && (lzs == 0)) icc = (k+lys+1)*nxzt + (i+lxs+1)*nzt-1;
	  else icc = (k+lys+1)*nxzt + (i+lxs)*nzt+j+lzs -1;
          icol = gindices[icc];
          val = lval[15];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
          if((k == Nz_local-1) && (lzs == 0))  icc = (k+lys+1)*nxzt + (i+lxs)*nzt;
          else icc = (k+lys+1)*nxzt + (i+lxs)*nzt+j+lzs+1;
          icol = gindices[icc];
          val = lval[17];
          MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
	  if(xbds == 0) {
	    icc = (k+lys+1)*nxzt + (i+lxs-1)*nzt+j+lzs;
            icol = gindices[icc];
            val = lval[14];
            MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
	  }
	  if(xbde == 0) {
	    icc = (k+lys+1)*nxzt + (i+lxs+1)*nzt+j+lzs;
            icol = gindices[icc];
            val = lval[18];
            MatSetValues(MatA,1,&irow,1,&icol,&val,INSERT_VALUES);
	  }
	}
      }    
    }
  }
  // Assemble Matrix
  MatAssemblyBegin( MatA, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatA, MAT_FINAL_ASSEMBLY );
}

void LaplacePetsc3DAmg::generateMatrixP(int kflag) {

  TRACE("LaplacePetsc3DAmg::generateMatrixP(int)");
  
  // Set (fine-level) matrix entries

  // Set (fine-level) matrix entries
  Coordinates *coords = mesh->coordinates();
  int i,j,k,i2,k2,j2,i2m,i2p,j2m,j2p,k2p,k2m,icc,irow,icol,nn,dz,*dzz,oz,*ozz;
  BoutReal ddx_C,ddz_C,ddy_C,ddx,ddy,ddz,dxdz,dxdy,dydz,dxd,dzd,dyd,ddJ,area,volm;
  PetscScalar lval[7],val,dhz,dhx,dhy,gt11,gt12,gt13,gt22,gt23,gt33;
  nn = Nx_local*Nz_local*Ny_local;
  int nxz = Nx_local*Nz_local;
  dz = 7;
  oz = 3; // For x- and y- directional decomposition only
  if(kflag >1) {
    int bcase,tdz,toz;
    dzz = new int[nn];
    ozz = new int[nn];
    for(k=0;k<Ny_local;k++) {
      for(i=0;i<Nx_local;i++) {
        for(j = 0;j<Nz_local;j++) {
	  icc = k*nxz + i*Nz_local + j;
	  bcase = 0;
	  if(k == 0) {
	    if(lys == 1) bcase += 1;
	    else if((yProcI == 0) && (ybdcon != 0)) bcase += 10;
	  }
	  else if(k == Ny_local -1) {
	    if(nyt == Ny_local+lys+1) bcase += 1;
	    else if((yProcI == yNP-1) && (ybdcon != 0)) bcase += 10;
	  }  
	  if(i == 0) {
	    if(lxs == 1) bcase += 1;
	    else if((xProcI == 0) && (xbdcon != 0)) bcase += 10;
	  }
	  else if(i == Nx_local -1) {
	    if(nxt == Nx_local + lxs + 1) bcase += 1;
	    else if((xProcI == xNP-1) && (xbdcon != 0)) bcase += 10;
	  }
	  if(j == 0) {
            if(lzs == 1) bcase += 1;
	    else if((zProcI == 0) && (zbdcon != 0)) bcase += 10;
	  }
	  else if(j == Nz_local - 1) {
	    if(lzs == Nz_local + lzs + 1) bcase += 1;
	    else if((zProcI == zNP-1) && (zbdcon != 0)) bcase += 10;
	  }
          int odz, tdz;
	  switch(bcase) {
	    case 0: tdz = 7;
	      odz = 0;
	      break;
	    case 1: tdz = 6;
	      odz = 1;
	      break; 
	    case 2: tdz = 5;
	      odz = 2;
	      break; 
	    case 3: tdz = 4;
	      odz = 3;
	      break; 
	    case 10: tdz = 6;
	      odz = 0;
	      break; 
	    case 11: tdz = 5;
	      odz = 1;
	      break; 
	    case 12: tdz = 4;
	      odz = 2;
	      break; 
	    case 20: tdz = 5;
	      odz = 0;
	      break; 
	    case 21: tdz = 4;
	      odz = 1;
	      break; 
	    case 30: tdz = 4;
	      odz = 0;
	      break; 
	    default: output<<"Error in BD at "<<icc<<":"<<bcase<<endl;
	      break;
	  }
	  dzz[icc] = tdz;
	  ozz[icc] = odz;
	}
      }
    }
    MatCreateAIJ(BoutComm::get(),nn,nn,PETSC_DETERMINE,PETSC_DETERMINE,dz,dzz,oz,ozz,&MatP );
    delete [] dzz;
    delete [] ozz;
  }
  else { 
    MatCreateAIJ(BoutComm::get(),nn,nn,PETSC_DETERMINE,PETSC_DETERMINE,dz,NULL,oz,NULL,&MatP );
  }
  MatSetFromOptions(MatP);
  
  dhz = coords->dz;
  int ybds,ybde,xbds,xbde;
  for (k = 0; k < Ny_local; k++) {
    k2 = k + mesh->ystart;
    k2p  = k2+1;            // How to handel when k2p = Ny_global and k2m = -1
    k2m  = k2-1; 
    for (i = 0; i < Nx_local; i++) {
      i2 = i+mesh->xstart;
      i2p  = i2+1;            // How to handel when i2p = Ny_global and i2m = -1
      i2m  = i2-1;
      dhx = coords->dx(i2,k2);
      dhy = coords->dy(i2,k2);
BOUT_OMP(parallel default(shared) private(j,j2,j2p,j2m))
BOUT_OMP(for)
      for (j = 0; j < Nz_local; j++) {
        j2 = j+mzstart;
        j2p  = (j2+1)%nzt;
        j2m  = (j2+Nz_global-1)%nzt;

	gt11 = coords->g11(i2,k2);
	gt12 = coords->g12(i2,k2);
	gt13 = coords->g13(i2,k2);
	gt22 = coords->g22(i2,k2) - 1.0/coords->g22(i2,k2);
	gt23 = coords->g23(i2,k2);
	gt33 = coords->g33(i2,k2);

	ddJ = (coords->J(i2,k2p)/coords->g22(i2,k2p) - coords->J(i2,k2m)/coords->g22(i2,k2m))
	    /2./dhy/coords->J(i2,k2);
        ddx_C = (C2(i2p,k2,j2)-C2(i2m,k2,j2))/2./dhx/C1(i2,k2,j2);
        ddy_C = (C2(i2,k2p,j2)-C2(i2,k2m,j2))/2./dhy/C1(i2,k2,j2);
        ddz_C = (C2(i2,k2,j2p)-C2(i2,k2,j2m))/2./dhz/C1(i2,k2,j2);
      
        ddx = D(i2,k2,j2)*gt11/dhx/dhx; 
        ddy = D(i2,k2,j2)*gt22/dhy/dhy;    
        ddz = D(i2,k2,j2)*gt33/dhz/dhz; 
      
        dxdy = 2.0*D(i2,k2,j2)*gt12/dhx/dhy; 
        dxdz = 2.0*D(i2,k2,j2)*gt13/dhx/dhz; 
        dydz = 2.0*D(i2,k2,j2)*gt23/dhy/dhz; 
       
        dxd = (D(i2,k2,j2)*coords->G1(i2,k2) + gt11*ddx_C + gt12*ddy_C + gt13*ddz_C)/dhx;
        dyd = (D(i2,k2,j2)*(coords->G2(i2,k2) - ddJ)+gt12*ddx_C + gt22*ddy_C + gt23*ddz_C)/dhy;
        dzd = (D(i2,k2,j2)*coords->G3(i2,k2) + gt33*ddz_C +gt13*ddx_C + gt23*ddy_C)/dhz;
        volm = dhx*dhy*dhz;
       
      // Put Matrix element with global numbering
        lval[0] = (ddy + dyd/2.0)*volm;

        lval[1] = (ddx - dxd/2.0)*volm;
        lval[2] = (ddz - dzd/2.0)*volm;
        lval[3] = (A(i2,k2,j2) - 2.0*(ddx+ddy+ddz))*volm;
        lval[4] = (ddz + dzd/2.0)*volm;
        lval[5] = (ddx+dxd/2.0)*volm;

        lval[6] = (ddy+dyd/2.0)*volm;
	ybds = 0;
	ybde = 0;
	xbds = 0;
	xbde = 0;
	if((yProcI == 0) && (k == 0)) {
	  if(ybdcon > 0) {
	    if(ybdcon%3 == 1) { // Neumann BD
	      lval[3] += lval[0];
	    }
	    else { // Dirichlet BD
	      lval[3] -= lval[0];
	    }
	    lval[0] = 0.0;
	    ybds = 1;
	  }
	}
        if((yProcI == yNP-1) && (k == Ny_local-1)) {
	  if(ybdcon > 0) {
	    if(ybdcon%5 == 1) { // Neumann BD
	      lval[3] += lval[6];
	    }
	    else { // Dirichlet BD
	      lval[3] -= lval[6];
	    }
	    lval[6] = 0.0;
	    ybde = 1;
	  }
	}
        if((xProcI == 0) && (i == 0)) {
	  if(xbdcon > 0) {
	    if(xbdcon%3 == 1) {
	      lval[3] += lval[1];
	    }
	    else {
	      lval[3] -= lval[1];
	    }
	    lval[1] = 0.0;
	    xbds = 1;
	  }
	}
        if((xProcI == xNP-1) && (i == Nx_local-1)) {
	  if(xbdcon > 0) {
	    if(xbdcon%5 == 1) {
	      lval[3] += lval[5];
	    }
	    else {
	      lval[3] -= lval[5];
	    }
	    lval[5] = 0.0;
	    xbde = 1;
	  }
	}
	  // Finish set BD condition
        icc = (k+lys)*nxzt + (i+lxs)*nzt+j+lzs;
        irow = gindices[icc];
        icol = irow;  
        val = lval[3];
        MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
        if((j == 0) && (lzs == 0)) icc = (k+lys)*nxzt + (i+lxs+1)*nzt-1;
	else icc = (k+lys)*nxzt + (i+lxs)*nzt+j+lzs -1;
        icol = gindices[icc];
        val = lval[2];
        MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
        if((k == Nz_local-1) && (lzs == 0))  icc = (k+lys)*nxzt + (i+lxs)*nzt;
        else icc = (k+lys)*nxzt + (i+lxs)*nzt+j+lzs+1;
        icol = gindices[icc];
        val = lval[4];
        MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
       
        if(xbds == 0) {
	  icc = (k+lys)*nxzt + (i+lxs-1)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[1];
          MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
	}
	if(xbde == 0) {
          icc = (k+lys)*nxzt + (i+lxs+1)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[5];
          MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
        }
        if(ybds == 0) {
	  icc = (k+lys-1)*nxzt + (i+lxs)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[0];
          MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
	}
        if(ybde == 0) {
	  icc = (k+lys+1)*nxzt + (i+lxs)*nzt+j+lzs;
          icol = gindices[icc];
          val = lval[6];
          MatSetValues(MatP,1,&irow,1,&icol,&val,INSERT_VALUES);
	}
      }
    }
  }
  // Assemble Matrix
  MatAssemblyBegin( MatP, MAT_FINAL_ASSEMBLY );
  MatAssemblyEnd( MatP, MAT_FINAL_ASSEMBLY );
}

