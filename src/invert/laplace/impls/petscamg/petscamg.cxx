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

//BoutReal amgsoltime=0.0,amgsettime=0.0;

LaplacePetscAmg::LaplacePetscAmg(Options *opt, const CELL_LOC loc, Mesh *mesh_in) :
  Laplacian(opt, loc, mesh_in),
  A(0.0), C1(1.0), C2(1.0), D(1.0) {

  TRACE("LaplacePetscAmg::LaplacePetscAmg(Options *opt)");

  // PetscInitialize(&argc,&args,(char*)0,help);
  // Get Options in Laplace Section
  //  if (!opt) opts = Options::getRoot()->getSection("petscamg");
  //  else opts=opt;
  opts = Options::getRoot()->getSection("petscamg");
  opts->get("rtol",rtol,pow(10.0,-11),true);
  opts->get("atol",atol,pow(10.0,-25),true);
  opts->get("dtol",dtol,pow(10.0,5),true);
  opts->get("maxits",maxits,200000,true);
  opts->get("rightpre",rightpre,0,true);
  opts->get("smtype",mgsm,1,true);
  opts->get("jacomega",omega,0.8,true);
  opts->get("checking",fcheck,1,true);
  opts->get("multigridlevel",mglevel,6,true);
  opts->get("solvertype",soltype,"gmres");

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
  
  commX = localmesh->getXcomm();    // Where to get This mesh //

  MPI_Comm_size(commX,&xNP);
  MPI_Comm_rank(commX,&xProcI); 
  Nx_local = localmesh->xend - localmesh->xstart + 1; // excluding guard cells
  Nx_global = localmesh->GlobalNx - 2*localmesh->xstart; // excluding guard cells
  mxstart = localmesh->xstart;
  
  if (mgcount == 0) {
    output <<"Nx="<<Nx_global<<"("<<Nx_local<<")"<<endl;
  }
  zNP = 1;
  zProcI = 0;
  Nz_global = localmesh->GlobalNz;
  Nz_local = Nz_global;
  mzstart = 0; //  
  // No parallelization in z-direction (for now)
  // 
  //else {
  //  Nz_local = localmesh->zend - localmesh->zstart + 1; // excluding guard cells
  //  Nz_global = localmesh->GlobalNz - 2*localmesh->zstart; // excluding guard cells
  // }
  if (mgcount==0) {
    output <<"Nz="<<Nz_global<<"("<<Nz_local<<")"<<endl;
  }

  int ig,jg,ll,lg;
  // Periodic boundary condition for z-direction
  // Nz_global = 0
  //
  
  nxt = Nx_local;
  nzt = Nz_local;
  xgstart = xProcI*nxt;
  zgstart = zProcI*nzt;
  lxs = 0;
  lzs = 0;
  if(xNP > 1) {
    if(xProcI > 0) {
      lxs += 1;
      nxt += 1;
    }
    if(xProcI < xNP -1) nxt += 1;
  }
  if(zNP > 1) {
    nzt += 2;
    lzs += 1;
  }
  Nlocal = nxt*nzt;
  Nglobal = Nz_global*Nx_global;

  if (mgcount == 0) {
    output <<"NP="<<xNP<<"N="<<Nlocal<<"("<<Nglobal<<")"<<endl;
  }
  gindices = new int[Nlocal];

  for(ig = 0;ig < Nx_local;ig++) {
    for(jg = 0;jg < Nz_local;jg++) {
      ll = (ig+lxs)*nzt+jg+lzs;
      lg = (ig+xgstart)*Nz_global + zgstart + jg;
      gindices[ll] = lg;
    }
  }
  if(zNP > 1) { // lzs = 1 
    if(zProcI == 0) { 
      for(ig = 0;ig < Nx_local;ig++) {
        ll = (ig+lxs)*nzt;
        lg = (ig+xgstart+1)*Nz_global - 1;
        gindices[ll] = lg;
      }
    }
    else {
      for(ig = 0;ig < Nx_local;ig++) {
        ll = (ig+lxs)*nzt;
        lg = (ig+xgstart)*Nz_global + zgstart - 1;
        gindices[ll] = lg;
      }
    }
    if(zProcI == zNP -1) {
      for(ig = 0;ig < Nx_local;ig++) {
        ll = (ig+lxs+1)*nzt - 1;
        lg = (ig+xgstart)*Nz_global;
        gindices[ll] = lg;
      }
    }
    else {
      for(ig = 0;ig < Nx_local;ig++) {
        ll = (ig+lxs+1)*nzt - 1;
        lg = (ig+xgstart)*Nz_global + zgstart + Nz_local;
        gindices[ll] = lg;
      }
    }
  }
  if(xNP > 1) { // For lzx = 1 
    if(xProcI > 0) { 
      for(jg = 0;jg < Nz_local;jg++) {
        ll = jg + lzs;
        lg = (xgstart-1)*Nz_global + jg + zgstart;
        gindices[ll] = lg;
      }
      if(zNP > 1) {
	gindices[0] = (xgstart-1)*Nz_global + zgstart - 1;
	gindices[nxt-1] = (xgstart-1)*Nz_global + zgstart+Nz_local;
      }
    }
    if(xProcI < xNP - 1) {
      for(jg = 0;jg < Nz_local;jg++) {
        ll = (nxt-1)*nzt + jg + lzs;
        lg = (Nx_local+xgstart)*Nz_global + zgstart+jg;
        gindices[ll] = lg;
      }
      if(zNP > 1) {
	gindices[(nxt-1)*nzt] = (Nx_local+xgstart)*Nz_global + zgstart - 1;
	gindices[nxt*nzt-1] = (Nx_local+xgstart)*Nz_global + zgstart+Nz_local;
      }
    }
  }
  ISLocalToGlobalMappingCreate(commX,1,Nlocal,gindices,PETSC_COPY_VALUES,&mgmapping);

  VecCreateMPI(commX,Nx_local*Nz_local,PETSC_DETERMINE,&xs);
  VecSetLocalToGlobalMapping(xs,mgmapping);
  VecSetFromOptions(xs);
  VecDuplicate(xs,&bs);
  opts->get("diffpre",diffpre,0,true);
  opts->get("elemf",elemf,0,true);     
}

FieldPerp LaplacePetscAmg::multiplyAx(const FieldPerp &x) {

  PetscErrorCode ierr;
  int i, k, i2, ind;
  BoutReal val;

  yindex = x.getIndex();

  for (i=0; i < Nx_local; i++) {
    i2 = i + mxstart;
    for (k=0; k < Nz_local; k++) {
      ind = gindices[(i+lxs)*nzt+k+lzs];
      val = x(i2, k);
      VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
    }
  }

  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);

  generateMatrixA(0);

  ierr = MatMult(MatA, xs, bs);
  if (ierr) throw BoutException("multiplyAx: Petsc error %i", ierr);

  FieldPerp result;
  result.allocate();
  result.setIndex(yindex);
  for(i = 0;i < Nx_local;i++) {
    for(k= 0;k < Nz_local;k++) {
      ind = gindices[(i+lxs)*nzt+k+lzs];
      VecGetValues(bs, 1, &ind, &val );
      result(i+mxstart,k+mzstart) = val;
    }
  }

  return result;
}
