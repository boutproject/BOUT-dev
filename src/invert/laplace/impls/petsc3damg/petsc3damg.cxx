/**************************************************************************
 * Perpendicular Laplacian inversion. 
 *                           Using Algebraic multigrid Solver
 *                              with PETSc library
 *
 * Equation solved is:
 *  d*\nabla^2_\perp x + (1/c1)\nabla_perp c2\cdot\nabla_\perp x + a x = b
 * for 3d formulation and solve by using PETSc
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
 * Ask how to define boundary conditions
 * Total Communictor
 * Matrix generations
 *
 **************************************************************************/
 
#include "bout/mesh.hxx"
#include "boutcomm.hxx"

#include "petsc3damg.hxx"

//BoutReal amgsoltime=0.0,amgsettime=0.0;

LaplacePetsc3DAmg::LaplacePetsc3DAmg(Options *opt) :
  Laplacian(opt),
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

  // For boundary condition for y-direction?? all Dirichlet
  ybdcon = 2;
  yNP = localmesh->getNYPE();
  yProcI = localmesh->getYProcIndex();
  Ny_local = localmesh->yend - localmesh->ystart + 1; // excluding guard cells
  Ny_global = localmesh->GlobalNy - 2*localmesh->ystart; // excluding guard cells
  mystart = localmesh->ystart;
  
  if (mgcount==0) {
    output <<"Ny="<<Ny_global<<"("<<Ny_local<<")"<<endl;
  }
  
  // For boundary condition for x-direction?? 0-Neumann 1-Dirichlet
  xbdcon = 4;
  xNP = localmesh->getNXPE();
  xProcI = localmesh->getXProcIndex();
  Nx_local = localmesh->xend - localmesh->xstart + 1; // excluding guard cells
  Nx_global = localmesh->GlobalNx - 2*localmesh->xstart; // excluding guard cells
  mxstart = localmesh->xstart;
  
  if (mgcount == 0) {
    output <<"Nx="<<Nx_global<<"("<<Nx_local<<")"<<endl;
  }
  zbdcon = 0;
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
  
  // Periodic boundary condition for z-direction
  // Nz_global = 0
  //
  
  nyt = Ny_local;
  nxt = Nx_local;
  nzt = Nz_local;
  ygstart = yProcI*nyt;
  xgstart = xProcI*nxt;
  zgstart = zProcI*nzt;
  lys = 0;
  lxs = 0;
  lzs = 0;
  
  if(yNP > 1) {
    if(ybdcon == 0) {
      lys += 1;
      nyt += 2;
    }
    else {
      if(yProcI > 0) {
        lys += 1;
        nyt += 1;
      }
      if(yProcI < yNP -1) nyt += 1;
    }
  }
  
  if(xNP > 1) {
    if(xbdcon == 0) {
      lxs += 1;
      nxt += 2;
    }
    else {
      if(xProcI > 0) {
        lxs += 1;
        nxt += 1;
      }
      if(xProcI < xNP -1) nxt += 1;
    }
  }
  
  if(zNP > 1) { // periodic Z all
    if(zbdcon == 0) {
      nzt += 2;
      lzs += 1;
    }
    else {
      if(zProcI > 0) {
        lzs += 1;
        nzt += 1;
      }
      if(zProcI < zNP -1) nzt += 1;
    }
  }
  Nlocal = nyt*nxt*nzt;
  Nglobal = Nz_global*Nx_global*Ny_global;

  if (mgcount == 0) {
    output <<"NP="<<tNP<<",("<<yNP<<","<<xNP<<") N="<<Nlocal<<"("<<Nglobal<<")"<<endl;
  }
  int ig,jg,kg,ll,lg,NxzG;
  nxzt = nzt*nxt;
  NxzG = Nz_global*Nx_global;

  gindices = new int[Nlocal];
  for(kg = 0;kg < Ny_local;kg++) {
    for(ig = 0;ig < Nx_local;ig++) {
      for(jg = 0;jg < Nz_local;jg++) {
        ll = (kg+lys)*nxzt + (ig+lxs)*nzt+jg+lzs;
        lg = (kg+ygstart)*NxzG + (ig+xgstart)*Nz_global + zgstart + jg;
        gindices[ll] = lg;
      }
    }
  }
  // Most inner loop (z-direction)
  if(zNP > 1) { // lzs = 1 always periodic
    // for jg = -1 = -lzs
    if(zProcI == 0) {
      for(kg = 0;kg <Ny_local;kg++) {
        for(ig = 0;ig < Nx_local;ig++) {
          ll = (kg+lys)*nxzt + (ig+lxs)*nzt;
          lg = (kg+ygstart)*NxzG + (ig+xgstart+1)*Nz_global - 1;
          gindices[ll] = lg;
        }
      }
    }
    else {
      for(kg = 0;kg < Ny_local;kg++) {
        for(ig = 0;ig < Nx_local;ig++) {
          ll = (kg+lys)*nxzt + (ig+lxs)*nzt;
          lg = (kg+ygstart)*NxzG + (ig+xgstart)*Nz_global + zgstart - 1;
          gindices[ll] = lg;
        }
      }
    }
    // for jg = Nz_local + lzs = nzt -1
    if(zProcI == zNP -1) {
      for(kg = 0;kg < Ny_local;kg++) {
        for(ig = 0;ig < Nx_local;ig++) {
          ll = (kg+lys)*nxzt + (ig+lxs+1)*nzt - 1;
          lg = (kg+ygstart)*NxzG + (ig+xgstart)*Nz_global;
          gindices[ll] = lg;
        }
      }
    }
    else {
      for(kg = 0;kg < Ny_local;kg++) {
        for(ig = 0;ig < Nx_local;ig++) {
          ll = (kg+lys)*nxzt + (ig+lxs+1)*nzt - 1;
          lg = (kg+ygstart)*NxzG + (ig+xgstart)*Nz_global + zgstart + Nz_local;
          gindices[ll] = lg;
	}
      }
    }
  }
  // Ends most inner loop
  
  if(xNP > 1) {
    // For middle loop (x-direction)
    // ig = -1 = -lxs 
    if(xProcI > 0) {
      for(kg = 0;kg < Ny_local;kg++) {
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (kg+lys)*nxzt + jg + lzs;
          lg = (kg+ygstart)*NxzG + (xgstart-1)*Nz_global + jg + zgstart;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) { // Periodic z-direction
	if(zProcI == 0) { // jg = -1 
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt;
            lg = (kg+ygstart)*NxzG + xgstart*Nz_global - 1;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt;
            lg = (kg+ygstart)*NxzG + (xgstart-1)*Nz_global + zgstart - 1;
	    gindices[ll] = lg;
	  }
	}
	if(zProcI == zNP - 1) { // jg = Nz_local
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt + nxt-1;
            lg = (kg+ygstart)*NxzG + (xgstart-1)*Nz_global;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt + nxt-1;
            lg = (kg+ygstart)*NxzG + (xgstart-1)*Nz_global + zgstart + Nz_local;
	    gindices[ll] = lg;
	  }
	}
      }
    }
    else if(lxs == 1) {
      // For lxs = 1 on xProcI = 0 (x-directional periodic BD condition)
      for(kg = 0;kg < Ny_local;kg++) {
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (kg+lys)*nxzt + jg + lzs;
          lg = (kg+ygstart)*NxzG + (Nx_global-1)*Nz_global + jg + zgstart;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) { // Periodic z-direction
	if(zProcI == 0) { // jg = -1
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt;
            lg = (kg+ygstart)*NxzG + Nx_global*Nz_global - 1;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt;
            lg = (kg+ygstart)*NxzG + (Nx_global-1)*Nz_global + zgstart - 1;
	    gindices[ll] = lg;
	  }
	}
	if(zProcI == zNP - 1) { // jg = Nz_local
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt + nxt-1;
            lg = (kg+ygstart)*NxzG + (Nx_global-1)*Nz_global;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
            ll = (kg+lys)*nxzt + nxt-1;
            lg = (kg+ygstart)*NxzG + Nx_global*Nz_global + zgstart;
	    gindices[ll] = lg;
	  }
	}
      }
    }
    
    // ig = Nx_local + lxs = nxt -1
    if(xProcI < xNP - 1) {
      for(kg = 0;kg < Ny_local;kg++) {
	// ig = nxt -1 = Nx_local + 1
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (kg+lys)*nxzt + (nxt-1)*nzt + jg + lzs;
          lg = (kg+lys)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart+jg;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) {  // Periodic z-direction
	if(zProcI == 0) { // jg = -1
          for(kg = 0;kg < Ny_local;kg++) {
  	    ll = (nxt-1)*nxt;
	    lg = (kg+lys)*NxzG + (Nx_local+xgstart+1)*Nz_global - 1;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
  	    ll = (nxt-1)*nxt;
	    lg = (kg+lys)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart - 1;
	    gindices[ll] = lg;
	  }
	}
	if(zProcI == zNP - 1) { // jg = Nz_local
          for(kg = 0;kg < Ny_local;kg++) {
	    ll = nxzt -1;
	    lg = (kg+lys)*NxzG + (Nx_local+xgstart)*Nz_global;
	    gindices[ll] = lg;
	  }
        }
	else {
          for(kg = 0;kg < Ny_local;kg++) {
	    ll = nxzt -1;
	    lg = (kg+lys)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart+Nz_local;
	    gindices[ll] = lg;
	  }
	}
      }
    }
    else if(nxt == Nx_local+lxs+1) {
      // For nxt = Nx_local + lxs + 1 on xProcI = xNP-1 
      // Periodic condition for x-direction
      for(kg = 0;kg < Ny_local;kg++) {
	// ig = nxt -1 = Nx_local + 1
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (kg+lys)*nxzt + (nxt-1)*nzt + jg + lzs;
          lg = (kg+lys)*NxzG + zgstart+jg;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) {  // Periodic z-direction
	if(zProcI == 0) { // jg = -1
          for(kg = 0;kg < Ny_local;kg++) {
	    ll = (nxt-1)*nxt;
	    lg = (kg+lys)*NxzG + Nz_global - 1;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
	    ll = (nxt-1)*nxt;
	    lg = (kg+lys)*NxzG + zgstart - 1;
	    gindices[ll] = lg;
	  }
	}
        if(zProcI == zNP - 1) { // jg = Nz_local
          for(kg = 0;kg < Ny_local;kg++) {
  	    ll = nxzt -1;
	    lg = (kg+lys)*NxzG;
	    gindices[ll] = lg;
	  }
	}
	else {
          for(kg = 0;kg < Ny_local;kg++) {
  	    ll = nxzt -1;
	    lg = (kg+lys)*NxzG + zgstart+Nz_local;
	    gindices[ll] = lg;
	  }
	}
      }      
    }
  }
  // End middle loop
  
  // For outer loop (y-dirction) 
  if(yNP > 1) { // For lys = 1 
    if(yProcI > 0) {
      // kg = -lys = -1
      for(ig = 0;ig < Nx_local;ig++) {
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (ig+lxs)*nzt + jg + lzs;
          lg = (ygstart-1)*NxzG + (ig+xgstart)*Nz_global + jg + zgstart;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) {  // Periodic z-direction
	  // jg = -1 = -lzs
	if(zProcI == 0) { // jg = -1
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (ig+lxs)*nzt;
	    lg = (ygstart-1)*NxzG + (ig+xgstart+1)*Nz_global - 1;
  	    gindices[ll] = lg;
	  }
	}
	else {
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (ig+lxs)*nzt;
	    lg = (ygstart-1)*NxzG + (ig+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	  }
	}
        if(zProcI == zNP - 1) { // jg = Nz_local
	  // jg = Nz_local = nzt - 1;
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (ig+lxs+1)*nzt - 1;
	    lg = (ygstart-1)*NxzG + (ig+xgstart)*Nz_global;	  
	    gindices[ll] = lg;
	  }
        }
	else {
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (ig+lxs+1)*nzt - 1;
	    lg = (ygstart-1)*NxzG + (ig+xgstart)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
	  }
	}
      }
      //Inner middle loop
      if(xNP > 1) {
	if(xProcI > 0) {
	  // ig = -lxs = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = jg + lzs;
            lg = (ygstart-1)*NxzG + (xgstart-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = 0;
    	    if(zProcI == 0) lg = (ygstart-1)*NxzG + xgstart*Nz_global - 1;
	    else lg = (ygstart-1)*NxzG + (xgstart-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart-1)*NxzG + (xgstart-1)*Nz_global;
	    else  lg = (ygstart-1)*NxzG + (xgstart-1)*Nz_global + zgstart + Nz_local;
	    gindices[ll] = lg;
          }
	}
	else if(lxs == 1) { // Periodic condition for x-direction on xProcI = 0
	  // kg = -1 and ig = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = jg + lzs;
            lg = (ygstart-1)*NxzG + (Nx_global-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = 0;
	    if(zProcI == 0) lg = (ygstart-1)*NxzG + Nx_global*Nz_global - 1;
	    else lg = (ygstart-1)*NxzG + (Nx_global-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart-1)*NxzG + (Nx_global-1)*Nz_global;	  
	    else lg = (ygstart-1)*NxzG + (Nx_global-1)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
	// ig = Nx_local + lxs, kg = -1
        if(xProcI < xNP - 1) {
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nxt-1)*nzt + jg + lzs;
            lg = (ygstart-1)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart+jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = 0;
	    if(zProcI == 0) lg = (ygstart-1)*NxzG + Nx_global*Nz_global - 1;
	    else lg =  (ygstart-1)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart-1)*NxzG + (Nx_global-1)*Nz_global;	  
	    else lg = (ygstart-1)*NxzG + Nx_global*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
	else if(nxt == Nx_local+lxs+1) {
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nxt-1)*nzt + jg + lzs;
            lg = (ygstart-1)*NxzG + zgstart + jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = 0;
	    if(zProcI == 0) lg = (ygstart-1)*NxzG + Nz_global - 1;
	    else lg =  (ygstart-1)*NxzG + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart-1)*NxzG;	  
	    else lg = (ygstart-1)*NxzG + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
      }      
    }
    else if(lys == 1) { // yProcI == 0 and periodic BD on y-direction
      // kg = -lys = -1
      for(ig = 0;ig < Nx_local;ig++) {
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (ig+lxs)*nzt + jg + lzs;
          lg = (Ny_global-1)*NxzG + (ig+xgstart)*Nz_global + jg + zgstart;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) {  // Periodic z-direction
	  // jg = -1 = -lzs
	if(zProcI == 0) { // jg = -1
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (ig+lxs)*nzt;
	    lg = (Ny_global-1)*NxzG + (ig+xgstart+1)*Nz_global - 1;
  	    gindices[ll] = lg;
	  }
	}
	else {
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (ig+lxs)*nzt;
	    lg = (Ny_global-1)*NxzG + (ig+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	  }
	}
        if(zProcI == zNP - 1) { // jg = Nz_local
	  // jg = Nz_local = nzt - 1;
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (ig+lxs+1)*nzt - 1;
	    lg = (Ny_global-1)*NxzG + (ig+xgstart)*Nz_global;	  
	    gindices[ll] = lg;
	  }
        }
	else {
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (ig+lxs+1)*nzt - 1;
	    lg = (Ny_global-1)*NxzG + (ig+xgstart)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
	  }
	}
      }
      
      if(xNP > 1) {
	if(xProcI > 0) {
	  // ig = -lxs = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = jg + lzs;
            lg = (Ny_global-1)*NxzG + (xgstart-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = 0;
    	    if(zProcI == 0) lg = (Ny_global-1)*NxzG + xgstart*Nz_global - 1;
	    else lg = (Ny_global-1)*NxzG + (xgstart-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (Ny_global-1)*NxzG + (xgstart-1)*Nz_global;
	    else  lg = (Ny_global-1)*NxzG + (xgstart-1)*Nz_global + zgstart + Nz_local;
	    gindices[ll] = lg;
          }
	}
	else if(lxs == 1) { // Periodic condition for x-direction on xProcI = 0
	  // kg = -1 and ig = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = jg + lzs;
            lg = (Ny_global-1)*NxzG + (Nx_global-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = 0;
	    if(zProcI == 0) lg = (Ny_global-1)*NxzG + Nx_global*Nz_global - 1;
	    else lg = (Ny_global-1)*NxzG + (Nx_global-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (Ny_global-1)*NxzG + (Nx_global-1)*Nz_global;	  
	    else lg = (Ny_global-1)*NxzG + (Nx_global-1)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
	// ig = Nx_local + lxs, kg = -1
        if(xProcI < xNP - 1) {
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nxt-1)*nzt + jg + lzs;
            lg = (Ny_global-1)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart+jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = 0;
	    if(zProcI == 0) lg = (Ny_global-1)*NxzG + Nx_global*Nz_global - 1;
	    else lg =  (Ny_global-1)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (Ny_global-1)*NxzG + (Nx_global-1)*Nz_global;	  
	    else lg = (Ny_global-1)*NxzG + Nx_global*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
	else if(nxt == Nx_local+lxs+1) {
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nxt-1)*nzt + jg + lzs;
            lg = (Ny_global-1)*NxzG + zgstart + jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = 0;
	    if(zProcI == 0) lg = (Ny_global-1)*NxzG + Nz_global - 1;
	    else lg =  (Ny_global-1)*NxzG + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nzt - 1;
	    if(zProcI == zNP - 1) lg = (Ny_global-1)*NxzG;	  
	    else lg = (Ny_global-1)*NxzG + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
      }
    }

    if(yProcI < yNP - 1) {
      // kg = Ny_local + lys  
      for(ig = 0;ig < Nx_local;ig++) {
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (nyt-1)*nxzt + (ig+lxs)*nzt + jg + lzs;
          lg = (ygstart+Ny_local)*NxzG + (ig+xgstart)*Nz_global + jg + zgstart;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) {  // Periodic z-direction
	  // jg = -1 = -lzs
	if(zProcI == 0) { // jg = -1
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (nyt-1)*nxzt + (ig+lxs)*nzt;
	    lg = (ygstart+Ny_local)*NxzG + (ig+xgstart+1)*Nz_global - 1;
  	    gindices[ll] = lg;
	  }
	}
	else {
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (nyt-1)*nxzt + (ig+lxs)*nzt;
	    lg = (ygstart+Ny_local)*NxzG + (ig+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	  }
	}
        if(zProcI == zNP - 1) { // jg = Nz_local
	  // jg = Nz_local = nzt - 1;
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (nyt-1)*nxzt + (ig+lxs+1)*nzt - 1;
	    lg = (ygstart+Ny_local)*NxzG + (ig+xgstart)*Nz_global;	  
	    gindices[ll] = lg;
	  }
        }
	else {
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (nyt-1)*nxzt + (ig+lxs+1)*nzt - 1;
	    lg = (ygstart+Ny_local)*NxzG + (ig+xgstart)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
	  }
	}
      }
      if(xNP > 1) {
	if(xProcI > 0) {
	  // ig = -lxs = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt-1)*nxzt + jg + lzs;
            lg = (ygstart+Ny_local)*NxzG + (xgstart-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt;
    	    if(zProcI == 0) lg = (ygstart+Ny_local)*NxzG + xgstart*Nz_global - 1;
	    else lg = (ygstart+Ny_local)*NxzG + (xgstart-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = (nyt-1)*nxzt + nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart+Ny_local)*NxzG + (xgstart-1)*Nz_global;
	    else  lg = (ygstart+Ny_local)*NxzG + (xgstart-1)*Nz_global + zgstart + Nz_local;
	    gindices[ll] = lg;
          }
	}
	else if(lxs == 1) { // Periodic condition for x-direction on xProcI = 0
	  //  kg = Ny_local + lys and ig = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt -1)*nxzt + jg + lzs;
            lg = (ygstart+Ny_local)*NxzG + (Nx_global-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt;
	    if(zProcI == 0) lg = (ygstart+Ny_local)*NxzG + Nx_global*Nz_global - 1;
	    else lg = (ygstart+Ny_local)*NxzG + (Nx_global-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = (nyt-1)*nxzt + nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart+Ny_local)*NxzG + (Nx_global-1)*Nz_global;	  
	    else lg = (ygstart+Ny_local)*NxzG + (Nx_global-1)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
        // ig = Nx_local + lxs,  kg = Ny_local + lys 
        if(xProcI < xNP - 1) {
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt-1)*nxzt + (nxt-1)*nzt + jg + lzs;
            lg = (ygstart+Ny_local)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart+jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll =  (nyt-1)*nxzt + (nxt-1)*nzt;
	    if(zProcI == 0) lg = (ygstart+Ny_local)*NxzG + Nx_global*Nz_global - 1;
	    else lg =  (ygstart+Ny_local)*NxzG + (Nx_local+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1
	    ll =  nyt*nxzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart+Ny_local)*NxzG + (Nx_global-1)*Nz_global;	  
	    else lg = (ygstart+Ny_local)*NxzG + (Nx_global+xgstart)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
	else if(nxt == Nx_local+lxs+1) { // xProcI = xNP-1 Periodic BD x-direction 
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt-1)*nxzt + (nxt-1)*nzt + jg + lzs;
            lg = (ygstart+Ny_local)*NxzG + zgstart + jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt + (nxt-1)*nzt;
	    if(zProcI == 0) lg = (ygstart+Ny_local)*NxzG + Nz_global - 1;
	    else lg =  (ygstart+Ny_local)*NxzG + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1
	    ll = nyt*nxzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart+Ny_local)*NxzG;	  
	    else lg = (ygstart+Ny_local)*NxzG + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
      }
    }
    else if(nyt == Ny_local + 2) { // yProcI = yNP -1  ////////////
      // kg = Ny_local + lys  
      for(ig = 0;ig < Nx_local;ig++) {
        for(jg = 0;jg < Nz_local;jg++) {
          ll = (nyt-1)*nxzt + (ig+lxs)*nzt + jg + lzs;
          lg = (ig+xgstart)*Nz_global + jg + zgstart;
          gindices[ll] = lg;
        }
      }
      if(zNP > 1) {  // Periodic z-direction
	  // jg = -1 = -lzs
	if(zProcI == 0) { // jg = -1
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (nyt-1)*nxzt + (ig+lxs)*nzt;
	    lg = (ig+xgstart+1)*Nz_global - 1;
  	    gindices[ll] = lg;
	  }
	}
	else {
          for(ig = 0;ig < Nx_local;ig++) {
  	    ll = (nyt-1)*nxzt + (ig+lxs)*nzt;
	    lg = (ig+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	  }
	}
        if(zProcI == zNP - 1) { // jg = Nz_local
	  // jg = Nz_local = nzt - 1;
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (nyt-1)*nxzt + (ig+lxs+1)*nzt - 1;
	    lg = (ig+xgstart)*Nz_global;	  
	    gindices[ll] = lg;
	  }
        }
	else {
          for(ig = 0;ig < Nx_local;ig++) {
	    ll = (nyt-1)*nxzt + (ig+lxs+1)*nzt - 1;
	    lg = (ig+xgstart)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
	  }
	}
      }
      
      if(xNP > 1) {
	if(xProcI > 0) {
	  // ig = -lxs = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt-1)*nxzt + jg + lzs;
            lg = (xgstart-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt;
    	    if(zProcI == 0) lg = (ygstart+Ny_local)*NxzG + xgstart*Nz_global - 1;
	    else lg = (xgstart-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = (nyt-1)*nxzt + nzt - 1;
	    if(zProcI == zNP - 1) lg = (ygstart+Ny_local)*NxzG + (xgstart-1)*Nz_global;
	    else  lg = (xgstart-1)*Nz_global + zgstart + Nz_local;
	    gindices[ll] = lg;
          }
	}
	else if(lxs == 1) { // Periodic condition for x-direction on xProcI = 0
	  //  kg = Ny_local + lys and ig = -1
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt -1)*nxzt + jg + lzs;
            lg = (Nx_global-1)*Nz_global + jg + zgstart;
            gindices[ll] = lg;
          }
          if(zNP > 1) {  // Periodic z-direction
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt;
	    if(zProcI == 0) lg = Nx_global*Nz_global - 1;
	    else lg = (Nx_global-1)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = (nyt-1)*nxzt + nzt - 1;
	    if(zProcI == zNP - 1) lg = (Nx_global-1)*Nz_global;	  
	    else lg = (Nx_global-1)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
        // ig = Nx_local + lxs,  kg = Ny_local + lys 
        if(xProcI < xNP - 1) {
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt-1)*nxzt + (nxt-1)*nzt + jg + lzs;
            lg = (Nx_local+xgstart)*Nz_global + zgstart+jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt + (nxt-1)*nzt;
	    if(zProcI == 0) lg = Nx_global*Nz_global - 1;
	    else lg = (Nx_local+xgstart)*Nz_global + zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll =  nyt*nxzt - 1;
	    if(zProcI == zNP - 1) lg = (Nx_global-1)*Nz_global;	  
	    else lg = (Nx_global+xgstart)*Nz_global + zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
	else if(nxt == Nx_local+lxs+1) { // xProcI = xNP-1 Periodic BD x-direction 
          for(jg = 0;jg < Nz_local;jg++) {
            ll = (nyt-1)*nxzt + (nxt-1)*nzt + jg + lzs;
            lg = zgstart + jg;
            gindices[ll] = lg;
          }
          if(zNP > 1) {
	    // jg = -1 = -lzs
	    ll = (nyt-1)*nxzt + (nxt-1)*nzt;
	    if(zProcI == 0) lg = Nz_global - 1;
	    else lg = zgstart - 1;
  	    gindices[ll] = lg;
	    // jg = Nz_local = nzt - 1;
	    ll = nyt*nxzt - 1;
	    if(zProcI == zNP - 1) lg = 0;	  
	    else lg = zgstart + Nz_local;	  
	    gindices[ll] = lg;
          }
	}
      }
    }
  }
  
  ISLocalToGlobalMappingCreate(BoutComm::get(),1,Nlocal,gindices,PETSC_COPY_VALUES,&mgmapping);

  VecCreateMPI(BoutComm::get(),Nlocal,PETSC_DETERMINE,&xs);
  VecSetLocalToGlobalMapping(xs,mgmapping);
  VecSetFromOptions(xs);
  VecDuplicate(xs,&bs);
  opts->get("diffpre",diffpre,0,true);
  opts->get("elemf",elemf,0,true);     
}

Field3D LaplacePetsc3DAmg::multiplyAx(const Field3D &x) {

  PetscErrorCode ierr;
  int i, k, j,  i2, k2, ind;
  BoutReal val;

  for (k=0; k < Ny_local; k++) {
    k2 = k + mystart;
    for (i=0; i < Nx_local; i++) {
      i2 = i + mxstart;
      for (j=0; j < Nz_local; j++) {
        ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
        val = x(i2, k2, j);
        VecSetValues( xs, 1, &ind, &val, INSERT_VALUES );
      }
    }
  }

  VecAssemblyBegin(xs);
  VecAssemblyEnd(xs);

  generateMatrixA(0);

  ierr = MatMult(MatA, xs, bs);
  if (ierr) throw BoutException("multiplyAx: Petsc error %i", ierr);

  Field3D result;
  result.allocate();
  for(k = 0;k<Ny_local;k++) {
    for(i = 0;i < Nx_local;i++) {
      for(j= 0;j < Nz_local;j++) {
        ind = gindices[(k+lys)*nxzt + (i+lxs)*nzt+j+lzs];
        VecGetValues(bs, 1, &ind, &val );
        result(i+mxstart,k+mystart,j+mzstart) = val;
      }
    }
  }

  return result;
}
