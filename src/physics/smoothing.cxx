/**************************************************************
 * Smoothing operators
 *
 *
 * 2014-10-29 Ben Dudson <bd512@york.ac.uk>
 *    * Moving averaging routines here from Mesh
 * 
 * 2010-05-17 Ben Dudson <bd512@york.ac.uk>
 *    * Added nonlinear filter
 * 
 **************************************************************
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
 **************************************************************/

#include <cmath>

#include "bout/build_config.hxx"

#include <bout/mesh.hxx>
#include <globals.hxx>
#include <smoothing.hxx>
#include <bout_types.hxx>
#include <msg_stack.hxx>

#include <bout/constants.hxx>
#include <bout/rvec.hxx>
#include <utils.hxx>

// Smooth using simple 1-2-1 filter
const Field3D smooth_x(const Field3D &f) {
  TRACE("smooth_x");
  Mesh *mesh = f.getMesh();
  Field3D result{emptyFrom(f)};
  
  // Copy boundary region
  for(int jy=0;jy<mesh->LocalNy;jy++)
    for(int jz=0;jz<mesh->LocalNz;jz++) {
      result(0,jy,jz) = f(0,jy,jz);
      result(mesh->LocalNx-1,jy,jz) = f(mesh->LocalNx-1,jy,jz);
    }

  // Smooth using simple 1-2-1 filter

  for(int jx=1;jx<mesh->LocalNx-1;jx++)
    for(int jy=0;jy<mesh->LocalNy;jy++)
      for(int jz=0;jz<mesh->LocalNz;jz++) {
	result(jx,jy,jz) = 0.5*f(jx,jy,jz) + 0.25*( f(jx-1,jy,jz) + f(jx+1,jy,jz) );
      }

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}


const Field3D smooth_y(const Field3D &f) {
  TRACE("smooth_y");
  Mesh *mesh = f.getMesh();
  Field3D result{emptyFrom(f)};
  
  // Copy boundary region
  for(int jx=0;jx<mesh->LocalNx;jx++)
    for(int jz=0;jz<mesh->LocalNz;jz++) {
      result(jx,0,jz) = f(jx,0,jz);
      result(jx,mesh->LocalNy-1,jz) = f(jx,mesh->LocalNy-1,jz);
    }
  
  // Smooth using simple 1-2-1 filter

  for(int jx=0;jx<mesh->LocalNx;jx++)
    for(int jy=1;jy<mesh->LocalNy-1;jy++)
      for(int jz=0;jz<mesh->LocalNz;jz++) {
	result(jx,jy,jz) = 0.5*f(jx,jy,jz) + 0.25*( f.ydown()(jx,jy-1,jz) + f.yup()(jx,jy+1,jz) );
      }

  // Need to communicate boundaries
  mesh->communicate(result);
  
  return result;
}

/*!

  Issues
  ======
  
  Assumes every processor has the same domain shape

  Will only work if X communicator is constant in Y
  so no processor/branch cuts in X
 */
const Field2D averageX(const Field2D &f) {
  TRACE("averageX(Field2D)");
  Mesh *mesh = f.getMesh();

  int ngx = mesh->LocalNx;
  int ngy = mesh->LocalNy;

  Array<BoutReal> input(ngy), result(ngy);
  
  // Average on this processor
  for(int y=0;y<ngy;y++) {
    input[y] = 0.;
    // Sum values, not including boundaries
    for(int x=mesh->xstart;x<=mesh->xend;x++) {
      input[y] += f(x,y);
    }
    input[y] /= (mesh->xend - mesh->xstart + 1);
  }

  Field2D r{emptyFrom(f)};

  MPI_Comm comm_x = mesh->getXcomm();

  int np;
  MPI_Comm_size(comm_x, &np);
  
  if(np == 1) {
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        r(x,y) = input[y];
  }else {
    MPI_Allreduce(input.begin(), result.begin(), ngy, MPI_DOUBLE, MPI_SUM, comm_x);
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        r(x,y) = result[y] / static_cast<BoutReal>(np);
  }
  
  return r;
}

/*!

  Issues
  ======

  Creates static arrays
  
  Not thread safe
  
  Assumes every processor has the same domain shape
  
  Will only work if X communicator is constant in Y
  so no processor/branch cuts in X
  
 */
const Field3D averageX(const Field3D &f) {
  TRACE("averageX(Field3D)");

  Mesh *mesh = f.getMesh();

  int ngx = mesh->LocalNx;
  int ngy = mesh->LocalNy;
  int ngz = mesh->LocalNz;

  auto input = Matrix<BoutReal>(ngy, ngz);
  auto result = Matrix<BoutReal>(ngy, ngz);

  // Average on this processor
  for(int y=0;y<ngy;y++)
    for(int z=0;z<ngz;z++) {
      input(y, z) = 0.;
      // Sum values, not including boundaries
      for(int x=mesh->xstart;x<=mesh->xend;x++) {
        input(y, z) += f(x, y, z);
      }
      input(y, z) /= (mesh->xend - mesh->xstart + 1);
    }

  Field3D r{emptyFrom(f)};
  
  MPI_Comm comm_x = mesh->getXcomm();
 
  int np;
  MPI_Comm_size(comm_x, &np);
  if(np > 1) {
    MPI_Allreduce(std::begin(input), std::begin(result), ngy*ngz, MPI_DOUBLE, MPI_SUM, comm_x);
    
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          r(x, y, z) = result(y, z) / static_cast<BoutReal>(np);
        }
  }else {
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          r(x, y, z) = input(y, z);
        }
  }
  
  return r;
}

const Field2D averageY(const Field2D &f) {
  TRACE("averageY(Field2D)");

  Mesh *mesh = f.getMesh();
  int ngx = mesh->LocalNx;
  int ngy = mesh->LocalNy;

  Array<BoutReal> input(ngx), result(ngx);
  
  // Average on this processor
  for(int x=0;x<ngx;x++) {
    input[x] = 0.;
    // Sum values, not including boundaries
    for(int y=mesh->ystart;y<=mesh->yend;y++) {
      input[x] += f(x,y);
    }
    input[x] /= (mesh->yend - mesh->ystart + 1);
  }

  Field2D r{emptyFrom(f)};

  /// NOTE: This only works if there are no branch-cuts
  MPI_Comm comm_inner = mesh->getYcomm(0);
  
  int np;
  MPI_Comm_size(comm_inner, &np);
  
  if(np == 1) {
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        r(x,y) = input[x];
  }else {
    MPI_Allreduce(input.begin(), result.begin(), ngx, MPI_DOUBLE, MPI_SUM, comm_inner);
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        r(x,y) = result[x] / static_cast<BoutReal>(np);
  }

  return r;
}

const Field3D averageY(const Field3D &f) {
  TRACE("averageY(Field3D)");

  Mesh *mesh = f.getMesh();

  int ngx = mesh->LocalNx;
  int ngy = mesh->LocalNy;
  int ngz = mesh->LocalNz;

  auto input = Matrix<BoutReal>(ngx, ngz);
  auto result = Matrix<BoutReal>(ngx, ngz);

  // Average on this processor
  for(int x=0;x<ngx;x++)
    for(int z=0;z<ngz;z++) {
      input(x, z) = 0.;
      // Sum values, not including boundaries
      for(int y=mesh->ystart;y<=mesh->yend;y++) {
        input(x, z) += f(x, y, z);
      }
      input(x, z) /= (mesh->yend - mesh->ystart + 1);
    }

  Field3D r{emptyFrom(f)};

  /// NOTE: This only works if there are no branch-cuts
  MPI_Comm comm_inner = mesh->getYcomm(0);
  
  int np;
  MPI_Comm_size(comm_inner, &np);
  if(np > 1) {
    MPI_Allreduce(std::begin(input), std::begin(result), ngx*ngz, MPI_DOUBLE, MPI_SUM, comm_inner);
    
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          r(x, y, z) = result(x, z) / static_cast<BoutReal>(np);
        }
  }else {
    for(int x=0;x<ngx;x++)
      for(int y=0;y<ngy;y++)
        for(int z=0;z<ngz;z++) {
          r(x, y, z) = input(x, z);
        }
  }
  
  return r;
}


BoutReal Average_XY(const Field2D &var) {
  Mesh *mesh = var.getMesh();
  BoutReal Vol_Loc, Vol_Glb;
  int i;
  Field2D result = averageY(var);

  Vol_Loc = 0.;
  Vol_Glb = 0.;

  for (i=mesh->xstart;i<=mesh->xend;i++)
    Vol_Loc +=  result(i,0);

  MPI_Comm comm_x = mesh->getXcomm();

  MPI_Allreduce(&Vol_Loc,&Vol_Glb,1,MPI_DOUBLE,MPI_SUM,comm_x);
  Vol_Glb /= static_cast<BoutReal>(mesh->GlobalNx-2*mesh->xstart);

  return Vol_Glb;
}

BoutReal Vol_Integral(const Field2D &var) {
#if BOUT_USE_METRIC_3D
  AUTO_TRACE();
  throw BoutException("Vol_Intregral currently incompatible with 3D metrics");
#else
  Mesh *mesh = var.getMesh();
  BoutReal Int_Glb;
  Coordinates *metric = var.getCoordinates();

  auto result = metric->J * var * metric->dx * metric->dy;

  Int_Glb = Average_XY(result);
  Int_Glb *= static_cast<BoutReal>(
               (mesh->GlobalNx-2*mesh->xstart)
               * (mesh->GlobalNy-mesh->numberOfYBoundaries()*2*mesh->ystart)) * PI * 2.;

  return Int_Glb;
#endif
}

const Field3D smoothXY(const Field3D &f) {
  Mesh *mesh = f.getMesh();
  Field3D result{emptyFrom(f)};

  for(int x=2;x<mesh->LocalNx-2;x++)
    for(int y=2;y<mesh->LocalNy-2;y++)
      for(int z=0;z<mesh->LocalNz;z++) {
        result(x,y,z) = 0.5*f(x,y,z) + 0.125*( 0.5*f(x+1,y,z) + 0.125*(f(x+2,y,z) + f(x,y,z) + f(x+1,y-1,z) + f(x+1,y+1,z)) +
                                                   0.5*f(x-1,y,z) + 0.125*(f(x,y,z) + f(x-2,y,z) + f(x-1,y-1,z) + f(x-1,y+1,z)) +
                                                   0.5*f(x,y-1,z) + 0.125*(f(x+1,y-1,z) + f(x-1,y-1,z) + f(x,y-2,z) + f(x,y,z)) +
                                                   0.5*f(x,y+1,z) + 0.125*(f(x+1,y+1,z) + f(x-1,y+1,z) + f(x,y,z) + f(x,y+2,z)));
      }
  
  return result;
}


void nl_filter(rvec &f, BoutReal w) {
  for(size_t i=1; i<f.size()-1; i++) {
    
    BoutReal dp = f[i+1] - f[i];
    BoutReal dm = f[i-1] - f[i];
    if(dp*dm > 0.) {
      // Local extrema - adjust
      BoutReal ep, em, e; // Amount to adjust by
      if(fabs(dp) > fabs(dm)) {
	ep = w*0.5*dp;
	em = w*dm;
	e = (fabs(ep) < fabs(em)) ? ep : em; // Pick smallest absolute
	// Adjust
	f[i+1] -= e;
	f[i] += e;
      }else {
	ep = w*0.5*dm;
	em = w*dp;
	e = (fabs(ep) < fabs(em)) ? ep : em; // Pick smallest absolute
	// Adjust
	f[i-1] -= e;
	f[i] += e;
      }
    }
  }
}

const Field3D nl_filter_x(const Field3D &f, BoutReal w) {
  TRACE("nl_filter_x( Field3D )");
  Mesh *mesh = f.getMesh();

  Field3D result{emptyFrom(f)};
  rvec v(mesh->LocalNx);
  
  for (int jy=0;jy<mesh->LocalNy;jy++) {
    for (int jz=0;jz<mesh->LocalNz;jz++) {
      for (int jx=0;jx<mesh->LocalNx;jx++) {
        v[jx] = f(jx,jy,jz);
      }
      nl_filter(v, w);
      for (int jx=0;jx<mesh->LocalNx;jx++) {
         result(jx,jy,jz) = v[jx];
      }
    }
  }
  
  return result;
}

const Field3D nl_filter_y(const Field3D &f, BoutReal w) {
  TRACE("nl_filter_x( Field3D )");

  Mesh *mesh = f.getMesh();

  Field3D result{emptyFrom(f)};

  rvec v(mesh->LocalNy); // Temporary array
  
  // Transform into field-aligned coordinates
  Field3D fs = toFieldAligned(f);

  for (int jx=0;jx<mesh->LocalNx;jx++) {
    for (int jz=0;jz<mesh->LocalNz;jz++) {
      for (int jy=0;jy<mesh->LocalNy;jy++) {
        v[jy] = fs(jx,jy,jz);
      }
      nl_filter(v, w);
      for (int jy=0;jy<mesh->LocalNy;jy++) {
        result(jx,jy,jz) = v[jy];
      }
    }
  }
  
  // Tranform the field back from field aligned coordinates
  return fromFieldAligned(result);
}

const Field3D nl_filter_z(const Field3D &fs, BoutReal w) {
  TRACE("nl_filter_z( Field3D )");

  Mesh *mesh = fs.getMesh();
  Field3D result{emptyFrom(fs)};
  
  rvec v(mesh->LocalNz);
  
  for (int jx=0;jx<mesh->LocalNx;jx++) {
    for (int jy=0;jy<mesh->LocalNy;jy++) {
      for (int jz=0;jz<mesh->LocalNz;jz++) {
        v[jz] = fs(jx,jy,jz);
      }
      nl_filter(v, w);
      for (int jz=0;jz<mesh->LocalNz;jz++) {
        result(jx,jy,jz) = v[jz];
      }
    }
  }
  
  return result;
}

const Field3D nl_filter(const Field3D &f, BoutReal w) {
  /// Perform filtering in Z, Y then X
  Field3D result = nl_filter_x(nl_filter_y(nl_filter_z(f, w), w), w);
  /// Communicate boundaries
  f.getMesh()->communicate(result);
  return result;
}
