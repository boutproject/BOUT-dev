/**************************************************************
 * Smoothing operators
 *
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

#include <math.h>

#include <globals.hxx>
#include <smoothing.hxx>
#include <bout_types.hxx>
#include <msg_stack.hxx>

// Smooth using simple 1-2-1 filter
const Field3D smooth_x(const Field3D &f, bool BoutRealspace) {
  Field3D fs, result;

  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f; 

  result.allocate();
  
  // Copy boundary region
  for(int jy=0;jy<mesh->ngy;jy++)
    for(int jz=0;jz<mesh->ngz;jz++) {
      result[0][jy][jz] = fs[0][jy][jz];
      result[mesh->ngx-1][jy][jz] = fs[mesh->ngx-1][jy][jz];
    }

  // Smooth using simple 1-2-1 filter

  for(int jx=1;jx<mesh->ngx-1;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	result[jx][jy][jz] = 0.5*fs[jx][jy][jz] + 0.25*( fs[jx-1][jy][jz] + fs[jx+1][jy][jz] );
      }

  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}


const Field3D smooth_y(const Field3D &f) {
  Field3D result;

  result.allocate();
  
  // Copy boundary region
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz;jz++) {
      result[jx][0][jz] = f[jx][0][jz];
      result[jx][mesh->ngy-1][jz] = f[jx][mesh->ngy-1][jz];
    }
  
  // Smooth using simple 1-2-1 filter

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=1;jy<mesh->ngy-1;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	result[jx][jy][jz] = 0.5*f[jx][jy][jz] + 0.25*( f[jx][jy-1][jz] + f[jx][jy+1][jz] );
      }

  // Need to communicate boundaries
  mesh->communicate(result);
  
  return result;
}

const Field2D averageY(const Field2D &f) {
  return mesh->averageY(f);
}

const Field3D averageY(const Field3D &f) {
  return mesh->averageY(f);
}

const Field3D smoothXY(const Field3D &f) {
  Field3D result;
  result.allocate();

  for(int x=2;x<mesh->ngx-2;x++)
    for(int y=2;y<mesh->ngy-2;y++)
      for(int z=0;z<mesh->ngz;z++) {
        result[x][y][z] = 0.5*f[x][y][z] + 0.125*( 0.5*f[x+1][y][z] + 0.125*(f[x+2][y][z] + f[x][y][z] + f[x+1][y-1][z] + f[x+1][y+1][z]) +
                                                   0.5*f[x-1][y][z] + 0.125*(f[x][y][z] + f[x-2][y][z] + f[x-1][y-1][z] + f[x-1][y+1][z]) +
                                                   0.5*f[x][y-1][z] + 0.125*(f[x+1][y-1][z] + f[x-1][y-1][z] + f[x][y-2][z] + f[x][y][z]) +
                                                   0.5*f[x][y+1][z] + 0.125*(f[x+1][y+1][z] + f[x-1][y+1][z] + f[x][y][z] + f[x][y+2][z]));
      }
  
  return result;
}

/// Nonlinear filtering to remove grid-scale noise
/*!
  From a paper:

  W.Shyy et. al. JCP 102 (1) September 1992 page 49

  "On the Suppression of Numerical Oscillations Using a Non-Linear Filter"
  
 */
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
#ifdef CHECK
  msg_stack.push("nl_filter_x( Field3D )");
#endif
  
  Field3D fs;
  fs = f.shiftZ(true); // Shift into BoutReal space
  Field3D result;
  rvec v;
  
  for(int jy=0;jy<mesh->ngy;jy++)
    for(int jz=0;jz<mesh->ngz-1;jz++) {
      fs.getXArray(jy, jz, v);
      nl_filter(v, w);
      result.setXArray(jy, jz, v);
    }
  
  result = result.shiftZ(false); // Shift back

#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D nl_filter_y(const Field3D &fs, BoutReal w) {
#ifdef CHECK
  msg_stack.push("nl_filter_x( Field3D )");
#endif
  
  Field3D result;
  rvec v;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz-1;jz++) {
      fs.getYArray(jx, jz, v);
      nl_filter(v, w);
      result.setYArray(jx, jz, v);
    }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D nl_filter_z(const Field3D &fs, BoutReal w) {
#ifdef CHECK
  msg_stack.push("nl_filter_z( Field3D )");
#endif
  
  Field3D result;
  rvec v;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++) {
      fs.getZArray(jx, jy, v);
      nl_filter(v, w);
      result.setZArray(jx, jy, v);
    }

#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D nl_filter(const Field3D &f, BoutReal w)
{
  Field3D result;
  /// Perform filtering in Z, Y then X
  result = nl_filter_x(nl_filter_y(nl_filter_z(f, w), w), w);
  /// Communicate boundaries
  mesh->communicate(result);
  return result;
}
