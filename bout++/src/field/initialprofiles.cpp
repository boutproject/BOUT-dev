/**************************************************************************
 * Sets initial profiles
 *
 * ChangeLog
 * =========
 *
 * 2010-05-12 Ben Dudson <bd512@york.ac.uk>
 *    
 *    * Changed random numbers to use a hash of the parameters
 *      so that the phase doesn't vary with number of processors or grid size
 *      User can vary phase to give a different random sequence
 * 
 **************************************************************************
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
 **************************************************************************/

#include "globals.h"
#include "initialprofiles.h"

#include <math.h>
#include <string.h>
#include <stdlib.h>

BoutReal Prof1D(BoutReal s, BoutReal s0, BoutReal sMin, BoutReal sMax, BoutReal sWidth, int nMode, BoutReal phase, int opt);

// Initial profile options

bool ShiftInitial; // Shift initial profile? (default true if shifting in X)
bool Ballooning;   // Use ballooning transform to enforce periodicity

// Get initial profile options
void get_profile_opts()
{
  static bool gotopts = false;

  if(gotopts)
    return;
  gotopts = true;
  
  output.write("Initial profile global options\n");
  options.setSection("");
  
  bool TwistShift;
  OPTION(TwistShift, false);
  OPTION(Ballooning, TwistShift); // Use if have TwistShift
  OPTION(ShiftInitial, mesh->ShiftXderivs && (!Ballooning));
  output.write("\n");
}

int initial_profile(const char *name, Field3D &var)
{
  BoutReal scale;
  int xs_opt, ys_opt, zs_opt;
  int xs_mode, ys_mode, zs_mode;
  BoutReal xs_phase, ys_phase, zs_phase;
  BoutReal xs_s0, ys_s0, zs_s0;
  BoutReal xs_wd, ys_wd, zs_wd;
  int jx, jy, jz;
  BoutReal cx, cy, cz;

  var = 0.0;

  /////////// get options //////////
  
  get_profile_opts();
      
  // Look in the section with the variable name
  // Otherwise look for global settings
  options.get(name, "All", "scale", scale, 1.0e-4);

  // What type of profile? 0 - constant, 1 - Gaussian, 2 - Sinusoidal  
  options.get(name, "All", "xs_opt", xs_opt, 1);

  options.get(name, "All", "ys_opt", ys_opt, 0);
  
  options.get(name, "All", "zs_opt", zs_opt, 2);

  // Mode number (for sinusoidal)
  
    options.get(name, "All", "xs_mode", xs_mode, 4);

    options.get(name, "All", "ys_mode", ys_mode, 4);

    options.get(name, "All", "zs_mode", zs_mode, 4);

  // Phase (for sinusoidal), in units of pi

    options.get(name, "All", "xs_phase", xs_phase, 0.);

    options.get(name, "All", "ys_phase", ys_phase, 0.);

    options.get(name, "All", "zs_phase", zs_phase, 0.);
  
  // Gaussian peak location
      
    options.get(name, "All", "xs_s0", xs_s0, 0.5);

    options.get(name, "All", "ys_s0", ys_s0, 0.5);
      
    options.get(name, "All", "zs_s0", zs_s0, 0.5);

  // Gaussian width

    options.get(name, "All", "xs_wd", xs_wd, 0.2);

    options.get(name, "All", "ys_wd", ys_wd, 0.2);

    options.get(name, "All", "zs_wd", zs_wd, 0.2);

  for (jx=0; jx < mesh->ngx; jx++) {
    BoutReal xcoord = mesh->GlobalX(jx);
    
    for (jz=0; jz < mesh->ngz; jz++) {
      for (jy=0; jy < mesh->ngy; jy++) {
	BoutReal ycoord = mesh->GlobalY(jy);
	
	cx=Prof1D(xcoord, xs_s0, 0., 1.0, xs_wd, xs_mode, xs_phase, xs_opt);
	cy=Prof1D(ycoord, ys_s0, 0., 1.0, ys_wd, ys_mode, ys_phase, ys_opt);
	cz=Prof1D((BoutReal) jz, zs_s0, 0., (BoutReal) (mesh->ngz-1), zs_wd, zs_mode, zs_phase, zs_opt);
	
	var[jx][jy][jz] = scale*cx*cy*cz;
	
	BoutReal ts; ///< Twist-shift angle
	if(mesh->surfaceClosed(jx, ts) && Ballooning) {
	  // Use a truncated Ballooning transform to enforce periodicity
	  
	  int ball_n = 3; // How many times around in each direction
	  
	  for(int i=1; i<= ball_n; i++) {
	    // y - i * nycore
	    cy=Prof1D(ycoord - i, ys_s0, 0., 1.0, ys_wd, ys_mode, ys_phase, ys_opt);
	    cz=Prof1D((BoutReal) jz + ((BoutReal) i)*ts/mesh->dz, zs_s0, 0., (BoutReal) (mesh->ngz-1), zs_wd, zs_mode, zs_phase, zs_opt);
	    var[jx][jy][jz] += scale*cx*cy*cz;
	    
	    // y + i * nycore
	    cy=Prof1D(ycoord + i, ys_s0, 0., 1., ys_wd, ys_mode, ys_phase, ys_opt);
	    cz=Prof1D((BoutReal) jz - ((BoutReal) i)*ts/mesh->dz, zs_s0, 0., (BoutReal) (mesh->ngz-1), zs_wd, zs_mode, zs_phase, zs_opt);
	    var[jx][jy][jz] += scale*cx*cy*cz;
	  }
	}
	
	/*
	  //NOTE: This code causes NANs if prof1D goes to zero
	if(!MYPE_IN_CORE) {
	  ly = YGLOBAL(jy);
	  if(ly <= jyseps1_1) {
	    // Inner lower
	    var[jx][jy][jz] *= Prof1D((BoutReal) (ly-jyseps1_1-1), ys_s0, 0., (BoutReal) nycore, ys_wd, ys_mode, ys_phase, ys_opt) / Prof1D(0.0, ys_s0, 0., (BoutReal) nycore, ys_wd, ys_mode, ys_phase, ys_opt);
	  }else if(ly > jyseps2_2) {
	    // Outer lower
	    var[jx][jy][jz] *= Prof1D((BoutReal) (ly-jyseps2_2+nycore-1), ys_s0, 0., (BoutReal) nycore, ys_wd, ys_mode, ys_phase, ys_opt) / Prof1D((BoutReal) (nycore-1), ys_s0, 0., (BoutReal) nycore, ys_wd, ys_mode, ys_phase, ys_opt);
	  }
	}
	*/
      }
    }
  }

  if(ShiftInitial)
    var = var.shiftZ(false);

  return(0);
}

// For 2D variables almost identical, just no z dependence
int initial_profile(const char *name, Field2D &var)
{
  BoutReal scale;
  int xs_opt, ys_opt;
  int xs_mode, ys_mode;
  BoutReal xs_phase, ys_phase;
  BoutReal xs_s0, ys_s0;
  BoutReal xs_wd, ys_wd;
  int jx, jy;
  BoutReal cx, cy;

  var = 0.0;

  /////////// get options //////////
  
  get_profile_opts();

  // Size of the initial perturbation
  // Look in the section with the variable name
  // Otherwise look for global settings
  options.get(name, "All", "scale", scale, 1e-4);

  // What type of profile? 0 - constant, 1 - Gaussian, 2 - Sinusoidal

    options.get(name, "All", "xs_opt", xs_opt, 1);

    options.get(name, "All", "ys_opt", ys_opt, 0);

  // Mode number (for sinusoidal)

    options.get(name, "All", "xs_mode", xs_mode, 4);

    options.get(name, "All", "ys_mode", ys_mode, 4);

  // Phase (for sinusoidal), in units of pi

    options.get(name, "All", "xs_phase", xs_phase, 0.);

    options.get(name, "All", "ys_phase", ys_phase, 0.);
  
  // Gaussian peak location

    options.get(name, "All", "xs_s0", xs_s0, 0.5);

    options.get(name, "All", "ys_s0", ys_s0, 0.5);

  // Gaussian width

    options.get(name, "All", "xs_wd", xs_wd, 0.2);

    options.get(name, "All", "ys_wd", ys_wd, 0.2);
  

  for (jx=0; jx < mesh->ngx; jx++) {
    BoutReal xcoord = mesh->GlobalX(jx);
    for (jy=0; jy < mesh->ngy; jy++) {
      BoutReal ycoord = mesh->GlobalY(jy);
      cx=Prof1D(xcoord, xs_s0, 0., 1., xs_wd, xs_mode, xs_phase, xs_opt);
      cy=Prof1D(ycoord, ys_s0, 0., 1., ys_wd, ys_mode, ys_phase, ys_opt);
      
      var[jx][jy] = scale*cx*cy;
    }
  }
  
  return(0);
}

int initial_profile(const char *name, Vector2D &var)
{
  char *s;

  s = (char*) malloc(strlen(name)+3);

  if(var.covariant) {
    sprintf(s, "%s_x", name);
    initial_profile(s, var.x);
    
    sprintf(s, "%s_y", name);
    initial_profile(s, var.y);

    sprintf(s, "%s_z", name);
    initial_profile(s, var.z);
  }else {
    sprintf(s, "%sx", name);
    initial_profile(s, var.x);
    
    sprintf(s, "%sy", name);
    initial_profile(s, var.y);

    sprintf(s, "%sz", name);
    initial_profile(s, var.z);
  }

  free(s);
  
  return 0;
}

int initial_profile(const char *name, Vector3D &var)
{
  char *s;

  s = (char*) malloc(strlen(name)+3);

  if(var.covariant) {
    sprintf(s, "%s_x", name);
    initial_profile(s, var.x);
    
    sprintf(s, "%s_y", name);
    initial_profile(s, var.y);

    sprintf(s, "%s_z", name);
    initial_profile(s, var.z);
  }else {
    sprintf(s, "%sx", name);
    initial_profile(s, var.x);
    
    sprintf(s, "%sy", name);
    initial_profile(s, var.y);

    sprintf(s, "%sz", name);
    initial_profile(s, var.z);
  }

  free(s);
  
  return 0;
}

/// Hash BoutReal values to produce a number between -1 and 1
BoutReal hash_BoutReals(BoutReal a, BoutReal b, BoutReal c)
{
  int hash = 0;
  
  unsigned char *p = (unsigned char*) &a;
  for(size_t i=0;i<sizeof(BoutReal);i++) {
    hash += (int) p[i];
    hash += hash << 10;
    hash ^= hash >> 6;
  }
  p = (unsigned char*) &b;
  for(size_t i=0;i<sizeof(BoutReal);i++) {
    hash += (int) p[i];
    hash += hash << 10;
    hash ^= hash >> 6;
  }
  p = (unsigned char*) &c;
  for(size_t i=0;i<sizeof(BoutReal);i++) {
    hash += (int) p[i];
    hash += hash << 10;
    hash ^= hash >> 6;
  }

  hash += hash << 3;
  hash ^= hash >> 11;
  hash += hash << 15;

  if(hash <= 0) hash += 2147483647;
  return ((BoutReal) hash)/2147483646.0;
}

BoutReal Prof1D(BoutReal s, BoutReal s0, BoutReal sMin, BoutReal sMax, BoutReal sWidth, int nMode, BoutReal phase, int opt)
/*Calculate initial perturbation for given seed parameters*/
{
  
  BoutReal res;
  BoutReal sNorm=s/(sMax-sMin);

  phase *= PI;

  switch (opt) 
    {
    case 1: 
      res=exp(-pow((sNorm-s0)/sWidth,2.));
      break;

    case 2:
      res=sin(nMode*sNorm*TWOPI + phase);
      break;

    case 3:
      res=(1./4./4.)*cos(1.*sNorm*TWOPI  + 0.01*hash_BoutReals(sNorm, phase, 1.)*PI)
	+(1./3./3.)*cos(2*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 2.)*PI)
	+(1./2./2.)*cos(3*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 3.)*PI)
	+(1./1.)*cos(4*sNorm*TWOPI       + 0.01*hash_BoutReals(sNorm, phase, 4.)*PI)
	+(1./2./2.)*cos(5*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 5.)*PI)
	+(1./3./3.)*cos(6*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 6.)*PI)
	+(1./4./4.)*cos(7*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 7.)*PI)
	+(1./5./5.)*cos(8*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 8.)*PI)
	+(1./6./6.)*cos(9*sNorm*TWOPI    + 0.01*hash_BoutReals(sNorm, phase, 9.)*PI)
	+(1./7./7.)*cos(10*sNorm*TWOPI   + 0.01*hash_BoutReals(sNorm, phase, 10.)*PI)
	+(1./8./8.)*cos(11*sNorm*TWOPI   + 0.01*hash_BoutReals(sNorm, phase, 11.)*PI)
	+(1./9./9.)*cos(12*sNorm*TWOPI   + 0.01*hash_BoutReals(sNorm, phase, 12.)*PI)
	+(1./10./10.)*cos(13*sNorm*TWOPI + 0.01*hash_BoutReals(sNorm, phase, 13.)*PI)
	+(1./11./11.)*cos(14*sNorm*TWOPI + 0.01*hash_BoutReals(sNorm, phase, 14.)*PI);
      break;
    case 4:
      res = sin(nMode*sNorm*TWOPI)
        + sin(2.*nMode*sNorm*TWOPI + hash_BoutReals(sNorm, phase, nMode)*TWOPI)
        + sin(3.*nMode*sNorm*TWOPI + hash_BoutReals(sNorm, phase, 2.*nMode)*TWOPI)
        + sin(4.*nMode*sNorm*TWOPI + hash_BoutReals(sNorm, phase, 3.*nMode)*TWOPI);
      break;

    default: res=1.0;
    }

  return res;
}

/**************************************************************************
 * Routines to generate profiles 
 **************************************************************************/

/// Generate a 3D field with a given Z oscillation
/*!
  @param[in] n      Mode number. Note that this is mode-number in the domain
  @param[in] phase  Phase shift in units of pi
*/
const Field3D genZMode(int n, BoutReal phase)
{
  Field3D result;

  result.allocate();
  BoutReal ***d = result.getData();
  
  for(int jz=0;jz<mesh->ngz;jz++) {
    BoutReal val = sin(phase*PI +  TWOPI * ((BoutReal) jz)/ ((BoutReal) mesh->ngz-1) );
    for(int jx=0;jx<mesh->ngx;jx++)
      for(int jy=0;jy<mesh->ngy;jy++)
	d[jx][jy][jz] = val;
  }
  return result;
}
