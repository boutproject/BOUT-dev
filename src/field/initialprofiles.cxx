/**************************************************************************
 * Sets initial profiles
 *
 * ChangeLog
 * =========
 * 
 * 2011-02-12 Ben Dudson <bd512@york.ac.uk>
 *    * Changed to use new options system. For now the structure of the
 *      options is the same, but this could be modified more easily in future
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

#include <globals.hxx>
#include <initialprofiles.hxx>
#include <boutexception.hxx>
#include <field_factory.hxx>

#include <math.h>
#include <string.h>
#include <stdlib.h>

BoutReal Prof1D(BoutReal s, BoutReal s0, BoutReal sMin, BoutReal sMax, BoutReal sWidth, int nMode, BoutReal phase, int opt);

// Initial profile options

bool ShiftInitial; // Shift initial profile? (default true if shifting in X)
bool Ballooning;   // Use ballooning transform to enforce periodicity

// Get initial profile options
void get_profile_opts(Options *opt)
{
  static bool gotopts = false;

  if(gotopts)
    return;
  gotopts = true;
  
  output.write("Initial profile global options\n");
  
  bool TwistShift;
  OPTION(opt, TwistShift, false);
  OPTION(opt, Ballooning, TwistShift); // Use if have TwistShift
  OPTION(opt, ShiftInitial, mesh->ShiftXderivs && (!Ballooning));
  output.write("\n");
}

// Macro which looks first in one section then another for an option
#define FIND_OPT(opt1, opt2, name, var, def) {  \
    if(opt1->isSet(name)) {                     \
      opt1->get(name, var, def);                \
    }else                                       \
      opt2->get(name, var, def);}

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
  
  Options *globalOpts = Options::getRoot();
  get_profile_opts(globalOpts);
  
  output.write("Setting initial value of %s\n", name);
  
  // Get the section for all variables
  Options *allOpts = globalOpts->getSection("All");
  
  // Get the section for this specific variable 
  Options *varOpts = globalOpts->getSection(name);
  
  FIND_OPT(varOpts, allOpts, "scale", scale, 1.0e-4);

  if(varOpts->isSet("function")) {
    // Form of perturbation specified as string. Use FieldFactory to generate values
    
    FieldFactory f(mesh);
    string s;
    varOpts->get("function", s, "");
    var = scale*f.create3D(s);
  }else {
    // Backwards-compatible method
    
    // What type of profile? 0 - constant, 1 - Gaussian, 2 - Sinusoidal
    
    FIND_OPT(varOpts, allOpts, "xs_opt", xs_opt, 1);
    FIND_OPT(varOpts, allOpts, "ys_opt", ys_opt, 0);
    FIND_OPT(varOpts, allOpts, "zs_opt", zs_opt, 2);
    
    // Mode number (for sinusoidal)
    
    FIND_OPT(varOpts, allOpts, "xs_mode", xs_mode, 4);
    FIND_OPT(varOpts, allOpts, "ys_mode", ys_mode, 4);
    FIND_OPT(varOpts, allOpts, "zs_mode", zs_mode, 4);
    
    // Phase (for sinusoidal), in units of pi
    
    FIND_OPT(varOpts, allOpts, "xs_phase", xs_phase, 0.);
    FIND_OPT(varOpts, allOpts, "ys_phase", ys_phase, 0.);
    FIND_OPT(varOpts, allOpts, "zs_phase", zs_phase, 0.);
    
    // Gaussian peak location
    
    FIND_OPT(varOpts, allOpts, "xs_s0", xs_s0, 0.5);
    FIND_OPT(varOpts, allOpts, "ys_s0", ys_s0, 0.5);
    FIND_OPT(varOpts, allOpts, "zs_s0", zs_s0, 0.5);
    
    // Gaussian width
    
    FIND_OPT(varOpts, allOpts, "xs_wd", xs_wd, 0.2);
    FIND_OPT(varOpts, allOpts, "ys_wd", ys_wd, 0.2);
    FIND_OPT(varOpts, allOpts, "zs_wd", zs_wd, 0.2);

    for (jx=0; jx < mesh->ngx; jx++) {
      BoutReal xcoord = mesh->GlobalX(jx);
    
      for (jz=0; jz < mesh->ngz; jz++) {
        for (jy=0; jy < mesh->ngy; jy++) {
          BoutReal ycoord = mesh->GlobalY(jy);
	
          cx=Prof1D(xcoord, xs_s0, 0., 1.0, xs_wd, xs_mode, xs_phase, xs_opt);
          cy=Prof1D(ycoord, ys_s0, 0., 1.0, ys_wd, ys_mode, ys_phase, ys_opt);
          cz=Prof1D((BoutReal) jz, zs_s0, 0., (BoutReal) (mesh->ngz-1), zs_wd, zs_mode, zs_phase, zs_opt);
	
          var[jx][jy][jz] = scale*cx*cy*cz;
          if(!finite(var[jx][jy][jz])) {
            output.write("%d, %d, %d -> %e, %e, %e, %e\n", jx, jy, jz, cx, cy, cz, var[jx][jy][jz]);
            output.write("%e, %e, %e, %e, %e, %e\n",
                         ycoord, ys_s0, ys_wd, ys_mode, ys_phase, ys_opt);
            throw BoutException("Invalid initial profiles for '%s' at (%d,%d,%d)\n", name, jx, jy, jz);
          }
          BoutReal ts; ///< Twist-shift angle
          if(mesh->periodicY(jx, ts) && Ballooning) {
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
          }else if(Ballooning) {
            // Open surfaces. Not sure what to do, so set to zero
            var[jx][jy][jz] = 0.;
          }
          
        }
      }
    }
  } // End of backwards-compatible method

  if(ShiftInitial)
    var = var.shiftZ(false);
  
  output << endl;

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
  
  Options *globalOpts = Options::getRoot();
  get_profile_opts(globalOpts);
  
  // Get the section for all variables
  Options *allOpts = globalOpts->getSection("All");
  
  // Get the section for this specific variable 
  Options *varOpts = globalOpts->getSection(name);
  
  // Size of the initial perturbation
  
  FIND_OPT(varOpts, allOpts, "scale", scale, 1.0e-4);

  if(varOpts->isSet("function")) {
    // Form of perturbation specified as string. Use FieldFactory to generate values
    
    FieldFactory f(mesh);
    string s;
    varOpts->get("function", s, "");
    var = scale*f.create2D(s);
  }else {
    // What type of profile? 0 - constant, 1 - Gaussian, 2 - Sinusoidal
    
    FIND_OPT(varOpts, allOpts, "xs_opt", xs_opt, 1);
    FIND_OPT(varOpts, allOpts, "ys_opt", ys_opt, 0);
    
    // Mode number (for sinusoidal)
  
    FIND_OPT(varOpts, allOpts, "xs_mode", xs_mode, 4);
    FIND_OPT(varOpts, allOpts, "ys_mode", ys_mode, 4);
  
    // Phase (for sinusoidal), in units of pi

    FIND_OPT(varOpts, allOpts, "xs_phase", xs_phase, 0.);
    FIND_OPT(varOpts, allOpts, "ys_phase", ys_phase, 0.);
  
    // Gaussian peak location

    FIND_OPT(varOpts, allOpts, "xs_s0", xs_s0, 0.5);
    FIND_OPT(varOpts, allOpts, "ys_s0", ys_s0, 0.5);

    // Gaussian width
  
    FIND_OPT(varOpts, allOpts, "xs_wd", xs_wd, 0.2);
    FIND_OPT(varOpts, allOpts, "ys_wd", ys_wd, 0.2);
  
    for (jx=0; jx < mesh->ngx; jx++) {
      BoutReal xcoord = mesh->GlobalX(jx);
      for (jy=0; jy < mesh->ngy; jy++) {
        BoutReal ycoord = mesh->GlobalY(jy);
        cx=Prof1D(xcoord, xs_s0, 0., 1., xs_wd, xs_mode, xs_phase, xs_opt);
        cy=Prof1D(ycoord, ys_s0, 0., 1., ys_wd, ys_mode, ys_phase, ys_opt);
      
        var[jx][jy] = scale*cx*cy;
      }
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
    case 5:
      res=cos(1.*sNorm*TWOPI  + hash_BoutReals(sNorm, phase, 1.)*PI)
	+cos(2*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 2.)*PI)
	+cos(3*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 3.)*PI)
	+cos(4*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 4.)*PI)
	+cos(5*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 5.)*PI)
	+cos(6*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 6.)*PI)
	+cos(7*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 7.)*PI)
	+cos(8*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 8.)*PI)
	+cos(9*sNorm*TWOPI    + hash_BoutReals(sNorm, phase, 9.)*PI)
	+cos(10*sNorm*TWOPI   + hash_BoutReals(sNorm, phase, 10.)*PI)
	+cos(11*sNorm*TWOPI   + hash_BoutReals(sNorm, phase, 11.)*PI)
	+cos(12*sNorm*TWOPI   + hash_BoutReals(sNorm, phase, 12.)*PI)
	+cos(13*sNorm*TWOPI + hash_BoutReals(sNorm, phase, 13.)*PI)
	+cos(14*sNorm*TWOPI + hash_BoutReals(sNorm, phase, 14.)*PI)
	+cos(15*sNorm*TWOPI + hash_BoutReals(sNorm, phase, 14.)*PI);
      break;
      //random start, note that using rand() to randomize phase will cause problems, its possible that different 
      // cpus will produce a different random value resulting in a sharp phase discontinuity leading to mode blowup 
      // above roughly k> N/3
    case 6:
      res=cos(1.*sNorm*TWOPI +  0*(rand()/(100.0 * RAND_MAX + 1.0) * 1.0 )*PI)
	+cos(2*sNorm*TWOPI)
	+cos(3*sNorm*TWOPI)
	+cos(4*sNorm*TWOPI)
	+cos(5*sNorm*TWOPI)
	+cos(6*sNorm*TWOPI +  0*(rand()/(100.0 * RAND_MAX + 1.0) * 1.0 )*PI)
	+cos(7*sNorm*TWOPI)
	+cos(8*sNorm*TWOPI)
	+cos(9*sNorm*TWOPI)
	+cos(10*sNorm*TWOPI)
	+cos(11*sNorm*TWOPI)
	+cos(12*sNorm*TWOPI)
	+cos(13*sNorm*TWOPI)
	+cos(14*sNorm*TWOPI);
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
