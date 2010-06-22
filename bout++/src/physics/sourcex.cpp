/**************************************************************
 * radial source and mask operators
 **************************************************************/

#include "globals.h"
#include <math.h>

#include "sourcex.h"

/* hyperbolic tangent */
//template<typename TYPE> inline TYPE TanH(TYPE &a)

real TanH(real a)
{
  real temp = exp(a);
  //  TanH = (temp - 1 / temp) / (temp + 1 / temp);
  return ((temp - 1.0 / temp) / (temp + 1.0 / temp));
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_tanhx(const Field2D &f,real swidth,real slength)
{
  Field2D  fs, result;

  result.allocate();
  
  // create a radial buffer zone to set jpar zero near radial boundary

  //  real slength = 0.5;
  //  real width = 20.0;
  real length  = slength;
  real width   = swidth;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      {
	real lx = mesh->GlobalX(jx) - length;
	real dampl = TanH(lx/width);
	result[jx][jy] = 0.5*(1.0 - dampl);
      }
  
  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_expx2(const Field2D &f,real swidth,real slength)
{
  Field2D  fs, result;

  result.allocate();
  
  // create a radial buffer zone to set jpar zero near radial boundary

  //  real slength = 0.5;
  //  real width = 20.0;

  //	    output.write("source, swidth=%e, width=%e, slength=%e, length=%e, nx=%d, ny=%d, MXG=%d, MYG=%d\n", swidth, width, slength, length, nx, ny, MXG, MYG);

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      {
	real lx = mesh->GlobalX(jx) - slength;
	real dampl = exp(-lx*lx/swidth/swidth);
	result[jx][jy] = dampl;
      }
  
  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhx(const Field2D &f0, const Field3D &f,real swidth,real slength, bool realspace)
//const Field3D sink_tanhx(const Field2D &f0, const Field3D &f, bool realspace)
{
  Field3D fs, result;
  Field2D fs0;
 
  if(realspace) {
    fs = f.shiftZ(true); // Shift into real space
  }else
    fs = f;
 
  result.allocate();
 
// create a radial buffer zone to set jpar zero near radial boundary
 
  //  real slength = 0.15;
  //  real width = 20.0;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
        real rlx = 1. - mesh->GlobalX(jx) - slength;
	real dampr = TanH(rlx/swidth);
	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]);
	//	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]-fs0[jx][jy]);
      }
 
  if(realspace)
    result = result.shiftZ(false); // Shift back
 
  // Need to communicate boundaries
  mesh->communicate(result);
  
  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D mask_x(const Field3D &f, bool realspace)
{
  Field3D fs, result;

  if(realspace) {
    fs = f.shiftZ(true); // Shift into real space
  }else
    fs = f; 

  result.allocate();
  
// create a radial buffer zone to set jpar zero near radial boundary

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	real lx = mesh->GlobalX(jx);
        real dampl = TanH(lx/40.0);
        real dampr = TanH((1. - lx)/40.0);
	
	result[jx][jy][jz] = (1.0 - dampl*dampr)*fs[jx][jy][jz];
	//	result[jx][jy][jz] = dampl*fs[jx][jy][jz]*dampr;
      }
  
  if(realspace)
    result = result.shiftZ(false); // Shift back

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxl(const Field2D &f0, const Field3D &f,real swidth,real slength, bool realspace)
{
  Field3D fs, result;
  Field2D fs0;
 
  if(realspace) {
    fs = f.shiftZ(true); // Shift into real space
  }else
    fs = f;
 
  result.allocate();
 
// create a radial buffer zone to set jpar zero near radial boundary
 
  //  real slength = 0.15;
  //  real width = 20.0;

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {

	real lx = mesh->GlobalX(jx) - slength;
	real dampl = TanH(lx/swidth);

	result[jx][jy][jz] = 0.5*(1.0 - dampl)*(fs[jx][jy][jz]);
	//	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]-fs0[jx][jy]);
      }
 
  if(realspace)
    result = result.shiftZ(false); // Shift back
 
  // Need to communicate boundaries
  mesh->communicate(result);
 
  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxr(const Field2D &f0, const Field3D &f,real swidth,real slength, bool realspace)
{
  Field3D fs, result;
  Field2D fs0;
 
  if(realspace) {
    fs = f.shiftZ(true); // Shift into real space
  }else
    fs = f;
 
  result.allocate();
 
// create a radial buffer zone to set jpar zero near radial boundary
 
  //  real slength = 0.15;
  //  real width = 20.0;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
        real rlx = 1. - mesh->GlobalX(jx) - slength;
	real dampr = TanH(rlx/swidth);

	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]);
	//	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]-fs0[jx][jy]);
      }
 
  if(realspace)
    result = result.shiftZ(false); // Shift back
 
  // Need to communicate boundaries
  mesh->communicate(result);
 
  return result;
}

// create radial buffer zones to damp Psi to zero near radial boundaries
const Field3D buff_x(const Field3D &f, bool realspace)
{
  Field3D fs, result;

  if(realspace) {
    fs = f.shiftZ(true); // Shift into real space
  }else
    fs = f; 

  result.allocate();
  
  // create a radial buffer zone to set jpar zero near radial boundary

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	real lx = mesh->GlobalX(jx);
	//        real rlx = nx-(XGLOBAL(jx)+1);
        real rlx = 1. - lx;
        real dampl = 1.e0;
        real dampr = 1.e0;
	real deltal = 0.05;
	real deltar = 0.05;
	
	result[jx][jy][jz] = (dampl*exp(- (real) (lx*lx)/(deltal*deltal))
			      +dampr*exp(-(real) ((rlx*rlx))/(deltar*deltar)))*fs[jx][jy][jz];
      }
  
  // Need to communicate boundaries
  mesh->communicate(result);

  if(realspace)
    result = result.shiftZ(false); // Shift back

  return result;
}
