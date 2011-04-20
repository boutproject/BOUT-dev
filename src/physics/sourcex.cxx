/**************************************************************
 * radial source and mask operators
 **************************************************************/

#include <globals.hxx>
#include <math.h>

#include <sourcex.hxx>

/* hyperbolic tangent */
//template<typename TYPE> inline TYPE TanH(TYPE &a)

BoutReal TanH(BoutReal a)
{
  BoutReal temp = exp(a);
  //  TanH = (temp - 1 / temp) / (temp + 1 / temp);
  return ((temp - 1.0 / temp) / (temp + 1.0 / temp));
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_tanhx(const Field2D &f,BoutReal swidth,BoutReal slength)
{
  Field2D  fs, result;

  result.allocate();
  
  // create a radial buffer zone to set jpar zero near radial boundary

  //  BoutReal slength = 0.5;
  //  BoutReal width = 20.0;
  BoutReal length  = slength;
  BoutReal width   = swidth;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      {
	BoutReal lx = mesh->GlobalX(jx) - length;
	BoutReal dampl = TanH(lx/width);
	result[jx][jy] = 0.5*(1.0 - dampl);
      }
  
  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field2D source_expx2(const Field2D &f,BoutReal swidth,BoutReal slength)
{
  Field2D  fs, result;

  result.allocate();
  
  // create a radial buffer zone to set jpar zero near radial boundary

  //  BoutReal slength = 0.5;
  //  BoutReal width = 20.0;

  //	    output.write("source, swidth=%e, width=%e, slength=%e, length=%e, nx=%d, ny=%d, MXG=%d, MYG=%d\n", swidth, width, slength, length, nx, ny, MXG, MYG);

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      {
	BoutReal lx = mesh->GlobalX(jx) - slength;
	BoutReal dampl = exp(-lx*lx/swidth/swidth);
	result[jx][jy] = dampl;
      }
  
  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhx(const Field2D &f0, const Field3D &f,BoutReal swidth,BoutReal slength, bool BoutRealspace)
//const Field3D sink_tanhx(const Field2D &f0, const Field3D &f, bool BoutRealspace)
{
  Field3D fs, result;
  Field2D fs0;
 
  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f;
 
  result.allocate();
 
// create a radial buffer zone to set jpar zero near radial boundary
 
  //  BoutReal slength = 0.15;
  //  BoutReal width = 20.0;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
        BoutReal rlx = 1. - mesh->GlobalX(jx) - slength;
	BoutReal dampr = TanH(rlx/swidth);
	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]);
	//	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]-fs0[jx][jy]);
      }
 
  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back
 
  // Need to communicate boundaries
  mesh->communicate(result);
  
  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D mask_x(const Field3D &f, bool BoutRealspace)
{
  Field3D fs, result;

  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f; 

  result.allocate();
  
// create a radial buffer zone to set jpar zero near radial boundary

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	BoutReal lx = mesh->GlobalX(jx);
        BoutReal dampl = TanH(lx/40.0);
        BoutReal dampr = TanH((1. - lx)/40.0);
	
	result[jx][jy][jz] = (1.0 - dampl*dampr)*fs[jx][jy][jz];
	//	result[jx][jy][jz] = dampl*fs[jx][jy][jz]*dampr;
      }
  
  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxl(const Field2D &f0, const Field3D &f,BoutReal swidth,BoutReal slength, bool BoutRealspace)
{
  Field3D fs, result;
  Field2D fs0;
 
  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f;
 
  result.allocate();
 
// create a radial buffer zone to set jpar zero near radial boundary
 
  //  BoutReal slength = 0.15;
  //  BoutReal width = 20.0;

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {

	BoutReal lx = mesh->GlobalX(jx) - slength;
	BoutReal dampl = TanH(lx/swidth);

	result[jx][jy][jz] = 0.5*(1.0 - dampl)*(fs[jx][jy][jz]);
	//	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]-fs0[jx][jy]);
      }
 
  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back
 
  // Need to communicate boundaries
  mesh->communicate(result);
 
  return result;
}

// create radial buffer zones to set jpar zero near radial boundaries
const Field3D sink_tanhxr(const Field2D &f0, const Field3D &f,BoutReal swidth,BoutReal slength, bool BoutRealspace)
{
  Field3D fs, result;
  Field2D fs0;
 
  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f;
 
  result.allocate();
 
// create a radial buffer zone to set jpar zero near radial boundary
 
  //  BoutReal slength = 0.15;
  //  BoutReal width = 20.0;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
        BoutReal rlx = 1. - mesh->GlobalX(jx) - slength;
	BoutReal dampr = TanH(rlx/swidth);

	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]);
	//	result[jx][jy][jz] = 0.5*(1.0 - dampr)*(fs[jx][jy][jz]-fs0[jx][jy]);
      }
 
  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back
 
  // Need to communicate boundaries
  mesh->communicate(result);
 
  return result;
}

// create radial buffer zones to damp Psi to zero near radial boundaries
const Field3D buff_x(const Field3D &f, bool BoutRealspace)
{
  Field3D fs, result;

  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f; 

  result.allocate();
  
  // create a radial buffer zone to set jpar zero near radial boundary

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	BoutReal lx = mesh->GlobalX(jx);
	//        BoutReal rlx = nx-(XGLOBAL(jx)+1);
        BoutReal rlx = 1. - lx;
        BoutReal dampl = 1.e0;
        BoutReal dampr = 1.e0;
	BoutReal deltal = 0.05;
	BoutReal deltar = 0.05;
	
	result[jx][jy][jz] = (dampl*exp(- (BoutReal) (lx*lx)/(deltal*deltal))
			      +dampr*exp(-(BoutReal) ((rlx*rlx))/(deltar*deltar)))*fs[jx][jy][jz];
      }
  
  // Need to communicate boundaries
  mesh->communicate(result);

  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back

  return result;
}
