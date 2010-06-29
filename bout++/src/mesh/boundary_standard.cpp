
#include "boundary_standard.h"
#include "globals.h"
#include "invert_laplace.h"
#include "fft.h"

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryDirichlet::clone(BoundaryRegion *region)
{
  return new BoundaryDirichlet(region);
}

void BoundaryDirichlet::apply(Field2D &f)
{
  // Just loop over all elements and set to zero
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = 0.;
}

void BoundaryDirichlet::apply(Field3D &f)
{
  // Just loop over all elements and set to zero
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f[bndry->x][bndry->y][z] = 0.;
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryNeumann::clone(BoundaryRegion *region)
{
  return new BoundaryNeumann(region);
}

void BoundaryNeumann::apply(Field2D &f)
{
  // Loop over all elements and set equal to the next point in
  for(bndry->first(); !bndry->isDone(); bndry->next())
    f[bndry->x][bndry->y] = f[bndry->x - bndry->bx][bndry->y - bndry->by];
}

void BoundaryNeumann::apply(Field3D &f)
{
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++)
      f[bndry->x][bndry->y][z] = f[bndry->x - bndry->bx][bndry->y - bndry->by][z];
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryZeroLaplace::clone(BoundaryRegion *region)
{
  return new BoundaryZeroLaplace(region);
}

void BoundaryZeroLaplace::apply(Field2D &f)
{
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    bout_error("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }

  // Constant X derivative
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    BoutReal g = (f[x-bx][y] - f[x-2*bx][y]) / mesh->dx[x-bx][y];
    // Loop in X towards edge of domain
    do {
      f[x][y] = f[x-bx][y] + g*mesh->dx[x][y];
      bndry->next();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

void BoundaryZeroLaplace::apply(Field3D &f)
{
  static dcomplex *c0 = (dcomplex*) NULL, *c1;
  int ncz = mesh->ngz-1;
  
  if(c0 == (dcomplex*) NULL) {
    // allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
  }
  
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    bout_error("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    // bndry->(x,y) is the first point in the boundary
    // bndry->(x-bx,y) is the last "real" point in the domain
    
    int x = bndry->x;
    int y = bndry->y;
    
    // Take FFT of last 2 points in domain
    ZFFT(f[x-bx][y], mesh->zShift[x-bx][y], c0);
    ZFFT(f[x-2*bx][y], mesh->zShift[x-2*bx][y], c1);
    c1[0] = c0[0] - c1[0]; // Only need gradient
    
    // Solve  mesh->g11*d2f/dx2 - mesh->g33*kz^2f = 0
    // Assume mesh->g11, mesh->g33 constant -> exponential growth or decay
    
    // Loop in X towards edge of domain
    do {
      // kz = 0 solution
      c0[0] += c1[0];  // Straight line
      
      // kz != 0 solution
      BoutReal coef = -1.0*sqrt(mesh->g33[x][y] / mesh->g11[x][y])*mesh->dx[x][y];
      for(int jz=1;jz<=ncz/2;jz++) {
	BoutReal kwave=jz*2.0*PI/mesh->zlength; // wavenumber in [rad^-1]
	c0[jz] *= exp(coef*kwave); // The decaying solution only
      }
      // Reverse FFT
      ZFFT_rev(c0, mesh->zShift[x][y], f[x][y]);
      
      bndry->next();
      x = bndry->x; y = bndry->y;
    }while(!bndry->isDone());
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryConstLaplace::clone(BoundaryRegion *region)
{
  return new BoundaryConstLaplace(region);
}

void BoundaryConstLaplace::apply(Field2D &f)
{
  
}

void BoundaryConstLaplace::apply(Field3D &f)
{
  if((bndry->location != BNDRY_XIN) && (bndry->location != BNDRY_XOUT)) {
    // Can't apply this boundary condition to non-X boundaries
    bout_error("ERROR: Can't apply Zero Laplace condition to non-X boundaries\n");
  }
  
  static dcomplex *c0 = (dcomplex*) NULL, *c1, *c2;
  int ncz = mesh->ngz-1;
  if(c0 == (dcomplex*) NULL) {
    //allocate memory
    c0 = new dcomplex[ncz/2 + 1];
    c1 = new dcomplex[ncz/2 + 1];
    c2 = new dcomplex[ncz/2 + 1];
  }
  
  int bx = bndry->bx;
  // Loop over the Y dimension
  for(bndry->first(); !bndry->isDone(); bndry->nextY()) {
    int x = bndry->x;
    int y = bndry->y;
    
    // Take FFT of last 3 points in domain
    ZFFT(f[x-bx][y], mesh->zShift[x-bx][y], c0);
    ZFFT(f[x-2*bx][y], mesh->zShift[x-2*bx][y], c1);
    ZFFT(f[x-3*bx][y], mesh->zShift[x-3*bx][y], c1);
    dcomplex k0lin = (c1[0] - c0[0])/mesh->dx[x-bx][y]; // for kz=0 solution
    
    // Calculate Delp2 on point MXG+1 (and put into c1)
    for(int jz=0;jz<=ncz/2;jz++) {
      dcomplex d,e,f;
      laplace_tridag_coefs(x-2*bx, y, jz, d, e, f);
      c1[jz] = d*c0[jz] + e*c1[jz] + f*c2[jz];
    }
  }
}

///////////////////////////////////////////////////////////////

BoundaryOp* BoundaryRelax::clone(BoundaryOp *operation)
{
  BoundaryRelax* result = new BoundaryRelax(r);
  result->op = operation;
  result->bndry = operation->bndry;
  
  return result;
}
  
void BoundaryRelax::apply(Field2D &f)
{
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply(Field3D &f)
{
  // Just apply the original boundary condition to f
  op->apply(f);
}

void BoundaryRelax::apply_ddt(Field2D &f)
{
#ifdef CHECK
  msg_stack.push("BoundaryRelax::apply_ddt(Field2D)");
#endif

  // Make a copy of f
  Field2D g = f;
  // Apply the boundary to g
  op->apply(g);
  
  bndry->first();
  exit(0);
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next())
    ddt(f)[bndry->x][bndry->y] = ddt(f)[bndry->x - bndry->bx][bndry->y - bndry->by] 
      + r * (g[bndry->x][bndry->y] - f[bndry->x][bndry->y]);

#ifdef CHECK
  msg_stack.pop();
#endif
}

void BoundaryRelax::apply_ddt(Field3D &f)
{
#ifdef CHECK
  msg_stack.push("BoundaryRelax::apply_ddt(Field2D)");
#endif
  
  // Make a copy of f
  Field3D g = f; // NOTE: This is not very efficient... copying entire field
  // Apply the boundary to g
  op->apply(g);
  // Set time-derivatives
  for(bndry->first(); !bndry->isDone(); bndry->next())
    for(int z=0;z<mesh->ngz;z++) {
      ddt(f)[bndry->x][bndry->y][z] = ddt(f)[bndry->x - bndry->bx][bndry->y - bndry->by][z]
	+ r * (g[bndry->x][bndry->y][z] - f[bndry->x][bndry->y][z]);
    }

#ifdef CHECK
  msg_stack.pop();
#endif
}

