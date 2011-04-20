
#include <globals.hxx>
#include <invert_laplace_gmres.hxx>
#include <invert_laplace.hxx>
#include <difops.hxx>

const Field3D LaplaceGMRES::invert(const Field3D &b, const Field3D &start, int inv_flags, bool precon, Field3D *a, Field3D *c)
{
  flags = inv_flags;

  enable_a = (a != NULL);
  if(enable_a)
    a3d = *a;
  enable_c = (c != NULL);
  if(enable_c)
    c3d = *c;
 
  Field3D rhs;

  if(precon) {
    /// Get DC components for preconditioner
    aptr = cptr = NULL;
    if(enable_a) {
      a2d = a3d.DC();
      aptr = &a2d;
    }
    if(enable_c) {
      c2d = c3d.DC();
      cptr = &c2d;
    }
    
    // Need to apply preconditioner to rhs
    rhs = invert_laplace(b, flags, aptr, cptr);
  }else
    rhs = b;
  
  int restart=10;
  int itmax=100;
  BoutReal tol=1.e-7;
  
  // Call the solver
  Field3D result;
  result = start;
  solve(rhs, result,
		      flags, 
		      restart, itmax, tol);
  return result;
}

const FieldPerp LaplaceGMRES::function(const FieldPerp &x)
{
  FieldPerp result = Delp2(x);
  
  return x;
}
