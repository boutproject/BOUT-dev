

#ifndef __COORDINATES_H__
#define __COORDINATES_H__

#include "mesh.hxx"
#include "datafile.hxx"

class Coordinates {
public:
  /// Constructor
  Coordinates(Mesh *mesh);
  
  ~Coordinates();
  
  /// Must be a better way so that Coordinates doesn't depend on Datafile
  void outputVars(Datafile &file);

  // Mesh spacing
  
  Field2D dx, dy; 
  BoutReal zlength, dz;
  
  Field2D d1_dx, d1_dy;  // 2nd-order correction for non-uniform meshes d/di(1/dx) and d/di(1/dy)
  
  Field2D J; // Jacobian

  Field2D Bxy; // Magnitude of B = nabla z times nabla x
  
  // Contravariant metric tensor (g^{ij})
  Field2D g11, g22, g33, g12, g13, g23; // These are read in grid.cpp
  
  // Covariant metric tensor
  Field2D g_11, g_22, g_33, g_12, g_13, g_23;
  
  // Christoffel symbol of the second kind (connection coefficients)
  Field2D G1_11, G1_22, G1_33, G1_12, G1_13;
  Field2D G2_11, G2_22, G2_33, G2_12, G2_23;
  Field2D G3_11, G3_22, G3_33, G3_13, G3_23;
  
  Field2D G1, G2, G3;

  /// Calculate differential geometry quantities from the metric tensor
  int geometry();
  int calcCovariant(); ///< Inverts contravatiant metric to get covariant
  int calcContravariant(); ///< Invert covariant metric to get contravariant
  int jacobian(); // Calculate J and Bxy
private:
  int gaussj(BoutReal **a, int n);
  int *indxc, *indxr, *ipiv, ilen;
};

/*
/// Standard coordinate system for tokamak simulations
class TokamakCoordinates : public Coordinates {
public:
  TokamakCoordinates(Mesh *mesh) : Coordinates(mesh) {
    
  }
private:
  
};
*/

#endif // __COORDINATES_H__
