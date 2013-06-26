
class SurfaceIter;
class DistribSurfaceIter;

#ifndef __SURFACEITER_H__
#define __SURFACEITER_H__

#include "mesh.hxx"

/// Iterates over Y-Z surfaces, optionally distributing work between processors
class SurfaceIter {
 public:
  SurfaceIter(Mesh *mesh) : m(mesh) {}
  
  int xpos; // X position
  int ySize(); // Return the size of the current surface
  
  bool closed(); // Test if the current surface is closed
  bool closed(BoutReal &ts);
  
  MPI_Comm communicator(); // Communicator for this surface
  
  int yGlobal(int yloc); // Return global y index of given local index

  bool firstY(); ///< Is this processor at the lower end?
  bool lastY();  ///< Is this processor at the upper end?

  void first();
  void next();
  bool isDone();
private:
  Mesh *m;
};

#endif // __SURFACEITER_H__
