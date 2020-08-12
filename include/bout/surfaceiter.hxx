/// \file surfaceiter.hxx
/// Defines a class for iterating over flux surfaces (surfaces of constant x)
/// 

class SurfaceIter;

#ifndef __SURFACEITER_H__
#define __SURFACEITER_H__

#include "mesh.hxx"

/// Iterates over Y-Z surfaces, optionally distributing work between processors
///
/// Example
/// -------
///
/// SurfaceIter si(mesh);
///
/// for( si.first(); !si.isDone(); si.next() ) {
///   // Perform operation at x = si.xpos
///   if(si.closed()) {
///     // A closed flux surface (no boundaries in Y)
///   }else {
///     // Open, so boundaries in Y
///     if(si.firstY()) {
///       // Boundary at lower Y on this processor
///     }
///     if(si.lastY()) {
///       // Boundary at upper Y on this processor
///     }
///   }
/// } 
class SurfaceIter {
 public:
  /// Constructor, needs a mesh to iterate over
  /// @param[in] mesh   The mesh to iterate over
  SurfaceIter(Mesh *mesh, bool include_guards=false)
    : m(mesh), firstpos(include_guards ? 0 : mesh->xstart),
      lastpos(include_guards ? mesh->LocalNx - 1 : mesh->xend) {}
  
  int xpos;    ///< X position where iteration is currently at
  int ySize(); ///< Return the length of the current surface in Y
  
  bool closed(); ///< Test if the current surface is closed

  /// Test if the current surface (x = xpos) is closed
  /// @param[out] ts   The twist-shift angle by which points are shifted in Z between the end and beginning of Y
  bool closed(BoutReal &ts);
  
  MPI_Comm communicator(); ///< Communicator for this surface
  
  int yGlobal(int yloc); ///< Return global y index of given local index \p yloc

  bool firstY(); ///< Is this processor at the lower end?
  bool lastY();  ///< Is this processor at the upper end?

  void first(); ///< Begin iteration
  void next();  ///< Move to next flux surface
  bool isDone(); ///< Are we done iterating?
private:
  Mesh *m; ///< The mesh being iterated over
  const int firstpos;
  const int lastpos;
};

#endif // __SURFACEITER_H__
