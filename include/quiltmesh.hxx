
#ifndef __QUILTMESH_H__
#define __QUILTMESH_H__

#include "mesh.hxx"

/*!
 * Mesh composed of a patchwork of regions. Generalised
 * code which extends the topology of BoutMesh
 */
class QuiltMesh : public Mesh {
 public:
  ~QuiltMesh();
  
  /// Read in the mesh from data sources
  int load();

  /////////////////////////////////////////////
  // Get data
  
  int get(int &ival, const char *name);
  int get(BoutReal &rval, const char *name);
  
  int get(Field2D &var, const char *name, BoutReal def=0.0);
  int get(Field2D &var, const string &name, BoutReal def=0.0);
  int get(Field3D &var, const char *name);
  int get(Field3D &var, const string &name);
  
  int communicate(FieldGroup &g); // Returns error code
  comm_handle send(FieldGroup &g);  // Return handle
  int wait(comm_handle handle); // Wait for the handle, return error code
  
  // X communications
  bool firstX();
  bool lastX();
  int sendXOut(BoutReal *buffer, int size, int tag);
  int sendXIn(BoutReal *buffer, int size, int tag);
  comm_handle irecvXOut(BoutReal *buffer, int size, int tag);
  comm_handle irecvXIn(BoutReal *buffer, int size, int tag);

  // Y-Z surface gather/scatter operations
  SurfaceIter* iterateSurfaces();
  const Field2D averageY(const Field2D &f);
  bool surfaceClosed(int jx, BoutReal &ts); ///< Test if a surface is closed, and if so get the twist-shift angle
  
  // Boundary region iteration
  RangeIter* iterateBndryLowerY();
  RangeIter* iterateBndryUpperY();
  
  // Boundary regions
  vector<BoundaryRegion*> getBoundaries();

  BoutReal GlobalX(int jx); ///< Continuous X index between 0 and 1
  BoutReal GlobalY(int jy); ///< Continuous Y index (0 -> 1)

  void outputVars(Datafile &file); ///< Add mesh vars to file
 private:
  struct MeshRegion {
    int id;     // Region number
    int x0, y0; // Lower-left corner in grid data
    int nx, ny; // Size of the region
    
  };
  struct MeshDomain {
    int region;  // Region number
    int x0, y0;  // Lower left corner in region
    int nx, ny;  // Number of x and y points
    
    int proc;    // The processor for this domain
    
    Domain *xin, *xout;
    int yup_xsplit;
    
  };
  vector<MeshDomain> domain;
  Domain *mydomain;
  
  const vector<int> readInts(const string &name);
};

#endif // __QUILTMESH_H__

