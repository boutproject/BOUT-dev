
#ifndef __QUILTMESH_H__
#define __QUILTMESH_H__

#include "mesh.h"

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
 private:
  struct Domain {
    int region;  // Region number
    int x0, y0;  // Lower left corner in region
    int nx, ny;  // Number of x and y points
    
    int proc;    // The processor for this domain
    
    Domain *xin, *xout;
    int yup_xsplit;
    
  };
  vector<Domain> domain;
  Domain *mydomain;
  
  const vector<int> readInts(const string &name);
};

#endif // __QUILTMESH_H__

