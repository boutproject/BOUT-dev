/*
 *
 *
 */

class GlobalField;
class GlobalField2D;
//class GlobalField3D;

#ifndef __GLOBALFIELD_H__
#define __GLOBALFIELD_H__

#include "mesh.hxx"

class GlobalField {
public:
  GlobalField(Mesh *m, int xsize, int ysize, int zsize) : mesh(m), nx(xsize), ny(ysize), nz(zsize);
  
  virtual bool valid() const = 0;  ///< Is the data valid on any processor?
  virtual bool dataIsLocal() const = 0; ///< Data is on this processor

  // Data access by index
  BoutReal& operator()(int jx, int jy, int jz) {return data[jz + nz*jy + nz*ny*jx];}
  const BoutReal& operator()(int jx, int jy, int jz) const {return data[jz + nz*jy + nz*ny*jx];}

  // Direct data access
  BoutReal* getData() {return data;}
protected:
  Mesh *mesh;
  
  int nx, ny, nz;
  BoutReal *data;
  
  MPI_Comm comm;
  int npes, mype;
};

class GlobalField2D : public GlobalField {
public:
  GlobalField2D(Mesh *m) : GlobalField(m) {}
  
  bool valid() const;
  bool dataIsLocal() const;
  
  void gather(const Field2D &f, int proc=0); ///< Gather all data onto one processor
  const Field2D scatter() const; ///< Scatter data back from one to many processors
  
  /// Assignment from a 2D field. Shorthand for a gather, and must be called on all processors
  /// The scatter assignment operator needs to be a member of Field2D.
  GlobalField& operator=(const Field2D &rhs) {
    gather(rhs);
    return *this;
  }
  
  // Data access by index
  BoutReal& operator()(int jx, int jy) {return (*this)(jx, jy, 0);}
  const BoutReal& operator()(int jx, int jy) const {return (*this)(jx, jy, 0);}
  
protected:
  

private:
  
  
};


/*
class GlobalField3D : public GlobalField {
public:
  
};
*/

#endif // __GLOBALFIELD_H__
