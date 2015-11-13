/**************************************************************************
 * Interface for mesh classes. Contains standard variables and useful
 * routines.
 * 
 * Changelog
 * =========
 *
 * 2014-12 Ben Dudson <bd512@york.ac.uk>
 *     * Removing coordinate system into separate
 *       Coordinates class
 *     * Adding index derivative functions from derivs.cxx
 * 
 * 2010-06 Ben Dudson, Sean Farley
 *     * Initial version, adapted from GridData class
 *     * Incorporates code from topology.cpp and Communicator
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

class Mesh;

#ifndef __MESH_H__
#define __MESH_H__

#include "mpi.h"

#include "field_data.hxx"
#include "bout_types.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "datafile.hxx"
#include "options.hxx"
#include "dcomplex.hxx" // For poloidal lowpass filter

#include "fieldgroup.hxx"

#include "boundary_region.hxx"
#include "parallel_boundary_region.hxx"

#include "sys/range.hxx" // RangeIterator

#include "bout/deprecated.hxx"

#include <bout/griddata.hxx>

#include "coordinates.hxx"    // Coordinates class

#include "paralleltransform.hxx" // ParallelTransform class

#include <list>
#include <memory>

typedef void* comm_handle;

class Mesh {
 public:
  
  Mesh(GridDataSource *s, Options *options);
  virtual ~Mesh();
  
  static Mesh* create(GridDataSource *source, Options *opt = NULL); ///< Create a Mesh object
  static Mesh* create(Options *opt = NULL);
  
  // Currently need to create and load mesh in separate calls. Will be removed
  virtual int load() {return 1;}
  virtual void outputVars(Datafile &file) {} ///< Output variables to a data file

  
  // Get routines to request data from mesh file
  int get(int &ival, const string &name); ///< Get an integer
  int get(BoutReal &rval, const string &name); ///< Get a BoutReal number
  
  int get(Field2D &var, const string &name, BoutReal def=0.0);
  int get(Field3D &var, const string &name, BoutReal def=0.0);
  
  int get(Vector2D &var, const string &name);
  int get(Vector3D &var, const string &name);
  
  // Communications
  /*!
   * Communicate a list of FieldData objects
   * Uses a variadic template (C++11) to pack all
   * arguments into a FieldGroup
   */
  template <typename... Ts>
  void communicate(Ts&... ts) {
    FieldGroup g(ts...);
    communicate(g);
  }

  /*!
   * Communicate a group of fields
   */
  void communicate(FieldGroup &g);

  /*!
   * Communicate an X-Z field
   */
  void communicate(FieldPerp &f); 

  /*!
   * Send a list of FieldData objects
   * Packs arguments into a FieldGroup and passes
   * to send(FieldGroup&).
   */
  template <typename... Ts>
  comm_handle send(Ts&... ts) {
    FieldGroup g(ts...);
    return send(g);
  }
  
  virtual comm_handle send(FieldGroup &g) = 0;  // Return handle
  virtual int wait(comm_handle handle) = 0; // Wait for the handle, return error code

  // non-local communications
  virtual MPI_Request sendToProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) = 0;
  virtual int getNXPE() = 0;
  virtual int getNYPE() = 0;
  virtual int getXProcIndex() = 0;
  virtual int getYProcIndex() = 0;
  
  // X communications
  virtual bool firstX() = 0;
  virtual bool lastX() = 0;
  bool periodicX; // Domain is periodic in X?

  int NXPE, PE_XIND; ///< Number of processors in X, and X processor index
  virtual int sendXOut(BoutReal *buffer, int size, int tag) = 0;
  virtual int sendXIn(BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle irecvXOut(BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle irecvXIn(BoutReal *buffer, int size, int tag) = 0;

  DEPRECATED(MPI_Comm getXcomm()) {return getXcomm(0);}
  virtual MPI_Comm getXcomm(int jy) const = 0; // Return X communicator
  virtual MPI_Comm getYcomm(int jx) const = 0; // Return Y communicator

  //virtual bool periodicX() const = 0;                     ///< Test if a surface is periodic in X
  virtual bool periodicY(int jx) const;                   ///< Test if a surface is closed in Y
  virtual bool periodicY(int jx, BoutReal &ts) const = 0; ///< Also get the twist-shift angle
  
  virtual int ySize(int jx) const;

  // Y communications
  virtual bool firstY() = 0;
  virtual bool lastY() = 0;
  virtual int UpXSplitIndex() = 0;
  virtual int DownXSplitIndex() = 0;
  virtual int sendYOutIndest(BoutReal *buffer, int size, int tag) = 0;
  virtual int sendYOutOutdest(BoutReal *buffer, int size, int tag) = 0;
  virtual int sendYInIndest(BoutReal *buffer, int size, int tag) = 0;
  virtual int sendYInOutdest(BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle irecvYOutIndest(BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle irecvYOutOutdest(BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle irecvYInIndest(BoutReal *buffer, int size, int tag) = 0;
  virtual comm_handle irecvYInOutdest(BoutReal *buffer, int size, int tag) = 0;
  
  // Boundary region iteration
  virtual const RangeIterator iterateBndryLowerY() const = 0;
  virtual const RangeIterator iterateBndryUpperY() const = 0;
  
  bool hasBndryLowerY();
  bool hasBndryUpperY();

  // Boundary regions
  virtual vector<BoundaryRegion*> getBoundaries() = 0;
  virtual void addBoundary(BoundaryRegion* bndry) {}
  virtual vector<BoundaryRegionPar*> getBoundariesPar() = 0;
  virtual void addBoundaryPar(BoundaryRegionPar* bndry) {}
  
  // Branch-cut special handling (experimental)
  virtual const Field3D smoothSeparatrix(const Field3D &f) {return f;}
  
  virtual BoutReal GlobalX(int jx) const = 0; ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(int jy) const = 0; ///< Continuous Y index (0 -> 1)
  virtual BoutReal GlobalX(BoutReal jx) const = 0; ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(BoutReal jy) const = 0; ///< Continuous Y index (0 -> 1)
  
  //////////////////////////////////////////////////////////
  
  int GlobalNx, GlobalNy, GlobalNz; // Size of the global arrays. Note: can have holes
  int OffsetX, OffsetY, OffsetZ;    // Offset of this mesh within the global array
                                    // so startx on this processor is OffsetX in global
  
  /// Global locator functions
  virtual int XGLOBAL(int xloc) const = 0;
  virtual int YGLOBAL(int yloc) const = 0;

  // poloidal lowpass filter for n=0 mode
  virtual void slice_r_y(BoutReal *, BoutReal *, int , int)=0;
  virtual void get_ri( dcomplex *, int, BoutReal *, BoutReal *)=0;
  virtual void set_ri( dcomplex *, int, BoutReal *, BoutReal *)=0;
  virtual const Field2D lowPass_poloidal(const Field2D &,int)=0;
  
  /// volume integral
  virtual const Field3D Switch_YZ(const Field3D &var) = 0;
  virtual const Field3D Switch_XZ(const Field3D &var) = 0;

  /// Size of the mesh on this processor including guard/boundary cells
  int ngx, ngy, ngz;
  
  /// Local ranges of data (inclusive), excluding guard cells
  int xstart, xend, ystart, yend;
  
  bool ShiftXderivs; // Use shifted X derivatives
  int  ShiftOrder;   // Order of shifted X derivative interpolation
  Field2D zShift; // Z shift for each point (radians)
  
  bool FCI; ///< Using Flux Coordinate Independent (FCI) method?
  
  bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).
  
  bool IncIntShear; // Include integrated shear (if shifting X)
  
  /// Coordinate system
  Coordinates* coordinates();
  
  bool freeboundary_xin, freeboundary_xout, freeboundary_ydown, freeboundary_yup;

  // First derivatives in index space
  // Implemented in src/mesh/index_derivs.hxx
  const Field3D indexDDX(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field2D indexDDX(const Field2D &f);
  
  const Field3D indexDDY(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field2D indexDDY(const Field2D &f);

  const Field3D indexDDZ(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method, bool inc_xbndry);
  const Field2D indexDDZ(const Field2D &f);

  // Second derivatives in index space
  // Implemented in src/mesh/index_derivs.hxx
  const Field3D indexD2DX2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field2D indexD2DX2(const Field2D &f);
  const Field3D indexD2DY2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field2D indexD2DY2(const Field2D &f);
  const Field3D indexD2DZ2(const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  
  // Fourth derivatives in index space
  const Field3D indexD4DX4(const Field3D &f);
  const Field2D indexD4DX4(const Field2D &f);
  const Field3D indexD4DY4(const Field3D &f);
  const Field2D indexD4DY4(const Field2D &f);
  const Field3D indexD4DZ4(const Field3D &f);
  
  // Advection schemes
  const Field2D indexVDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field3D indexVDDX(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field2D indexVDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field3D indexVDDY(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field3D indexVDDZ(const Field &v, const Field &f, CELL_LOC outloc, DIFF_METHOD method);

  const Field2D indexFDDX(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field3D indexFDDX(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field2D indexFDDY(const Field2D &v, const Field2D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field3D indexFDDY(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  const Field3D indexFDDZ(const Field3D &v, const Field3D &f, CELL_LOC outloc, DIFF_METHOD method);
  
  typedef BoutReal (*deriv_func)(stencil &); // f
  typedef BoutReal (*upwind_func)(stencil &, stencil &); // v, f

  typedef struct {
    BoutReal inner;
    BoutReal outer;
  } boundary_derivs_pair;
  // More types for forward/backward differences to calculate derivatives in boundary guard cells for free boundary conditions
  typedef boundary_derivs_pair (*inner_boundary_deriv_func)(forward_stencil &); // f
  typedef boundary_derivs_pair (*outer_boundary_deriv_func)(backward_stencil &); // f
  typedef boundary_derivs_pair (*inner_boundary_upwind_func)(forward_stencil &); // v,f
  typedef boundary_derivs_pair (*outer_boundary_upwind_func)(backward_stencil &); // v,f

  /// Transform a field into field-aligned coordinates
  const Field3D toFieldAligned(const Field3D &f) {
    return getParallelTransform().toFieldAligned(f);
  }
  /// Convert back into standard form
  const Field3D fromFieldAligned(const Field3D &f) {
    return getParallelTransform().fromFieldAligned(f);
  }

  /*!
   * Unique pointer to ParallelTransform object
   */
  typedef std::unique_ptr<ParallelTransform> PTptr;
  
  /*!
   * Set the parallel (y) transform for this mesh.
   * Unique pointer used so that ParallelTransform will be deleted
   */
  void setParallelTransform(PTptr pt) {
    transform = std::move(pt);
  }
  
  
 protected:
  
  GridDataSource *source; ///< Source for grid data
  
  Coordinates *coords;    ///< Coordinate system. Initialised to Null

  Options *options; ///< Mesh options section
  
  /*!
   * 
   */
  ParallelTransform& getParallelTransform();
  
  PTptr transform; ///< Handles calculation of yup and ydown

  /// Read a 1D array of integers
  const vector<int> readInts(const string &name, int n);
  
  /// Calculates the size of a message for a given x and y range
  int msg_len(const vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt);
  
  // Initialise derivatives
  void derivs_init(Options* options);
  
  // Loop over mesh, applying a stencil in the X direction
  const Field2D applyXdiff(const Field2D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);
  const Field3D applyXdiff(const Field3D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);
  
  const Field2D applyYdiff(const Field2D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);
  const Field3D applyYdiff(const Field3D &var, deriv_func func, inner_boundary_deriv_func func_in, outer_boundary_deriv_func func_out, CELL_LOC loc = CELL_DEFAULT);

  const Field3D applyZdiff(const Field3D &var, Mesh::deriv_func func, CELL_LOC loc = CELL_DEFAULT);
  
private:
  
};

#endif // __MESH_H__

