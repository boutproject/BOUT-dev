/**************************************************************************
 * Interface for mesh classes. Contains standard variables and useful
 * routines.
 * 
 * Changelog
 * =========
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
#include "sys/range.hxx" // RangeIterator

#include <bout/griddata.hxx>

#include <list>

typedef void* comm_handle;

class Mesh {
 public:
  
  Mesh(GridDataSource *s);
  virtual ~Mesh();
  
  static Mesh* create(GridDataSource *source, Options *opt = NULL); ///< Create a Mesh object
  static Mesh* create(Options *opt = NULL);
  
  // Currently need to create and load mesh in separate calls. Will be removed
  virtual int load() {return 1;}
  virtual void outputVars(Datafile &file) {} ///< Output variables to a data file

  // Get routines to request data from mesh file
  virtual int get(int &ival, const string &name); ///< Get an integer
  virtual int get(BoutReal &rval, const string &name); ///< Get a BoutReal number
  
  virtual int get(Field2D &var, const string &name, BoutReal def=0.0);
  virtual int get(Field3D &var, const string &name);
  
  int get(Vector2D &var, const string &name);
  int get(Vector3D &var, const string &name);
  
  // Communications
  
  int communicate(FieldData &f);  // Returns error code
  int communicate(FieldData &f1, FieldData &f2);
  int communicate(FieldData &f1, FieldData &f2, FieldData &f3);
  int communicate(FieldData &f1, FieldData &f2, FieldData &f3, FieldData &f4);
  virtual int communicate(FieldGroup &g) = 0; // Returns error code
  int communicate(FieldPerp &f); // Communicate an X-Z field
  
  virtual comm_handle send(FieldGroup &g) = 0;  // Return handle
  comm_handle send(FieldData &f);   // Send a single field
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

  virtual MPI_Comm getXcomm() const = 0; // Return X communicator
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
  
  
  virtual const Field2D averageY(const Field2D &f) = 0;
  virtual const Field3D averageY(const Field3D &f);
  
  // Boundary region iteration
  virtual const RangeIterator iterateBndryLowerY() const = 0;
  virtual const RangeIterator iterateBndryUpperY() const = 0;
  
  bool hasBndryLowerY();
  bool hasBndryUpperY();

  // Boundary regions
  virtual vector<BoundaryRegion*> getBoundaries() = 0;
  
  // Branch-cut special handling (experimental)
  virtual const Field3D smoothSeparatrix(const Field3D &f) {return f;}
  
  virtual BoutReal GlobalX(int jx) const = 0; ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(int jy) const = 0; ///< Continuous Y index (0 -> 1)
  
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
  virtual BoutReal Average_XY(const Field2D &var)=0;
  virtual const Field3D Switch_YZ(const Field3D &var)=0;
  virtual BoutReal Vol_Integral(const Field2D &var)=0;

  /// Size of the mesh on this processor including guard/boundary cells
  int ngx, ngy, ngz;
  
  /// Local ranges of data (inclusive), excluding guard cells
  int xstart, xend, ystart, yend;

  // These used for differential operators 
  Field2D dx, dy;      // Read in grid.cpp
  Field2D d1_dx, d1_dy;  // 2nd-order correction for non-uniform meshes d/di(1/dx) and d/di(1/dy)
  
  
  BoutReal zlength, dz;    // Derived from options in grid.cpp (in radians)
  
  bool ShiftXderivs; // Use shifted X derivatives
  int  ShiftOrder;   // Order of shifted X derivative interpolation
  Field2D zShift; // Z shift for each point (radians)
  
  int  TwistOrder;   // Order of twist-shift interpolation
  bool BoundaryOnCell; // NB: DOESN'T REALLY BELONG HERE
  bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).
  
  Field2D ShiftTorsion; // d <pitch angle> / dx. Needed for vector differentials (Curl)
  Field2D IntShiftTorsion; // Integrated shear (I in BOUT notation)
  bool IncIntShear; // Include integrated shear (if shifting X)
  
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
  
  bool non_uniform; // Use corrections for non-uniform meshes
  
 protected:
  
  GridDataSource *source; ///< Source for grid data

  /// Read a 1D array of integers
  const vector<int> readInts(const string &name, int n);
  
  /// Calculates the size of a message for a given x and y range
  int msg_len(const vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt);
  
  /// Read a 2D array of data from a data source (e.g. file)
  void read2Dvar(GridDataSource *s, const string &name, 
                 int xs, int ys,  // Source origin
                 int xd, int yd,  // Destination origin
                 int nx, int ny,  // Size of the domain to copy
                 BoutReal **data);
  
 private:
  int gaussj(BoutReal **a, int n);
  int *indxc, *indxr, *ipiv, ilen;
};

#endif // __MESH_H__

