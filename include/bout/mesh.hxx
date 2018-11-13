/*!************************************************************************
 * \file mesh.hxx
 *
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

#include <bout/deriv_store.hxx>
#include "field_data.hxx"
#include "bout_types.hxx"
#include "field2d.hxx"
#include "field3d.hxx"
#include "datafile.hxx"
#include "options.hxx"

#include "fieldgroup.hxx"

#include "boundary_region.hxx"
#include "parallel_boundary_region.hxx"

#include "sys/range.hxx" // RangeIterator

#include <bout/griddata.hxx>

#include "coordinates.hxx"    // Coordinates class

#include "paralleltransform.hxx" // ParallelTransform class

#include "unused.hxx"

#include <bout/region.hxx>

#include <list>
#include <memory>
#include <map>

/// Type used to return pointers to handles
typedef void* comm_handle;

class Mesh {
 public:

  /// Constructor for a "bare", uninitialised Mesh
  /// Only useful for testing
  Mesh() : source(nullptr), options(nullptr) {}

  /// Constructor
  /// @param[in] s  The source to be used for loading variables
  /// @param[in] options  The options section for settings
  Mesh(GridDataSource *s, Options *options);

  /// Destructor
  virtual ~Mesh();

  /// Create a Mesh object
  ///
  /// @param[in] source  The data source to use for loading variables
  /// @param[in] opt     The option section. By default this is "mesh"
  static Mesh *create(GridDataSource *source, Options *opt = nullptr);

  /// Create a Mesh object
  ///
  /// The source is determined by
  ///  1) If "file" is set in the options, read that
  ///  2) If "grid" is set in global options, read that
  ///  3) Use options as data source
  ///
  /// @param[in] opt  Input options. Default is "mesh" section
  static Mesh *create(Options *opt = nullptr);

  /// Loads the mesh values
  /// 
  /// Currently need to create and load mesh in separate calls
  /// because creating Fields uses the global "mesh" pointer
  /// which isn't created until Mesh is constructed
  virtual int load() {return 1;}
  
  /// Add output variables to a data file
  /// These are used for post-processing
  virtual void outputVars(Datafile &UNUSED(file)) {} 
  
  // Get routines to request data from mesh file
  
  /// Get an integer from the input source
  /// 
  /// @param[out] ival  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  ///
  /// @returns zero if successful, non-zero on failure
  int get(int &ival, const string &name);

  /// Get a BoutReal from the input source
  /// 
  /// @param[out] rval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  ///
  /// @returns zero if successful, non-zero on failure
  int get(BoutReal &rval, const string &name); 

  /// Get a Field2D from the input source
  /// including communicating guard cells
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(Field2D &var, const string &name, BoutReal def=0.0);

  /// Get a Field3D from the input source
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero if successful, non-zero on failure
  int get(Field3D &var, const string &name, BoutReal def=0.0, bool communicate=true);

  /// Get a Vector2D from the input source.
  /// If \p var is covariant then this gets three
  /// Field2D variables with "_x", "_y", "_z" appended to \p name
  /// If \p var is contravariant, then "x", "y", "z" are appended to \p name
  ///
  /// By default all fields revert to zero
  ///
  /// @param[in] var  This will be set to the value read
  /// @param[in] name  The name of the vector. Individual fields are read based on this name by appending. See above
  ///
  /// @returns zero always. 
  int get(Vector2D &var, const string &name);

  /// Get a Vector3D from the input source.
  /// If \p var is covariant then this gets three
  /// Field3D variables with "_x", "_y", "_z" appended to \p name
  /// If \p var is contravariant, then "x", "y", "z" are appended to \p name
  ///
  /// By default all fields revert to zero
  ///
  /// @param[in] var  This will be set to the value read
  /// @param[in] name  The name of the vector. Individual fields are read based on this name by appending. See above
  ///
  /// @returns zero always. 
  int get(Vector3D &var, const string &name);

  /// Wrapper for GridDataSource::hasVar
  bool sourceHasVar(const string &name);
  
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

  template <typename... Ts>
  void communicateXZ(Ts&... ts) {
    FieldGroup g(ts...);
    communicateXZ(g);
  }

  /*!
   * Communicate a group of fields
   */
  void communicate(FieldGroup &g);

  /// Communcate guard cells in XZ only
  /// i.e. no Y communication
  ///
  /// @param g  The group of fields to communicate. Guard cells will be modified
  void communicateXZ(FieldGroup &g);

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

  /// Perform communications without waiting for them
  /// to finish. Requires a call to wait() afterwards.
  ///
  /// \param g Group of fields to communicate
  /// \returns handle to be used as input to wait()
  virtual comm_handle send(FieldGroup &g) = 0;  
  virtual int wait(comm_handle handle) = 0; ///< Wait for the handle, return error code

  // non-local communications

  /// Low-level communication routine
  /// Send a buffer of data from this processor to another
  /// This must be matched by a corresponding call to
  /// receiveFromProc on the receiving processor
  ///
  /// @param[in] xproc  X index of processor to send to
  /// @param[in] yproc  Y index of processor to send to
  /// @param[in] buffer A buffer of data to send
  /// @param[in] size   The length of \p buffer
  /// @param[in] tag    A label, must be the same at receive
  virtual MPI_Request sendToProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) = 0;

  /// Low-level communication routine
  /// Receive a buffer of data from another processor
  /// Must be matched by corresponding sendToProc call
  /// on the sending processor
  ///
  /// @param[in] xproc X index of sending processor
  /// @param[in] yproc Y index of sending processor
  /// @param[inout] buffer  The buffer to fill with data. Must already be allocated of length \p size
  /// @param[in] size  The length of \p buffer
  /// @param[in] tag   A label, must be the same as send
  virtual comm_handle receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) = 0;
  
  virtual int getNXPE() = 0; ///< The number of processors in the X direction
  virtual int getNYPE() = 0; ///< The number of processors in the Y direction
  virtual int getXProcIndex() = 0; ///< This processor's index in X direction
  virtual int getYProcIndex() = 0; ///< This processor's index in Y direction
  
  // X communications
  virtual bool firstX() = 0;  ///< Is this processor first in X? i.e. is there a boundary to the left in X?
  virtual bool lastX() = 0; ///< Is this processor last in X? i.e. is there a boundary to the right in X?
  bool periodicX; ///< Domain is periodic in X?

  int NXPE, PE_XIND; ///< Number of processors in X, and X processor index

  /// Send a buffer of data to processor at X index +1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  virtual int sendXOut(BoutReal *buffer, int size, int tag) = 0;

  /// Send a buffer of data to processor at X index -1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  virtual int sendXIn(BoutReal *buffer, int size, int tag) = 0;

  /// Receive a buffer of data from X index +1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  virtual comm_handle irecvXOut(BoutReal *buffer, int size, int tag) = 0;

  /// Receive a buffer of data from X index -1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  virtual comm_handle irecvXIn(BoutReal *buffer, int size, int tag) = 0;

  MPI_Comm getXcomm() {return getXcomm(0);} ///< Return communicator containing all processors in X
  virtual MPI_Comm getXcomm(int jy) const = 0; ///< Return X communicator
  virtual MPI_Comm getYcomm(int jx) const = 0; ///< Return Y communicator
  
  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  virtual bool periodicY(int jx) const;

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  /// \param[out] ts  The Twist-Shift angle if periodic
  virtual bool periodicY(int jx, BoutReal &ts) const = 0;
  
  virtual int ySize(int jx) const; ///< The number of points in Y at fixed X index \p jx

  // Y communications
  virtual bool firstY() const = 0; ///< Is this processor first in Y? i.e. is there a boundary at lower Y?
  virtual bool lastY() const = 0; ///< Is this processor last in Y? i.e. is there a boundary at upper Y?
  virtual bool firstY(int xpos) const = 0; ///< Is this processor first in Y? i.e. is there a boundary at lower Y?
  virtual bool lastY(int xpos) const = 0; ///< Is this processor last in Y? i.e. is there a boundary at upper Y?
  virtual int UpXSplitIndex() = 0;  ///< If the upper Y guard cells are split in two, return the X index where the split occurs
  virtual int DownXSplitIndex() = 0; ///< If the lower Y guard cells are split in two, return the X index where the split occurs

  /// Send data
  virtual int sendYOutIndest(BoutReal *buffer, int size, int tag) = 0;

  /// 
  virtual int sendYOutOutdest(BoutReal *buffer, int size, int tag) = 0;

  ///
  virtual int sendYInIndest(BoutReal *buffer, int size, int tag) = 0;

  ///
  virtual int sendYInOutdest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  virtual comm_handle irecvYOutIndest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  virtual comm_handle irecvYOutOutdest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  virtual comm_handle irecvYInIndest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  virtual comm_handle irecvYInOutdest(BoutReal *buffer, int size, int tag) = 0;
  
  // Boundary region iteration

  /// Iterate over the lower Y boundary
  virtual const RangeIterator iterateBndryLowerY() const = 0;

  /// Iterate over the upper Y boundary
  virtual const RangeIterator iterateBndryUpperY() const = 0;
  virtual const RangeIterator iterateBndryLowerOuterY() const = 0;
  virtual const RangeIterator iterateBndryLowerInnerY() const = 0;
  virtual const RangeIterator iterateBndryUpperOuterY() const = 0;
  virtual const RangeIterator iterateBndryUpperInnerY() const = 0;
  
  bool hasBndryLowerY(); ///< Is there a boundary on the lower guard cells in Y?
  bool hasBndryUpperY(); ///< Is there a boundary on the upper guard cells in Y?

  // Boundary regions

  /// Return a vector containing all the boundary regions on this processor
  virtual vector<BoundaryRegion*> getBoundaries() = 0;

  /// Add a boundary region to this processor
  virtual void addBoundary(BoundaryRegion* UNUSED(bndry)) {}

  /// Get all the parallel (Y) boundaries on this processor 
  virtual vector<BoundaryRegionPar*> getBoundariesPar() = 0;

  /// Add a parallel(Y) boundary to this processor 
  virtual void addBoundaryPar(BoundaryRegionPar* UNUSED(bndry)) {}
  
  /// Branch-cut special handling (experimental)
  virtual const Field3D smoothSeparatrix(const Field3D &f) {return f;}
  
  virtual BoutReal GlobalX(int jx) const = 0; ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(int jy) const = 0; ///< Continuous Y index (0 -> 1)
  virtual BoutReal GlobalX(BoutReal jx) const = 0; ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(BoutReal jy) const = 0; ///< Continuous Y index (0 -> 1)
  
  //////////////////////////////////////////////////////////
  
  int GlobalNx, GlobalNy, GlobalNz; ///< Size of the global arrays. Note: can have holes
  int OffsetX, OffsetY, OffsetZ;    ///< Offset of this mesh within the global array
                                    ///< so startx on this processor is OffsetX in global
  
  /// Returns the global X index given a local indexs
  /// If the local index includes the boundary cells, then so does the global.
  virtual int XGLOBAL(int xloc) const = 0;
  /// Returns the global Y index given a local index
  /// The local index must include the boundary, the global index does not.
  virtual int YGLOBAL(int yloc) const = 0;

  /// Size of the mesh on this processor including guard/boundary cells
  int LocalNx, LocalNy, LocalNz;
  
  /// Local ranges of data (inclusive), excluding guard cells
  int xstart, xend, ystart, yend;
  
  bool StaggerGrids;    ///< Enable staggered grids (Centre, Lower). Otherwise all vars are cell centred (default).
  
  bool IncIntShear; ///< Include integrated shear (if shifting X)

  /// Coordinate system
  Coordinates *getCoordinates(const CELL_LOC location = CELL_CENTRE) {
    ASSERT1(location != CELL_DEFAULT);
    ASSERT1(location != CELL_VSHIFT);

    if (coords_map.count(location)) { // True branch most common, returns immediately
      return coords_map[location].get();
    } else {
      // No coordinate system set. Create default
      // Note that this can't be allocated here due to incomplete type
      // (circular dependency between Mesh and Coordinates)
      coords_map.emplace(location, createDefaultCoordinates(location));
      return coords_map[location].get();
    }
  }

  /// Returns the non-CELL_CENTRE location
  /// allowed as a staggered location
  template<DIRECTION direction>
  CELL_LOC getAllowedStaggerLoc() const;

  /// Returns the number of grid points in the
  /// particular direction
  template<DIRECTION direction>
  int getNpoints() const;

  /// Returns the number of guard points in the
  /// particular direction
  template<DIRECTION direction>
  int getNguard() const;
  
  ///////////////////////////////////////////////////////////
  // INDEX DERIVATIVE OPERATORS
  ///////////////////////////////////////////////////////////

  ////// Utilties and parameters
  
  BoutReal fft_derivs_filter; ///< Fraction of modes to filter. This is set in derivs_init from option "ddz:fft_filter"

  /// Determines the resultant output stagger location in derivatives
  /// given the input and output location. Also checks that the
  /// combination of locations is allowed
  STAGGER getStagger(const CELL_LOC inloc, const CELL_LOC outloc, const CELL_LOC allowedloc) const;
  
  /// Determines the resultant output stagger location in derivatives
  /// given the input and output location. Also checks that the
  /// combination of locations is allowed. This overload also checks
  /// the location of a second input field (velocity) is consistent.
  STAGGER getStagger(const CELL_LOC vloc, const CELL_LOC inloc, const CELL_LOC outloc, const CELL_LOC allowedloc) const;  
  
  ////// STANDARD OPERATORS

  ////////////// X DERIVATIVE /////////////////
  template<typename T>
  T indexDDX(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexStandardDerivative<T, DIRECTION::X, 1>(f, outloc, method, region);
  }

  template<typename T>
  T indexD2DX2(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexStandardDerivative<T, DIRECTION::X, 2>(f, outloc, method, region);
  }

  template<typename T>
  T indexD4DX4(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexStandardDerivative<T, DIRECTION::X, 4>(f, outloc, method, region);
  }

  ////////////// Y DERIVATIVE /////////////////

  template<typename T>
  T indexDDY(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    if (std::is_base_of<Field3D, T>::value && f.hasYupYdown() &&
      ((&f.yup() != &f) || (&f.ydown() != &f))) {    
      return indexStandardDerivative<T, DIRECTION::YOrthogonal, 1>(f, outloc, method, region);
    } else {
      const T f_aligned = f.getMesh()->toFieldAligned(f);
      T result = indexStandardDerivative<T, DIRECTION::Y, 1>(f, outloc, method, region);
      return f.getMesh()->fromFieldAligned(result); 
    }
  }

  template<typename T>
  T indexD2DY2(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    if (std::is_base_of<Field3D, T>::value && f.hasYupYdown() &&
      ((&f.yup() != &f) || (&f.ydown() != &f))) {    
      return indexStandardDerivative<T, DIRECTION::YOrthogonal, 2>(f, outloc, method, region);
    } else {
      const T f_aligned = f.getMesh()->toFieldAligned(f);
      T result = indexStandardDerivative<T, DIRECTION::Y, 2>(f, outloc, method, region);
      return f.getMesh()->fromFieldAligned(result); 
    }
  }

  template<typename T>
  T indexD4DY4(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    if (std::is_base_of<Field3D, T>::value && f.hasYupYdown() &&
      ((&f.yup() != &f) || (&f.ydown() != &f))) {    
      return indexStandardDerivative<T, DIRECTION::YOrthogonal, 4>(f, outloc, method, region);  
    } else {
      const T f_aligned = f.getMesh()->toFieldAligned(f);
      T result = indexStandardDerivative<T, DIRECTION::Y, 4>(f, outloc, method, region);
      return f.getMesh()->fromFieldAligned(result); 
    }
  }

  ////////////// Z DERIVATIVE /////////////////
  template<typename T>
  T indexDDZ(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexStandardDerivative<T, DIRECTION::Z, 1>(f, outloc, method, region);
  }

  template<typename T>
  T indexD2DZ2(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexStandardDerivative<T, DIRECTION::Z, 2>(f, outloc, method, region);
  }

  template<typename T>
  T indexD4DZ4(const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexStandardDerivative<T, DIRECTION::Z, 4>(f, outloc, method, region);
  }

  ////// ADVECTION AND FLUX OPERATORS

  /// Advection operator in index space in [] direction
  ///
  /// \f[
  ///   v \frac{d}{di} f
  /// \f]
  ///
  /// @param[in] v  The velocity in the Y direction
  /// @param[in] f  The field being advected
  /// @param[in] outloc The cell location where the result is desired. The default is the same as \p f
  /// @param[in] method  The differencing method to use
  /// @param[in] region  The region of the grid for which the result is calculated.

  ////////////// X DERIVATIVE /////////////////

  template<typename T>
  T indexVDDX(const T& vel, const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexFlowDerivative<T, DIRECTION::X, DERIV::Upwind>(vel, f, outloc, method, region);
  }

  template<typename T>
  T indexFDDX(const T& vel, const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexFlowDerivative<T, DIRECTION::X, DERIV::Flux>(vel, f, outloc, method, region);
  }

  ////////////// Y DERIVATIVE /////////////////

  template<typename T>
  T indexVDDY(const T& vel, const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    bool fHasParallelSlices = (std::is_base_of<Field3D, T>::value && f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f)));
    bool velHasParallelSlices = (std::is_base_of<Field3D, T>::value && vel.hasYupYdown() && ((&vel.yup() != &vel) || (&vel.ydown() != &vel)));
    if (fHasParallelSlices && velHasParallelSlices) {
      return indexFlowDerivative<T, DIRECTION::YOrthogonal, DERIV::Upwind>(vel, f, outloc, method, region);
    } else {
      const T f_aligned = f.getMesh()->toFieldAligned(f);
      const T vel_aligned = vel.getMesh()->toFieldAligned(vel);      
      T result = indexFlowDerivative<T, DIRECTION::Y, DERIV::Upwind>(vel_aligned, f_aligned, outloc, method, region);
      return f.getMesh()->fromFieldAligned(result); 
    }
  }

  template<typename T>
  T indexFDDY(const T& vel, const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    bool fHasParallelSlices = (std::is_base_of<Field3D, T>::value && f.hasYupYdown() && ((&f.yup() != &f) || (&f.ydown() != &f)));
    bool velHasParallelSlices = (std::is_base_of<Field3D, T>::value && vel.hasYupYdown() && ((&vel.yup() != &vel) || (&vel.ydown() != &vel)));
    if (fHasParallelSlices && velHasParallelSlices) {
      return indexFlowDerivative<T, DIRECTION::YOrthogonal, DERIV::Flux>(vel, f, outloc, method, region);
    } else {
      const T f_aligned = f.getMesh()->toFieldAligned(f);
      const T vel_aligned = vel.getMesh()->toFieldAligned(vel);      
      T result = indexFlowDerivative<T, DIRECTION::Y, DERIV::Flux>(vel_aligned, f_aligned, outloc, method, region);
      return f.getMesh()->fromFieldAligned(result); 
    }
  }

  ////////////// Z DERIVATIVE /////////////////

  template<typename T>
  T indexVDDZ(const T& vel, const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexFlowDerivative<T, DIRECTION::Z, DERIV::Upwind>(vel, f, outloc, method, region);
  }

  template<typename T>
  T indexFDDZ(const T& vel, const T &f, CELL_LOC outloc = CELL_DEFAULT, DIFF_METHOD method = DIFF_DEFAULT, REGION region = RGN_NOBNDRY) const {
    return indexFlowDerivative<T, DIRECTION::Z, DERIV::Flux>(vel, f, outloc, method, region);
  }
  
  ///////////////////////////////////////////////////////////
  // PARALLEL TRANSFORMS
  ///////////////////////////////////////////////////////////

  /// Transform a field into field-aligned coordinates
  const Field3D toFieldAligned(const Field3D &f) {
    return getParallelTransform().toFieldAligned(f);
  }
  const Field2D toFieldAligned(const Field2D &f) {
    return f;
  }
  
  /// Convert back into standard form
  const Field3D fromFieldAligned(const Field3D &f) {
    return getParallelTransform().fromFieldAligned(f);
  }
  const Field2D fromFieldAligned(const Field2D &f) {
    return f;
  }

  bool canToFromFieldAligned() {
    return getParallelTransform().canToFromFieldAligned();
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
  /*!
   * Set the parallel (y) transform from the options file
   */
  void setParallelTransform();

  /*!
   * Return the parallel transform, setting it if need be
   */
  ParallelTransform& getParallelTransform();
  
  ///////////////////////////////////////////////////////////
  // REGION RELATED ROUTINES
  ///////////////////////////////////////////////////////////

  // The maxregionblocksize to use when creating the default regions.
  // Can be set in the input file and the global default is set by,
  // MAXREGIONBLOCKSIZE in include/bout/region.hxx
  int maxregionblocksize;
  
  /// Get the named region from the region_map for the data iterator
  ///
  /// Throws if region_name not found
  const Region<> &getRegion(const std::string &region_name) const{
    return getRegion3D(region_name);
  }
  const Region<Ind3D> &getRegion3D(const std::string &region_name) const;
  const Region<Ind2D> &getRegion2D(const std::string &region_name) const;
  const Region<IndPerp> &getRegionPerp(const std::string &region_name) const;

  /// Add a new region to the region_map for the data iterator
  ///
  /// Outputs an error message if region_name already exists
  void addRegion(const std::string &region_name, const Region<> &region) {
    return addRegion3D(region_name, region);
  }
  void addRegion(const std::string &region_name, const Region<Ind2D> &region) {
    return addRegion2D(region_name, region);
  }
  void addRegion(const std::string &region_name, const Region<IndPerp> &region) {
    return addRegionPerp(region_name, region);
  }
  void addRegion3D(const std::string &region_name, const Region<Ind3D> &region);
  void addRegion2D(const std::string &region_name, const Region<Ind2D> &region);
  void addRegionPerp(const std::string &region_name, const Region<IndPerp> &region);

  /// Converts an Ind2D to an Ind3D using calculation
  Ind3D ind2Dto3D(const Ind2D &ind2D, int jz = 0) { return {ind2D.ind * LocalNz + jz, LocalNy, LocalNz}; }

  /// Converts an Ind3D to an Ind2D using calculation
  Ind2D ind3Dto2D(const Ind3D &ind3D) { return {ind3D.ind / LocalNz, LocalNy, 1}; }

  /// Converts an Ind3D to an IndPerp using calculation
  IndPerp ind3DtoPerp(const Ind3D &ind3D) { return {ind3D.x() * LocalNz + ind3D.z(), 1, LocalNz}; }

  /// Converts an IndPerp to an Ind3D using calculation
  Ind3D indPerpto3D(const IndPerp &indPerp, int jy = 0) {
    int jz = indPerp.z();
    return { (indPerp.ind - jz) * LocalNy + LocalNz * jy + jz , LocalNy, LocalNz};
  }
  
  /// Converts an Ind3D to an Ind2D representing a 2D index using a lookup -- to be used with care
  Ind2D map3Dto2D(const Ind3D &ind3D){
    return {indexLookup3Dto2D[ind3D.ind], LocalNy, 1};
  }
  
  /// Create the default regions for the data iterator
  ///
  /// Creates RGN_{ALL,NOBNDRY,NOX,NOY}
  void createDefaultRegions();
    
 protected:
  
  GridDataSource *source; ///< Source for grid data
  
  std::map<CELL_LOC, std::shared_ptr<Coordinates> > coords_map; ///< Coordinate systems at different CELL_LOCs

  Options *options; ///< Mesh options section
  

  PTptr transform; ///< Handles calculation of yup and ydown

  /// Read a 1D array of integers
  const vector<int> readInts(const string &name, int n);
  
  /// Calculates the size of a message for a given x and y range
  int msg_len(const vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt);
  
  /// Initialise derivatives
  void derivs_init(Options* options);
  
private:
  
  /// The main kernel used for all standard derivatives
  template<typename T, DIRECTION direction, int order>
  T indexStandardDerivative(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) const;

  /// The main kernel used for all upwind and flux derivatives
  template<typename T, DIRECTION direction, DERIV derivType>
  T indexFlowDerivative(const T &vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) const;

  /// Allocates default Coordinates objects
  std::shared_ptr<Coordinates> createDefaultCoordinates(const CELL_LOC location);

  //Internal region related information
  std::map<std::string, Region<Ind3D>> regionMap3D;
  std::map<std::string, Region<Ind2D>> regionMap2D;
  std::map<std::string, Region<IndPerp>> regionMapPerp;
  Array<int> indexLookup3Dto2D;
};

template<typename T, DIRECTION direction, DERIV derivType>
T Mesh::indexFlowDerivative(const T &vel, const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) const {
  TRACE("%s", __thefunc__);
  
  // Checks
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "indexDDX only works on Field2D or Field3D input");

  // Special handling for SPLIT method
  if ((derivType == DERIV::Flux) && (method == DIFF_SPLIT)) {
    // Split into an upwind and a central differencing part
    // d/dx(v*f) = v*d/dx(f) + f*d/dx(v)
    auto tmp = indexFlowDerivative<T, direction, DERIV::Upwind>(vel, f, outloc, DIFF_DEFAULT, region);
    tmp += indexStandardDerivative<T, direction, 1>(vel, outloc, DIFF_DEFAULT, region) * interp_to(f, tmp.getLocation());
    return tmp;      
  }
  
  // Check that the mesh is correct
  ASSERT1(this == f.getMesh());
  ASSERT1(this == vel.getMesh());  
  // Check that the input variable has data
  ASSERT1(f.isAllocated());
  ASSERT1(vel.isAllocated());  
  
  //Check the input data is valid
  {
    TRACE("Checking inputs");
    checkData(f);
    checkData(vel);    
  }

  // Define properties of this approach
  const CELL_LOC allowedStaggerLoc = getAllowedStaggerLoc<direction>();

  // Handle the staggering
  const CELL_LOC inloc = f.getLocation(); // Input locations
  const CELL_LOC vloc = vel.getLocation();
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  const STAGGER stagger = getStagger(vloc, inloc, outloc, allowedStaggerLoc);

  // Check for early exit
  const int nPoint = getNpoints<direction>();

  if (nPoint == 1) {
    auto tmp = T(0., f.getMesh());
    tmp.setLocation(outloc);
    return tmp;
  }
  
  // Lookup the method
  auto& derivativeStore = DerivativeStore<T>::getInstance();
  typename DerivativeStore<T>::upwindFunc derivativeMethod;
  if (derivType == DERIV::Upwind) {
    derivativeMethod = derivativeStore.getUpwindDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else if (derivType == DERIV::Flux) {
    derivativeMethod = derivativeStore.getFluxDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else {
    throw BoutException("Invalid derivative type in call to indexFlowDerivative.");
  }
  
  // Create the result field
  T result(f.getMesh());
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  // Apply method
  derivativeMethod(vel, f, result, region);

  //Check the result is valid
  {
    TRACE("Checking result");
    checkData(result);
  }
  
  return result;
}

template<typename T, DIRECTION direction, int order>
T Mesh::indexStandardDerivative(const T &f, CELL_LOC outloc, DIFF_METHOD method, REGION region) const {
  TRACE("%s", __thefunc__);
  
  // Checks
  static_assert(std::is_base_of<Field2D, T>::value || std::is_base_of<Field3D, T>::value,
                "indexDDX only works on Field2D or Field3D input");
  // Check that the mesh is correct
  ASSERT1(this == f.getMesh());
  // Check that the input variable has data
  ASSERT1(f.isAllocated());

  //Check the input data is valid
  {
    TRACE("Checking input");
    checkData(f);
  }
  
  // Define properties of this approach
  const CELL_LOC allowedStaggerLoc = getAllowedStaggerLoc<direction>();

  // Handle the staggering
  const CELL_LOC inloc = f.getLocation(); // Input location
  if (outloc == CELL_DEFAULT)
    outloc = inloc;
  const STAGGER stagger = getStagger(inloc, outloc, allowedStaggerLoc);

  // Check for early exit
  const int nPoint = getNpoints<direction>();

  if (nPoint == 1) {
    auto tmp = T(0., f.getMesh());
    tmp.setLocation(outloc);
    return tmp;
  }
  
  // Lookup the method
  auto& derivativeStore = DerivativeStore<T>::getInstance();
  typename DerivativeStore<T>::standardFunc derivativeMethod;
  
  if (order == 1) {
    derivativeMethod = derivativeStore.getStandardDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else if (order == 2) {
    derivativeMethod = derivativeStore.getStandard2ndDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else if (order == 4) {
    derivativeMethod = derivativeStore.getStandard4thDerivative(DIFF_METHOD_STRING(method), direction, stagger);
  } else {
    throw BoutException("Invalid order used in indexStandardDerivative.");
  }
  
  // Create the result field
  T result(f.getMesh());
  result.allocate(); // Make sure data allocated
  result.setLocation(outloc);

  // Apply method
  derivativeMethod(f, result, region);

  //Check the result is valid
  {
    TRACE("Checking result");
    checkData(result);
  }
  
  return result;
}

/*******************************************************************************
 * Helper routines
 *******************************************************************************/

template<DIRECTION direction>
CELL_LOC Mesh::getAllowedStaggerLoc() const {
  switch(direction) {
  case(DIRECTION::X):
    return CELL_XLOW;
  case(DIRECTION::Y):
  case(DIRECTION::YOrthogonal):
  case(DIRECTION::YAligned):    
    return CELL_YLOW;
  case(DIRECTION::Z):
    return CELL_ZLOW;
  }
};


template<DIRECTION direction>
int Mesh::getNpoints() const {
  switch(direction) {
  case(DIRECTION::X):
    return LocalNx;
  case(DIRECTION::Y):
  case(DIRECTION::YOrthogonal):
  case(DIRECTION::YAligned):    
    return LocalNy;
  case(DIRECTION::Z):
    return LocalNz;
  }
};

template<DIRECTION direction>
int Mesh::getNguard() const {
  switch(direction) {
  case(DIRECTION::X):
    return xstart;
  case(DIRECTION::Y):
  case(DIRECTION::YOrthogonal):
  case(DIRECTION::YAligned):    
    return ystart;
  case(DIRECTION::Z):
    return 2;
  }
};

#endif // __MESH_H__

