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

#include <bout/deprecated.hxx>
#include <bout/deriv_store.hxx>
#include <bout/index_derivs_interface.hxx>
#include <bout/mpi_wrapper.hxx>

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

#include "unused.hxx"

#include <bout/region.hxx>
#include "bout/generic_factory.hxx"

#include <list>
#include <memory>
#include <map>

class MeshFactory : public Factory<
  Mesh, MeshFactory,
  std::function<std::unique_ptr<Mesh>(GridDataSource*, Options*)>> {
public:
  static constexpr auto type_name = "Mesh";
  static constexpr auto section_name = "mesh";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "bout";

  ReturnType create(Options* options = nullptr, GridDataSource* source = nullptr);
};

template <class DerivedType>
class RegisterMesh {
public:
  RegisterMesh(const std::string& name) {
    MeshFactory::getInstance().add(
        name, [](GridDataSource* source, Options* options) -> std::unique_ptr<Mesh> {
          return std::make_unique<DerivedType>(source, options);
        });
  }
};

/// Type used to return pointers to handles
using comm_handle = void*;

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
  
  /// Get a string from the input source
  /// 
  /// @param[out] sval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(std::string& sval, const std::string& name, const std::string& def="");

  /// Get an integer from the input source
  /// 
  /// @param[out] ival  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(int &ival, const std::string &name, int def=0);

  /// Get a BoutReal from the input source
  /// 
  /// @param[out] rval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(BoutReal& rval, const std::string& name, BoutReal def=0.0);

  /// Get a bool from the input source
  ///
  /// @param[out] bval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(bool &bval, const std::string &name, bool def=false);

  /// Get a Field2D from the input source
  /// including communicating guard cells
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(Field2D &var, const std::string &name, BoutReal def=0.0);

  /// Get a Field3D from the input source
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero if successful, non-zero on failure
  int get(Field3D &var, const std::string &name, BoutReal def=0.0, bool communicate=true);

  /// Get a FieldPerp from the input source
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero if successful, non-zero on failure
  int get(FieldPerp &var, const std::string &name, BoutReal def=0.0, bool communicate=true);

  /// Get a Vector2D from the input source.
  /// If \p var is covariant then this gets three
  /// Field2D variables with "_x", "_y", "_z" appended to \p name
  /// If \p var is contravariant, then "x", "y", "z" are appended to \p name
  ///
  /// By default all fields revert to zero
  ///
  /// @param[in] var  This will be set to the value read
  /// @param[in] name  The name of the vector. Individual fields are read based on this name by appending. See above
  /// @param[in] def   The default value if not found (used for all the components)
  ///
  /// @returns zero always. 
  int get(Vector2D &var, const std::string &name, BoutReal def=0.0);

  /// Get a Vector3D from the input source.
  /// If \p var is covariant then this gets three
  /// Field3D variables with "_x", "_y", "_z" appended to \p name
  /// If \p var is contravariant, then "x", "y", "z" are appended to \p name
  ///
  /// By default all fields revert to zero
  ///
  /// @param[in] var  This will be set to the value read
  /// @param[in] name  The name of the vector. Individual fields are read based on this name by appending. See above
  /// @param[in] def    The default value if not found (used for all the components)
  ///
  /// @returns zero always. 
  int get(Vector3D &var, const std::string &name, BoutReal def=0.0);

  /// Test if input source was a grid file
  bool isDataSourceGridFile() const;

  /// Wrapper for GridDataSource::hasVar
  bool sourceHasVar(const std::string &name);

  /// Wrapper for GridDataSource::hasXBoundaryGuards
  bool sourceHasXBoundaryGuards();

  /// Wrapper for GridDataSource::hasYBoundaryGuards
  bool sourceHasYBoundaryGuards();
  
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
  virtual void communicate(FieldPerp& f);

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
  [[deprecated("This experimental functionality will be removed in 5.0")]]
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
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual comm_handle receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag) = 0;
  
  virtual int getNXPE() = 0; ///< The number of processors in the X direction
  virtual int getNYPE() = 0; ///< The number of processors in the Y direction
  virtual int getXProcIndex() = 0; ///< This processor's index in X direction
  virtual int getYProcIndex() = 0; ///< This processor's index in Y direction
  
  // X communications
  virtual bool firstX() = 0;  ///< Is this processor first in X? i.e. is there a boundary to the left in X?
  virtual bool lastX() = 0; ///< Is this processor last in X? i.e. is there a boundary to the right in X?

  /// Domain is periodic in X?
  bool periodicX{false};

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
  MpiWrapper& getMpi() { return *mpi; } ///< Return pointer to the mesh's MPI Wrapper object

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  virtual bool periodicY(int jx) const;

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  /// \param[out] ts  The Twist-Shift angle if periodic
  virtual bool periodicY(int jx, BoutReal &ts) const = 0;

  /// Get number of boundaries in the y-direction, i.e. locations where there are boundary
  /// cells in the global grid
  virtual int numberOfYBoundaries() const = 0;

  /// Is there a branch cut at this processor's lower y-boundary?
  ///
  /// @param[in] jx             The local (on this processor) index in X
  /// @returns pair<bool, BoutReal> - bool is true if there is a branch cut,
  ///                                 BoutReal gives the total zShift for a 2pi
  ///                                 poloidal circuit if there is a branch cut
  virtual std::pair<bool, BoutReal> hasBranchCutLower(int jx) const = 0;

  /// Is there a branch cut at this processor's upper y-boundary?
  ///
  /// @param[in] jx             The local (on this processor) index in X
  /// @returns pair<bool, BoutReal> - bool is true if there is a branch cut,
  ///                                 BoutReal gives the total zShift for a 2pi
  ///                                 poloidal circuit if there is a branch cut
  virtual std::pair<bool, BoutReal> hasBranchCutUpper(int jx) const = 0;
  
  virtual int ySize(int jx) const; ///< The number of points in Y at fixed X index \p jx

  // Y communications
  virtual bool firstY() const = 0; ///< Is this processor first in Y? i.e. is there a boundary at lower Y?
  virtual bool lastY() const = 0; ///< Is this processor last in Y? i.e. is there a boundary at upper Y?
  virtual bool firstY(int xpos) const = 0; ///< Is this processor first in Y? i.e. is there a boundary at lower Y?
  virtual bool lastY(int xpos) const = 0; ///< Is this processor last in Y? i.e. is there a boundary at upper Y?
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual int UpXSplitIndex() = 0;  ///< If the upper Y guard cells are split in two, return the X index where the split occurs
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual int DownXSplitIndex() = 0; ///< If the lower Y guard cells are split in two, return the X index where the split occurs

  /// Send data
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual int sendYOutIndest(BoutReal *buffer, int size, int tag) = 0;

  ///
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual int sendYOutOutdest(BoutReal *buffer, int size, int tag) = 0;

  ///
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual int sendYInIndest(BoutReal *buffer, int size, int tag) = 0;

  ///
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual int sendYInOutdest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual comm_handle irecvYOutIndest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual comm_handle irecvYOutOutdest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  [[deprecated("This experimental functionality will be removed in 5.0")]]
  virtual comm_handle irecvYInIndest(BoutReal *buffer, int size, int tag) = 0;

  /// Non-blocking receive. Must be followed by a call to wait()
  ///
  /// @param[out] buffer  A buffer of length \p size which must already be allocated
  /// @param[in] size The number of BoutReals expected
  /// @param[in] tag  The tag number of the expected message
  [[deprecated("This experimental functionality will be removed in 5.0")]]
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
  virtual std::vector<BoundaryRegion*> getBoundaries() = 0;

  /// Add a boundary region to this processor
  virtual void addBoundary(BoundaryRegion* UNUSED(bndry)) {}

  /// Get all the parallel (Y) boundaries on this processor 
  virtual std::vector<BoundaryRegionPar*> getBoundariesPar() = 0;

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
  
  /// Returns the global X index given a local index
  /// If the local index includes the boundary cells, then so does the global.
  [[deprecated("Use getGlobalXIndex instead")]]
  int XGLOBAL(int xloc) const { return getGlobalXIndex(xloc); }
  /// Returns the global Y index given a local index
  /// The local index must include the boundary, the global index does not.
  [[deprecated("Use getGlobalYIndex or getGlobalYIndexNoBoundaries instead")]]
  virtual int YGLOBAL(int yloc) const { return getGlobalYIndexNoBoundaries(yloc); }

  /// Returns the local X index given a global index
  /// If the global index includes the boundary cells, then so does the local.
  virtual int XLOCAL(int xglo) const = 0;
  /// Returns the local Y index given a global index
  /// If the global index includes the boundary cells, then so does the local.
  virtual int YLOCAL(int yglo) const = 0;

  /// Returns a global X index given a local index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalXIndex(int xlocal) const = 0;

  /// Returns a global X index given a local index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalXIndexNoBoundaries(int xlocal) const = 0;

  /// Returns a global Y index given a local index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalYIndex(int ylocal) const = 0;

  /// Returns a global Y index given a local index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalYIndexNoBoundaries(int ylocal) const = 0;

  /// Returns a global Z index given a local index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalZIndex(int zlocal) const = 0;

  /// Returns a global Z index given a local index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalZIndexNoBoundaries(int zlocal) const = 0;

  /// Size of the mesh on this processor including guard/boundary cells
  int LocalNx, LocalNy, LocalNz;
  
  /// Local ranges of data (inclusive), excluding guard cells
  int xstart, xend, ystart, yend, zstart, zend;
  
  /// Enable staggered grids (Centre, Lower). Otherwise all vars are
  /// cell centred (default).
  bool StaggerGrids{false};
  
  /// Include integrated shear (if shifting X)
  bool IncIntShear{false};

  int numberOfXPoints{0};

  /// Coordinate system
  Coordinates *getCoordinates(const CELL_LOC location = CELL_CENTRE) {
    return getCoordinatesSmart(location).get();
  };

  std::shared_ptr<Coordinates>
  getCoordinatesSmart(const CELL_LOC location = CELL_CENTRE) {
    ASSERT1(location != CELL_DEFAULT);
    ASSERT1(location != CELL_VSHIFT);

    auto found = coords_map.find(location);
    if (found != coords_map.end()) {
      // True branch most common, returns immediately
      return found->second;
    }

    // No coordinate system set. Create default
    // Note that this can't be allocated here due to incomplete type
    // (circular dependency between Mesh and Coordinates)
    auto inserted = coords_map.emplace(location, nullptr);
    inserted.first->second = createDefaultCoordinates(location);
    return inserted.first->second;
  }

  /// Returns the non-CELL_CENTRE location
  /// allowed as a staggered location
  CELL_LOC getAllowedStaggerLoc(DIRECTION direction) const {
    AUTO_TRACE();
    switch (direction) {
    case (DIRECTION::X):
      return CELL_XLOW;
    case (DIRECTION::Y):
    case (DIRECTION::YOrthogonal):
    case (DIRECTION::YAligned):
      return CELL_YLOW;
    case (DIRECTION::Z):
      return CELL_ZLOW;
    default:
      throw BoutException("Unhandled direction encountered in getAllowedStaggerLoc");
    }
  };

  /// Returns the number of grid points in the
  /// particular direction
  int getNpoints(DIRECTION direction) const {
    AUTO_TRACE();
    switch (direction) {
    case (DIRECTION::X):
      return LocalNx;
    case (DIRECTION::Y):
    case (DIRECTION::YOrthogonal):
    case (DIRECTION::YAligned):
      return LocalNy;
    case (DIRECTION::Z):
      return LocalNz;
    default:
      throw BoutException("Unhandled direction encountered in getNpoints");
    }
  };

  /// Returns the number of guard points in the
  /// particular direction
  int getNguard(DIRECTION direction) const {
    AUTO_TRACE();
    switch (direction) {
    case (DIRECTION::X):
      return xstart;
    case (DIRECTION::Y):
    case (DIRECTION::YOrthogonal):
    case (DIRECTION::YAligned):
      return ystart;
    case (DIRECTION::Z):
      return 2;
    default:
      throw BoutException("Unhandled direction encountered in getNguard");
    }
  };

  /// Re-calculate staggered Coordinates, useful if CELL_CENTRE Coordinates are changed
  void recalculateStaggeredCoordinates();

  ///////////////////////////////////////////////////////////
  // INDEX DERIVATIVE OPERATORS
  ///////////////////////////////////////////////////////////

  ////// Utilties and parameters
  
  /// Fraction of modes to filter. This is set in derivs_init from option "ddz:fft_filter"
  BoutReal fft_derivs_filter{0.0};

  /// Determines the resultant output stagger location in derivatives
  /// given the input and output location. Also checks that the
  /// combination of locations is allowed
  STAGGER getStagger(const CELL_LOC inloc, const CELL_LOC outloc,
                     const CELL_LOC allowedloc) const;

  /// Determines the resultant output stagger location in derivatives
  /// given the input and output location. Also checks that the
  /// combination of locations is allowed. This overload also checks
  /// the location of a second input field (velocity) is consistent.
  STAGGER getStagger(const CELL_LOC vloc, const CELL_LOC inloc, const CELL_LOC outloc,
                     const CELL_LOC allowedloc) const;

  // All of these derivative routines should probably be moved out of mesh to become
  // free functions. As an intermediate step the member routines could just call the
  // free functions.

  ////// STANDARD OPERATORS

  ////////////// X DERIVATIVE /////////////////
  template <typename T>
  DEPRECATED(T indexDDX(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::DDX(f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexD2DX2(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                          const std::string& method = "DEFAULT",
                          REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::D2DX2(f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexD4DX4(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                          const std::string& method = "DEFAULT",
                          REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::D4DX4(f, outloc, method, region);
  }

  ////////////// Y DERIVATIVE /////////////////

  template <typename T>
  DEPRECATED(T indexDDY(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::DDY(f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexD2DY2(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                          const std::string& method = "DEFAULT",
                          REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::D2DY2(f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexD4DY4(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                          const std::string& method = "DEFAULT",
                          REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::D4DY4(f, outloc, method, region);
  }

  ////////////// Z DERIVATIVE /////////////////
  template <typename T>
  DEPRECATED(T indexDDZ(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                        const std::string& method = "DEFAULT",
                        REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::DDZ(f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexD2DZ2(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                          const std::string& method = "DEFAULT",
                          REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::D2DZ2(f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexD4DZ4(const T& f, CELL_LOC outloc = CELL_DEFAULT,
                          const std::string& method = "DEFAULT",
                          REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::D4DZ4(f, outloc, method, region);
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
  /// @param[in] outloc The cell location where the result is desired. The default is the
  /// same as \p f
  /// @param[in] method  The differencing method to use
  /// @param[in] region  The region of the grid for which the result is calculated.

  ////////////// X DERIVATIVE /////////////////

  template <typename T>
  DEPRECATED(T indexVDDX(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT",
                         REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::VDDX(vel, f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexFDDX(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT",
                         REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::FDDX(vel, f, outloc, method, region);
  }

  ////////////// Y DERIVATIVE /////////////////

  template <typename T>
  DEPRECATED(T indexVDDY(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT",
                         REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::VDDY(vel, f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexFDDY(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT",
                         REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::FDDY(vel, f, outloc, method, region);
  }

  ////////////// Z DERIVATIVE /////////////////

  template <typename T>
  DEPRECATED(T indexVDDZ(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT",
                         REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::VDDZ(vel, f, outloc, method, region);
  }

  template <typename T>
  DEPRECATED(T indexFDDZ(const T& vel, const T& f, CELL_LOC outloc = CELL_DEFAULT,
                         const std::string& method = "DEFAULT",
                         REGION region = RGN_NOBNDRY) const) {
    AUTO_TRACE();
    return bout::derivatives::index::FDDZ(vel, f, outloc, method, region);
  }

  [[deprecated("Please use free function toFieldAligned instead")]]
  const Field3D toFieldAligned(const Field3D &f, const REGION region = RGN_ALL) {
    return ::toFieldAligned(f, toString(region));
  }

  [[deprecated("Please use free function fromFieldAligned instead")]]
  const Field3D fromFieldAligned(const Field3D &f, const REGION region = RGN_ALL) {
    return ::fromFieldAligned(f, toString(region));
  }

  [[deprecated("Please use free function toFieldAligned instead")]]
  const Field2D toFieldAligned(const Field2D &f, const REGION region = RGN_ALL) {
    return ::toFieldAligned(f, toString(region));
  }

  [[deprecated("Please use free function fromFieldAligned instead")]]
  const Field2D fromFieldAligned(const Field2D &f, const REGION region = RGN_ALL) {
    return ::fromFieldAligned(f, toString(region));
  }

  [[deprecated("Please use "
      "Coordinates::getParallelTransform().canToFromFieldAligned instead")]]
  bool canToFromFieldAligned() {
    return getCoordinates()->getParallelTransform().canToFromFieldAligned();
  }

  [[deprecated("Please use Coordinates::setParallelTransform instead")]]
  void setParallelTransform(std::unique_ptr<ParallelTransform> pt) {
    getCoordinates()->setParallelTransform(std::move(pt));
  }

  [[deprecated("This call is now unnecessary")]]
  void setParallelTransform() {
    // The ParallelTransform is set from options in the Coordinates
    // constructor, so this method doesn't need to do anything
  }

  [[deprecated("Please use Coordinates::getParallelTransform instead")]]
  ParallelTransform& getParallelTransform() {
    return getCoordinates()->getParallelTransform();
  }


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

  /// Indicate if named region has already been defined
  bool hasRegion3D(const std::string& region_name) const;
  bool hasRegion2D(const std::string& region_name) const;
  bool hasRegionPerp(const std::string& region_name) const;

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

  /// Source for grid data
  GridDataSource* source{nullptr};

  /// Coordinate systems at different CELL_LOCs
  std::map<CELL_LOC, std::shared_ptr<Coordinates>> coords_map;

  /// Mesh options section
  Options *options{nullptr};

  /// Set whether to call calcParallelSlices on all communicated fields (true) or not (false)
  bool calcParallelSlices_on_communicate{true};

  /// Read a 1D array of integers
  const std::vector<int> readInts(const std::string &name, int n);
  
  /// Calculates the size of a message for a given x and y range
  int msg_len(const std::vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt);
  
  /// Initialise derivatives
  void derivs_init(Options* options);

  /// Pointer to the global MPI wrapper, for convenience
  MpiWrapper* mpi = nullptr;

private:

  /// Allocates default Coordinates objects
  /// By default attempts to read staggered Coordinates from grid data source,
  /// interpolating from CELL_CENTRE if not present. Set
  /// force_interpolate_from_centre argument to true to always interpolate
  /// (useful if CELL_CENTRE Coordinates have been changed, so reading from file
  /// would not be correct).
  std::shared_ptr<Coordinates> createDefaultCoordinates(const CELL_LOC location,
      bool force_interpolate_from_centre=false);

  //Internal region related information
  std::map<std::string, Region<Ind3D>> regionMap3D;
  std::map<std::string, Region<Ind2D>> regionMap2D;
  std::map<std::string, Region<IndPerp>> regionMapPerp;
  Array<int> indexLookup3Dto2D;

  int localNumCells3D = -1, localNumCells2D = -1, localNumCellsPerp = -1;
};

#endif // __MESH_H__
