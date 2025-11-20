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
 * Copyright 2010-2025 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_MESH_H
#define BOUT_MESH_H

#include "bout/bout_enum_class.hxx"
#include "bout/bout_types.hxx"
#include "bout/coordinates.hxx" // Coordinates class
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/field_data.hxx"
#include "bout/fieldgroup.hxx"
#include "bout/generic_factory.hxx"
#include "bout/index_derivs_interface.hxx"
#include "bout/mpi_wrapper.hxx"
#include "bout/options.hxx"
#include "bout/region.hxx"
#include "bout/sys/range.hxx" // RangeIterator
#include "bout/unused.hxx"

#include "mpi.h"

#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>

class BoundaryRegion;
class BoundaryRegionPar;
class GridDataSource;

class MeshFactory : public Factory<Mesh, MeshFactory, GridDataSource*, Options*> {
public:
  static constexpr auto type_name = "Mesh";
  static constexpr auto section_name = "mesh";
  static constexpr auto option_name = "type";
  static constexpr auto default_type = "bout";

  ReturnType create(const std::string& type, Options* options = nullptr,
                    GridDataSource* source = nullptr) const;
  ReturnType create(Options* options = nullptr, GridDataSource* source = nullptr) const;
};

BOUT_ENUM_CLASS(BoundaryParType, all, xin, xout, fwd, bwd, xin_fwd, xout_fwd, xin_bwd,
                xout_bwd, SIZE);

template <class DerivedType>
using RegisterMesh = MeshFactory::RegisterInFactory<DerivedType>;

/// Type used to return pointers to handles
using comm_handle = void*;

class Mesh {
public:
  /// Constructor for a "bare", uninitialised Mesh
  /// Only useful for testing
  Mesh() : source(nullptr), options(nullptr), include_corner_cells(true) {}

  /// Constructor
  /// @param[in] s  The source to be used for loading variables
  /// @param[in] options  The options section for settings
  Mesh(GridDataSource* s, Options* options);

  /// Destructor
  virtual ~Mesh();

  /// Create a Mesh object
  ///
  /// @param[in] source  The data source to use for loading variables
  /// @param[in] opt     The option section. By default this is "mesh"
  static Mesh* create(GridDataSource* source, Options* opt = nullptr);

  /// Create a Mesh object
  ///
  /// The source is determined by
  ///  1) If "file" is set in the options, read that
  ///  2) If "grid" is set in global options, read that
  ///  3) Use options as data source
  ///
  /// @param[in] opt  Input options. Default is "mesh" section
  static Mesh* create(Options* opt = nullptr);

  /// Loads the mesh values
  ///
  /// Currently need to create and load mesh in separate calls
  /// because creating Fields uses the global "mesh" pointer
  /// which isn't created until Mesh is constructed
  virtual int load() { return 1; }

  /// Add output variables to \p output_options
  /// These are used for post-processing
  virtual void outputVars([[maybe_unused]] Options& output_options) {}

  // Get routines to request data from mesh file

  /// Get a string from the input source
  ///
  /// @param[out] sval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(std::string& sval, const std::string& name, const std::string& def = "");

  /// Get an integer from the input source
  ///
  /// @param[out] ival  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(int& ival, const std::string& name, int def = 0);

  /// Get a BoutReal from the input source
  ///
  /// @param[out] rval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(BoutReal& rval, const std::string& name, BoutReal def = 0.0);

  /// Get a bool from the input source
  ///
  /// @param[out] bval  The value will be put into this variable
  /// @param[in] name   The name of the variable to read
  /// @param[in] def    The default value if not found
  ///
  /// @returns zero if successful, non-zero on failure
  int get(bool& bval, const std::string& name, bool def = false);

  /// Get a Field2D from the input source
  /// including communicating guard cells
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero if successful, non-zero on failure
  int get(Field2D& var, const std::string& name, BoutReal def = 0.0,
          bool communicate = true, CELL_LOC location = CELL_DEFAULT);

  /// Get a Field3D from the input source
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero if successful, non-zero on failure
  int get(Field3D& var, const std::string& name, BoutReal def = 0.0,
          bool communicate = true, CELL_LOC location = CELL_DEFAULT);

  /// Get a FieldPerp from the input source
  ///
  /// @param[out] var   This will be set to the value. Will be allocated if needed
  /// @param[in] name   Name of the variable to read
  /// @param[in] def    The default value if not found
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero if successful, non-zero on failure
  int get(FieldPerp& var, const std::string& name, BoutReal def = 0.0,
          bool communicate = true, CELL_LOC location = CELL_DEFAULT);

  /// Get a Vector2D from the input source.
  /// If \p var is covariant then this gets three
  /// Field2D variables with "_x", "_y", "_z" appended to \p name
  /// If \p var is contravariant, then "x", "y", "z" are appended to \p name
  ///
  /// By default all fields revert to zero
  ///
  /// @param[in] var  This will be set to the value read
  /// @param[in] name  The name of the vector. Individual fields are read based on this
  /// name by appending. See above
  /// @param[in] def   The default value if not found (used for all the components)
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero always.
  int get(Vector2D& var, const std::string& name, BoutReal def = 0.0,
          bool communicate = true);

  /// Get a Vector3D from the input source.
  /// If \p var is covariant then this gets three
  /// Field3D variables with "_x", "_y", "_z" appended to \p name
  /// If \p var is contravariant, then "x", "y", "z" are appended to \p name
  ///
  /// By default all fields revert to zero
  ///
  /// @param[in] var  This will be set to the value read
  /// @param[in] name  The name of the vector. Individual fields are read based on this
  /// name by appending. See above
  /// @param[in] def    The default value if not found (used for all the components)
  /// @param[in] communicate  Should the field be communicated to fill guard cells?
  ///
  /// @returns zero always.
  int get(Vector3D& var, const std::string& name, BoutReal def = 0.0,
          bool communicate = true);

  /// Test if input source was a grid file
  bool isDataSourceGridFile() const;

  /// Wrapper for GridDataSource::hasVar
  bool sourceHasVar(const std::string& name);

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

  template <typename... Ts>
  void communicateYZ(Ts&... ts) {
    FieldGroup g(ts...);
    communicateYZ(g);
  }

  /*!
   * Communicate a group of fields
   */
  void communicate(FieldGroup& g);

  /// Communcate guard cells in XZ only
  /// i.e. no Y communication
  ///
  /// @param g  The group of fields to communicate. Guard cells will be modified
  void communicateXZ(FieldGroup& g);

  /// Communcate guard cells in YZ only
  /// i.e. no X communication
  ///
  /// @param g  The group of fields to communicate. Guard cells will be modified
  void communicateYZ(FieldGroup& g);

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

  /// Send guard cells from a list of FieldData objects in the x-direction
  /// Packs arguments into a FieldGroup and passes to send(FieldGroup&).
  /// Perform communications without waiting for them to finish.
  /// Requires a call to wait() afterwards.
  template <typename... Ts>
  comm_handle sendX(Ts&... ts) {
    FieldGroup g(ts...);
    return sendX(g);
  }

  /// Send guard cells from a list of FieldData objects in the y-direction
  /// Packs arguments into a FieldGroup and passes to send(FieldGroup&).
  template <typename... Ts>
  comm_handle sendY(Ts&... ts) {
    FieldGroup g(ts...);
    return sendY(g);
  }

  /// Perform communications without waiting for them
  /// to finish. Requires a call to wait() afterwards.
  ///
  /// \param g Group of fields to communicate
  /// \returns handle to be used as input to wait()
  virtual comm_handle send(FieldGroup& g) = 0;

  /// Send only the x-guard cells
  virtual comm_handle sendX(FieldGroup& g, comm_handle handle = nullptr,
                            bool disable_corners = false) = 0;

  /// Send only the y-guard cells
  virtual comm_handle sendY(FieldGroup& g, comm_handle handle = nullptr) = 0;

  /// Wait for the handle, return error code
  virtual int wait(comm_handle handle) = 0; ///< Wait for the handle, return error code

  // non-local communications

  virtual int getNXPE() = 0;       ///< The number of processors in the X direction
  virtual int getNYPE() = 0;       ///< The number of processors in the Y direction
  virtual int getXProcIndex() = 0; ///< This processor's index in X direction
  virtual int getYProcIndex() = 0; ///< This processor's index in Y direction

  // X communications
  virtual bool firstX()
      const = 0; ///< Is this processor first in X? i.e. is there a boundary to the left in X?
  virtual bool lastX()
      const = 0; ///< Is this processor last in X? i.e. is there a boundary to the right in X?

  /// Domain is periodic in X?
  bool periodicX{false};

  int NXPE, PE_XIND; ///< Number of processors in X, and X processor index

  /// Send a buffer of data to processor at X index +1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  virtual int sendXOut(BoutReal* buffer, int size, int tag) = 0;

  /// Send a buffer of data to processor at X index -1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  virtual int sendXIn(BoutReal* buffer, int size, int tag) = 0;

  /// Receive a buffer of data from X index +1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  virtual comm_handle irecvXOut(BoutReal* buffer, int size, int tag) = 0;

  /// Receive a buffer of data from X index -1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  virtual comm_handle irecvXIn(BoutReal* buffer, int size, int tag) = 0;

  MPI_Comm getXcomm() {
    return getXcomm(0);
  } ///< Return communicator containing all processors in X
  virtual MPI_Comm getXcomm(int jy) const = 0; ///< Return X communicator
  virtual MPI_Comm getYcomm(int jx) const = 0; ///< Return Y communicator

  /// Return pointer to the mesh's MPI Wrapper object
  MpiWrapper& getMpi() { return *mpi; }

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  virtual bool periodicY(int jx) const;

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  /// \param[out] ts  The Twist-Shift angle if periodic
  virtual bool periodicY(int jx, BoutReal& ts) const = 0;

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

  /// Is this processor first in Y?
  /// Note: First on the global grid, not necessarily at a boundary
  [[deprecated("Please use firstY(xind) instead")]] virtual bool firstY() const = 0;

  /// Is this processor last in Y?
  /// Note: Last on the global grid, not necessarily at a boundary
  [[deprecated("Please use lastY(xind) instead")]] virtual bool lastY() const = 0;

  /// Is this processor first in Y?
  /// Note: Not necessarily at a boundary, but first in the Y communicator
  ///       for the flux surface through local X index xpos
  virtual bool firstY(int xpos) const = 0;

  /// Is this processor last in Y?
  /// Note: Not necessarily at a boundary, but last in the Y communicator
  ///       for the flux surface through local X index xpos
  virtual bool lastY(int xpos) const = 0;

  // Boundary region iteration

  /// Iterate over the lower Y boundary
  virtual RangeIterator iterateBndryLowerY() const = 0;

  /// Iterate over the upper Y boundary
  virtual RangeIterator iterateBndryUpperY() const = 0;
  virtual RangeIterator iterateBndryLowerOuterY() const = 0;
  virtual RangeIterator iterateBndryLowerInnerY() const = 0;
  virtual RangeIterator iterateBndryUpperOuterY() const = 0;
  virtual RangeIterator iterateBndryUpperInnerY() const = 0;

  /// Is there a boundary on the lower guard cells in Y
  /// on any processor along the X direction?
  bool hasBndryLowerY();

  /// Is there a boundary on the upper guard cells in Y
  /// on any processor along the X direction?
  bool hasBndryUpperY();
  // Boundary regions

  /// Return a vector containing all the boundary regions on this processor
  virtual std::vector<BoundaryRegion*> getBoundaries() = 0;

  /// Get the set of all possible boundaries in this configuration
  virtual std::set<std::string> getPossibleBoundaries() const { return {}; }

  /// Add a boundary region to this processor
  virtual void addBoundary(BoundaryRegion* UNUSED(bndry)) {}

  /// Get the list of parallel boundary regions. The option specifies with
  /// region to get. Default is to get all regions. All possible options are
  /// listed at the top of this file, see BoundaryParType.
  /// For example:
  /// get all regions:
  /// mesh->getBoundariesPar(Mesh::BoundaryParType::all)
  /// get only xout:
  /// mesh->getBoundariesPar(Mesh::BoundaryParType::xout)
  virtual std::vector<std::shared_ptr<BoundaryRegionPar>>
  getBoundariesPar(BoundaryParType type = BoundaryParType::all) = 0;

  /// Add a parallel(Y) boundary to this processor
  virtual void addBoundaryPar(std::shared_ptr<BoundaryRegionPar> UNUSED(bndry),
                              BoundaryParType UNUSED(type)) {}

  /// Branch-cut special handling (experimental)
  virtual Field3D smoothSeparatrix(const Field3D& f) { return f; }

  virtual BoutReal GlobalX(int jx) const = 0;      ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(int jy) const = 0;      ///< Continuous Y index (0 -> 1)
  virtual BoutReal GlobalX(BoutReal jx) const = 0; ///< Continuous X index between 0 and 1
  virtual BoutReal GlobalY(BoutReal jy) const = 0; ///< Continuous Y index (0 -> 1)

  //////////////////////////////////////////////////////////

  int GlobalNx, GlobalNy, GlobalNz; ///< Size of the global arrays. Note: can have holes
  /// Size of the global arrays excluding boundary points.
  int GlobalNxNoBoundaries, GlobalNyNoBoundaries, GlobalNzNoBoundaries;

  /// Note: These offsets only correct if Y guards are not included in the global array
  ///       and are corrected in gridfromfile.cxx
  int OffsetX, OffsetY, OffsetZ; ///< Offset of this mesh within the global array
                                 ///< so startx on this processor is OffsetX in global

  /// Map between local and global indices
  /// (MapGlobalX, MapGlobalY, MapGlobalZ) in the global index space maps to (MapLocalX, MapLocalY, MapLocalZ) locally.
  /// Note that boundary cells are included in the global index space, but communication
  /// guard cells are not.
  int MapGlobalX, MapGlobalY, MapGlobalZ; ///< Start global indices
  int MapLocalX, MapLocalY, MapLocalZ;    ///< Start local indices
  int MapCountX, MapCountY, MapCountZ;    ///< Size of the mapped region

  /// Returns the number of unique cells (i.e., ones not used for
  /// communication) on this processor for 3D fields. Boundaries
  /// are only included to a depth of 1.
  virtual int localSize3D();
  /// Returns the number of unique cells (i.e., ones not used for
  /// communication) on this processor for 2D fields. Boundaries
  /// are only included to a depth of 1.
  virtual int localSize2D();
  /// Returns the number of unique cells (i.e., ones not used for
  /// communication) on this processor for perpendicular fields.
  /// Boundaries are only included to a depth of 1.
  virtual int localSizePerp();

  /// Get the value of the first global 3D index on this processor.
  virtual int globalStartIndex3D();
  /// Get the value of the first global 2D index on this processor.
  virtual int globalStartIndex2D();
  /// Get the value of the first global perpendicular index on this processor.
  virtual int globalStartIndexPerp();

  /// Returns a global X index given a local index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalXIndex(int xlocal) const = 0;

  /// Returns a global X index given a local index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalXIndexNoBoundaries(int xlocal) const = 0;

  /// Returns a local X index given a global index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getLocalXIndex(int xglobal) const = 0;

  /// Returns a local X index given a global index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getLocalXIndexNoBoundaries(int xglobal) const = 0;

  /// Returns a global Y index given a local index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalYIndex(int ylocal) const = 0;

  /// Returns a global Y index given a local index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalYIndexNoBoundaries(int ylocal) const = 0;

  /// Returns a local Y index given a global index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getLocalYIndex(int yglobal) const = 0;

  /// Returns a local Y index given a global index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getLocalYIndexNoBoundaries(int yglobal) const = 0;

  /// Returns a global Z index given a local index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalZIndex(int zlocal) const = 0;

  /// Returns a global Z index given a local index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getGlobalZIndexNoBoundaries(int zlocal) const = 0;

  /// Returns a local Z index given a global index.
  /// Global index includes boundary cells, local index includes boundary or guard cells.
  virtual int getLocalZIndex(int zglobal) const = 0;

  /// Returns a local Z index given a global index.
  /// Global index excludes boundary cells, local index includes boundary or guard cells.
  virtual int getLocalZIndexNoBoundaries(int zglobal) const = 0;

  /// Size of the mesh on this processor including guard/boundary cells
  int LocalNx, LocalNy, LocalNz;

  /// Local ranges of data (inclusive), excluding guard cells
  int xstart, xend, ystart, yend, zstart, zend;

  /// Include integrated shear (if shifting X)
  bool IncIntShear{false};

  int numberOfXPoints{0};

  /// Coordinate system
  Coordinates* getCoordinates(const CELL_LOC location = CELL_CENTRE) {
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
    inserted.first->second->geometry(false);
    return inserted.first->second;
  }

  std::shared_ptr<Coordinates>
  getCoordinatesConst(const CELL_LOC location = CELL_CENTRE) const {
    ASSERT1(location != CELL_DEFAULT);
    ASSERT1(location != CELL_VSHIFT);

    auto found = coords_map.find(location);
    if (found != coords_map.end()) {
      // True branch most common, returns immediately
      return found->second;
    }
    throw BoutException("Coordinates not yet set. Use non-const version!");
  }

  /// Returns the non-CELL_CENTRE location
  /// allowed as a staggered location
  static CELL_LOC getAllowedStaggerLoc(DIRECTION direction) {
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

  ///////////////////////////////////////////////////////////
  // REGION RELATED ROUTINES
  ///////////////////////////////////////////////////////////

  /// Get the named region from the region_map for the data iterator
  ///
  /// Throws if region_name not found
  template <class T>
  const Region<typename T::ind_type>& getRegion(const std::string& region_name) const;

  const Region<>& getRegion(const std::string& region_name) const {
    return getRegion3D(region_name);
  }
  const Region<Ind3D>& getRegion3D(const std::string& region_name) const;
  const Region<Ind2D>& getRegion2D(const std::string& region_name) const;
  const Region<IndPerp>& getRegionPerp(const std::string& region_name) const;

  /// Indicate if named region has already been defined
  bool hasRegion3D(const std::string& region_name) const;
  bool hasRegion2D(const std::string& region_name) const;
  bool hasRegionPerp(const std::string& region_name) const;

  /// Add a new region to the region_map for the data iterator
  ///
  /// Outputs an error message if region_name already exists
  void addRegion(const std::string& region_name, const Region<>& region) {
    return addRegion3D(region_name, region);
  }
  void addRegion(const std::string& region_name, const Region<Ind2D>& region) {
    return addRegion2D(region_name, region);
  }
  void addRegion(const std::string& region_name, const Region<IndPerp>& region) {
    return addRegionPerp(region_name, region);
  }
  void addRegion3D(const std::string& region_name, const Region<Ind3D>& region);
  void addRegion2D(const std::string& region_name, const Region<Ind2D>& region);
  void addRegionPerp(const std::string& region_name, const Region<IndPerp>& region);

  /// Converts an Ind2D to an Ind3D using calculation
  Ind3D ind2Dto3D(const Ind2D& ind2D, int jz = 0) {
    return {ind2D.ind * LocalNz + jz, LocalNy, LocalNz};
  }

  /// Converts an Ind3D to an Ind2D using calculation
  Ind2D ind3Dto2D(const Ind3D& ind3D) { return {ind3D.ind / LocalNz, LocalNy, 1}; }

  /// Converts an Ind3D to an IndPerp using calculation
  IndPerp ind3DtoPerp(const Ind3D& ind3D) {
    return {ind3D.x() * LocalNz + ind3D.z(), 1, LocalNz};
  }

  /// Converts an IndPerp to an Ind3D using calculation
  Ind3D indPerpto3D(const IndPerp& indPerp, int jy = 0) {
    int jz = indPerp.z();
    return {(indPerp.ind - jz) * LocalNy + LocalNz * jy + jz, LocalNy, LocalNz};
  }

  /// Converts an Ind3D to an Ind2D representing a 2D index using a lookup -- to be used with care
  Ind2D map3Dto2D(const Ind3D& ind3D) {
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
  Options* options{nullptr};

  /// Set whether to call calcParallelSlices on all communicated fields (true) or not (false)
  bool calcParallelSlices_on_communicate{true};

  /// Read a 1D array of integers
  const std::vector<int> readInts(const std::string& name, int n);

  /// Calculates the size of a message for a given x and y range
  int msg_len(const std::vector<FieldData*>& var_list, int xge, int xlt, int yge,
              int ylt);

  /// Initialise derivatives
  void derivs_init(Options* options);

  /// Pointer to the global MPI wrapper, for convenience
  MpiWrapper* mpi = nullptr;

public:
  // The maxregionblocksize to use when creating the default regions.
  // Can be set in the input file and the global default is set by,
  // MAXREGIONBLOCKSIZE in include/bout/region.hxx
  int maxregionblocksize{MAXREGIONBLOCKSIZE};

  /// Enable staggered grids (Centre, Lower). Otherwise all vars are
  /// cell centred (default).
  bool StaggerGrids{false};

  // Switch for communication of corner guard and boundary cells
  const bool include_corner_cells;

  std::optional<size_t> getCommonRegion(std::optional<size_t>, std::optional<size_t>);
  size_t getRegionID(const std::string& region) const;
  const Region<Ind3D>& getRegion(size_t RegionID) const { return region3D[RegionID]; }
  const Region<Ind3D>& getRegion(std::optional<size_t> RegionID) const {
    ASSERT1(RegionID.has_value());
    return region3D[RegionID.value()];
  }
  bool isFci() const {
    const auto coords = this->getCoordinatesConst();
    if (coords == nullptr) {
      return false;
    }
    if (not coords->hasParallelTransform()) {
      return false;
    }
    return not coords->getParallelTransform().canToFromFieldAligned();
  }

private:
  /// Allocates default Coordinates objects
  /// By default attempts to read staggered Coordinates from grid data source,
  /// interpolating from CELL_CENTRE if not present. Set
  /// force_interpolate_from_centre argument to true to always interpolate
  /// (useful if CELL_CENTRE Coordinates have been changed, so reading from file
  /// would not be correct).
  std::shared_ptr<Coordinates>
  createDefaultCoordinates(CELL_LOC location, bool force_interpolate_from_centre = false);

  //Internal region related information
  std::map<std::string, size_t> regionMap3D;
  std::vector<Region<Ind3D>> region3D;
  std::vector<std::optional<size_t>> region3Dintersect;
  std::map<std::string, Region<Ind2D>> regionMap2D;
  std::map<std::string, Region<IndPerp>> regionMapPerp;
  Array<int> indexLookup3Dto2D;

  int localNumCells3D = -1, localNumCells2D = -1, localNumCellsPerp = -1;
};

template <>
inline const Region<Ind3D>&
Mesh::getRegion<Field3D>(const std::string& region_name) const {
  return getRegion3D(region_name);
}
template <>
inline const Region<Ind2D>&
Mesh::getRegion<Field2D>(const std::string& region_name) const {
  return getRegion2D(region_name);
}
template <>
inline const Region<IndPerp>&
Mesh::getRegion<FieldPerp>(const std::string& region_name) const {
  return getRegionPerp(region_name);
}

#endif // BOUT_MESH_H
