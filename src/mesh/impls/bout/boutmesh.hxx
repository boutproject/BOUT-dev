
#ifndef __BOUTMESH_H__
#define __BOUTMESH_H__

#include "mpi.h"

#include <bout/mesh.hxx>
#include "unused.hxx"

#include <list>
#include <vector>
#include <cmath>

/// Implementation of Mesh (mostly) compatible with BOUT
///
/// Topology and communications compatible with BOUT
/// conventions.
class BoutMesh : public Mesh {
 public:
  BoutMesh(GridDataSource *s, Options *options = nullptr);
  ~BoutMesh() override;

  /// Read in the mesh from data sources
  int load() override;

  /////////////////////////////////////////////
  // Communicate variables

  /// Send data between processors
  /// Does not wait for communications
  /// to complete, so wait() must be called
  /// before guard cell values are used
  ///
  /// @param[in] g  A group of fields to communicate
  ///
  /// Example
  /// -------
  ///
  /// comm_handle handle = mesh->send(group);
  /// ...
  /// mesh->wait(handle);
  ///
  comm_handle send(FieldGroup& g) override;

  /// Send only in the x-direction
  comm_handle sendX(FieldGroup& g, comm_handle handle = nullptr,
                    bool disable_corners = false) override;

  /// Send only in the y-direction
  comm_handle sendY(FieldGroup& g, comm_handle handle = nullptr) override;

  /// Wait for a send operation to complete
  /// @param[in] handle  The handle returned by send()
  int wait(comm_handle handle) override;

  /////////////////////////////////////////////
  // non-local communications

  MPI_Request sendToProc(int xproc, int yproc, BoutReal* buffer, int size,
                         int tag) override;
  comm_handle receiveFromProc(int xproc, int yproc, BoutReal* buffer, int size,
                              int tag) override;

  int getNXPE() override;       ///< The number of processors in the X direction
  int getNYPE() override;       ///< The number of processors in the Y direction
  int getXProcIndex() override; ///< This processor's index in X direction
  int getYProcIndex() override; ///< This processor's index in Y direction

  /////////////////////////////////////////////
  // X communications

  bool firstX() const override; ///< Is this processor the first in X? i.e. is there a boundary
                          ///< to the left in X?
  bool lastX() const override; ///< Is this processor last in X? i.e. is there a boundary to the
                         ///< right in X?

  /// Send a buffer of data to processor at X index +1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  int sendXOut(BoutReal* buffer, int size, int tag) override;

  /// Send a buffer of data to processor at X index -1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  int sendXIn(BoutReal* buffer, int size, int tag) override;

  /// Receive a buffer of data from X index +1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  comm_handle irecvXOut(BoutReal* buffer, int size, int tag) override;

  /// Receive a buffer of data from X index -1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  comm_handle irecvXIn(BoutReal* buffer, int size, int tag) override;

  /// Return communicator containing all processors in X
  MPI_Comm getXcomm(int UNUSED(jy)) const override { return comm_x; }
  /// Return communicator containing all processors in Y
  MPI_Comm getYcomm(int jx) const override;

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  /// \param[out] ts  The Twist-Shift angle if periodic
  bool periodicY(int jx, BoutReal& ts) const override;

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  bool periodicY(int jx) const override;

  /// Get number of boundaries in the y-direction, i.e. locations where there are boundary
  /// cells in the global grid
  int numberOfYBoundaries() const override;

  /// Is there a branch cut at this processor's lower boundary?
  ///
  /// @param[in] jx             The local (on this processor) index in X
  /// @returns pair<bool, BoutReal> - bool is true if there is a branch cut,
  ///                                 BoutReal gives the total zShift for a 2pi
  ///                                 poloidal circuit if there is a branch cut
  std::pair<bool, BoutReal> hasBranchCutLower(int jx) const override;

  /// Is there a branch cut at this processor's upper boundary?
  ///
  /// @param[in] jx             The local (on this processor) index in X
  /// @returns pair<bool, BoutReal> - bool is true if there is a branch cut,
  ///                                 BoutReal gives the total zShift for a 2pi
  ///                                 poloidal circuit if there is a branch cut
  std::pair<bool, BoutReal> hasBranchCutUpper(int jx) const override;

  int ySize(int jx) const override; ///< The number of points in Y at fixed X index \p jx

  /////////////////////////////////////////////
  // Y communications

  bool firstY() const override;
  bool lastY() const override;
  bool firstY(int xpos) const override;
  bool lastY(int xpos) const override;
  int UpXSplitIndex() override;
  int DownXSplitIndex() override;
  int sendYOutIndest(BoutReal* buffer, int size, int tag) override;
  int sendYOutOutdest(BoutReal* buffer, int size, int tag) override;
  int sendYInIndest(BoutReal* buffer, int size, int tag) override;
  int sendYInOutdest(BoutReal* buffer, int size, int tag) override;
  comm_handle irecvYOutIndest(BoutReal* buffer, int size, int tag) override;
  comm_handle irecvYOutOutdest(BoutReal* buffer, int size, int tag) override;
  comm_handle irecvYInIndest(BoutReal* buffer, int size, int tag) override;
  comm_handle irecvYInOutdest(BoutReal* buffer, int size, int tag) override;

  // Boundary iteration
  const RangeIterator iterateBndryLowerY() const override;
  const RangeIterator iterateBndryUpperY() const override;
  const RangeIterator iterateBndryLowerInnerY() const override;
  const RangeIterator iterateBndryLowerOuterY() const override;
  const RangeIterator iterateBndryUpperInnerY() const override;
  const RangeIterator iterateBndryUpperOuterY() const override;

  // Boundary regions
  std::vector<BoundaryRegion*> getBoundaries() override;
  std::vector<BoundaryRegionPar*> getBoundariesPar() override;
  void addBoundaryPar(BoundaryRegionPar* bndry) override;

  const Field3D smoothSeparatrix(const Field3D& f) override;

  int getNx() const { return nx; }
  int getNy() const { return ny; }

  BoutReal GlobalX(int jx) const override;
  BoutReal GlobalY(int jy) const override;
  BoutReal GlobalX(BoutReal jx) const override;
  BoutReal GlobalY(BoutReal jy) const override;

  BoutReal getIxseps1() const { return ixseps1; }
  BoutReal getIxseps2() const { return ixseps2; }

  void outputVars(Datafile& file) override;

  int getGlobalXIndex(int xlocal) const override;
  int getGlobalXIndexNoBoundaries(int xlocal) const override;
  int getLocalXIndex(int xglobal) const override;
  int getLocalXIndexNoBoundaries(int xglobal) const override;
  int getGlobalYIndex(int ylocal) const override;
  int getGlobalYIndexNoBoundaries(int ylocal) const override;
  int getLocalYIndex(int yglobal) const override;
  int getLocalYIndexNoBoundaries(int yglobal) const override;
  int getGlobalZIndex(int zlocal) const override;
  int getGlobalZIndexNoBoundaries(int zlocal) const override;
  int getLocalZIndex(int zglobal) const override;
  int getLocalZIndexNoBoundaries(int zglobal) const override;

protected:
  /// A constructor used when making fake meshes for testing. This
  /// will make a mesh which thinks it corresponds to the subdomain on
  /// one processor, even though it's actually being run in serial.
  ///
  /// Pass \p create_topology = false to not set up topology, regions etc.
  BoutMesh(int input_nx, int input_ny, int input_nz, int mxg, int myg, int nxpe, int nype,
           int pe_xind, int pe_yind, bool create_topology = true);

  /// Very basic initialisation, only suitable for testing
  BoutMesh(int input_nx, int input_ny, int input_nz, int mxg, int myg, int input_npes)
      : nx(input_nx), ny(input_ny), nz(input_nz), NPES(input_npes), ixseps1(nx),
        ixseps2(nx), jyseps1_1(-1), jyseps2_1(ny / 2), jyseps1_2(jyseps2_1),
        jyseps2_2(ny - 1), ny_inner(jyseps2_1), MXG(mxg), MYG(myg), MZG(0) {}

  /// For debugging purposes (when creating fake parallel meshes), make
  /// the send and receive buffers share memory. This allows for
  /// communications to be faked between meshes as though they were on
  /// different processors.
  void overlapHandleMemory(BoutMesh* yup, BoutMesh* ydown, BoutMesh* xin,
			   BoutMesh* xout);

  /// Set the various y-decomposition indices, enforcing the following invariants:
  ///
  /// - `jyseps1_1` >= -1
  /// - `jyseps2_1` >= `jyseps1_1`
  /// - `jyseps1_2` >= `jyseps2_1`
  /// - `jyseps2_2` >= `jyseps1_2`
  /// - `jyseps2_2` < `ny`
  ///
  /// Inputs inconsistent with these invariants will be set to the
  /// minimum acceptable value.
  ///
  /// Also sets `numberOfXPoints` consistently (0, 1, or 2).
  void setYDecompositionIndices(int jyseps1_1_, int jyseps2_1_, int jyseps1_2_,
                                int jyseps2_2_, int ny_inner_);

  /// Structure for `setYDecompositionIndices` input/output values
  struct YDecompositionIndices {
    int jyseps1_1;
    int jyseps2_1;
    int jyseps1_2;
    int jyseps2_2;
    int ny_inner;
  };

  /// Version of `setYDecompositionindices` that returns the values
  /// used, useful for testing
  YDecompositionIndices setYDecompositionIndices(const YDecompositionIndices& indices);

  /// Choose NXPE (or NYPE) based on user input
  void chooseProcessorSplit(Options& options);
  /// Find a value for NXPE
  void findProcessorSplit();

  struct XDecompositionIndices {
    int ixseps1;
    int ixseps2;
  };

  /// Set the two x-decomposition indices. No invariants are enforced
  void setXDecompositionIndices(const XDecompositionIndices& indices);

  /// Set the derived grid sizes like `LocalN*`, `GlobalN*`, `Offset*`, etc
  ///
  /// Requires the following to be set first:
  ///
  /// - nx, ny, nz
  /// - MXG, MYG, MZG
  /// - NXPE, NYPE, NZPE
  /// - PE_XIND, PE_YIND
  /// - jyseps1_2, jyseps2_1
  void setDerivedGridSizes();

  /// Create the boundary regions in X
  void createXBoundaries();
  /// Create the boundary regions in Y
  void createYBoundaries();
private:
  std::string gridname;
  int nx, ny, nz; ///< Size of the grid in the input file
  int MX, MY, MZ; ///< size of the grid excluding boundary regions

  int MYSUB, MXSUB, MZSUB; ///< Size of the grid on this processor

  int NPES; ///< Number of processors
  int MYPE; ///< Rank of this processor

  int PE_YIND; ///< Y index of this processor
  int NYPE;    // Number of processors in the Y direction

  int NZPE;

  int MYPE_IN_CORE; // 1 if processor in core

  using Mesh::YGLOBAL;
  int XGLOBAL(BoutReal xloc, BoutReal& xglo) const;
  int YGLOBAL(BoutReal yloc, BoutReal& yglo) const;

  // Topology
  int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;
  int ixseps_inner, ixseps_outer, ixseps_upper, ixseps_lower;
  int ny_inner;

  std::vector<BoutReal> ShiftAngle; ///< Angle for twist-shift location

protected:
  // These are protected so we can make them public in the test suite
  // for testing

  // Processor number, local <-> global translation
  /// Returns the processor number, given X (\p xind) and Y (\p yind)
  /// processor indices. Returns -1 if out of range (no processor)
  int PROC_NUM(int xind, int yind) const;
  int YGLOBAL(int yloc, int yproc) const;
  int YLOCAL(int yglo, int yproc) const;
  /// Return the Y processor number given a global Y index
  int YPROC(int yind) const;
  /// Return the X processor number given a global X index
  int XPROC(int xind) const;

  /// Communication parameters calculated by topology
  struct ConnectionInfo {
    bool TS_up_in, TS_up_out, TS_down_in, TS_down_out;
    int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
    int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
    int IDATA_DEST, ODATA_DEST; // X inner and outer destinations
  };

  /// Return the communication parameters as calculated by `topology`
  ConnectionInfo getConnectionInfo() const {
    return {TS_up_in,      TS_up_out,     TS_down_in,   TS_down_out,
            UDATA_INDEST,  UDATA_OUTDEST, UDATA_XSPLIT, DDATA_INDEST,
            DDATA_OUTDEST, DDATA_XSPLIT,  IDATA_DEST,   ODATA_DEST};
  }

private:
  // Twist-shift switches
  bool TS_up_in, TS_up_out, TS_down_in, TS_down_out;

  // Communication parameters calculated by topology
  int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
  int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
  int IDATA_DEST, ODATA_DEST; // X inner and outer destinations

  // Settings
  bool TwistShift; // Use a twist-shift condition in core?

  bool symmetricGlobalX; ///< Use a symmetric definition in GlobalX() function
  bool symmetricGlobalY;

  int zperiod;
  BoutReal ZMIN, ZMAX; // Range of the Z domain (in fractions of 2pi)

  int MXG, MYG, MZG; // Boundary sizes

  // Grid file provenance tracking info
  std::string grid_id = "";
  std::string hypnotoad_version = "";
  std::string hypnotoad_git_hash = "";
  std::string hypnotoad_git_diff = "";
  std::string hypnotoad_geqdsk_filename = "";

protected:
  // These are protected so we can make them public in the test suite
  // for testing

  /// Connection initialisation: Set processors in a simple 2D grid
  void default_connections();
  /// Add a topology connection
  ///
  /// Set \p ypos1 and \p ypos2 to be neighbours in the range \p xge <= x < \p xlt.
  /// Optional argument \p ts sets whether to use twist-shift condition
  void set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts = false);
  /// Add a divertor target or limiter
  ///
  /// \p ypos is the y index which will become an upper target.
  /// `ypos+1` will become a lower target.
  /// Target created in the range \p xge <= x < \p xlt.
  void add_target(int ypos, int xge, int xlt);
  /// Create the communication connections between processors to set
  /// the layout of the global grid. See the manual section on "BOUT++
  /// Topology" for more information
  void topology();

  /// Adds 2D and 3D regions for boundaries
  void addBoundaryRegions();

private:
  std::vector<BoundaryRegion*> boundary;        // Vector of boundary regions
  std::vector<BoundaryRegionPar*> par_boundary; // Vector of parallel boundary regions

  //////////////////////////////////////////////////
  // Communications

  bool async_send; ///< Switch to asyncronous sends (ISend, not Send)

  /// Communication handle
  /// Used to keep track of communications between send and receive
  struct CommHandle {
    /// Array of receive requests. One for each possible neighbour; one each way in X, two
    /// each way in Y
    MPI_Request request[6];
    /// Array of send requests (for non-blocking send). One for each possible neighbour;
    /// one each way in X, two each way in Y
    MPI_Request sendreq[6];
    /// Length of the buffers used to send/receive (in BoutReals)
    int xbufflen, ybufflen;
    /// Sending buffers
    Array<BoutReal> umsg_sendbuff, dmsg_sendbuff, imsg_sendbuff, omsg_sendbuff;
    /// Receiving buffers
    Array<BoutReal> umsg_recvbuff, dmsg_recvbuff, imsg_recvbuff, omsg_recvbuff;
    /// Is the communication still going?
    bool in_progress;
    /// Are corner cells included in x-communication?
    bool include_x_corners;
    /// Is there a y-communication
    bool has_y_communication;
    /// List of fields being communicated
    FieldGroup var_list;
  };
  void free_handle(CommHandle* h);
  CommHandle* get_handle(int xlen, int ylen);
  void clear_handles();
  std::list<CommHandle*> comm_list; // List of allocated communication handles

  //////////////////////////////////////////////////
  // X communicator

  /// Communicator containing all processors in X
  MPI_Comm comm_x{MPI_COMM_NULL};

  //////////////////////////////////////////////////
  // Surface communications

  /// Communicator in Y inside both separatrices
  MPI_Comm comm_inner{MPI_COMM_NULL};
  /// Communicator in Y between separatrices
  MPI_Comm comm_middle{MPI_COMM_NULL};
  /// Communicator in Y outside both separatrices
  MPI_Comm comm_outer{MPI_COMM_NULL};

  //////////////////////////////////////////////////
  // Communication routines

  /// Create the MPI requests to receive data in the x-direction. Non-blocking call.
  void post_receiveX(CommHandle& ch);

  /// Create the MPI requests to receive data in the y-direction. Non-blocking call.
  void post_receiveY(CommHandle& ch);

  /// Take data from objects and put into a buffer
  int pack_data(const std::vector<FieldData*>& var_list, int xge, int xlt, int yge,
                int ylt, BoutReal* buffer);
  /// Copy data from a buffer back into the fields

  int unpack_data(const std::vector<FieldData*>& var_list, int xge, int xlt, int yge,
                  int ylt, BoutReal* buffer);
};

namespace {
RegisterMesh<BoutMesh> registermeshbout{"bout"};
}

namespace bout {
/// Result of `bout::checkMesh`: bool + reason for failure
///
/// This would be nicer as std::optional
struct CheckMeshResult {
  bool success;
  std::string reason;
};

/// Check that \p total_processors can be decomposed into \p
/// num_y_processors in Y for the given `BoutMesh` topology parameters
CheckMeshResult checkBoutMeshYDecomposition(int total_processors, int num_y_processors,
                                            int ny, int num_y_guards, int jyseps1_1,
                                            int jyseps2_1, int jyseps1_2, int jyseps2_2,
                                            int ny_inner);
}

#endif // __BOUTMESH_H__
