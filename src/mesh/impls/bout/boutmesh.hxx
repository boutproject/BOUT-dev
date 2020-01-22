
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
  int getGlobalYIndex(int ylocal) const override;
  int getGlobalYIndexNoBoundaries(int ylocal) const override;
  int getGlobalZIndex(int zlocal) const override;
  int getGlobalZIndexNoBoundaries(int zlocal) const override;

  int XLOCAL(int xglo) const override;
  int YLOCAL(int yglo) const override;

  // Switch for communication of corner guard and boundary cells
  const bool include_corner_cells;

protected:
  BoutMesh(int input_nx, int input_ny, int input_nz, int mxg, int myg, int nxpe, int nype,
           int pe_xind, int pe_yind, bool include_corners=true);
  /// For debugging purposes (when creating fake parallel meshes), make
  /// the send and receive buffers share memory. This allows for
  /// communications to be faked between meshes as though they were on
  /// different processors.
  void overlapHandleMemory(BoutMesh* yup, BoutMesh* ydown, BoutMesh* xin,
			   BoutMesh* xout, BoutMesh* xinyup, BoutMesh* xinydown,
			   BoutMesh* xoutyup, BoutMesh* xoutydown,
			   bool xinyupSendsInner, bool xinydownSendsInner,
			   bool xoutyupSendsInner, bool xoutydownSendsInner);

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

  // Processor number, local <-> global translation
  int PROC_NUM(int xind, int yind); // (PE_XIND, PE_YIND) -> MYPE
  int YGLOBAL(int yloc, int yproc) const;
  int YLOCAL(int yglo, int yproc) const;
  int YPROC(int yind);
  int XPROC(int xind);

  // Twist-shift switches
  bool TS_up_in, TS_up_out, TS_down_in, TS_down_out;

  // Communication parameters calculated by topology
  int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
  int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
  int IDATA_DEST, ODATA_DEST; // X inner and outer destinations
  int lower_inner_corner_dest, upper_inner_corner_dest, lower_outer_corner_dest,
      upper_outer_corner_dest; // destinations for the corner cells
  int lower_inner_corner_orig, upper_inner_corner_orig, lower_outer_corner_orig,
      upper_outer_corner_orig; // origins for the corner guard cells
  // y-limits of buffers communicated in x-direction. Include y-boundary cells but not
  // y-guard cells. Need different variables for sending and receiving because y-boundary
  // might be present on sending proc but not receiving proc or vice versa
  int IDATA_buff_lowerY_send, IDATA_buff_upperY_send, ODATA_buff_lowerY_send,
      ODATA_buff_upperY_send;
  int IDATA_buff_lowerY_recv, IDATA_buff_upperY_recv, ODATA_buff_lowerY_recv,
      ODATA_buff_upperY_recv;
  // x-limits of buffers communicated in y-direction. Include x-boundary cells but not
  // x-guard cells.
  int YDATA_buff_innerX, YDATA_buff_outerX;

  // Settings
  bool TwistShift; // Use a twist-shift condition in core?

  bool symmetricGlobalX; ///< Use a symmetric definition in GlobalX() function
  bool symmetricGlobalY;

  int zperiod;
  BoutReal ZMIN, ZMAX; // Range of the Z domain (in fractions of 2pi)

  int MXG, MYG, MZG; // Boundary sizes

  void default_connections();
  void set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts = false);
  void add_target(int ypos, int xge, int xlt);
  void topology();

  void addBoundaryRegions(); ///< Adds 2D and 3D regions for boundaries

  std::vector<BoundaryRegion*> boundary;        // Vector of boundary regions
  std::vector<BoundaryRegionPar*> par_boundary; // Vector of parallel boundary regions

  //////////////////////////////////////////////////
  // Communications

  bool async_send; ///< Switch to asyncronous sends (ISend, not Send)

  /// Communication handle
  /// Used to keep track of communications between send and receive
  struct CommHandle {
    /// Array of receive requests. One for each possible neighbour; one each way in X, two
    /// each way in Y. Plus four for the corners.
    MPI_Request request[10] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL};
    /// Array of send requests (for non-blocking send). One for each possible neighbour;
    /// one each way in X, two each way in Y. Plus four for the corners.
    MPI_Request sendreq[10] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                               MPI_REQUEST_NULL};
    /// Length of the buffers used to send/receive (in BoutReals)
    int xbufflen, ybufflen, cornerbufflen;
    /// Sending buffers
    Array<BoutReal> umsg_sendbuff, dmsg_sendbuff, imsg_sendbuff, omsg_sendbuff,
                    lowin_corner_sendbuff, upin_corner_sendbuff, lowout_corner_sendbuff,
                    upout_corner_sendbuff;
    /// Receiving buffers
    Array<BoutReal> umsg_recvbuff, dmsg_recvbuff, imsg_recvbuff, omsg_recvbuff,
                    lowin_corner_recvbuff, upin_corner_recvbuff, lowout_corner_recvbuff,
                    upout_corner_recvbuff;
    /// Is the communication still going?
    bool in_progress;
    /// List of fields being communicated
    FieldGroup var_list;
  };
  void free_handle(CommHandle* h);
  CommHandle* get_handle(int xlen, int ylen, int cornerlen = 0);
  void clear_handles();
  std::list<CommHandle*> comm_list; // List of allocated communication handles

  //////////////////////////////////////////////////
  // X communicator

  MPI_Comm comm_x; ///< Communicator containing all processors in X

  //////////////////////////////////////////////////
  // Surface communications

  /// Communicators in Y. Inside both separatrices; between separatrices;
  /// and outside both separatrices
  MPI_Comm comm_inner, comm_middle, comm_outer;

  //////////////////////////////////////////////////
  // Communication routines

  /// Create the MPI requests to receive data. Non-blocking call.
  void post_receive(CommHandle& ch);

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

#endif // __BOUTMESH_H__
