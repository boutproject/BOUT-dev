
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
  ~BoutMesh();

  /// Read in the mesh from data sources
  int load();

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
  comm_handle send(FieldGroup &g);

  /// Wait for a send operation to complete
  /// @param[in] handle  The handle returned by send()
  int wait(comm_handle handle);

  /////////////////////////////////////////////
  // non-local communications

  MPI_Request sendToProc(int xproc, int yproc, BoutReal *buffer, int size, int tag);
  comm_handle receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag);

  int getNXPE(); ///< The number of processors in the X direction
  int getNYPE(); ///< The number of processors in the Y direction
  int getXProcIndex();  ///< This processor's index in X direction
  int getYProcIndex();  ///< This processor's index in Y direction

  /////////////////////////////////////////////
  // X communications

  bool firstX(); ///< Is this processor the first in X? i.e. is there a boundary to the left in X?
  bool lastX();  ///< Is this processor last in X? i.e. is there a boundary to the right in X?

  /// Send a buffer of data to processor at X index +1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  int sendXOut(BoutReal *buffer, int size, int tag);

  /// Send a buffer of data to processor at X index -1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of BoutReals to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  int sendXIn(BoutReal *buffer, int size, int tag);

  /// Receive a buffer of data from X index +1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  comm_handle irecvXOut(BoutReal *buffer, int size, int tag);

  /// Receive a buffer of data from X index -1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of BoutReals to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  comm_handle irecvXIn(BoutReal *buffer, int size, int tag);

  /// Send a buffer of complex data to processor at X index +1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of dcomplex to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  comm_handle isendXOut(dcomplex *buffer, int size, int tag);

  /// Send a buffer of complex data to processor at X index -1
  ///
  /// @param[in] buffer  The data to send. Must be at least length \p size
  /// @param[in] size    The number of dcomplex to send
  /// @param[in] tag     A label for the communication. Must be the same at receive
  comm_handle isendXIn(dcomplex *buffer, int size, int tag);

  /// Receive a buffer of complex data from X index +1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of dcomplex to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  comm_handle irecvXOut(dcomplex *buffer, int size, int tag);

  /// Receive a buffer of complex data from X index -1
  ///
  /// @param[in] buffer  A buffer to put the data in. Must already be allocated of length \p size
  /// @param[in] size    The number of dcomplex to receive and put in \p buffer
  /// @param[in] tag     A label for the communication. Must be the same as sent
  comm_handle irecvXIn(dcomplex *buffer, int size, int tag);

  MPI_Comm getXcomm(int UNUSED(jy)) const {return comm_x; } ///< Return communicator containing all processors in X
  MPI_Comm getYcomm(int jx) const; ///< Return communicator containing all processors in Y

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  /// \param[out] ts  The Twist-Shift angle if periodic
  bool periodicY(int jx, BoutReal &ts) const;

  /// Is local X index \p jx periodic in Y?
  ///
  /// \param[in] jx   The local (on this processor) index in X
  bool periodicY(int jx) const;

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

  int ySize(int jx) const; ///< The number of points in Y at fixed X index \p jx

  /////////////////////////////////////////////
  // Y communications

  bool firstY() const;
  bool lastY() const;
  bool firstY(int xpos) const;
  bool lastY(int xpos) const;
  int UpXSplitIndex();
  int DownXSplitIndex();
  int sendYOutIndest(BoutReal *buffer, int size, int tag);
  int sendYOutOutdest(BoutReal *buffer, int size, int tag);
  int sendYInIndest(BoutReal *buffer, int size, int tag);
  int sendYInOutdest(BoutReal *buffer, int size, int tag);
  comm_handle irecvYOutIndest(BoutReal *buffer, int size, int tag);
  comm_handle irecvYOutOutdest(BoutReal *buffer, int size, int tag);
  comm_handle irecvYInIndest(BoutReal *buffer, int size, int tag);
  comm_handle irecvYInOutdest(BoutReal *buffer, int size, int tag);

  // Boundary iteration
  const RangeIterator iterateBndryLowerY() const;
  const RangeIterator iterateBndryUpperY() const;
  const RangeIterator iterateBndryLowerInnerY() const;
  const RangeIterator iterateBndryLowerOuterY() const;
  const RangeIterator iterateBndryUpperInnerY() const;
  const RangeIterator iterateBndryUpperOuterY() const;


  // Boundary regions
  std::vector<BoundaryRegion*> getBoundaries();
  std::vector<BoundaryRegionPar*> getBoundariesPar();
  void addBoundaryPar(BoundaryRegionPar* bndry);

  const Field3D smoothSeparatrix(const Field3D &f);

  int getNx() const {return nx;}
  int getNy() const {return ny;}

  BoutReal GlobalX(int jx) const;
  BoutReal GlobalY(int jy) const;
  BoutReal GlobalX(BoutReal jx) const;
  BoutReal GlobalY(BoutReal jy) const;

  BoutReal getIxseps1() const {return ixseps1;}
  BoutReal getIxseps2() const {return ixseps2;}

  void outputVars(Datafile &file);

  int XGLOBAL(int xloc) const;
  int YGLOBAL(int yloc) const;
  int XGLOBAL(BoutReal xloc, BoutReal &xglo) const;
  int YGLOBAL(BoutReal yloc, BoutReal &yglo) const;

  int XLOCAL(int xglo) const;
  int YLOCAL(int yglo) const;

 private:
  std::string gridname;
  int nx, ny, nz; ///< Size of the grid in the input file
  int MX, MY, MZ; ///< size of the grid excluding boundary regions

  int MYSUB, MXSUB, MZSUB; ///< Size of the grid on this processor

  int NPES; ///< Number of processors
  int MYPE; ///< Rank of this processor

  int PE_YIND; ///< Y index of this processor
  int NYPE; // Number of processors in the Y direction

  int NZPE;

  int MYPE_IN_CORE;  // 1 if processor in core

  // Topology
  int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;
  int ixseps_inner, ixseps_outer, ixseps_upper, ixseps_lower;
  int ny_inner;

  std::vector<BoutReal> ShiftAngle;  ///< Angle for twist-shift location

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

  // Settings
  bool TwistShift;   // Use a twist-shift condition in core?

  bool symmetricGlobalX; ///< Use a symmetric definition in GlobalX() function
  bool symmetricGlobalY;

  int  zperiod;
  BoutReal ZMIN, ZMAX;   // Range of the Z domain (in fractions of 2pi)

  int MXG, MYG, MZG; // Boundary sizes

  void default_connections();
  void set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts = false);
  void add_target(int ypos, int xge, int xlt);
  void topology();
  
  void addBoundaryRegions(); ///< Adds 2D and 3D regions for boundaries
  
  std::vector<BoundaryRegion*> boundary; // Vector of boundary regions
  std::vector<BoundaryRegionPar*> par_boundary; // Vector of parallel boundary regions

  //////////////////////////////////////////////////
  // Communications

  bool async_send;   ///< Switch to asyncronous sends (ISend, not Send)

  /// Communication handle
  /// Used to keep track of communications between send and receive
  struct CommHandle {
    /// Array of receive requests. One for each possible neighbour; one each way in X, two each way in Y
    MPI_Request request[6];
    /// Array of send requests (for non-blocking send). One for each possible neighbour; one each way in X, two each way in Y
    MPI_Request sendreq[6];
    int xbufflen, ybufflen;  ///< Length of the buffers used to send/receive (in BoutReals)
    Array<BoutReal> umsg_sendbuff, dmsg_sendbuff, imsg_sendbuff, omsg_sendbuff; ///< Sending buffers
    Array<BoutReal> umsg_recvbuff, dmsg_recvbuff, imsg_recvbuff, omsg_recvbuff; ///< Receiving buffers
    bool in_progress; ///< Is the communication still going?

    /// List of fields being communicated
    FieldGroup var_list;
  };
  void free_handle(CommHandle *h);
  CommHandle* get_handle(int xlen, int ylen);
  void clear_handles();
  std::list<CommHandle*> comm_list; // List of allocated communication handles

  //////////////////////////////////////////////////
  // X communicator

  MPI_Comm comm_x; ///< Communicator containing all processors in X

  //////////////////////////////////////////////////
  // Surface communications

  MPI_Comm comm_inner, comm_middle, comm_outer; ///< Communicators in Y. Inside both separatrices; between separatrices; and outside both separatrices

  //////////////////////////////////////////////////
  // Communication routines

  /// Create the MPI requests to receive data. Non-blocking call.
  void post_receive(CommHandle &ch);

  /// Take data from objects and put into a buffer
  int pack_data(const std::vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, BoutReal *buffer);
  /// Copy data from a buffer back into the fields
  int unpack_data(const std::vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, BoutReal *buffer);
};

#endif // __BOUTMESH_H__
