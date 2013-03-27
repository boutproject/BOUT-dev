
#ifndef __BOUTMESH_H__
#define __BOUTMESH_H__

#include "mpi.h"

#include <bout/mesh.hxx>

#include <list>
#include <vector>
#include <cmath>

using std::list;
using std::vector;

class BoutMesh : public Mesh {
 public:
  BoutMesh(GridDataSource *s, Options *options = NULL);
  ~BoutMesh();
  
  /// Read in the mesh from data sources
  int load();
  
  /////////////////////////////////////////////
  // Get data
  
  int get(Field2D &var, const char *name, BoutReal def=0.0);
  int get(Field2D &var, const string &name, BoutReal def=0.0);
  int get(Field3D &var, const char *name);
  int get(Field3D &var, const string &name);
  
  /////////////////////////////////////////////
  // Communicate variables
  
  int communicate(FieldGroup &g);
  comm_handle send(FieldGroup &g);
  int wait(comm_handle handle);
  
  /////////////////////////////////////////////
  // non-local communications
  MPI_Request sendToProc(int xproc, int yproc, BoutReal *buffer, int size, int tag);
  comm_handle receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size, int tag);
  int getNXPE();
  int getNYPE();
  int getXProcIndex();
  int getYProcIndex();
  
  /////////////////////////////////////////////
  // X communications
  
  bool firstX();
  bool lastX();
  int sendXOut(BoutReal *buffer, int size, int tag);
  int sendXIn(BoutReal *buffer, int size, int tag);
  comm_handle irecvXOut(BoutReal *buffer, int size, int tag);
  comm_handle irecvXIn(BoutReal *buffer, int size, int tag);
  
  MPI_Comm getXcomm() const {return comm_x; }
  MPI_Comm getYcomm(int jx) const;
  
  bool periodicY(int jx, BoutReal &ts) const;
  bool periodicY(int jx) const;

  int ySize(int jx) const;

  /////////////////////////////////////////////
  // Y communications
  
  bool firstY();
  bool lastY();
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
  
  /////////////////////////////////////////////
  // Y-Z communications

  const Field2D averageY(const Field2D&);
  const Field3D averageY(const Field3D &f);

  // Boundary iteration
  const RangeIterator iterateBndryLowerY() const;
  const RangeIterator iterateBndryUpperY() const;

  // Boundary regions
  vector<BoundaryRegion*> getBoundaries();

  const Field3D smoothSeparatrix(const Field3D &f);

  BoutReal GlobalX(int jx) const;
  BoutReal GlobalY(int jy) const;

  void outputVars(Datafile &file);

  int XGLOBAL(int xloc) const;
  int YGLOBAL(int yloc) const;

  // poloidal lowpass filtering for n=0 mode
  void slice_r_y(BoutReal *fori, BoutReal * fxy, int ystart, int ncy);
  void get_ri( dcomplex * ayn, int ncy, BoutReal * ayn_Real, BoutReal * ayn_Imag);
  void set_ri( dcomplex * ayn, int ncy, BoutReal * ayn_Real, BoutReal * ayn_Imag);
  const Field2D lowPass_poloidal(const Field2D &var,int mmax);

  //added for volume average and integral
  const Field3D Switch_YZ(const Field3D &var);
  BoutReal Average_XY(const Field2D &var);
  BoutReal Vol_Integral(const Field2D &var);

 private:
  int nx, ny;        ///< Size of the grid in the input file
  int MX, MY;        ///< size of the grid excluding boundary regions
  
  int MYSUB, MXSUB;  ///< Size of the grid on this processor
  
  int NPES; ///< Number of processors
  int MYPE; ///< Rank of this processor

  int PE_YIND; ///< Y index of this processor
  int NYPE; // Number of processors in the Y direction
  
  int MYPE_IN_CORE;  // 1 if processor in core
  
  // Topology
  int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;
  int ixseps_inner, ixseps_outer, ixseps_upper, ixseps_lower;
  int ny_inner;
  
  BoutReal *ShiftAngle;  ///< Angle for twist-shift location
  
  // Processor number, local <-> global translation
  int PROC_NUM(int xind, int yind); // (PE_XIND, PE_YIND) -> MYPE
  bool IS_MYPROC(int xind, int yind);
  int XLOCAL(int xglo) const;
  int YGLOBAL(int yloc, int yproc) const;
  int YLOCAL(int yglo) const;
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
  //int  TwistOrder;   // Order of twist-shift interpolation
  
  bool symmetricGlobalX; ///< Use a symmetric definition in GlobalX() function

  int  zperiod; 
  BoutReal ZMIN, ZMAX;   // Range of the Z domain (in fractions of 2pi)
  
  int  MXG, MYG;     // Boundary sizes

  void default_connections();
  void set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts = false);
  void add_target(int ypos, int xge, int xlt);
  void topology();

  vector<BoundaryRegion*> boundary; // Vector of boundary regions
  
  //////////////////////////////////////////////////
  // Communications
  
  bool async_send;   ///< Switch to asyncronous sends (ISend, not Send)
  
  struct CommHandle {
    MPI_Request request[6];
    /// Array of send requests (for non-blocking send)
    MPI_Request sendreq[6];
    int xbufflen, ybufflen;  ///< Length of the buffers used to send/receive (in BoutReals)
    BoutReal *umsg_sendbuff, *dmsg_sendbuff, *imsg_sendbuff, *omsg_sendbuff;
    BoutReal *umsg_recvbuff, *dmsg_recvbuff, *imsg_recvbuff, *omsg_recvbuff;
    bool in_progress;
    
    /// List of fields being communicated
    vector<FieldData*> var_list;
  };
  void free_handle(CommHandle *h);
  CommHandle* get_handle(int xlen, int ylen);
  void clear_handles();
  list<CommHandle*> comm_list; // List of allocated communication handles

  //////////////////////////////////////////////////
  // X communicator
  
  MPI_Comm comm_x;
  
  //////////////////////////////////////////////////
  // Surface communications
  
  MPI_Comm comm_inner, comm_middle, comm_outer;
  
  //////////////////////////////////////////////////
  // Data reading
  
  /// Read in a portion of the X-Y domain
  int readgrid_3dvar(GridDataSource *s, const char *name, 
	             int yread, int ydest, int ysize, 
                     int xge, int xlt, BoutReal ***var);
  
  /// Copy a section of a 3D variable
  void cpy_3d_data(int yfrom, int yto, int xge, int xlt, BoutReal ***var);

  int readgrid_2dvar(GridDataSource *s, const char *varname, 
                     int yread, int ydest, int ysize, 
                     int xge, int xlt, BoutReal **var);
  void cpy_2d_data(int yfrom, int yto, int xge, int xlt, BoutReal **var);

  //////////////////////////////////////////////////
  // Communication routines
  
  void post_receive(CommHandle &ch);

  /// Take data from objects and put into a buffer
  int pack_data(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, BoutReal *buffer);
  /// Copy data from a buffer back into the fields
  int unpack_data(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, BoutReal *buffer);
};

#endif // __BOUTMESH_H__
