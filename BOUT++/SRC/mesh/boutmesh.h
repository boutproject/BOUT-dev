
#ifndef __BOUTMESH_H__
#define __BOUTMESH_H__

#include "mpi.h"

#include "mesh.h"
#include "field3d.h"

#include <list>
#include <vector>
#include <cmath>

using std::list;
using std::vector;

class BoutMesh : public Mesh {
 public: 
  /// Read in the mesh from data sources
  int load();
  
  /////////////////////////////////////////////
  // Get data
  
  int get(int &ival, const char *name);
  int get(real &rval, const char *name);
  
  int get(Field2D &var, const char *name, real def=0.0);
  int get(Field2D &var, const string &name, real def=0.0);
  int get(Field3D &var, const char *name);
  int get(Field3D &var, const string &name);
  
  /////////////////////////////////////////////
  // Communicate variables
  
  int communicate(FieldGroup &g);
  comm_handle send(FieldGroup &g);
  int wait(comm_handle handle);
  
  /////////////////////////////////////////////
  // X communications
  
  bool first_x();
  bool last_x();
  int send_xout(real *buffer, int size, int tag);
  int send_xin(real *buffer, int size, int tag);
  comm_handle irecv_xout(real *buffer, int size, int tag);
  comm_handle irecv_xin(real *buffer, int size, int tag);
  
  /////////////////////////////////////////////
  // Y-Z communications
  
  SurfaceIter* iterateSurfaces();
  friend class BoutSurfaceIter;
  const Field2D average_y(const Field2D&);
  bool surfaceClosed(int jx, real &ts);

  // Boundary iteration
  RangeIter* iterateBndryLowerY();
  RangeIter* iterateBndryUpperY();
  friend class BoutRangeIter;

  real GlobalX(int jx);
  real GlobalY(int jy);

  void outputVars(Datafile &file);
 private:
  int nx, ny;        ///< Size of the grid in the input file
  int MX, MY;        ///< size of the grid excluding boundary regions
  
  int MYSUB, MXSUB;  ///< Size of the grid on this processor
  
  int NPES; ///< Number of processors
  int MYPE; ///< Rank of this processor

  int PE_YIND; ///< Y index of this processor
  int NYPE; // Number of processors in the Y direction
  
  // Topology
  int ixseps1, ixseps2, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;
  int ixseps_inner, ixseps_outer, ixseps_upper, ixseps_lower;
  int ny_inner;
  
  real *ShiftAngle;  ///< Angle for twist-shift location
  
  // Processor number, local <-> global translation
  int PROC_NUM(int xind, int yind); // (PE_XIND, PE_YIND) -> MYPE
  bool IS_MYPROC(int xind, int yind);
  int XGLOBAL(int xloc);
  int XLOCAL(int xglo);
  int YGLOBAL(int yloc);
  int YGLOBAL(int yloc, int yproc);
  int YLOCAL(int yglo);
  int YLOCAL(int yglo, int yproc);
  int YPROC(int yind);
  int XPROC(int xind);
  
  // Twist-shift switches
  bool TS_up_in, TS_up_out, TS_down_in, TS_down_out;
  
  // Communication parameters calculated by topology
  int UDATA_INDEST, UDATA_OUTDEST, UDATA_XSPLIT;
  int DDATA_INDEST, DDATA_OUTDEST, DDATA_XSPLIT;
  int IDATA_DEST, ODATA_DEST; // X inner and outer destinations
  
  //int MYPE_IN_CORE; // 1 if processor in core (topology.cpp)
  
  // Settings
  bool TwistShift;   // Use a twist-shift condition in core?
  //int  TwistOrder;   // Order of twist-shift interpolation
  
  int  zperiod; 
  real ZMIN, ZMAX;   // Range of the Z domain (in fractions of 2pi)
  
  int  MXG, MYG;     // Boundary sizes

  void default_connections();
  void set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts = false);
  void topology();

  //////////////////////////////////////////////////
  // Communications
  
  bool async_send;   ///< Switch to asyncronous sends (ISend, not Send)
  
  struct CommHandle {
    MPI_Request request[6];
    /// Array of send requests (for non-blocking send)
    MPI_Request sendreq[6];
    int xbufflen, ybufflen;  ///< Length of the buffers used to send/receive (in reals)
    real *umsg_sendbuff, *dmsg_sendbuff, *imsg_sendbuff, *omsg_sendbuff;
    real *umsg_recvbuff, *dmsg_recvbuff, *imsg_recvbuff, *omsg_recvbuff;
    bool in_progress;
    
    /// List of fields being communicated
    vector<FieldData*> var_list;
  };
  void free_handle(CommHandle *h);
  CommHandle* get_handle(int xlen, int ylen);
  void clear_handles();
  list<CommHandle*> comm_list; // List of allocated communication handles

  //////////////////////////////////////////////////
  // Surface communications
  
  MPI_Comm comm_inner, comm_middle, comm_outer;
  
  //////////////////////////////////////////////////
  // Data reading
  
  /// Read in a portion of the X-Y domain
  int readgrid_3dvar(GridDataSource *s, const char *name, 
	             int yread, int ydest, int ysize, 
                     int xge, int xlt, real ***var);
  
  /// Copy a section of a 3D variable
  void cpy_3d_data(int yfrom, int yto, int xge, int xlt, real ***var);

  int readgrid_2dvar(GridDataSource *s, const char *varname, 
                     int yread, int ydest, int ysize, 
                     int xge, int xlt, real **var);
  void cpy_2d_data(int yfrom, int yto, int xge, int xlt, real **var);

  //////////////////////////////////////////////////
  // Communication routines
  
  void post_receive(CommHandle &ch);

  /// Take data from objects and put into a buffer
  int pack_data(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, real *buffer);
  /// Copy data from a buffer back into the fields
  int unpack_data(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt, real *buffer);
  /// Calculates the size of a message for a given x and y range
  int msg_len(vector<FieldData*> &var_list, int xge, int xlt, int yge, int ylt);
};

/*
class BoutSurfaceIter : public SurfaceIter {
 public:
  BoutSurfaceIter(BoutMesh* mi);
  int ysize(); // Return the size of the current surface
  bool closed(real &ts);

  void first();
  void next();
  bool isDone();
  
  int gather(const Field2D &f, real *data);
  int gather(const Field3D &f, real **data);
  
  int scatter(real *data, Field2D &f);
  int scatter(real **data, Field3D &f);
 private:
  BoutMesh* m;
};
*/

class BoutRangeIter : public RangeIter {
 public:
  BoutRangeIter(int start, int end);
  void first();
  void next();
  bool isDone();
  
  int ind; // The index
  
 private:
  int s,e;
};

#endif // __BOUTMESH_H__
