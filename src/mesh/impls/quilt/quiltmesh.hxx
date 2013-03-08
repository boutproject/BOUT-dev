
#ifndef __QUILTMESH_H__
#define __QUILTMESH_H__

#include <bout/mesh.hxx>
#include <options.hxx>
#include <bout/sys/range.hxx>

#include "quiltdomain.hxx"

/*!
 * Mesh composed of a patchwork of regions. Generalised
 * code which extends the topology of BoutMesh
 */
class QuiltMesh : public Mesh {
 public:
  QuiltMesh(GridDataSource *s, Options *opt);
  ~QuiltMesh();
  
  int load();

  /// Read in the mesh from data sources
  int load(MPI_Comm comm);
  
  /////////////////////////////////////////////
  // Communicate variables
  
  int communicate(FieldGroup &g); // Returns error code
  comm_handle send(FieldGroup &g);  // Return handle
  int wait(comm_handle handle); // Wait for the handle, return error code
  
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
  
  MPI_Comm getXcomm() const;
  MPI_Comm getYcomm(int jx) const;
  
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
  // Y-Z surface gather/scatter operations
  const Field2D averageY(const Field2D &f);
  const Field3D averageY(const Field3D &f);
  
  bool periodicY(int jx, BoutReal &ts) const;
  
  // Boundary region iteration
  const RangeIterator iterateBndryLowerY() const;
  const RangeIterator iterateBndryUpperY() const;
  
  // Boundary regions
  vector<BoundaryRegion*> getBoundaries();

  BoutReal GlobalX(int jx) const; ///< Continuous X index between 0 and 1
  BoutReal GlobalY(int jy) const; ///< Continuous Y index (0 -> 1)

  void outputVars(Datafile &file); ///< Add mesh vars to file
  
  /// Global locator functions
  int XGLOBAL(int xloc) const;
  int YGLOBAL(int yloc) const;
  
 private:
  Options *options;  // Configuration options to use
  
  // Settings
  bool async_send;
  int MXG, MYG;
  
  // Describes regions of the mesh and connections between them
  QuiltDomain *mydomain;
  
  /// Describes guard cell regions
  struct GuardRange {
    int xmin, xmax;
    int ymin, ymax;
    
    int proc; // Neighbour
  };
  vector<GuardRange*> all_guards;
  vector<GuardRange*> xlow_guards;
  vector<GuardRange*> xhigh_guards;
  
  RangeIterator xlow_range;  // X lower boundary 
  RangeIterator xhigh_range; // X upper boundary
  RangeIterator ylow_range;  // Lower Y boundary
  RangeIterator yhigh_range; // Upper Y boundary
  
  vector<BoundaryRegion*> boundaries; // All boundaries

  RangeIterator yperiodic;  // Range of X for which Y is periodic
  vector<BoutReal> twistshift; // Value of twist-shift for each X location
  
  /// A single MPI request with metadata
  struct QMRequest {
    MPI_Request request;
    vector<BoutReal> buffer; // Buffer to store the incoming data
    GuardRange* guard; ///< Where is the data going?
  };

  /// Handle for communications
  struct QMCommHandle {
    bool inProgress;   ///< Denotes if communication is in progress
    vector<FieldData*> var_list; ///< List of fields being communicated
    
    // May have several different requests
    vector<QMRequest*> request;
    vector<MPI_Request> mpi_rq;  ///< All requests, so can use MPI_Waitall
  };
  
  QMCommHandle *getCommHandle();
  void freeCommHandle(QMCommHandle *h);
  QMRequest* getRequestHandle(int n);
  void freeRequestHandle(QMRequest* r);
  
  list<QMCommHandle*> comm_list;

  

  void packData(const vector<FieldData*> &vars, GuardRange* range, vector<BoutReal> &data);
  void unpackData(vector<BoutReal> &data, GuardRange* range, vector<FieldData*> &vars);
  
  void read2Dvar(GridDataSource *s, const char *name, 
                 int xs, int ys,
                 int xd, int yd,
                 int nx, int ny,
                 BoutReal **data);
                 
};

#endif // __QUILTMESH_H__

