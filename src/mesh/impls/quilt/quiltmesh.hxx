
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
  QuiltMesh(Options *opt);
  ~QuiltMesh();
  
  int load();

  /// Read in the mesh from data sources
  int load(MPI_Comm comm);

  /////////////////////////////////////////////
  // Get data
  
  int get(Field2D &var, const char *name, BoutReal def=0.0);
  int get(Field2D &var, const string &name, BoutReal def=0.0);
  int get(Field3D &var, const char *name);
  int get(Field3D &var, const string &name);
  
  /////////////////////////////////////////////
  // Communicate variables
  
  int communicate(FieldGroup &g); // Returns error code
  comm_handle send(FieldGroup &g);  // Return handle
  int wait(comm_handle handle); // Wait for the handle, return error code
  
  /////////////////////////////////////////////
  // X communications
  
  bool firstX();
  bool lastX();
  int sendXOut(BoutReal *buffer, int size, int tag);
  int sendXIn(BoutReal *buffer, int size, int tag);
  comm_handle irecvXOut(BoutReal *buffer, int size, int tag);
  comm_handle irecvXIn(BoutReal *buffer, int size, int tag);
  
  MPI_Comm getXcomm() const;
  
  /////////////////////////////////////////////
  // Y-Z surface gather/scatter operations
  SurfaceIter* iterateSurfaces();
  const Field2D averageY(const Field2D &f);
  const Field3D averageY(const Field3D &f);
  
  bool surfaceClosed(int jx, BoutReal &ts); ///< Test if a surface is closed, and if so get the twist-shift angle
  
  // Boundary region iteration
  const RangeIterator iterateBndryLowerY() const;
  const RangeIterator iterateBndryUpperY() const;
  
  // Boundary regions
  vector<BoundaryRegion*> getBoundaries();

  BoutReal GlobalX(int jx); ///< Continuous X index between 0 and 1
  BoutReal GlobalY(int jy); ///< Continuous Y index (0 -> 1)

  void outputVars(Datafile &file); ///< Add mesh vars to file
  
  /// Global locator functions
  int XGLOBAL(int xloc);
  int YGLOBAL(int yloc);
  
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

