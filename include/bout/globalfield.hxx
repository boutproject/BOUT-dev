/*!
 * Provides global gather/scatter operations for fields
 *
 */

class GlobalField;
class GlobalField2D;

#ifndef __GLOBALFIELD_H__
#define __GLOBALFIELD_H__

#include "mesh.hxx"

/*!
 * This provides a method for gathering and scattering a field
 * which takes into account the local and global indices
 * 
 * This is a base class which is inherited by GlobalField2D and GlobalField3D
 */ 
class GlobalField {
public:
  GlobalField() = delete;
  virtual ~GlobalField() = default;
  virtual bool valid() const = 0;  ///< Is the data valid on any processor?
  bool dataIsLocal() const {return valid() && (data_on_proc == mype);} ///< Data is on this processor

  /*!
   * Data access by index. This doesn't perform any checks,
   * so the user should first test if the data is available
   * on this processor by calling dataIsLocal()
   */
  BoutReal& operator()(int jx, int jy, int jz) {return data[jz + nz*jy + nz*ny*jx];}

  /// Const data access by index
  const BoutReal& operator()(int jx, int jy, int jz) const {return data[jz + nz*jy + nz*ny*jx];}
  
  /// Size of the field in X
  int xSize() const {return nx;}
  /// Size of the field in Y
  int ySize() const {return ny;}
  /// Size of the field in Z
  int zSize() const {return nz;}
  
  /*
   * Direct data access
   *
   * Stored as a 1D array [x*ny*nz + y*nz + z]
   */ 
  Array<BoutReal>& getData() {return data;}
protected:
  GlobalField(Mesh *m, int proc, int xsize, int ysize, int zsize);
  
  Mesh *mesh; ///< The mesh we're gathering/scattering over

  int data_on_proc; ///< Which processor is this data on?  
  int nx, ny, nz; ///< Global field sizes
  Array<BoutReal> data; ///< The global data, if on this processor

  MPI_Comm comm; ///< Communicator for all mesh
  int npes, mype; ///< Number of MPI processes, this processor index

  void proc_local_origin(int proc, int *x, int *y, int *z = nullptr) const;
  /// Return the global origin of processor proc
  void proc_origin(int proc, int *x, int *y, int *z = nullptr) const;
  /// Return the array size of processor proc
  void proc_size(int proc, int *lx, int *ly, int *lz = nullptr) const;
};

/*!
 * Gather and scatter a Field2D
 *
 * Example
 * -------
 *
 * To create a GlobalField2D, pass a mesh pointer
 * 
 *     GlobalField2D g2d(mesh); 
 *
 * By default data is gathered and scattered to/from processor 0. 
 * To change this, pass the processor number as a second argument:
 *
 *     GlobalField3D g2d(mesh, 1); // Gather onto processor 1
 *
 * Gather and scatter methods operate on Field2D objects:
 *
 *     Field2D localdata;
 *     
 *     g2d.gather(localdata); // Gather onto one processsor
 *
 * To scatter data back, use the scatter method:
 *
 *     localdata = g2d.scatter();
 *
 * Note that both gather and scatter are collective operations, 
 * which must be performed by all processors.
 *
 * To test if the data is available on a processor, use:
 *
 *     if(g2d.dataIsLocal()) {
 *       // g2d data on this processor
 *     }
 * 
 * The data in a GlobalField2D can be accessed using (x,y,z) indexing,
 * with the index ranges given by xSize, ySize, zSize methods. 
 * 
 *     for ( int x=0; x<g2d.xSize(); x++)
 *       for ( int y=0; y<g2d.ySize(); y++)
 *         output.write(" Value at ({:d} ,{:d}) is {:e}\n" ,
 *                        x, y,
 *                        g2d(x, y));
 *
 */ 
class GlobalField2D : public GlobalField {
public:
  /// Can't be constructed without args
  GlobalField2D() = delete;
  /// Construct, giving a mesh and an optional processor
  ///
  /// @param[in] mesh   The mesh to gather over
  /// @param[in] proc   The processor index where everything will be gathered/scattered to/from
  GlobalField2D(Mesh *mesh, int proc = 0);

  /// Destructor
  ~GlobalField2D() override;

  /// Is the data valid and on this processor?
  bool valid() const override { return data_valid; }

  /// Gather all data onto one processor
  void gather(const Field2D &f);
  /// Scatter data back from one to many processors
  const Field2D scatter() const;
  
  /// Assignment from a 2D field. Shorthand for a gather, and must be called on all processors
  /// The scatter assignment operator needs to be a member of Field2D.
  GlobalField2D& operator=(const Field2D &rhs) {
    gather(rhs);
    return *this;
  }
  
  /// Data access by global index
  BoutReal& operator()(int jx, int jy) {return GlobalField::operator()(jx, jy, 0);}
  const BoutReal& operator()(int jx, int jy) const {return GlobalField::operator()(jx, jy, 0);}
  
protected:
  
private:

  /// Buffer for sending and receiving. First index
  /// is the processor index, and second is the data
  BoutReal** buffer;
  
  /// The length of message (in BoutReals) to
  /// be sent to or from processor \p proc
  ///
  /// @param[in] proc  MPI processor index
  int msg_len(int proc) const;

  /// Is the data valid and on this processor?
  bool data_valid;
};

/*!
 * Gather and scatter a Field3D to/from one processor
 *
 * Example
 * -------
 *
 * To create a GlobalField3D, pass a mesh pointer
 * 
 *     GlobalField3D g3d(mesh); 
 *
 * By default data is gathered and scattered to/from processor 0. 
 * To change this, pass the processor number as a second argument:
 *
 *     GlobalField3D g3d(mesh, 1); // Gather onto processor 1
 *
 * Gather and scatter methods operate on Field3D objects:
 *
 *     Field3D localdata;
 *     
 *     g3d.gather(localdata); // Gather onto one processsor
 *
 * To scatter data back, use the scatter method:
 *
 *     localdata = g3d.scatter();
 *
 * Note that both gather and scatter are collective operations, 
 * which must be performed by all processors.
 *
 * To test if the data is available on a processor, use:
 *
 *     if(g3d.dataIsLocal()) {
 *       // g3d data on this processor
 *     }
 * 
 * The data in a GlobalField3D can be accessed using (x,y,z) indexing,
 * with the index ranges given by xSize, ySize, zSize methods. 
 * 
 *     for ( int x=0; x<g3d.xSize(); x++)
 *       for ( int y=0; y<g3d.ySize(); y++)
 *         for ( int z=0; z<g3d.zSize(); z++)
 *           output.write(" Value at ({:d} ,{:d} ,{:d}) is {:e}\n" ,
 *                        x, y, z,
 *                        g3d(x, y, z));
 *
 */
class GlobalField3D : public GlobalField {
public:
  /// Can't be constructed without args
  GlobalField3D() = delete;

  /// Construct, giving a mesh and an optional processor
  ///
  /// @param[in] mesh   The mesh to gather over
  /// @param[in] proc   The processor index where everything will be gathered/scattered to/from
  GlobalField3D(Mesh *mesh, int proc = 0);

  /// Destructor
  ~GlobalField3D() override;

  /// Test if the data is valid i.e. has been allocated
  bool valid() const override { return data_valid; }

  /// Gather all data onto one processor
  void gather(const Field3D &f);
  /// Scatter data back from one to many processors
  const Field3D scatter() const;
  
  /// Assignment from a 2D field. Shorthand for a gather, and must be called on all processors
  /// The scatter assignment operator needs to be a member of Field2D.
  GlobalField3D& operator=(const Field3D &rhs) {
    gather(rhs);
    return *this;
  }
  
protected:
  
private:

  /// Buffer for sending and receiving. First index
  /// is the processor index, and second is the data
  BoutReal** buffer; 

  /// The length of message (in BoutReals) to
  /// be sent to or from processor \p proc
  ///
  /// @param[in] proc  MPI processor index
  int msg_len(int proc) const;

  /// Is the data valid and on this processor?
  bool data_valid;
};

#endif // __GLOBALFIELD_H__
