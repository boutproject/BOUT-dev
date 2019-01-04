#ifndef __MESH_FACTORY_H__
#define __MESH_FACTORY_H__

#include "bout/sys/uncopyable.hxx"

class GridDataSource;
class Mesh;
class Options;

class MeshFactory : private Uncopyable {
 public:
  /// Return a pointer to the only instance
  static MeshFactory* getInstance();

  Mesh *createMesh(GridDataSource *source, Options *options = nullptr);

private:
  MeshFactory() {} // Prevent instantiation of this class
  static MeshFactory* instance; ///< The only instance of this class (Singleton)
};

#endif // __LAPLACE_FACTORY_H__

