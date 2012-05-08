
class MeshFactory;

#ifndef __MESH_FACTORY_H__
#define __MESH_FACTORY_H__

#include <bout/mesh.hxx>
#include <uncopyable.hxx>

class MeshFactory : private Uncopyable {
 public:
  /// Return a pointer to the only instance
  static MeshFactory* getInstance();
  
  Mesh* createMesh(Options *options = NULL);
  
private:
  MeshFactory() {} // Prevent instantiation of this class
  static MeshFactory* instance; ///< The only instance of this class (Singleton)
};

#endif // __LAPLACE_FACTORY_H__

