#ifndef BOUT_GLOBALINDEXER_H
#define BOUT_GLOBALINDEXER_H

#include <bout/mesh.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/region.hxx>
#include <bout_types.hxx>
#include <boutcomm.hxx>
#include "bout/petsclib.hxx"

class GlobalIndexer;
using IndexerPtr = std::shared_ptr<GlobalIndexer>;
using InterpolationWeights = std::vector<ParallelTransform::PositionsAndWeights>;

/*!
 * A singleton which accepts index objects produced by iterating over
 * fields and returns a global index. This index can be used when
 * constructing PETSc arrays. Guard regions used for communication
 * between processes will have the indices of the part of the interior
 * region they are mirroring. Boundaries will have unique indices, but
 * are only included to a depth of 1.
 */
class GlobalIndexer {
public:
  /// If \p localmesh is the same as the global one, return a pointer
  /// to the global instance. Otherwise create a new one.
  static IndexerPtr getInstance(Mesh* localmesh);
  /// Call this immediately after construction when running unit tests.
  void initialiseTest();
  /// Finish setting up the indexer, communicating indices across processes.
  void initialise();
  Mesh* getMesh() { return fieldmesh; }

  /// Convert the local index object to a global index which can be
  /// used in PETSc vectors and matrices.
  int getGlobal(Ind2D ind) const;
  int getGlobal(Ind3D ind) const;
  int getGlobal(IndPerp ind) const;

  // Mark the globalIndexer instance to be recreated next time it is
  // requested (useful for unit tests)
  static void recreateGlobalInstance();

  // Delete the global instance; necessary to also cleanup the
  // PetscLib instance
  static void cleanup();

protected:
  GlobalIndexer(Mesh* localmesh);
  Field3D& getIndices3D() { return indices3D; }
  Field2D& getIndices2D() { return indices2D; }
  FieldPerp& getIndicesPerp() { return indicesPerp; }

private:
  /// This gets called by initialiseTest and is used to register
  /// fields with fake parallel meshes.
  virtual void registerFieldForTest(FieldData& f);
  virtual void registerFieldForTest(FieldPerp& f);

#ifdef BOUT_HAS_PETSC
  // Not sure this is the best approach: maybe hold PetscLib elsewhere?
  PetscLib lib{};
#endif
  Mesh* fieldmesh;

  /// Fields containing the indices for each element (as reals)
  Field3D indices3D;
  Field2D indices2D;
  FieldPerp indicesPerp;

  /// The only instance of this class acting on the global Mesh
  static IndexerPtr globalInstance;
  static bool initialisedGlobal;
};

#endif // BOUT_GLOBALINDEXER_H
