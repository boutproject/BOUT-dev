#ifndef BOUT_GLOBALINDEXER_H
#define BOUT_GLOBALINDEXER_H

#include <bout/mesh.hxx>
#include <bout/operatorstencil.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/region.hxx>
#include <bout/traits.hxx>
#include <bout_types.hxx>
#include <boutcomm.hxx>

template <class T>
class GlobalIndexer;

template <class T>
using IndexerPtr = std::shared_ptr<GlobalIndexer<T>>;

using InterpolationWeights = std::vector<ParallelTransform::PositionsAndWeights>;

/*!
 * An object which accepts index objects produced by iterating over
 * fields and returns a global index. This index can be used when
 * constructing PETSc arrays. Guard regions used for communication
 * between processes will have the indices of the part of the interior
 * region they are mirroring. Boundaries required by the stencil will
 * have unique indices. If no stencil is provided then only the
 * interior region will be assigned global indices. By default, the
 * indexer is fully initialised so that guard cells are communicated
 * (ensuring they hold the appropriate global indices). However, by
 * passing `autoInitialise = false` this behaviour can be turned off
 * and the user can then manually call the `initialise()` method
 * later. This can be useful for mocking/faking the class when
 * testing.
 */
template <class T>
class GlobalIndexer {
public:
  static_assert(bout::utils::is_Field<T>::value, "GlobalIndexer only works with Fields");
  using ind_type = typename T::ind_type;

  GlobalIndexer() = default;

  GlobalIndexer(Mesh* localmesh,
                OperatorStencil<ind_type> stencil = OperatorStencil<ind_type>(),
                bool autoInitialise = true)
      : fieldmesh(localmesh), indices(-1., localmesh), stencils(std::move(stencil)) {

    Region<ind_type> allCandidate, bndryCandidate;
    if (stencils.getNumParts() > 0) {
      std::set<ind_type> allIndices(getRegionNobndry().getIndices().begin(),
                                    getRegionNobndry().getIndices().end()),
          newIndices;
      BOUT_FOR_SERIAL(i, getRegionNobndry()) {
        for (const IndexOffset<ind_type> j : stencils.getStencilPart(i)) {
          insertIndex(i + j, allIndices, newIndices);
        }
      }
      std::set<ind_type> candidateIndices = newIndices;
      while (candidateIndices.size() > 0) {
        newIndices.clear();
        for (const ind_type i : candidateIndices) {
          insertIndex(i, allIndices, newIndices);
        }
        candidateIndices = newIndices;
      }
      std::vector<ind_type> allIndicesVec(allIndices.begin(), allIndices.end());
      allCandidate = Region<ind_type>(allIndicesVec);
    }

    bndryCandidate = mask(allCandidate, getRegionNobndry());

    regionInnerX = getUnion(bndryCandidate, indices.getRegion("RGN_INNER_X"));
    regionOuterX = getUnion(bndryCandidate, indices.getRegion("RGN_OUTER_X"));
    if (std::is_same<T, FieldPerp>::value) {
      regionLowerY = Region<ind_type>({});
      regionUpperY = Region<ind_type>({});
    } else {
      regionLowerY = getUnion(bndryCandidate, indices.getRegion("RGN_LOWER_Y"));
      regionUpperY = getUnion(bndryCandidate, indices.getRegion("RGN_UPPER_Y"));
    }
    regionBndry = regionLowerY + regionInnerX + regionOuterX + regionUpperY;
    regionAll = getRegionNobndry() + regionBndry;
    regionBndry.sort();
    regionAll.sort();

    int localSize = size();
    MPI_Comm comm =
        std::is_same<T, FieldPerp>::value ? fieldmesh->getXcomm() : BoutComm::get();
    fieldmesh->getMpi().MPI_Scan(&localSize, &globalEnd, 1, MPI_INT, MPI_SUM, comm);
    globalEnd--;
    int counter = globalStart = globalEnd - size() + 1;

    BOUT_FOR_SERIAL(i, regionAll) { indices[i] = counter++; }

    if (autoInitialise) {
      initialise();
    }
  }

  virtual ~GlobalIndexer() {}
  
  /// Call this immediately after construction when running unit tests.
  void initialiseTest() {}

  /// Finish setting up the indexer, communicating indices across
  /// processes and, if possible, calculating the sparsity pattern of
  /// any matrices.
  void initialise() {
    fieldmesh->communicate(indices);
  }

  Mesh* getMesh() const { return fieldmesh; }

  /// Convert the local index object to a global index which can be
  /// used in PETSc vectors and matrices.
  int getGlobal(const ind_type& ind) const {
    return static_cast<int>(std::round(indices[ind]));
  }

  /// Check whether the local index corresponds to an element which is
  /// stored locally.
  bool isLocal(const ind_type& ind) const {
    if (ind.ind < 0) {
      return false;
    }
    int index = getGlobal(ind);
    return (globalStart <= index) && (index <= globalEnd);
  }

  int getGlobalStart() const { return globalStart; }

  const Region<ind_type>& getRegionAll() const { return regionAll; }
  const Region<ind_type>& getRegionNobndry() const {
    return indices.getRegion("RGN_NOBNDRY");
  }
  const Region<ind_type>& getRegionBndry() const { return regionBndry; }
  const Region<ind_type>& getRegionLowerY() const { return regionLowerY; }
  const Region<ind_type>& getRegionUpperY() const { return regionUpperY; }
  const Region<ind_type>& getRegionInnerX() const { return regionInnerX; }
  const Region<ind_type>& getRegionOuterX() const { return regionOuterX; }

  bool sparsityPatternAvailable() const { return stencils.getNumParts() > 0; }

  const std::vector<int>& getNumDiagonal() const {
    ASSERT2(sparsityPatternAvailable());
    if (!sparsityCalculated) {
      calculateSparsity();
    }
    return numDiagonal;
  }

  const std::vector<int>& getNumOffDiagonal() const {
    ASSERT2(sparsityPatternAvailable());
    if (!sparsityCalculated) {
      calculateSparsity();
    }
    return numOffDiagonal;
  }

  int size() const { return regionAll.size(); }

protected:
  // Must not be const as the index field needs to be mutable in order
  // to fake parallel communication in the unit tests.
  T& getIndices() { return indices; }

private:
  static void insertIndex(const ind_type i, std::set<ind_type>& allInds,
                          std::set<ind_type>& newInds) {
    auto result = allInds.insert(i);
    if (result.second) {
      newInds.insert(i);
    }
  }

  /// This gets called by initialiseTest and is used to register
  /// fields with fake parallel meshes.
  virtual void registerFieldForTest(T& UNUSED(f)) {
    // This is a place-holder which does nothing. It can be overridden
    // by descendent classes if necessary to set up testing.
    return;
  }

  void calculateSparsity() const {
    numDiagonal = std::vector<int>(size());
    numOffDiagonal = std::vector<int>(size(), 0);

    // Set initial guess for number of on-diagonal elements
    BOUT_FOR(i, regionAll) {
      numDiagonal[getGlobal(i) - globalStart] = stencils.getStencilSize(i);
    }

    BOUT_FOR_SERIAL(i, indices.getRegion("RGN_GUARDS")) {
      if (getGlobal(i) >= 0 && !isLocal(i)) {
        for (const auto& j : stencils.getIndicesWithStencilIncluding(i)) {
          if (isLocal(j)) {
            const int n = getGlobal(j) - globalStart;
            numDiagonal[n] -= 1;
            numOffDiagonal[n] += 1;
          }
        }
      }
    }

    sparsityCalculated = true;
  }

  Mesh* fieldmesh;

  /// Fields containing the indices for each element (as reals)
  T indices;
  /// The first and last global index on this processor (inclusive in
  /// both cases)
  int globalStart, globalEnd;

  /// Stencil for which this indexer has been configured
  OperatorStencil<ind_type> stencils;

  /// Regions containing the elements for which there are global indices
  Region<ind_type> regionAll, regionLowerY, regionUpperY, regionInnerX, regionOuterX,
      regionBndry;

  mutable bool sparsityCalculated = false;
  mutable std::vector<int> numDiagonal, numOffDiagonal;
};

#endif // BOUT_GLOBALINDEXER_H
