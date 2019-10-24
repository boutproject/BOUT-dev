/*!***********************************************************************
 * Classes to wrap PETSc matrices and vectors, providing a convenient
 * interface to them. In particular, they will internally convert
 * between BOUT++ indices and PETSc ones, making it far easier to set
 * up a linear system.
 *
 **************************************************************************
 * Copyright 2013 J. Buchanan, J.Omotani
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/
#include "mpi.h"
#include <memory>
#include <vector>

#include <bout/mesh.hxx>
#include <bout/paralleltransform.hxx>
#include <bout/petsc_interface.hxx>
#include <bout/petsclib.hxx>
#include <bout/region.hxx>
#include <bout_types.hxx>
#include <boutcomm.hxx>

#ifdef BOUT_HAS_PETSC

// GlobalIndexer implementation

bool GlobalIndexer::initialisedGlobal = false;
IndexerPtr GlobalIndexer::globalInstance;
Mesh* GlobalIndexer::globalmesh = nullptr;

IndexerPtr GlobalIndexer::getInstance(Mesh* localmesh) {
  // Check that the identity of bout::globals::mesh hasn't changed
  // since last call. (Needed for unit tests)
  if (bout::globals::mesh != globalmesh) {
    globalmesh = bout::globals::mesh;
    initialisedGlobal = false;
  }
  if (localmesh == globalmesh) {
    if (!initialisedGlobal) {
      globalInstance = std::shared_ptr<GlobalIndexer>(new GlobalIndexer(localmesh));
      initialisedGlobal = true;
    }
    return globalInstance;
  } else {
    return std::shared_ptr<GlobalIndexer>(new GlobalIndexer(localmesh));
  }
}

void GlobalIndexer::initialiseTest() {
  registerFieldForTest(indices3D);
  registerFieldForTest(indices2D);
  registerFieldForTest(indicesPerp);
}

void GlobalIndexer::initialise() {
  fieldmesh->communicate(indices3D, indices2D);
  fieldmesh->communicate(indicesPerp);
  // Communicate a second time to get any corner values
  fieldmesh->communicate(indices3D, indices2D);
  fieldmesh->communicate(indicesPerp);
}

Mesh* GlobalIndexer::getMesh() { return fieldmesh; }

PetscInt GlobalIndexer::getGlobal(Ind2D ind) {
  return static_cast<PetscInt>(indices2D[ind] + 0.5);
}

PetscInt GlobalIndexer::getGlobal(Ind3D ind) {
  return static_cast<PetscInt>(indices3D[ind] + 0.5);
}

PetscInt GlobalIndexer::getGlobal(IndPerp ind) {
  return static_cast<PetscInt>(indicesPerp[ind] + 0.5);
}

void GlobalIndexer::registerFieldForTest(FieldData& UNUSED(f)) {
  // This is a place-holder which does nothing. It can be overridden
  // by descendent classes if necessary to set up testing.
  return;
}

void GlobalIndexer::registerFieldForTest(FieldPerp& UNUSED(f)) {
  // This is a place-holder which does nothing. It can be overridden
  // by descendent classes if necessary to set up testing.
  return;
}

GlobalIndexer::GlobalIndexer(Mesh* localmesh)
    : fieldmesh(localmesh), indices3D(-2., localmesh), indices2D(-2., localmesh),
      indicesPerp(-2., localmesh) {
  // Set up the 3D indices
  int counter = localmesh->globalStartIndex3D();
  for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    for (int z = 0; z < localmesh->LocalNz; z++) {
      if (it.ind == localmesh->xstart)
        indices3D(it.ind - 1, localmesh->ystart - 1, z) = counter++;
      if (it.ind == localmesh->xend)
        indices3D(it.ind + 1, localmesh->ystart - 1, z) = counter++;
      indices3D(it.ind, localmesh->ystart - 1, z) = counter++;
    }
  }
  for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    for (int z = 0; z < localmesh->LocalNz; z++) {
      if (it.ind == localmesh->xstart)
        indices3D(it.ind - 1, localmesh->yend + 1, z) = counter++;
      if (it.ind == localmesh->xend)
        indices3D(it.ind + 1, localmesh->yend + 1, z) = counter++;
      indices3D(it.ind, localmesh->yend + 1, z) = counter++;
    }
  }
  BOUT_FOR_SERIAL(i, localmesh->getRegion3D("RGN_NOY")) {
    if ((i.x() >= localmesh->xstart && i.x() <= localmesh->xend)
        || (i.x() == localmesh->xstart - 1 && localmesh->firstX())
        || (i.x() == localmesh->xend + 1 && localmesh->lastX())) {
      indices3D[i] = counter++;
    }
  }

  // Set up the 2D indices
  counter = localmesh->globalStartIndex2D();
  for (RangeIterator it = localmesh->iterateBndryLowerY(); !it.isDone(); it++) {
    if (it.ind == localmesh->xstart)
      indices2D(it.ind - 1, localmesh->ystart - 1) = counter++;
    if (it.ind == localmesh->xend)
      indices2D(it.ind + 1, localmesh->ystart - 1) = counter++;
    indices2D(it.ind, localmesh->ystart - 1) = counter++;
  }
  for (RangeIterator it = localmesh->iterateBndryUpperY(); !it.isDone(); it++) {
    if (it.ind == localmesh->xstart)
      indices2D(it.ind - 1, localmesh->yend + 1) = counter++;
    if (it.ind == localmesh->xend)
      indices2D(it.ind + 1, localmesh->ystart + 1) = counter++;
    indices2D(it.ind, localmesh->yend + 1) = counter++;
  }
  BOUT_FOR_SERIAL(i, localmesh->getRegion2D("RGN_NOY")) {
    if ((i.x() >= localmesh->xstart && i.x() <= localmesh->xend)
        || (i.x() == localmesh->xstart - 1 && localmesh->firstX())
        || (i.x() == localmesh->xend + 1 && localmesh->lastX())) {
      indices2D[i] = counter++;
    }
  }

  // Set up the Perp indices; will these work in general or will
  // different ones be needed for each value of y?
  counter = localmesh->globalStartIndexPerp();
  BOUT_FOR_SERIAL(i, localmesh->getRegionPerp("RGN_NOY")) {
    if ((i.x() >= localmesh->xstart && i.x() <= localmesh->xend)
        || (i.x() == localmesh->xstart - 1 && localmesh->firstX())
        || (i.x() == localmesh->xend + 1 && localmesh->lastX())) {
      indicesPerp[i] = counter++;
    }
  }
}

#endif // BOUT_HAS_PETSC
