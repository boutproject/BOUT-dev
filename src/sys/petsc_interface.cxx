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

IndexerPtr GlobalIndexer::getInstance(Mesh* localmesh) {
  if (localmesh == bout::globals::mesh) {
    if (!initialisedGlobal) {
      globalInstance = std::shared_ptr<GlobalIndexer>(new GlobalIndexer(localmesh));
      globalInstance->initialise();
      initialisedGlobal = true;
    }
    return globalInstance;
  } else {
    IndexerPtr indexer(new GlobalIndexer(localmesh));
    indexer->initialise();
    return indexer;
  }
}

void GlobalIndexer::initialise() {
  fieldmesh->communicate(indices3D, indices2D);
  fieldmesh->communicate(indicesPerp);
  // Communicate a second time to get any corner values
  fieldmesh->communicate(indices3D, indices2D);
  fieldmesh->communicate(indicesPerp);
}

PetscInt GlobalIndexer::getGlobal(Ind2D ind) {
  return static_cast<PetscInt>(std::round(indices2D[ind]));
}

PetscInt GlobalIndexer::getGlobal(Ind3D ind) {
  return static_cast<PetscInt>(std::round(indices3D[ind]));
}

PetscInt GlobalIndexer::getGlobal(IndPerp ind) {
  return static_cast<PetscInt>(std::round(indicesPerp[ind]));
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
    : fieldmesh(localmesh), indices3D(-1., localmesh), indices2D(-1., localmesh),
      indicesPerp(-1., localmesh) {
  // Set up the 3D indices
  if (!localmesh->hasRegion3D("RGN_ALL_THIN")) {
    Region<Ind3D> bndry3d = localmesh->getRegion3D("RGN_LOWER_Y_THIN")
                            + localmesh->getRegion3D("RGN_UPPER_Y_THIN")
                            + localmesh->getRegion3D("RGN_INNER_X_THIN")
                            + localmesh->getRegion3D("RGN_NOBNDRY")
                            + localmesh->getRegion3D("RGN_OUTER_X_THIN");
    bndry3d.unique();
    localmesh->addRegion3D("RGN_ALL_THIN", bndry3d);
  }
  int counter = localmesh->globalStartIndex3D();
  BOUT_FOR_SERIAL(i, localmesh->getRegion3D("RGN_ALL_THIN")) { indices3D[i] = counter++; }

  // Set up the 2D indices
  if (!localmesh->hasRegion2D("RGN_ALL_THIN")) {
    Region<Ind2D> bndry2d = localmesh->getRegion2D("RGN_LOWER_Y_THIN")
                            + localmesh->getRegion2D("RGN_UPPER_Y_THIN")
                            + localmesh->getRegion2D("RGN_INNER_X_THIN")
                            + localmesh->getRegion2D("RGN_NOBNDRY")
                            + localmesh->getRegion2D("RGN_OUTER_X_THIN");
    bndry2d.unique();
    localmesh->addRegion2D("RGN_ALL_THIN", bndry2d);
  }
  counter = localmesh->globalStartIndex2D();
  BOUT_FOR_SERIAL(i, localmesh->getRegion2D("RGN_ALL_THIN")) { indices2D[i] = counter++; }

  // Set up the Perp indices; will these work in general or will
  // different ones be needed for each value of y?
  if (!localmesh->hasRegionPerp("RGN_ALL_THIN")) {
    Region<IndPerp> bndryPerp = localmesh->getRegionPerp("RGN_INNER_X_THIN")
                                + localmesh->getRegionPerp("RGN_NOBNDRY")
                                + localmesh->getRegionPerp("RGN_OUTER_X_THIN");
    bndryPerp.unique();
    localmesh->addRegionPerp("RGN_ALL_THIN", bndryPerp);
  }
  counter = localmesh->globalStartIndexPerp();
  BOUT_FOR_SERIAL(i, localmesh->getRegionPerp("RGN_ALL_THIN")) {
    indicesPerp[i] = counter++;
  }
}

void GlobalIndexer::recreateGlobalInstance() { initialisedGlobal = false; }

#endif // BOUT_HAS_PETSC
