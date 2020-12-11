/**************************************************************************
 * Create initial profiles for testing 3d inversion of perpendicular Laplacian
 *
 **************************************************************************
 * Copyright 2019 J.T. Omotani, C. MacMackin
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

#include "bout/physicsmodel.hxx"
#include "initialprofiles.hxx"

class CreateInitialProfiles : public PhysicsModel {
  int init(bool) {
    auto& opt = Options::root();
    parallel = opt["parallel"].withDefault(false);
    perpendicular = opt["perpendicular"].withDefault(false);

    SOLVE_FOR(initial);

    Field3D input1, input2, input3, input4, input5, input6;
    initial_profile("input1", input1);
    initial_profile("input2", input2);
    initial_profile("input3", input3);
    initial_profile("input4", input4);
    initial_profile("input5", input5);
    initial_profile("input6", input6);

    int jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2;
    mesh->get(jyseps1_1, "jyseps1_1");
    mesh->get(jyseps2_1, "jyseps2_1");
    mesh->get(jyseps1_2, "jyseps1_2");
    mesh->get(jyseps2_2, "jyseps2_2");

    // outboard midplane
    int yind = mesh->getLocalYIndexNoBoundaries((jyseps1_2 + jyseps2_2)/2);
    if (yind >= mesh->ystart and yind <= mesh->yend) {
      initial = sliceXZ(input1, yind);
    }

    // inboard midplane
    yind = mesh->getLocalYIndexNoBoundaries((jyseps1_1 + jyseps2_1)/2);
    if (yind >= mesh->ystart and yind <= mesh->yend) {
      initial = sliceXZ(input2, yind);
    }

    if (mesh->firstY()) {
      // lower, inner divertor
      for (auto it = mesh->iterateBndryLowerY(); not it.isDone(); it.next()) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          initial(it.ind, mesh->ystart, z) = input3(it.ind, mesh->ystart, z);
        }
      }
    } else {
      // upper, outer divertor, if there is one
      for (auto it = mesh->iterateBndryLowerY(); not it.isDone(); it.next()) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          initial(it.ind, mesh->ystart, z) = input5(it.ind, mesh->ystart, z);
        }
      }
    }

    if (mesh->lastY()) {
      // lower, outer divertor
      for (auto it = mesh->iterateBndryUpperY(); not it.isDone(); it.next()) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          initial(it.ind, mesh->yend, z) = input4(it.ind, mesh->yend, z);
        }
      }
    } else {
      // upper, inner divertor, if there is one
      for (auto it = mesh->iterateBndryUpperY(); not it.isDone(); it.next()) {
        for (int z = mesh->zstart; z <= mesh->zend; z++) {
          initial(it.ind, mesh->yend, z) = input6(it.ind, mesh->yend, z);
        }
      }
    }

    return 0;
  }

  int rhs(BoutReal) {
    mesh->communicate(initial);
    initial.applyParallelBoundary();

    ddt(initial) = 0.;
    if (parallel) {
      ddt(initial) += Grad2_par2(initial);
    }
    if (perpendicular) {
      ddt(initial) += Delp2(initial, CELL_DEFAULT, false);
    }

    return 0;
  }

  Field3D initial;
  bool parallel, perpendicular;
};

BOUTMAIN(CreateInitialProfiles);
