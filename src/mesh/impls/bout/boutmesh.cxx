/**************************************************************************
 * Implementation of the Mesh class, handling input files compatible with
 * BOUT / BOUT-06.
 *
 * Changelog
 * ---------
 *
 * 2015-01 Ben Dudson <benjamin.dudson@york.ac.uk>
 *      *
 *
 * 2010-05 Ben Dudson <bd512@york.ac.uk>
 *      * Initial version, adapted from grid.cpp and topology.cpp
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
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

#include "boutmesh.hxx"

#include <bout/constants.hxx>
#include <bout/sys/timer.hxx>
#include <boutcomm.hxx>
#include <boutexception.hxx>
#include <dcomplex.hxx>
#include <derivs.hxx>
#include <fft.hxx>
#include <msg_stack.hxx>
#include <options.hxx>
#include <output.hxx>
#include <utils.hxx>

/// MPI type of BoutReal for communications
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

BoutMesh::BoutMesh(GridDataSource *s, Options *opt) : Mesh(s, opt) {
  OPTION(options, symmetricGlobalX, true);
  if (!options->isSet("symmetricGlobalY")) {
    std::string optionfile;
    OPTION(Options::getRoot(), optionfile, "");
    output_warn << "WARNING: The default of this option has changed in release 4.1.\n\
If you want the old setting, you have to specify mesh:symmetricGlobalY=false in "
                << optionfile << "\n";
  }
  OPTION(options, symmetricGlobalY, true);

  comm_x = MPI_COMM_NULL;
  comm_inner = MPI_COMM_NULL;
  comm_middle = MPI_COMM_NULL;
  comm_outer = MPI_COMM_NULL;
}

BoutMesh::~BoutMesh() {
  // Delete the communication handles
  clear_handles();

  // Delete the boundary regions
  for (const auto &bndry : boundary)
    delete bndry;
  for (const auto &bndry : par_boundary)
    delete bndry;

  if (comm_x != MPI_COMM_NULL)
    MPI_Comm_free(&comm_x);
  if (comm_inner != MPI_COMM_NULL)
    MPI_Comm_free(&comm_inner);
  if (comm_outer != MPI_COMM_NULL)
    MPI_Comm_free(&comm_outer);
}

int BoutMesh::load() {
  TRACE("BoutMesh::load()");

  output_progress << "Loading mesh" << endl;

  // Use root level options
  Options *options = Options::getRoot();

  //////////////
  // Number of processors

  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);

  //////////////
  // Grid sizes

  if (Mesh::get(nx, "nx"))
    throw BoutException("Mesh must contain nx");

  if (Mesh::get(ny, "ny"))
    throw BoutException("Mesh must contain ny");

  int MZ;

  if (Mesh::get(MZ, "nz")) {
    // No "nz" variable in the grid file. Instead read MZ from options

    OPTION(options, MZ, 64);
    if (!is_pow2(MZ)) {
      // Should be a power of 2 for efficient FFTs
      output_warn.write("WARNING: Number of toroidal points should be 2^n\n", MZ);
    }
  } else {
    output_info.write("\tRead nz from input grid file\n");
  }

  output_info << "\tGrid size: " << nx << " x " << ny << " x " << MZ << endl;

  // Get guard cell sizes
  // Try to read from grid file first, then if not found
  // get from options
  if (Mesh::get(MXG, "MXG")) {
    // Error code returned
    options->get("MXG", MXG, 2);
  }
  ASSERT0(MXG >= 0);

  if (Mesh::get(MYG, "MYG")) {
    options->get("MYG", MYG, 2);
  }
  ASSERT0(MYG >= 0);

  output_info << "\tGuard cells (x,y): " << MXG << ", " << MYG << std::endl;

  // Check that nx is large enough
  if (nx <= 2 * MXG) {
    throw BoutException("Error: nx must be greater than 2 times MXG (2 * %d)", MXG);
  }

  // Set global grid sizes
  GlobalNx = nx;
  GlobalNy = ny + 2 * MYG;
  GlobalNz = MZ;

  if (2 * MXG >= nx)
    throw BoutException("nx must be greater than 2*MXG");

  // separatrix location
  if (Mesh::get(ixseps1, "ixseps1")) {
    ixseps1 = GlobalNx;
    output_warn.write(
        "\tWARNING: Separatrix location 'ixseps1' not found. Setting to %d\n", ixseps1);
  }
  if (Mesh::get(ixseps2, "ixseps2")) {
    ixseps2 = GlobalNx;
    output_warn.write(
        "\tWARNING: Separatrix location 'ixseps2' not found. Setting to %d\n", ixseps2);
  }
  if (Mesh::get(jyseps1_1, "jyseps1_1")) {
    jyseps1_1 = -1;
    output_warn.write("\tWARNING: Branch-cut 'jyseps1_1' not found. Setting to %d\n",
                      jyseps1_1);
  }
  if (Mesh::get(jyseps1_2, "jyseps1_2")) {
    jyseps1_2 = ny / 2;
    output_warn.write("\tWARNING: Branch-cut 'jyseps1_2' not found. Setting to %d\n",
                      jyseps1_2);
  }
  if (Mesh::get(jyseps2_1, "jyseps2_1")) {
    jyseps2_1 = jyseps1_2;
    output_warn.write("\tWARNING: Branch-cut 'jyseps2_1' not found. Setting to %d\n",
                      jyseps2_1);
  }
  if (Mesh::get(jyseps2_2, "jyseps2_2")) {
    jyseps2_2 = ny - 1;
    output_warn.write("\tWARNING: Branch-cut 'jyseps2_2' not found. Setting to %d\n",
                      jyseps2_2);
  }

  if (Mesh::get(ny_inner, "ny_inner")) {
    ny_inner = jyseps2_1;
    output_warn.write(
        "\tWARNING: Number of inner y points 'ny_inner' not found. Setting to %d\n",
        ny_inner);
  }

  /// Check inputs
  if (jyseps1_1 < -1) {
    output_warn.write("\tWARNING: jyseps1_1 (%d) must be >= -1. Setting to -1\n",
                      jyseps1_1);
    jyseps1_1 = -1;
  }

  if (jyseps2_1 <= jyseps1_1) {
    output_warn.write(
        "\tWARNING: jyseps2_1 (%d) must be > jyseps1_1 (%d). Setting to %d\n", jyseps2_1,
        jyseps1_1, jyseps1_1 + 1);
    jyseps2_1 = jyseps1_1 + 1;
  }
  if (jyseps1_2 < jyseps2_1) {
    output_warn.write(
        "\tWARNING: jyseps1_2 (%d) must be >= jyseps2_1 (%d). Setting to %d\n", jyseps1_2,
        jyseps2_1, jyseps2_1);
    jyseps1_2 = jyseps2_1;
  }
  if (jyseps2_2 >= ny) {
    output_warn.write("\tWARNING: jyseps2_2 (%d) must be < ny (%d). Setting to %d\n",
                      jyseps2_2, ny, ny - 1);
    jyseps2_2 = ny - 1;
  }
  if (jyseps2_2 < jyseps1_2) {
    if (jyseps1_2 >= ny) {
      throw BoutException("jyseps1_2 (%d) must be < ny (%d).", jyseps1_2, ny);
    }
    output_warn.write("\tWARNING: jyseps2_2 (%d) must be >= jyseps1_2 (%d). Setting to %d\n",
                      jyseps2_2, jyseps1_2, jyseps1_2);
    jyseps2_2 = jyseps1_2;
  }

  if (options->isSet("NXPE")) {    // Specified NXPE
    options->get("NXPE", NXPE, 1); // Decomposition in the radial direction
    if ((NPES % NXPE) != 0) {
      throw BoutException(
          "Number of processors (%d) not divisible by NPs in x direction (%d)\n", NPES,
          NXPE);
    }

    NYPE = NPES / NXPE;

    int nyp = NPES / NXPE;
    int ysub = ny / NYPE;

    // Check size of Y mesh
    if (ysub < MYG) {
      throw BoutException("\t -> ny/NYPE (%d/%d = %d) must be >= MYG (%d)\n", ny, nyp,
                        ysub, MYG);
    }
    // Check branch cuts
    if ((jyseps1_1 + 1) % ysub != 0) {
      throw BoutException(
          "\t -> Leg region jyseps1_1+1 (%d) must be a multiple of MYSUB (%d)\n",
          jyseps1_1 + 1, ysub);
    }

    if (jyseps2_1 != jyseps1_2) {
      // Double Null

      if ((jyseps2_1 - jyseps1_1) % ysub != 0) {
        throw BoutException("\t -> Core region jyseps2_1-jyseps1_1 (%d-%d = %d) must "
                          "be a multiple of MYSUB (%d)\n",
                          jyseps2_1, jyseps1_1, jyseps2_1 - jyseps1_1, ysub);
      }

      if ((jyseps2_2 - jyseps1_2) % ysub != 0) {
        throw BoutException("\t -> Core region jyseps2_2-jyseps1_2 (%d-%d = %d) must "
                          "be a multiple of MYSUB (%d)\n",
                          jyseps2_2, jyseps1_2, jyseps2_2 - jyseps1_2, ysub);
      }

      // Check upper legs
      if ((ny_inner - jyseps2_1 - 1) % ysub != 0) {
        throw BoutException("\t -> leg region ny_inner-jyseps2_1-1 (%d-%d-1 = %d) must "
                          "be a multiple of MYSUB (%d)\n",
                          ny_inner, jyseps2_1, ny_inner - jyseps2_1 - 1, ysub);
      }
      if ((jyseps1_2 - ny_inner + 1) % ysub != 0) {
        throw BoutException("\t -> leg region jyseps1_2-ny_inner+1 (%d-%d+1 = %d) must "
                          "be a multiple of MYSUB (%d)\n",
                          jyseps1_2, ny_inner, jyseps1_2 - ny_inner + 1, ysub);
      }
    } else {
      // Single Null
      if ((jyseps2_2 - jyseps1_1) % ysub != 0) {
        throw BoutException("\t -> Core region jyseps2_2-jyseps1_1 (%d-%d = %d) must "
                          "be a multiple of MYSUB (%d)\n",
                          jyseps2_2, jyseps1_1, jyseps2_2 - jyseps1_1, ysub);
      }
    }

    if ((ny - jyseps2_2 - 1) % ysub != 0) {
      throw BoutException("\t -> leg region ny-jyseps2_2-1 (%d-%d-1 = %d) must be a "
                        "multiple of MYSUB (%d)\n",
                        ny, jyseps2_2, ny - jyseps2_2 - 1, ysub);
    }
  } else {
    // Choose NXPE

    MX = nx - 2 * MXG;

    NXPE = -1; // Best option
    
    BoutReal ideal = sqrt(MX * NPES / static_cast<BoutReal>(ny)); // Results in square domains

    output_info.write("Finding value for NXPE (ideal = %f)\n", ideal);

    for (int i = 1; i <= NPES; i++) { // Loop over all possibilities
      if ((NPES % i == 0) &&          // Processors divide equally
          (MX % i == 0) &&            // Mesh in X divides equally
          (ny % (NPES / i) == 0)) {   // Mesh in Y divides equally

        output_info.write("\tCandidate value: %d\n", i);

        int nyp = NPES / i;
        int ysub = ny / nyp;

        // Check size of Y mesh
        if (ysub < MYG) {
          output_info.write("\t -> ny/NYPE (%d/%d = %d) must be >= MYG (%d)\n", ny, nyp,
                            ysub, MYG);
          continue;
        }
        // Check branch cuts
        if ((jyseps1_1 + 1) % ysub != 0) {
          output_info.write(
              "\t -> Leg region jyseps1_1+1 (%d) must be a multiple of MYSUB (%d)\n",
              jyseps1_1 + 1, ysub);
          continue;
        }

        if (jyseps2_1 != jyseps1_2) {
          // Double Null

          if ((jyseps2_1 - jyseps1_1) % ysub != 0) {
            output_info.write("\t -> Core region jyseps2_1-jyseps1_1 (%d-%d = %d) must "
                              "be a multiple of MYSUB (%d)\n",
                              jyseps2_1, jyseps1_1, jyseps2_1 - jyseps1_1, ysub);
            continue;
          }

          if ((jyseps2_2 - jyseps1_2) % ysub != 0) {
            output_info.write("\t -> Core region jyseps2_2-jyseps1_2 (%d-%d = %d) must "
                              "be a multiple of MYSUB (%d)\n",
                              jyseps2_2, jyseps1_2, jyseps2_2 - jyseps1_2, ysub);
            continue;
          }

          // Check upper legs
          if ((ny_inner - jyseps2_1 - 1) % ysub != 0) {
            output_info.write("\t -> leg region ny_inner-jyseps2_1-1 (%d-%d-1 = %d) must "
                              "be a multiple of MYSUB (%d)\n",
                              ny_inner, jyseps2_1, ny_inner - jyseps2_1 - 1, ysub);
            continue;
          }
          if ((jyseps1_2 - ny_inner + 1) % ysub != 0) {
            output_info.write("\t -> leg region jyseps1_2-ny_inner+1 (%d-%d+1 = %d) must "
                              "be a multiple of MYSUB (%d)\n",
                              jyseps1_2, ny_inner, jyseps1_2 - ny_inner + 1, ysub);
            continue;
          }
        } else {
          // Single Null
          if ((jyseps2_2 - jyseps1_1) % ysub != 0) {
            output_info.write("\t -> Core region jyseps2_2-jyseps1_1 (%d-%d = %d) must "
                              "be a multiple of MYSUB (%d)\n",
                              jyseps2_2, jyseps1_1, jyseps2_2 - jyseps1_1, ysub);
            continue;
          }
        }

        if ((ny - jyseps2_2 - 1) % ysub != 0) {
          output_info.write("\t -> leg region ny-jyseps2_2-1 (%d-%d-1 = %d) must be a "
                            "multiple of MYSUB (%d)\n",
                            ny, jyseps2_2, ny - jyseps2_2 - 1, ysub);
          continue;
        }
        output_info.write("\t -> Good value\n");
        // Found an acceptable value
        if ((NXPE < 1) || (fabs(ideal - i) < fabs(ideal - NXPE)))
          NXPE = i; // Keep value nearest to the ideal
      }
    }

    if (NXPE < 1)
      throw BoutException(
          "Could not find a valid value for NXPE. Try a different number of processors.");

    NYPE = NPES / NXPE;

    output_progress.write(
        "\tDomain split (NXPE=%d, NYPE=%d) into domains (localNx=%d, localNy=%d)\n", NXPE,
        NYPE, MX / NXPE, ny / NYPE);
  }

  /// Get X and Y processor indices
  PE_YIND = MYPE / NXPE;
  PE_XIND = MYPE % NXPE;

  // Work out other grid size quantities

  /// MXG at each end needed for edge boundary regions
  MX = nx - 2 * MXG;

  /// Split MX points between NXPE processors
  MXSUB = MX / NXPE;
  if ((MX % NXPE) != 0) {
    throw BoutException("Cannot split %d X points equally between %d processors\n", MX,
                        NXPE);
  }

  /// NOTE: No grid data reserved for Y boundary cells - copy from neighbours
  MY = ny;
  MYSUB = MY / NYPE;
  if ((MY % NYPE) != 0) {
    throw BoutException(
        "\tERROR: Cannot split %d Y points equally between %d processors\n", MY, NYPE);
  }

  /// Get mesh options
  OPTION(options, IncIntShear, false);
  OPTION(options, periodicX, false); // Periodic in X

  OPTION(options, async_send, false); // Whether to use asyncronous sends

  // Set global offsets

  OffsetX = PE_XIND * MXSUB;
  OffsetY = PE_YIND * MYSUB;
  OffsetZ = 0;

  if (options->isSet("zperiod")) {
    OPTION(options, zperiod, 1);
    ZMIN = 0.0;
    ZMAX = 1.0 / static_cast<BoutReal>(zperiod);
  } else {
    OPTION(options, ZMIN, 0.0);
    OPTION(options, ZMAX, 1.0);

    zperiod = ROUND(1.0 / (ZMAX - ZMIN));
  }

  /// Number of grid cells is ng* = M*SUB + guard/boundary cells
  LocalNx = MXSUB + 2 * MXG;
  LocalNy = MYSUB + 2 * MYG;
  LocalNz = MZ;

  // Set local index ranges

  xstart = MXG;
  xend = MXG + MXSUB - 1;

  ystart = MYG;
  yend = MYG + MYSUB - 1;

  ///////////////////// TOPOLOGY //////////////////////////
  /// Call topology to set layout of grid
  topology();

  OPTION(options, TwistShift, false);

  if (TwistShift) {
    output_info.write("Applying Twist-Shift condition. Interpolation: FFT\n");

    // Try to read the shift angle from the grid file
    // NOTE: All processors should know the twist-shift angle (for invert_parderiv)

    ShiftAngle.resize(LocalNx);

    if (!source->get(this, ShiftAngle, "ShiftAngle", LocalNx, XGLOBAL(0))) {
      throw BoutException("ERROR: Twist-shift angle 'ShiftAngle' not found.");
    }
  }

  //////////////////////////////////////////////////////
  /// Communicator

  MPI_Group group_world;
  MPI_Comm_group(BoutComm::get(), &group_world); // Get the entire group

  //////////////////////////////////////////////////////
  /// Communicator in X

  MPI_Group group;
  MPI_Comm comm_tmp;

  int proc[3]; // Processor range

  for (int yp = 0; yp < NYPE; yp++) {
    proc[0] = PROC_NUM(0, yp);        // First
    proc[1] = PROC_NUM(NXPE - 1, yp); // Last
    proc[2] = 1;                      // stride

    output_debug << "XCOMM " << proc[0] << ", " << proc[1] << endl;

    if (MPI_Group_range_incl(group_world, 1, &proc, &group) != MPI_SUCCESS)
      throw BoutException(
          "Could not create X communication group for yp=%d (xind=%d,yind=%d)\n", yp,
          PE_XIND, PE_YIND);
    if (MPI_Comm_create(BoutComm::get(), group, &comm_tmp) != MPI_SUCCESS)
      throw BoutException("Could not create X communicator for yp=%d (xind=%d,yind=%d)\n",
                          yp, PE_XIND, PE_YIND);
    MPI_Group_free(&group);

    if (yp == PE_YIND) {
      // Should be in this group
      if (comm_tmp == MPI_COMM_NULL)
        throw BoutException("X communicator null");

      comm_x = comm_tmp;
    } else {
      if (comm_tmp != MPI_COMM_NULL)
        throw BoutException("X communicator should be null");
    }
  }

  //////////////////////////////////////////////////////
  /// Communicators for Y gather/scatter

  MPI_Group group_tmp1, group_tmp2;

  proc[2] = NXPE; // Stride in processor rank

  // Outer SOL regions
  if (jyseps1_2 == jyseps2_1) {
    // Single-null. All processors with same PE_XIND
    TRACE("Creating Outer SOL communicators for Single Null operation");

    for (int i = 0; i < NXPE; i++) {
      proc[0] = PROC_NUM(i, 0);
      proc[1] = PROC_NUM(i, NYPE - 1);

      output_debug << "Outer SOL " << proc[0] << ", " << proc[1] << endl;

      MPI_Group_range_incl(group_world, 1, &proc, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if (i == PE_XIND) {
        // Should be part of this communicator
        if (comm_tmp == MPI_COMM_NULL) {
          // error
          throw BoutException("Single null outer SOL not correct\n");
        }
        comm_outer = comm_tmp;
      } else if (comm_tmp != MPI_COMM_NULL) {
        // Not part of this communicator so should be NULL
        throw BoutException("Single null outer SOL not correct\n");
      }
      MPI_Group_free(&group);
    }
  } else {
    // Double null
    TRACE("Creating Outer SOL communicators for Double Null operation");

    for (int i = 0; i < NXPE; i++) {
      // Inner SOL
      proc[0] = PROC_NUM(i, 0);
      proc[1] = PROC_NUM(i, YPROC(ny_inner - 1));

      output_debug << "Double Null inner SOL " << proc[0] << ", " << proc[1] << endl;

      if (MPI_Group_range_incl(group_world, 1, &proc, &group) != MPI_SUCCESS)
        throw BoutException("MPI_Group_range_incl failed for xp = %d", NXPE);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if (comm_tmp != MPI_COMM_NULL)
        comm_outer = comm_tmp;
      MPI_Group_free(&group);

      // Outer SOL
      proc[0] = PROC_NUM(i, YPROC(ny_inner));
      proc[1] = PROC_NUM(i, NYPE - 1);

      output_debug << "Double Null outer SOL " << proc[0] << ", " << proc[1] << endl;

      MPI_Group_range_incl(group_world, 1, &proc, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if (comm_tmp != MPI_COMM_NULL) {
        comm_outer = comm_tmp;
      }
      MPI_Group_free(&group);
    }
  }

  for (int i = 0; i < NXPE; i++) {
    // Lower PF region

    if ((jyseps1_1 >= 0) || (jyseps2_2 + 1 < ny)) {
      // A lower PF region exists
      TRACE("Creating lower PF communicators for xp=%d", i);

      output_debug << "Creating lower PF communicators for xp = " << i << endl;

      if (jyseps1_1 >= 0) {
        proc[0] = PROC_NUM(i, 0);
        proc[1] = PROC_NUM(i, YPROC(jyseps1_1));

        output_debug << "PF1 " << proc[0] << ", " << proc[1] << endl;

        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      } else
        group_tmp1 = MPI_GROUP_EMPTY;

      if (jyseps2_2 + 1 < ny) {
        proc[0] = PROC_NUM(i, YPROC(jyseps2_2 + 1));
        proc[1] = PROC_NUM(i, NYPE - 1);

        output_debug << "PF2 " << proc[0] << ", " << proc[1] << endl;

        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      } else {
        group_tmp2 = MPI_GROUP_EMPTY;
      }

      MPI_Group_union(group_tmp1, group_tmp2, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if (comm_tmp != MPI_COMM_NULL) {
        comm_inner = comm_tmp;
        if (ixseps_lower == ixseps_outer) {
          // Between the separatrices is still in the PF region

          output_debug << "-> Inner and middle\n";

          comm_middle = comm_inner;
        } else {

          output_debug << "-> Outer and middle\n";

          comm_middle = comm_outer;
        }
      }

      output_debug << "Freeing\n";

      MPI_Group_free(&group);
      if (group_tmp1 != MPI_GROUP_EMPTY) {
        MPI_Group_free(&group_tmp1);
      }
      if (group_tmp2 != MPI_GROUP_EMPTY) {
        MPI_Group_free(&group_tmp2);
      }

      output_debug << "done lower PF\n";
    }

    if (jyseps2_1 != jyseps1_2) {
      // Upper PF region
      // Note need to order processors so that a continuous surface is formed
      TRACE("Creating upper PF communicators for xp=%d", i);

      output_debug << "Creating upper PF communicators for xp = " << i << endl;

      proc[0] = PROC_NUM(i, YPROC(ny_inner));
      proc[1] = PROC_NUM(i, YPROC(jyseps1_2));

      output_debug << "PF3 " << proc[0] << ", " << proc[1] << endl;

      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      proc[0] = PROC_NUM(i, YPROC(jyseps2_1 + 1));
      proc[1] = PROC_NUM(i, YPROC(ny_inner - 1));

      output_debug << "PF4 " << proc[0] << ", " << proc[1] << endl;

      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      MPI_Group_union(group_tmp1, group_tmp2, &group);
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if (comm_tmp != MPI_COMM_NULL) {
        comm_inner = comm_tmp;
        if (ixseps_upper == ixseps_outer) {

          output_debug << "-> Inner and middle\n";

          comm_middle = comm_inner;
        } else {

          output_debug << "-> Outer and middle\n";

          comm_middle = comm_outer;
          // MPI_Comm_dup(comm_outer, &comm_middle);
        }
      }

      output_debug << "Freeing\n";

      MPI_Group_free(&group);
      if (group_tmp1 != MPI_GROUP_EMPTY)
        MPI_Group_free(&group_tmp1);
      if (group_tmp2 != MPI_GROUP_EMPTY)
        MPI_Group_free(&group_tmp2);

      output_debug << "done upper PF\n";
    }

    // Core region
    TRACE("Creating core communicators");
    proc[0] = PROC_NUM(i, YPROC(jyseps1_1 + 1));
    proc[1] = PROC_NUM(i, YPROC(jyseps2_1));

    output_debug << "CORE1 " << proc[0] << ", " << proc[1] << endl;

    if ((proc[0] < 0) || (proc[1] < 0))
      throw BoutException("Invalid processor range for core processors");
    MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);

    proc[0] = PROC_NUM(i, YPROC(jyseps1_2 + 1));
    proc[1] = PROC_NUM(i, YPROC(jyseps2_2));

    output_debug << "CORE2 " << proc[0] << ", " << proc[1] << endl;

    if ((proc[0] < 0) || (proc[1] < 0)) {
      group_tmp2 = MPI_GROUP_EMPTY;
    } else {
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
    }

    MPI_Group_union(group_tmp1, group_tmp2, &group);
    MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
    if (comm_tmp != MPI_COMM_NULL) {
      comm_inner = comm_tmp;

      if (ixseps_inner == ixseps_outer)
        MPI_Comm_dup(comm_inner, &comm_middle);
    }

    if (group_tmp1 != MPI_GROUP_EMPTY)
      MPI_Group_free(&group_tmp1);
    if (group_tmp2 != MPI_GROUP_EMPTY)
      MPI_Group_free(&group_tmp2);
    MPI_Group_free(&group);
  }

  if (ixseps_inner == ixseps_outer) {
    // Balanced null, so no middle
    MPI_Comm_dup(comm_inner, &comm_middle);
  } else {
    // Need to handle unbalanced double-null case

    output_debug << "Unbalanced " << endl;

    if (ixseps_upper > ixseps_lower) {
      // middle is connected to the bottom
      TRACE("Creating unbalanced lower communicators");

      for (int i = 0; i < NXPE; i++) {
        proc[0] = PROC_NUM(i, 0);
        proc[1] = PROC_NUM(i, YPROC(jyseps2_1));
        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
        proc[0] = PROC_NUM(i, YPROC(jyseps1_2 + 1));
        proc[1] = PROC_NUM(i, NYPE - 1);
        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
        MPI_Group_union(group_tmp1, group_tmp2, &group);
        MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
        if (comm_tmp != MPI_COMM_NULL)
          comm_middle = comm_tmp;

        if (group_tmp1 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp1);
        if (group_tmp2 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp2);
        MPI_Group_free(&group);
      }
    } else {
      // middle is connected to the top
      TRACE("Creating unbalanced upper communicators");

      for (int i = 0; i < NXPE; i++) {
        proc[0] = PROC_NUM(i, YPROC(ny_inner));
        proc[1] = PROC_NUM(i, YPROC(jyseps2_2));
        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
        proc[0] = PROC_NUM(i, YPROC(jyseps1_1 + 1));
        proc[1] = PROC_NUM(i, YPROC(ny_inner - 1));
        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
        MPI_Group_union(group_tmp1, group_tmp2, &group);
        MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
        if (comm_tmp != MPI_COMM_NULL)
          comm_middle = comm_tmp;

        if (group_tmp1 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp1);
        if (group_tmp2 != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_tmp2);
        MPI_Group_free(&group);
      }
    }
  }
  MPI_Group_free(&group_world);
  // Now have communicators for all regions.

  output_debug << "Got communicators" << endl;

  //////////////////////////////////////////////////////
  // Boundary regions
  if (!periodicX && (MXG > 0)) {
    // Need boundaries in X if not periodic and have X guard cells
    if (PE_XIND == 0) {
      // Inner either core or PF

      int yg = YGLOBAL(MYG); // Get a global index in this processor

      if (((yg > jyseps1_1) && (yg <= jyseps2_1)) ||
          ((yg > jyseps1_2) && (yg <= jyseps2_2))) {
        // Core
        boundary.push_back(new BoundaryRegionXIn("core", ystart, yend, this));
      } else {
        // PF region
        boundary.push_back(new BoundaryRegionXIn("pf", ystart, yend, this));
      }
    }
    if (PE_XIND == (NXPE - 1)) {
      // Outer SOL
      boundary.push_back(new BoundaryRegionXOut("sol", ystart, yend, this));
    }
  }

  if (MYG > 0) {
    // Need boundaries in Y

    if ((UDATA_INDEST < 0) && (UDATA_XSPLIT > xstart))
      boundary.push_back(new BoundaryRegionYUp("upper_target", xstart, UDATA_XSPLIT - 1, this));
    if ((UDATA_OUTDEST < 0) && (UDATA_XSPLIT <= xend))
      boundary.push_back(new BoundaryRegionYUp("upper_target", UDATA_XSPLIT, xend, this));

    if ((DDATA_INDEST < 0) && (DDATA_XSPLIT > xstart))
      boundary.push_back(
          new BoundaryRegionYDown("lower_target", xstart, DDATA_XSPLIT - 1, this));
    if ((DDATA_OUTDEST < 0) && (DDATA_XSPLIT <= xend))
      boundary.push_back(new BoundaryRegionYDown("lower_target", DDATA_XSPLIT, xend, this));
  }

  if (!boundary.empty()) {
    output_info << "Boundary regions in this processor: ";
    for (const auto &bndry : boundary) {
      output_info << bndry->label << ", ";
    }
    output_info << endl;
  } else {
    output_info << "No boundary regions in this processor" << endl;
  }
  
  output_info << "Constructing default regions" << endl;
  createDefaultRegions();

  // Add boundary regions
  addBoundaryRegions();

  output_info.write("\tdone\n");

  return 0;
}

/****************************************************************
 *                 COMMUNICATIONS
 ****************************************************************/

const int IN_SENT_UP = 0;    ///< Data lower in X than branch-cut, at upper boundary in Y
const int OUT_SENT_UP = 1;   ///< Data higher in X than branch-cut, at upper boundary in Y
const int IN_SENT_DOWN = 2;  ///< Data lower in X than branch-cut, at lower boundary in Y
const int OUT_SENT_DOWN = 3; ///< Data higher in X than branch-cut, at lower boundary in Y
// X communication signals
const int IN_SENT_OUT = 4; ///< Data going in positive X direction (in to out)
const int OUT_SENT_IN = 5; ///< Data going in negative X direction (out to in)

void BoutMesh::post_receive(CommHandle &ch) {
  BoutReal *inbuff;
  int len;

  /// Post receive data from above (y+1)

  len = 0;
  if (UDATA_INDEST != -1) {
    len = msg_len(ch.var_list.get(), 0, UDATA_XSPLIT, 0, MYG);
    MPI_Irecv(std::begin(ch.umsg_recvbuff), len, PVEC_REAL_MPI_TYPE, UDATA_INDEST,
              IN_SENT_DOWN, BoutComm::get(), &ch.request[0]);
  }
  if (UDATA_OUTDEST != -1) {
    inbuff = &ch.umsg_recvbuff[len]; // pointer to second half of the buffer
    MPI_Irecv(inbuff, msg_len(ch.var_list.get(), UDATA_XSPLIT, LocalNx, 0, MYG),
              PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, OUT_SENT_DOWN, BoutComm::get(),
              &ch.request[1]);
  }

  /// Post receive data from below (y-1)

  len = 0;

  if (DDATA_INDEST != -1) { // If sending & recieving data from a processor
    len = msg_len(ch.var_list.get(), 0, DDATA_XSPLIT, 0, MYG);
    MPI_Irecv(std::begin(ch.dmsg_recvbuff), len, PVEC_REAL_MPI_TYPE, DDATA_INDEST,
              IN_SENT_UP, BoutComm::get(), &ch.request[2]);
  }
  if (DDATA_OUTDEST != -1) {
    inbuff = &ch.dmsg_recvbuff[len];
    MPI_Irecv(inbuff, msg_len(ch.var_list.get(), DDATA_XSPLIT, LocalNx, 0, MYG),
              PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, OUT_SENT_UP, BoutComm::get(),
              &ch.request[3]);
  }

  /// Post receive data from left (x-1)

  if (IDATA_DEST != -1) {
    MPI_Irecv(std::begin(ch.imsg_recvbuff), msg_len(ch.var_list.get(), 0, MXG, 0, MYSUB),
              PVEC_REAL_MPI_TYPE, IDATA_DEST, OUT_SENT_IN, BoutComm::get(),
              &ch.request[4]);
  }

  // Post receive data from right (x+1)

  if (ODATA_DEST != -1) {
    MPI_Irecv(std::begin(ch.omsg_recvbuff), msg_len(ch.var_list.get(), 0, MXG, 0, MYSUB),
              PVEC_REAL_MPI_TYPE, ODATA_DEST, IN_SENT_OUT, BoutComm::get(),
              &ch.request[5]);
  }
}

comm_handle BoutMesh::send(FieldGroup &g) {
  /// Start timer
  Timer timer("comms");

  /// Work out length of buffer needed
  int xlen = msg_len(g.get(), 0, MXG, 0, MYSUB);
  int ylen = msg_len(g.get(), 0, LocalNx, 0, MYG);

  /// Get a communications handle of (at least) the needed size
  CommHandle *ch = get_handle(xlen, ylen);
  ch->var_list = g; // Group of fields to send

  /// Post receives
  post_receive(*ch);

  //////////////////////////////////////////////////

  /// Send data going up (y+1)

  int len = 0;
  BoutReal *outbuff;

  if (UDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(ch->var_list.get(), 0, UDATA_XSPLIT, MYSUB, MYSUB + MYG,
                    std::begin(ch->umsg_sendbuff));
    // Send the data to processor UDATA_INDEST

    if (async_send) {
      MPI_Isend(std::begin(ch->umsg_sendbuff), // Buffer to send
                len,                           // Length of buffer in BoutReals
                PVEC_REAL_MPI_TYPE,            // Real variable type
                UDATA_INDEST,                  // Destination processor
                IN_SENT_UP,                    // Label (tag) for the message
                BoutComm::get(), &(ch->sendreq[0]));
    } else
      MPI_Send(std::begin(ch->umsg_sendbuff), len, PVEC_REAL_MPI_TYPE, UDATA_INDEST,
               IN_SENT_UP, BoutComm::get());
  }
  if (UDATA_OUTDEST != -1) {             // if destination for outer x data
    outbuff = &(ch->umsg_sendbuff[len]); // A pointer to the start of the second part
                                         // of the buffer
    len =
        pack_data(ch->var_list.get(), UDATA_XSPLIT, LocalNx, MYSUB, MYSUB + MYG, outbuff);
    // Send the data to processor UDATA_OUTDEST
    if (async_send) {
      MPI_Isend(outbuff, len, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, OUT_SENT_UP,
                BoutComm::get(), &(ch->sendreq[1]));
    } else
      MPI_Send(outbuff, len, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, OUT_SENT_UP,
               BoutComm::get());
  }

  /// Send data going down (y-1)

  len = 0;
  if (DDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(ch->var_list.get(), 0, DDATA_XSPLIT, MYG, 2 * MYG,
                    std::begin(ch->dmsg_sendbuff));
    // Send the data to processor DDATA_INDEST
    if (async_send) {
      MPI_Isend(std::begin(ch->dmsg_sendbuff), len, PVEC_REAL_MPI_TYPE, DDATA_INDEST,
                IN_SENT_DOWN, BoutComm::get(), &(ch->sendreq[2]));
    } else
      MPI_Send(std::begin(ch->dmsg_sendbuff), len, PVEC_REAL_MPI_TYPE, DDATA_INDEST,
               IN_SENT_DOWN, BoutComm::get());
  }
  if (DDATA_OUTDEST != -1) {             // if destination for outer x data
    outbuff = &(ch->dmsg_sendbuff[len]); // A pointer to the start of the second part
                                         // of the buffer
    len = pack_data(ch->var_list.get(), DDATA_XSPLIT, LocalNx, MYG, 2 * MYG, outbuff);
    // Send the data to processor DDATA_OUTDEST

    if (async_send) {
      MPI_Isend(outbuff, len, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, OUT_SENT_DOWN,
                BoutComm::get(), &(ch->sendreq[3]));
    } else
      MPI_Send(outbuff, len, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, OUT_SENT_DOWN,
               BoutComm::get());
  }

  /// Send to the left (x-1)

  if (IDATA_DEST != -1) {
    len = pack_data(ch->var_list.get(), MXG, 2 * MXG, MYG, MYG + MYSUB,
                    std::begin(ch->imsg_sendbuff));
    if (async_send) {
      MPI_Isend(std::begin(ch->imsg_sendbuff), len, PVEC_REAL_MPI_TYPE, IDATA_DEST,
                IN_SENT_OUT, BoutComm::get(), &(ch->sendreq[4]));
    } else
      MPI_Send(std::begin(ch->imsg_sendbuff), len, PVEC_REAL_MPI_TYPE, IDATA_DEST,
               IN_SENT_OUT, BoutComm::get());
  }

  /// Send to the right (x+1)

  if (ODATA_DEST != -1) {
    len = pack_data(ch->var_list.get(), MXSUB, MXSUB + MXG, MYG, MYG + MYSUB,
                    std::begin(ch->omsg_sendbuff));
    if (async_send) {
      MPI_Isend(std::begin(ch->omsg_sendbuff), len, PVEC_REAL_MPI_TYPE, ODATA_DEST,
                OUT_SENT_IN, BoutComm::get(), &(ch->sendreq[5]));
    } else
      MPI_Send(std::begin(ch->omsg_sendbuff), len, PVEC_REAL_MPI_TYPE, ODATA_DEST,
               OUT_SENT_IN, BoutComm::get());
  }

  /// Mark communication handle as in progress
  ch->in_progress = true;

  return static_cast<void *>(ch);
}

int BoutMesh::wait(comm_handle handle) {
  TRACE("BoutMesh::wait(comm_handle)");

  if (handle == nullptr)
    return 1;

  CommHandle *ch = static_cast<CommHandle *>(handle);

  if (!ch->in_progress)
    return 2;

  /// Start timer
  Timer timer("comms");

  ///////////// WAIT FOR DATA //////////////

  int ind, len;
  MPI_Status status;

  if (ch->var_list.size() == 0) {

    // Just waiting for a single MPI request
    MPI_Wait(ch->request, &status);
    free_handle(ch);

    return 0;
  }

  do {
    MPI_Waitany(6, ch->request, &ind, &status);
    switch (ind) {
    case 0: { // Up, inner
      unpack_data(ch->var_list.get(), 0, UDATA_XSPLIT, MYSUB + MYG, MYSUB + 2 * MYG,
                  std::begin(ch->umsg_recvbuff));
      break;
    }
    case 1: { // Up, outer
      len = msg_len(ch->var_list.get(), 0, UDATA_XSPLIT, 0, MYG);
      unpack_data(ch->var_list.get(), UDATA_XSPLIT, LocalNx, MYSUB + MYG, MYSUB + 2 * MYG,
                  &(ch->umsg_recvbuff[len]));
      break;
    }
    case 2: { // Down, inner
      unpack_data(ch->var_list.get(), 0, DDATA_XSPLIT, 0, MYG,
                  std::begin(ch->dmsg_recvbuff));
      break;
    }
    case 3: { // Down, outer
      len = msg_len(ch->var_list.get(), 0, DDATA_XSPLIT, 0, MYG);
      unpack_data(ch->var_list.get(), DDATA_XSPLIT, LocalNx, 0, MYG,
                  &(ch->dmsg_recvbuff[len]));
      break;
    }
    case 4: { // inner
      unpack_data(ch->var_list.get(), 0, MXG, MYG, MYG + MYSUB,
                  std::begin(ch->imsg_recvbuff));
      break;
    }
    case 5: { // outer
      unpack_data(ch->var_list.get(), MXSUB + MXG, MXSUB + 2 * MXG, MYG, MYG + MYSUB,
                  std::begin(ch->omsg_recvbuff));
      break;
    }
    }
    if (ind != MPI_UNDEFINED)
      ch->request[ind] = MPI_REQUEST_NULL;
  } while (ind != MPI_UNDEFINED);

  if (async_send) {
    /// Asyncronous sending: Need to check if sends have completed (frees MPI memory)
    MPI_Status async_status;

    if (UDATA_INDEST != -1)
      MPI_Wait(ch->sendreq, &async_status);
    if (UDATA_OUTDEST != -1)
      MPI_Wait(ch->sendreq + 1, &async_status);
    if (DDATA_INDEST != -1)
      MPI_Wait(ch->sendreq + 2, &async_status);
    if (DDATA_OUTDEST != -1)
      MPI_Wait(ch->sendreq + 3, &async_status);
    if (IDATA_DEST != -1)
      MPI_Wait(ch->sendreq + 4, &async_status);
    if (ODATA_DEST != -1)
      MPI_Wait(ch->sendreq + 5, &async_status);
  }

  // TWIST-SHIFT CONDITION
  if (TwistShift) {
    int jx, jy;

    // Perform Twist-shift using shifting method
    // Loop over 3D fields
    for (const auto &var : ch->var_list.field3d()) {
      // Lower boundary
      if (TS_down_in && (DDATA_INDEST != -1)) {
        for (jx = 0; jx < DDATA_XSPLIT; jx++)
          for (jy = 0; jy != MYG; jy++)
            shiftZ(*var, jx, jy, ShiftAngle[jx]);
      }
      if (TS_down_out && (DDATA_OUTDEST != -1)) {
        for (jx = DDATA_XSPLIT; jx < LocalNx; jx++)
          for (jy = 0; jy != MYG; jy++)
            shiftZ(*var, jx, jy, ShiftAngle[jx]);
      }

      // Upper boundary
      if (TS_up_in && (UDATA_INDEST != -1)) {
        for (jx = 0; jx < UDATA_XSPLIT; jx++)
          for (jy = LocalNy - MYG; jy != LocalNy; jy++)
            shiftZ(*var, jx, jy, -ShiftAngle[jx]);
      }
      if (TS_up_out && (UDATA_OUTDEST != -1)) {
        for (jx = UDATA_XSPLIT; jx < LocalNx; jx++)
          for (jy = LocalNy - MYG; jy != LocalNy; jy++)
            shiftZ(*var, jx, jy, -ShiftAngle[jx]);
      }
    }
  }

#if CHECK > 0
  // Keeping track of whether communications have been done
  for (const auto &var : ch->var_list)
    var->doneComms();
#endif

  free_handle(ch);

  return 0;
}

/***************************************************************
 *             Non-Local Communications
 ***************************************************************/

MPI_Request BoutMesh::sendToProc(int xproc, int yproc, BoutReal *buffer, int size,
                                 int tag) {
  Timer timer("comms");

  MPI_Request request;

  MPI_Isend(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(xproc, yproc), tag,
            BoutComm::get(), &request);

  return request;
}

comm_handle BoutMesh::receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size,
                                      int tag) {
  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(xproc, yproc), tag,
            BoutComm::get(), ch->request);

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

int BoutMesh::getNXPE() { return NXPE; }

int BoutMesh::getNYPE() { return NYPE; }

int BoutMesh::getXProcIndex() { return PE_XIND; }

int BoutMesh::getYProcIndex() { return PE_YIND; }

/****************************************************************
 *                 X COMMUNICATIONS
 *
 * Intended mainly to handle the perpendicular inversion operators
 ****************************************************************/

bool BoutMesh::firstX() { return PE_XIND == 0; }

bool BoutMesh::lastX() { return PE_XIND == NXPE - 1; }

int BoutMesh::sendXOut(BoutReal *buffer, int size, int tag) {
  if (PE_XIND == NXPE - 1)
    return 1;

  Timer timer("comms");

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND + 1, PE_YIND), tag,
           BoutComm::get());

  return 0;
}

int BoutMesh::sendXIn(BoutReal *buffer, int size, int tag) {
  if (PE_XIND == 0)
    return 1;

  Timer timer("comms");

  MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND - 1, PE_YIND), tag,
           BoutComm::get());

  return 0;
}

comm_handle BoutMesh::irecvXOut(BoutReal *buffer, int size, int tag) {
  if (PE_XIND == NXPE - 1)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND + 1, PE_YIND), tag,
            BoutComm::get(), ch->request);

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvXIn(BoutReal *buffer, int size, int tag) {
  if (PE_XIND == 0)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND - 1, PE_YIND), tag,
            BoutComm::get(), ch->request);

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

/****************************************************************
 *                 Y COMMUNICATIONS
 *
 * Intended mainly to handle the non-local heat flux integrations
 ****************************************************************/

bool BoutMesh::firstY() const { return PE_YIND == 0; }

bool BoutMesh::lastY() const { return PE_YIND == NYPE - 1; }

bool BoutMesh::firstY(int xpos) const {
  int xglobal = XGLOBAL(xpos);
  int rank;

  if (xglobal < ixseps_inner) {
    MPI_Comm_rank(comm_inner, &rank);
  } else if (xglobal < ixseps_outer) {
    MPI_Comm_rank(comm_middle, &rank);
  } else {
    MPI_Comm_rank(comm_outer, &rank);
  }
  return rank == 0;
}

bool BoutMesh::lastY(int xpos) const {
  int xglobal = XGLOBAL(xpos);
  int rank;
  int size;

  if (xglobal < ixseps_inner) {
    MPI_Comm_size(comm_inner, &size);
    MPI_Comm_rank(comm_inner, &rank);
  } else if (xglobal < ixseps_outer) {
    MPI_Comm_size(comm_middle, &size);
    MPI_Comm_rank(comm_middle, &rank);
  } else {
    MPI_Comm_size(comm_outer, &size);
    MPI_Comm_rank(comm_outer, &rank);
  }
  return rank == size - 1;
}

int BoutMesh::UpXSplitIndex() { return UDATA_XSPLIT; }

int BoutMesh::DownXSplitIndex() { return DDATA_XSPLIT; }

int BoutMesh::sendYOutIndest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1)
    return 1;

  Timer timer("comms");

  if (UDATA_INDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_INDEST, tag, BoutComm::get());
  else
    throw BoutException("Expected UDATA_INDEST to exist, but it does not.");
  return 0;
}

int BoutMesh::sendYOutOutdest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1)
    return 1;

  Timer timer("comms");

  if (UDATA_OUTDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, tag, BoutComm::get());
  else
    throw BoutException("Expected UDATA_OUTDEST to exist, but it does not.");

  return 0;
}

int BoutMesh::sendYInIndest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == 0)
    return 1;

  Timer timer("comms");

  if (DDATA_INDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_INDEST, tag, BoutComm::get());
  else
    throw BoutException("Expected DDATA_INDEST to exist, but it does not.");

  return 0;
}

int BoutMesh::sendYInOutdest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == 0)
    return 1;

  Timer timer("comms");

  if (DDATA_OUTDEST != -1)
    MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, tag, BoutComm::get());
  else
    throw BoutException("Expected DDATA_OUTDEST to exist, but it does not.");

  return 0;
}

comm_handle BoutMesh::irecvYOutIndest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  if (UDATA_INDEST != -1)
    MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_INDEST, tag, BoutComm::get(),
              ch->request);
  else
    throw BoutException("Expected UDATA_INDEST to exist, but it does not.");

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvYOutOutdest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  if (UDATA_OUTDEST != -1)
    MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, tag, BoutComm::get(),
              ch->request);
  else
    throw BoutException("Expected UDATA_OUTDEST to exist, but it does not.");

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvYInIndest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == 0)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  if (DDATA_INDEST != -1)
    MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_INDEST, tag, BoutComm::get(),
              ch->request);
  else
    throw BoutException("Expected DDATA_INDEST to exist, but it does not.");

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvYInOutdest(BoutReal *buffer, int size, int tag) {
  if (PE_YIND == 0)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  if (DDATA_OUTDEST != -1)
    MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, tag, BoutComm::get(),
              ch->request);
  else
    throw BoutException("Expected DDATA_OUTDEST to exist, but it does not.");

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

/****************************************************************
 *                 GRID INDEX ROUTINES
 *
 * These routines translate between local and global coordinates
 * This so that more flexible indexing (e.g. not equal on processors)
 * can be implemented later (maybe)
 ****************************************************************/

/// Returns the processor number, given X and Y processor indices.
/*!
 * If out of range returns -1 (no processor)
 */

int BoutMesh::PROC_NUM(int xind, int yind) {
  if ((xind >= NXPE) || (xind < 0))
    return -1;
  if ((yind >= NYPE) || (yind < 0))
    return -1;

  return yind * NXPE + xind;
}

/// Returns the global X index given a local index
int BoutMesh::XGLOBAL(int xloc) const { return xloc + PE_XIND * MXSUB; }

/// Returns the global X index given a local index
int BoutMesh::XGLOBAL(BoutReal xloc, BoutReal &xglo) const {
  xglo = xloc + PE_XIND * MXSUB;
  return static_cast<int>(xglo);
}

/// Returns a local X index given a global index
int BoutMesh::XLOCAL(int xglo) const { return xglo - PE_XIND * MXSUB; }

/// Returns the global Y index given a local index
int BoutMesh::YGLOBAL(int yloc) const { return yloc + PE_YIND * MYSUB - MYG; }

/// Returns the global Y index given a local index
int BoutMesh::YGLOBAL(BoutReal yloc, BoutReal &yglo) const {
  yglo = yloc + PE_YIND * MYSUB - MYG;
  return static_cast<int>(yglo);
}

/// Global Y index given local index and processor
int BoutMesh::YGLOBAL(int yloc, int yproc) const { return yloc + yproc * MYSUB - MYG; }

/// Returns a local Y index given a global index
int BoutMesh::YLOCAL(int yglo) const { return yglo - PE_YIND * MYSUB + MYG; }

int BoutMesh::YLOCAL(int yglo, int yproc) const { return yglo - yproc * MYSUB + MYG; }

/// Return the Y processor number given a global Y index
int BoutMesh::YPROC(int yind) {
  if ((yind < 0) || (yind > ny))
    return -1;
  return yind / MYSUB;
}

/// Return the X processor number given a global X index
int BoutMesh::XPROC(int xind) { return (xind >= MXG) ? (xind - MXG) / MXSUB : 0; }

/****************************************************************
 *                       CONNECTIONS
 ****************************************************************/

/// Connection initialisation: Set processors in a simple 2D grid
void BoutMesh::default_connections() {
  DDATA_XSPLIT = UDATA_XSPLIT = 0;  // everything by default outside (arb. choice)
  DDATA_INDEST = UDATA_INDEST = -1; // since nothing inside

  DDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND - 1);
  UDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND + 1);

  IDATA_DEST = PROC_NUM(PE_XIND - 1, PE_YIND);
  ODATA_DEST = PROC_NUM(PE_XIND + 1, PE_YIND);

  TS_up_in = TS_up_out = TS_down_in = TS_down_out = false; // No twist-shifts

  /// Check if X is periodic
  if (periodicX) {
    if (PE_XIND == (NXPE - 1))
      ODATA_DEST = PROC_NUM(0, PE_YIND);

    if (PE_XIND == 0)
      IDATA_DEST = PROC_NUM(NXPE - 1, PE_YIND);
  }
}

/// Add a topology connection
/*!
 * Set ypos1 and ypos2 to be neighbours in the range xge <= x < xlt.
 * Optional argument ts sets whether to use twist-shift condition
 */
void BoutMesh::set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts) {
  int ype1, ype2; // the two Y processor indices
  int ypeup, ypedown;
  int yind1, yind2;

  if (xlt <= xge)
    return;

  if ((ypos1 < 0) || (ypos1 >= MY)) {
    output_warn.write("WARNING adding connection: poloidal index %d out of range\n",
                      ypos1);
    return;
  }
  if ((ypos2 < 0) || (ypos2 >= MY)) {
    output_warn.write("WARNING adding connection: poloidal index %d out of range\n",
                      ypos2);
    return;
  }

  ype1 = YPROC(ypos1);
  ype2 = YPROC(ypos2);

  /* y index within processors */
  yind1 = YLOCAL(ypos1, ype1);
  yind2 = YLOCAL(ypos2, ype2);

  /* Check which boundary the connection is on */
  if ((yind1 == MYG) && (yind2 == MYSUB + MYG - 1)) {
    ypeup = ype2;   /* processor sending data up (+ve y) */
    ypedown = ype1; /* processor sending data down (-ve y) */
  } else if ((yind2 == MYG) && (yind1 == MYSUB + MYG - 1)) {
    ypeup = ype1;
    ypedown = ype2;
  } else {
    throw BoutException(
        "ERROR adding connection: y index %d or %d not on processor boundary\n", ypos1,
        ypos2);
  }

  /* check the x ranges are possible */
  if ((xge != 0) && (xlt != MX)) {
    throw BoutException(
        "ERROR adding connection(%d,%d,%d,%d): can only divide X domain in 2\n", ypos1,
        ypos2, xge, xlt);
  }

  output_info.write(
      "Connection between top of Y processor %d and bottom of %d in range %d <= x < %d\n",
      ypeup, ypedown, xge, xlt);

  // Convert X coordinates into local indices

  xge = XLOCAL(xge);
  xlt = XLOCAL(xlt);

  if ((xge >= LocalNx) || (xlt <= 0)) {
    return; // Not in this x domain
  }

  if (xge < 0)
    xge = 0;
  if (xlt > LocalNx)
    xlt = LocalNx;

  if (MYPE == PROC_NUM(PE_XIND, ypeup)) { /* PROCESSOR SENDING +VE Y */
    /* Set the branch cut x position */
    if (xge <= MXG) {
      /* Connect on the inside */
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = PROC_NUM(PE_XIND, ypedown);
      if (UDATA_XSPLIT == LocalNx)
        UDATA_OUTDEST = -1;

      TS_up_in = ts; // Twist-shift

      output_info.write("=> This processor sending in up\n");
    } else {
      /* Connect on the outside */
      if (UDATA_XSPLIT <= 0)
        UDATA_INDEST = UDATA_OUTDEST;
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = PROC_NUM(PE_XIND, ypedown);
      if (UDATA_XSPLIT <= 0)
        UDATA_INDEST = -1;

      TS_up_out = ts;
      output_info.write("=> This processor sending out up\n");
    }
  }

  if (MYPE == PROC_NUM(PE_XIND, ypedown)) { /* PROCESSOR SENDING -VE Y */
    /* Set the branch cut x position */
    if (xge <= MXG) {
      /* Connect on the inside */
      DDATA_XSPLIT = xlt;
      DDATA_INDEST = PROC_NUM(PE_XIND, ypeup);
      if (DDATA_XSPLIT == LocalNx)
        DDATA_OUTDEST = -1;

      TS_down_in = ts;

      output_info.write("=> This processor sending in down\n");
    } else {
      /* Connect on the outside */
      if (DDATA_XSPLIT <= 0)
        DDATA_INDEST = DDATA_OUTDEST;
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = PROC_NUM(PE_XIND, ypeup);
      if (DDATA_XSPLIT == 0)
        DDATA_INDEST = -1;

      TS_down_out = ts;

      output_info.write("=> This processor sending out down\n");
    }
  }
}

/// Add a divertor target or limiter
/*!
 * ypos is the y index which will become an upper target
 * ypos+1 will become a lower target.
 * Target created in the range xge <= x < xlt.
 */
void BoutMesh::add_target(int ypos, int xge, int xlt) {
  if (xlt <= xge)
    return;

  if ((ypos < 0) || (ypos >= MY)) {
    output_warn.write("WARNING adding target: poloidal index %d out of range\n", ypos);
    return;
  }

  int ypeup = YPROC(ypos);
  int ypedown = YPROC(ypos + 1);
  if (ypeup == ypedown) {
    throw BoutException("Adding target at y=%d in middle of processor %d\n", ypos, ypeup);
  }

  output_info.write(
      "Target at top of Y processor %d and bottom of %d in range %d <= x < %d\n", ypeup,
      ypedown, xge, xlt);

  // Convert X coordinates into local indices
  xge = XLOCAL(xge);
  xlt = XLOCAL(xlt);
  if ((xge >= LocalNx) || (xlt <= 0)) {
    return; // Not in this x domain
  }

  if (MYPE == PROC_NUM(PE_XIND, ypeup)) {
    // Target on upper processor boundary
    if (xge <= MXG) {
      // Target on inside
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = -1;
      if (xlt >= LocalNx)
        UDATA_OUTDEST = -1;
      output_info.write("=> This processor has target upper inner\n");
    } else {
      // Target on outside
      if (UDATA_XSPLIT <= 0)
        UDATA_INDEST = UDATA_OUTDEST;
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = -1;
      if (xge <= 0)
        UDATA_INDEST = -1;
      output_info.write("=> This processor has target upper outer\n");
    }
  }
  if (MYPE == PROC_NUM(PE_XIND, ypedown)) {
    // Target on upper processor boundary
    if (xge <= MXG) {
      // Target on inside
      DDATA_XSPLIT = xlt;
      DDATA_INDEST = -1;
      if (xlt >= LocalNx)
        DDATA_OUTDEST = -1;
      output_info.write("=> This processor has target lower inner\n");
    } else {
      // Target on outside
      if (DDATA_XSPLIT <= 0)
        DDATA_INDEST = DDATA_OUTDEST;
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = -1;
      if (xge <= 0)
        DDATA_INDEST = -1;
      output_info.write("=> This processor has target lower outer\n");
    }
  }
}

/****************************************************************
 *                MAIN TOPOLOGY FUNCTION
 ****************************************************************/

void BoutMesh::topology() {
  // Perform checks common to all topologies

  if (NPES != NXPE * NYPE) {
    throw BoutException("\tTopology error: npes=%d is not equal to NXPE*NYPE=%d\n", NPES,
                        NXPE * NYPE);
  }
  if (MYSUB * NYPE != MY) {
    throw BoutException("\tTopology error: MYSUB[%d] * NYPE[%d] != MY[%d]\n", MYSUB, NYPE,
                        MY);
  }
  if (MXSUB * NXPE != MX) {
    throw BoutException("\tTopology error: MXSUB[%d] * NXPE[%d] != MX[%d]\n", MXSUB, NXPE,
                        MX);
  }

  if ((NXPE > 1) && (MXSUB < MXG)) {
    throw BoutException("\tERROR: Grid X size must be >= guard cell size\n");
  }
  if (MYSUB < MYG) {
    throw BoutException("\tERROR: Grid Y size must be >= guard cell size\n");
  }

  if (jyseps2_1 == jyseps1_2) {
    /********* SINGLE NULL OPERATION *************/
    output_info.write("\tEQUILIBRIUM IS SINGLE NULL (SND) \n");

    /* Set separatrices - x location all the same */
    ixseps_inner = ixseps_outer = ixseps_upper = ixseps_lower = ixseps1;

    default_connections();
    set_connection(jyseps1_1 + 1, jyseps2_2, 0, ixseps1,
                   true);                                 // Twist-shift this connection
    set_connection(jyseps1_1, jyseps2_2 + 1, 0, ixseps1); // No twist-shift in PF region

  } else {
    /*************** DOUBLE NULL OPERATION *******************/
    /* UPPER LEGS: Do not have to be the same length as each
       other or lower legs, but do have to have an integer number
       of processors */
    if ((ny_inner - jyseps2_1 - 1) % MYSUB != 0) {
      throw BoutException("\tTopology error: Upper inner leg does not have integer "
                          "number of processors\n");
    }
    if ((jyseps1_2 - ny_inner + 1) % MYSUB != 0) {
      throw BoutException("\tTopology error: Upper outer leg does not have integer "
                          "number of processors\n");
    }

    if (ixseps1 == ixseps2) {
      /*************** CONNECTED (balanced) DOUBLE NULL ******************/
      output_info.write("\tEQUILIBRIUM IS CONNECTED DOUBLE NULL (CDND)\n");
      /* all separatrix indices the same */
      ixseps_inner = ixseps_outer = ixseps_upper = ixseps_lower = ixseps1;

    } else if (ixseps1 < ixseps2) {
      /*************** LOWER DOUBLE NULL **********************/
      output_info.write("\tEQUILIBRIUM IS LOWER DOUBLE NULL (LDND)\n");
      ixseps_inner = ixseps_lower = ixseps1;
      ixseps_outer = ixseps_upper = ixseps2;
    } else {
      /*************** UPPER DOUBLE NULL **********************/
      output_info.write("\tEQUILIBRIUM IS UPPER DOUBLE NULL (UDND)\n");
      ixseps_inner = ixseps_upper = ixseps2;
      ixseps_outer = ixseps_lower = ixseps1;
    }

    /* Following code works for any Double Null */

    /********* DND CONNECTIONS **********/
    default_connections();
    /* Lower x-point */
    set_connection(jyseps1_1 + 1, jyseps2_2, 0, ixseps_lower,
                   ixseps1 <= ixseps2);                        /* Core */
    set_connection(jyseps1_1, jyseps2_2 + 1, 0, ixseps_lower); /* PF   */
    /* Upper x-point */
    set_connection(jyseps2_1, jyseps1_2 + 1, 0, ixseps_upper,
                   ixseps1 > ixseps2);                         /* Core */
    set_connection(jyseps2_1 + 1, jyseps1_2, 0, ixseps_upper); /* PF   */

    // Add target plates at the top
    add_target(ny_inner - 1, 0, nx);
  }

  MYPE_IN_CORE = 0; // processor not in core
  if ((ixseps_inner > 0) &&
      (((PE_YIND * MYSUB > jyseps1_1) && (PE_YIND * MYSUB <= jyseps2_1)) ||
       ((PE_YIND * MYSUB > jyseps1_2) && (PE_YIND * MYSUB <= jyseps2_2)))) {
    MYPE_IN_CORE = 1; /* processor is in the core */
  }

  if (DDATA_XSPLIT > LocalNx)
    DDATA_XSPLIT = LocalNx;
  if (UDATA_XSPLIT > LocalNx)
    UDATA_XSPLIT = LocalNx;

  // Print out settings
  output_info.write("\tMYPE_IN_CORE = %d\n", MYPE_IN_CORE);
  output_info.write("\tDXS = %d, DIN = %d. DOUT = %d\n", DDATA_XSPLIT, DDATA_INDEST,
                    DDATA_OUTDEST);
  output_info.write("\tUXS = %d, UIN = %d. UOUT = %d\n", UDATA_XSPLIT, UDATA_INDEST,
                    UDATA_OUTDEST);
  output_info.write("\tXIN = %d, XOUT = %d\n", IDATA_DEST, ODATA_DEST);

  output_info.write("\tTwist-shift: ");
  if (TS_down_in)
    output_info.write("DI ");
  if (TS_down_out)
    output_info.write("DO ");
  if (TS_up_in)
    output_info.write("UI ");
  if (TS_up_out)
    output_info.write("UO ");
  output_info.write("\n");
}

/****************************************************************
 *                     Communication handles
 ****************************************************************/

BoutMesh::CommHandle *BoutMesh::get_handle(int xlen, int ylen) {
  if (comm_list.empty()) {
    // Allocate a new CommHandle

    auto *ch = new CommHandle;
    for (auto &i : ch->request)
      i = MPI_REQUEST_NULL;

    if (ylen > 0) {
      ch->umsg_sendbuff = Array<BoutReal>(ylen);
      ch->dmsg_sendbuff = Array<BoutReal>(ylen);
      ch->umsg_recvbuff = Array<BoutReal>(ylen);
      ch->dmsg_recvbuff = Array<BoutReal>(ylen);
    }

    if (xlen > 0) {
      ch->imsg_sendbuff = Array<BoutReal>(xlen);
      ch->omsg_sendbuff = Array<BoutReal>(xlen);
      ch->imsg_recvbuff = Array<BoutReal>(xlen);
      ch->omsg_recvbuff = Array<BoutReal>(xlen);
    }

    ch->xbufflen = xlen;
    ch->ybufflen = ylen;

    ch->in_progress = false;

    return ch;
  }

  // Pop first pointer off the list
  CommHandle *ch = comm_list.front();
  comm_list.pop_front();

  // Check that the buffers are big enough (NOTE: Could search list for bigger buffers)
  if (ch->ybufflen < ylen) {
    ch->umsg_sendbuff = Array<BoutReal>(ylen);
    ch->dmsg_sendbuff = Array<BoutReal>(ylen);
    ch->umsg_recvbuff = Array<BoutReal>(ylen);
    ch->dmsg_recvbuff = Array<BoutReal>(ylen);

    ch->ybufflen = ylen;
  }
  if (ch->xbufflen < xlen) {
    ch->imsg_sendbuff = Array<BoutReal>(xlen);
    ch->omsg_sendbuff = Array<BoutReal>(xlen);
    ch->imsg_recvbuff = Array<BoutReal>(xlen);
    ch->omsg_recvbuff = Array<BoutReal>(xlen);

    ch->xbufflen = xlen;
  }

  ch->in_progress = false;

  ch->var_list.clear();

  return ch;
}

void BoutMesh::free_handle(CommHandle *h) {
  h->var_list.clear();
  comm_list.push_front(h);
}

void BoutMesh::clear_handles() {
  while (!comm_list.empty()) {
    CommHandle *ch = comm_list.front();

    delete ch;

    comm_list.pop_front();
  }
}

/****************************************************************
 *                   Communication utilities
 ****************************************************************/

int BoutMesh::pack_data(const vector<FieldData *> &var_list, int xge, int xlt, int yge,
                        int ylt, BoutReal *buffer) {

  int len = 0;

  /// Loop over variables
  for (const auto &var : var_list) {
    if (var->is3D()) {
      // 3D variable
      ASSERT2(static_cast<Field3D *>(var)->isAllocated());
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++) {
          for (int jz = 0; jz < LocalNz; jz++, len++) {
            buffer[len] = (*static_cast<Field3D *>(var))(jx, jy, jz);
          }
        }
      }
    } else {
      // 2D variable
      ASSERT2(static_cast<Field2D *>(var)->isAllocated());
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++, len++) {
          buffer[len] = (*static_cast<Field2D *>(var))(jx, jy);
        }
      }
    }
  }

  return (len);
}

int BoutMesh::unpack_data(const vector<FieldData *> &var_list, int xge, int xlt, int yge,
                          int ylt, BoutReal *buffer) {

  int len = 0;

  /// Loop over variables
  for (const auto &var : var_list) {
    if (var->is3D()) {
      // 3D variable
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++) {
          for (int jz = 0; jz < LocalNz; jz++, len++) {
            (*static_cast<Field3D *>(var))(jx, jy, jz) = buffer[len];
          }
        }
      }
    } else {
      // 2D variable
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++, len++) {
          (*static_cast<Field2D *>(var))(jx, jy) = buffer[len];
        }
      }
    }
  }

  return (len);
}

/****************************************************************
 *                 SURFACE ITERATION
 ****************************************************************/

bool BoutMesh::periodicY(int jx) const {
  return (XGLOBAL(jx) < ixseps_inner) && MYPE_IN_CORE;
}

bool BoutMesh::periodicY(int jx, BoutReal &ts) const {
  ts = 0.;
  if ((XGLOBAL(jx) < ixseps_inner) && MYPE_IN_CORE) {
    if (TwistShift)
      ts = ShiftAngle[jx];
    return true;
  }
  return false;
}

int BoutMesh::ySize(int xpos) const {
  int xglobal = XGLOBAL(xpos);
  int yglobal = YGLOBAL(MYG);

  if ((xglobal < ixseps_lower) && ((yglobal <= jyseps1_1) || (yglobal > jyseps2_2))) {
    // Lower PF region
    return (jyseps1_1 + 1) + (ny - jyseps2_2);

  } else if ((xglobal < ixseps_upper) && (yglobal > jyseps2_1) &&
             (yglobal >= jyseps1_2)) {
    // Upper PF region
    return jyseps1_2 - jyseps2_1;

  } else if (xglobal < ixseps_inner) {
    // Core
    return (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

  } else if (jyseps2_1 == jyseps1_2) {
    // Single null, so in the SOL
    return ny;

  } else if ((xglobal >= ixseps_inner) && (xglobal < ixseps_outer)) {
    // Intermediate SOL in DND

    if (ixseps_lower < ixseps_upper) {
      // Connects to lower divertor
      return (jyseps2_1 + 1) + (ny - jyseps1_2);
    } else {
      // Connects to upper divertor
      return jyseps2_2 - jyseps1_1;
    }
  } else if (yglobal < ny_inner) {
    // Inner SOL
    return ny_inner;
  }
  // Outer SOL
  return ny - ny_inner;
}

MPI_Comm BoutMesh::getYcomm(int xpos) const {
  int xglobal = XGLOBAL(xpos);

  if (xglobal < ixseps_inner) {
    return comm_inner;
  } else if (xglobal < ixseps_outer) {
    return comm_middle;
  }
  return comm_outer;
}

/****************************************************************
 *                 Range iteration
 ****************************************************************/

void BoutMesh::addBoundaryRegions() {
  std::list<std::string> all_boundaries; ///< Keep track of all boundary regions
  
  // Lower Inner Y
  int xs = 0;
  int xe = LocalNx - 1;

  if (!firstY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
      xs = DDATA_XSPLIT;
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1))
      xe = DDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  }
  
  addRegion3D("RGN_LOWER_INNER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                                 LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_LOWER_INNER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                                 LocalNy, 1, maxregionblocksize));

  all_boundaries.emplace_back("RGN_LOWER_INNER_Y");

  // Lower Outer Y
  
  xs = 0;
  xe = LocalNx - 1;
  if (!firstY()) {
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
      xs = DDATA_XSPLIT;
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1))
      xe = DDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  } else {
    xs = -1;
    xe = -2;
  }

  addRegion3D("RGN_LOWER_OUTER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                                 LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_LOWER_OUTER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                                 LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_LOWER_OUTER_Y");
  
  // Lower Y

  xs = 0;
  xe = LocalNx - 1;
  if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
    xs = DDATA_XSPLIT;
  if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1))
    xe = DDATA_XSPLIT - 1;

  if (xs < xstart)
    xs = xstart;
  if (xe > xend)
    xe = xend;

  addRegion3D("RGN_LOWER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                           LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_LOWER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                           LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_LOWER_Y");
  
  // Upper Inner Y

  xs = 0;
  xe = LocalNx - 1;

  if (!lastY()) {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
      xs = UDATA_XSPLIT;
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1))
      xe = UDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  } else {
    xs = -1;
    xe = -2;
  }
  
  addRegion3D("RGN_UPPER_INNER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                                 LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_UPPER_INNER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                                 LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_UPPER_INNER_Y");

  // Upper Outer Y
  
  xs = 0;
  xe = LocalNx - 1;

  if (!lastY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
      xs = UDATA_XSPLIT;
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1))
      xe = UDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  }

  addRegion3D("RGN_UPPER_OUTER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                                 LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_UPPER_OUTER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                                 LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_UPPER_OUTER_Y");

  // Upper Y

  xs = 0;
  xe = LocalNx - 1;
  if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
    xs = UDATA_XSPLIT;
  if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1))
    xe = UDATA_XSPLIT - 1;

  if (xs < xstart)
    xs = xstart;
  if (xe > xend)
    xe = xend;

  addRegion3D("RGN_UPPER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                           LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_UPPER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                           LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_UPPER_Y");
  
  // Inner X
  if(mesh->firstX() && !mesh->periodicX) {
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, xstart-1, ystart, yend, 0, LocalNz-1,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_INNER_X", Region<Ind2D>(0, xstart-1, ystart, yend, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
    all_boundaries.emplace_back("RGN_INNER_X");
    
    output_info.write("\tBoundary region inner X\n");
  } else {
    // Empty region
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_INNER_X", Region<Ind2D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
  }

  // Outer X
  if(mesh->firstX() && !mesh->periodicX) {
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(xend+1, LocalNx-1, ystart, yend, 0, LocalNz-1,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_OUTER_X", Region<Ind2D>(xend+1, LocalNx-1, ystart, yend, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
    all_boundaries.emplace_back("RGN_OUTER_X");
    
    output_info.write("\tBoundary region outer X\n");
  } else {
    // Empty region
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_OUTER_X", Region<Ind2D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
  }

  // Join boundary regions together
  
  Region<Ind3D> bndry3d; // Empty
  for (const auto &region_name : all_boundaries) {
    bndry3d += getRegion3D(region_name);
  }
  bndry3d.unique(); // Ensure that the points are unique

  // Create a region which is all boundaries
  addRegion3D("RGN_BNDRY", bndry3d);

  Region<Ind2D> bndry2d; // Empty
  for (const auto &region_name : all_boundaries) {
    bndry2d += getRegion2D(region_name);
  }
  bndry2d.unique(); // Ensure that the points are unique

  // Create a region which is all boundaries
  addRegion2D("RGN_BNDRY", bndry2d);
}

const RangeIterator BoutMesh::iterateBndryLowerInnerY() const {

  int xs = 0;
  int xe = LocalNx - 1;

  if (!firstY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
      xs = DDATA_XSPLIT;
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1))
      xe = DDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  }
  return RangeIterator(xs, xe);
}

const RangeIterator BoutMesh::iterateBndryLowerOuterY() const {

  int xs = 0;
  int xe = LocalNx - 1;
  if (!firstY()) {
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
      xs = DDATA_XSPLIT;
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1))
      xe = DDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  } else {
    xs = -1;
    xe = -2;
  }
  return RangeIterator(xs, xe);
}

const RangeIterator BoutMesh::iterateBndryLowerY() const {
  int xs = 0;
  int xe = LocalNx - 1;
  if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart))
    xs = DDATA_XSPLIT;
  if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1))
    xe = DDATA_XSPLIT - 1;

  if (xs < xstart)
    xs = xstart;
  if (xe > xend)
    xe = xend;

  return RangeIterator(xs, xe);
}

const RangeIterator BoutMesh::iterateBndryUpperInnerY() const {
  int xs = 0;
  int xe = LocalNx - 1;

  if (!lastY()) {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
      xs = UDATA_XSPLIT;
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1))
      xe = UDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  } else {
    xs = -1;
    xe = -2;
  }
  return RangeIterator(xs, xe);
}

const RangeIterator BoutMesh::iterateBndryUpperOuterY() const {
  int xs = 0;
  int xe = LocalNx - 1;

  if (!lastY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
      xs = UDATA_XSPLIT;
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1))
      xe = UDATA_XSPLIT - 1;

    if (xs < xstart)
      xs = xstart;
    if (xe > xend)
      xe = xend;
  }
  return RangeIterator(xs, xe);
}

const RangeIterator BoutMesh::iterateBndryUpperY() const {
  int xs = 0;
  int xe = LocalNx - 1;
  if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart))
    xs = UDATA_XSPLIT;
  if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1))
    xe = UDATA_XSPLIT - 1;

  if (xs < xstart)
    xs = xstart;
  if (xe > xend)
    xe = xend;

  return RangeIterator(xs, xe);
}

vector<BoundaryRegion *> BoutMesh::getBoundaries() { return boundary; }

vector<BoundaryRegionPar *> BoutMesh::getBoundariesPar() { return par_boundary; }

void BoutMesh::addBoundaryPar(BoundaryRegionPar *bndry) {
  output_info << "Adding new parallel boundary: " << bndry->label << endl;
  par_boundary.push_back(bndry);
}

const Field3D BoutMesh::smoothSeparatrix(const Field3D &f) {
  Field3D result(f);
  if ((ixseps_inner > 0) && (ixseps_inner < nx - 1)) {
    result.allocate();
    if (XPROC(ixseps_inner) == PE_XIND) {
      int x = XLOCAL(ixseps_inner);
      for (int y = 0; y < LocalNy; y++)
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x - 1, y, z));
        }
    }
    if (XPROC(ixseps_inner - 1) == PE_XIND) {
      int x = XLOCAL(ixseps_inner - 1);
      for (int y = 0; y < LocalNy; y++)
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x + 1, y, z));
        }
    }
  }
  if ((ixseps_outer > 0) && (ixseps_outer < nx - 1) && (ixseps_outer != ixseps_inner)) {
    result.allocate();
    if (XPROC(ixseps_outer) == PE_XIND) {
      int x = XLOCAL(ixseps_outer);
      for (int y = 0; y < LocalNy; y++)
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x - 1, y, z));
        }
    }
    if (XPROC(ixseps_outer - 1) == PE_XIND) {
      int x = XLOCAL(ixseps_outer - 1);
      for (int y = 0; y < LocalNy; y++)
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x + 1, y, z));
        }
    }
  }
  return result;
}

BoutReal BoutMesh::GlobalX(int jx) const {
  if (symmetricGlobalX) {
    // With this definition the boundary sits dx/2 away form the first/last inner points
    return (0.5 + XGLOBAL(jx) - (nx - MX) * 0.5) / static_cast<BoutReal>(MX);
  }
  return static_cast<BoutReal>(XGLOBAL(jx)) / static_cast<BoutReal>(MX);
}

BoutReal BoutMesh::GlobalX(BoutReal jx) const {

  // Get global X index as a BoutReal
  BoutReal xglo;
  XGLOBAL(jx, xglo);

  if (symmetricGlobalX) {
    // With this definition the boundary sits dx/2 away form the first/last inner points
    return (0.5 + xglo - (nx - MX) * 0.5) / static_cast<BoutReal>(MX);
  }
  return xglo / static_cast<BoutReal>(MX);
}

BoutReal BoutMesh::GlobalY(int jy) const {
  if (symmetricGlobalY) {
    BoutReal yi = YGLOBAL(jy);
    int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

    if (yi < ny_inner) {
      yi -= jyseps1_1 + 0.5;
    } else {
      // Result in core between 0.5 and 1.0
      yi -= jyseps1_1 + 0.5 + (jyseps1_2 - jyseps2_1);
    }
    return yi / nycore;
  }

  int ly = YGLOBAL(jy); // global poloidal index across subdomains
  int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

  if (MYPE_IN_CORE) {
    // Turn ly into an index over the core cells only
    if (ly <= jyseps2_1) {
      ly -= jyseps1_1 + 1;
    } else
      ly -= jyseps1_1 + 1 + (jyseps1_2 - jyseps2_1);
  } else {
    // Not in core. Need to get the last "core" value
    if (ly <= jyseps1_1) {
      // Inner lower leg
      ly = 0;
    } else if ((ly > jyseps2_1) && (ly <= jyseps1_2)) {
      // Upper legs
      ly = jyseps2_1 - jyseps1_1;
    } else if (ly > jyseps2_2) {
      // Outer lower leg
      ly = nycore;
    }
  }

  return static_cast<BoutReal>(ly) / static_cast<BoutReal>(nycore);
}

BoutReal BoutMesh::GlobalY(BoutReal jy) const {

  // Get global Y index as a BoutReal
  BoutReal yglo;
  YGLOBAL(jy, yglo);

  if (symmetricGlobalY) {
    BoutReal yi = yglo;
    int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

    if (yi < ny_inner) {
      yi -= jyseps1_1 + 0.5;
    } else {
      // Result in core between 0.5 and 1.0
      yi -= jyseps1_1 + 0.5 + (jyseps1_2 - jyseps2_1);
    }
    return yi / nycore;
  }

  int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

  if (MYPE_IN_CORE) {
    // Turn yglo into an index over the core cells onyglo
    if (yglo <= jyseps2_1) {
      yglo -= jyseps1_1 + 1;
    } else {
      yglo -= jyseps1_1 + 1 + (jyseps1_2 - jyseps2_1);
    }
  } else {
    // Not in core. Need to get the last "core" value
    if (yglo <= jyseps1_1) {
      // Inner lower leg
      yglo = 0;
    } else if ((yglo > jyseps2_1) && (yglo <= jyseps1_2)) {
      // Upper legs
      yglo = jyseps2_1 - jyseps1_1;
    } else if (yglo > jyseps2_2) {
      // Outer lower leg
      yglo = nycore;
    }
  }

  return yglo / static_cast<BoutReal>(nycore);
}

void BoutMesh::outputVars(Datafile &file) {
  file.add(zperiod, "zperiod", false);
  file.add(MXSUB, "MXSUB", false);
  file.add(MYSUB, "MYSUB", false);
  file.add(MXG, "MXG", false);
  file.add(MYG, "MYG", false);
  file.add(nx, "nx", false);
  file.add(ny, "ny", false);
  file.add(LocalNz, "MZ", false);
  file.add(NXPE, "NXPE", false);
  file.add(NYPE, "NYPE", false);
  file.add(ZMAX, "ZMAX", false);
  file.add(ZMIN, "ZMIN", false);
  file.add(ixseps1, "ixseps1", false);
  file.add(ixseps2, "ixseps2", false);
  file.add(jyseps1_1, "jyseps1_1", false);
  file.add(jyseps1_2, "jyseps1_2", false);
  file.add(jyseps2_1, "jyseps2_1", false);
  file.add(jyseps2_2, "jyseps2_2", false);

  coordinates()->outputVars(file);
}
