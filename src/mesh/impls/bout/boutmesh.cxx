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

#include <algorithm>
#include <iterator>
#include <set>

/// MPI type of BoutReal for communications
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

BoutMesh::BoutMesh(GridDataSource *s, Options *opt) : Mesh(s, opt) {
  OPTION(options, symmetricGlobalX, true);
  if (!options->isSet("symmetricGlobalY")) {
    std::string optionfile = Options::root()["optionfile"].withDefault("");
    output_warn << "WARNING: The default of this option has changed in release 4.1.\n\
If you want the old setting, you have to specify mesh:symmetricGlobalY=false in "
                << optionfile << "\n";
  }
  OPTION(options, symmetricGlobalY, true);

  comm_x = MPI_COMM_NULL;
  comm_inner = MPI_COMM_NULL;
  comm_middle = MPI_COMM_NULL;
  comm_outer = MPI_COMM_NULL;

  mpi = bout::globals::mpi;
}

BoutMesh::~BoutMesh() {
  // Delete the communication handles
  clear_handles();

  // Delete the boundary regions
  for (const auto& bndry : boundary) {
    delete bndry;
  }
  for (const auto& bndry : par_boundary) {
    delete bndry;
  }

  if (comm_x != MPI_COMM_NULL) {
    MPI_Comm_free(&comm_x);
  }
  if (comm_inner != MPI_COMM_NULL) {
    MPI_Comm_free(&comm_inner);
  }
  if (comm_outer != MPI_COMM_NULL) {
    MPI_Comm_free(&comm_outer);
  }
}

BoutMesh::YDecompositionIndices
BoutMesh::setYDecompositionIndices(const BoutMesh::YDecompositionIndices& indices) {
  setYDecompositionIndices(indices.jyseps1_1, indices.jyseps2_1, indices.jyseps1_2,
                           indices.jyseps2_2, indices.ny_inner);

  return {jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2, ny_inner};
}

void BoutMesh::setYDecompositionIndices(int jyseps1_1_, int jyseps2_1_, int jyseps1_2_,
                                        int jyseps2_2_, int ny_inner_) {
  // Set member variables
  this->jyseps1_1 = jyseps1_1_;
  this->jyseps2_1 = jyseps2_1_;
  this->jyseps1_2 = jyseps1_2_;
  this->jyseps2_2 = jyseps2_2_;
  this->ny_inner = ny_inner_;

  // Check inputs
  if (jyseps1_1 < -1) {
    output_warn.write("\tWARNING: jyseps1_1 ({:d}) must be >= -1. Setting to -1\n",
                      jyseps1_1);
    jyseps1_1 = -1;
  }

  if (jyseps2_1 < jyseps1_1) {
    output_warn.write(
        "\tWARNING: jyseps2_1 ({:d}) must be >= jyseps1_1 ({:d}). Setting to {:d}\n",
        jyseps2_1, jyseps1_1, jyseps1_1 + 1);
    jyseps2_1 = jyseps1_1 + 1;
  }
  if (jyseps1_2 < jyseps2_1) {
    output_warn.write(
        "\tWARNING: jyseps1_2 ({:d}) must be >= jyseps2_1 ({:d}). Setting to {:d}\n",
        jyseps1_2, jyseps2_1, jyseps2_1);
    jyseps1_2 = jyseps2_1;
  }
  if (jyseps2_2 >= ny) {
    output_warn.write(
        "\tWARNING: jyseps2_2 ({:d}) must be < ny ({:d}). Setting to {:d}\n", jyseps2_2,
        ny, ny - 1);
    jyseps2_2 = ny - 1;
  }
  if (jyseps2_2 < jyseps1_2) {
    if (jyseps1_2 >= ny) {
      throw BoutException("jyseps1_2 ({:d}) must be < ny ({:d}).", jyseps1_2, ny);
    }
    output_warn.write(
        "\tWARNING: jyseps2_2 ({:d}) must be >= jyseps1_2 ({:d}). Setting to {:d}\n",
        jyseps2_2, jyseps1_2, jyseps1_2);
    jyseps2_2 = jyseps1_2;
  }

  if (jyseps1_1 < 0 and jyseps2_2 >= ny - 1) {
    numberOfXPoints = 0;
  } else if (jyseps2_1 == jyseps1_2) {
    numberOfXPoints = 1;
  } else {
    numberOfXPoints = 2;
  }
}

void BoutMesh::setXDecompositionIndices(const XDecompositionIndices& indices) {
  ixseps1 = indices.ixseps1;
  ixseps2 = indices.ixseps2;
}

namespace bout {
CheckMeshResult checkBoutMeshYDecomposition(int total_processors, int num_y_processors,
                                            int ny, int num_y_guards, int jyseps1_1,
                                            int jyseps2_1, int jyseps1_2, int jyseps2_2,
                                            int ny_inner) {

  const int num_local_y_points = ny / num_y_processors;

  // Check size of Y mesh if we've got multiple processors
  if (num_local_y_points < num_y_guards and total_processors != 1) {
    return {false,
            fmt::format(_("\t -> ny/NYPE ({:d}/{:d} = {:d}) must be >= MYG ({:d})\n"), ny,
                        num_y_processors, num_local_y_points, num_y_guards)};
  }
  // Check branch cuts
  if ((jyseps1_1 + 1) % num_local_y_points != 0) {
    return {false, fmt::format(_("\t -> Leg region jyseps1_1+1 ({:d}) must be a "
                                 "multiple of MYSUB ({:d})\n"),
                               jyseps1_1 + 1, num_local_y_points)};
  }

  if (jyseps2_1 != jyseps1_2) {
    // Double Null

    if ((jyseps2_1 - jyseps1_1) % num_local_y_points != 0) {
      return {
          false,
          fmt::format(_("\t -> Core region jyseps2_1-jyseps1_1 ({:d}-{:d} = {:d}) must "
                        "be a multiple of MYSUB ({:d})\n"),
                      jyseps2_1, jyseps1_1, jyseps2_1 - jyseps1_1, num_local_y_points)};
    }

    if ((jyseps2_2 - jyseps1_2) % num_local_y_points != 0) {
      return {
          false,
          fmt::format(_("\t -> Core region jyseps2_2-jyseps1_2 ({:d}-{:d} = {:d}) must "
                        "be a multiple of MYSUB ({:d})\n"),
                      jyseps2_2, jyseps1_2, jyseps2_2 - jyseps1_2, num_local_y_points)};
    }

    // Check upper legs
    if ((ny_inner - jyseps2_1 - 1) % num_local_y_points != 0) {
      return {
          false,
          fmt::format(_("\t -> leg region ny_inner-jyseps2_1-1 ({:d}-{:d}-1 = {:d}) must "
                        "be a multiple of MYSUB ({:d})\n"),
                      ny_inner, jyseps2_1, ny_inner - jyseps2_1 - 1, num_local_y_points)};
    }
    if ((jyseps1_2 - ny_inner + 1) % num_local_y_points != 0) {
      return {
          false,
          fmt::format(_("\t -> leg region jyseps1_2-ny_inner+1 ({:d}-{:d}+1 = {:d}) must "
                        "be a multiple of MYSUB ({:d})\n"),
                      jyseps1_2, ny_inner, jyseps1_2 - ny_inner + 1, num_local_y_points)};
    }
  } else {
    // Single Null
    if ((jyseps2_2 - jyseps1_1) % num_local_y_points != 0) {
      return {
          false,
          fmt::format(_("\t -> Core region jyseps2_2-jyseps1_1 ({:d}-{:d} = {:d}) must "
                        "be a multiple of MYSUB ({:d})\n"),
                      jyseps2_2, jyseps1_1, jyseps2_2 - jyseps1_1, num_local_y_points)};
    }
  }

  if ((ny - jyseps2_2 - 1) % num_local_y_points != 0) {
    return {false, fmt::format(
                       _("\t -> leg region ny-jyseps2_2-1 ({:d}-{:d}-1 = {:d}) must be a "
                         "multiple of MYSUB ({:d})\n"),
                       ny, jyseps2_2, ny - jyseps2_2 - 1, num_local_y_points)};
  }

  return {true, ""};
}
} // namespace bout

void BoutMesh::chooseProcessorSplit(Options& options) {
  // Possible issues:
  // - can set both NXPE and NYPE, but only NXPE is used
  // - can set NXPE > nx
  // - can set NYPE > ny (if only one processor)

  if (options.isSet("NXPE")) {
    NXPE = options["NXPE"]
               .doc("Decomposition in the radial direction. If not given then calculated "
                    "automatically.")
               .withDefault(1);
    if ((NPES % NXPE) != 0) {
      throw BoutException(
          _("Number of processors ({:d}) not divisible by NPs in x direction ({:d})\n"),
          NPES, NXPE);
    }

    NYPE = NPES / NXPE;
  } else {
    // NXPE not set, but NYPE is
    NYPE = options["NYPE"]
               .doc("Decomposition in the parallel direction. Can be given instead of "
                    "NXPE. If neither is given, then calculated automatically.")
               .withDefault(1);
    if ((NPES % NYPE) != 0) {
      throw BoutException(
          _("Number of processors ({:d}) not divisible by NPs in y direction ({:d})\n"),
          NPES, NYPE);
    }

    NXPE = NPES / NYPE;
  }

  auto result = bout::checkBoutMeshYDecomposition(
      NPES, NYPE, ny, MYG, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2, ny_inner);

  if (not result.success) {
    throw BoutException(result.reason);
  }
}

void BoutMesh::findProcessorSplit() {
  MX = nx - 2 * MXG;

  NXPE = -1; // Best option

  // Results in square domains
  const BoutReal ideal = sqrt(MX * NPES / static_cast<BoutReal>(ny));

  output_info.write(_("Finding value for NXPE (ideal = {:f})\n"), ideal);

  for (int i = 1; i <= NPES; i++) { // Loop over all possibilities
    if ((NPES % i == 0) &&          // Processors divide equally
        (MX % i == 0) &&            // Mesh in X divides equally
        (ny % (NPES / i) == 0)) {   // Mesh in Y divides equally

      output_info.write(_("\tCandidate value: {:d}\n"), i);

      const int nyp = NPES / i;

      auto result = bout::checkBoutMeshYDecomposition(
          NPES, nyp, ny, MYG, jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2, ny_inner);

      if (not result.success) {
        output_info.write(result.reason);
        continue;
      }

      output_info.write(_("\t -> Good value\n"));
      // Found an acceptable value
      if ((NXPE < 1) || (fabs(ideal - i) < fabs(ideal - NXPE))) {
        NXPE = i; // Keep value nearest to the ideal
      }
    }
  }

  if (NXPE < 1) {
    throw BoutException(_("Could not find a valid value for NXPE. Try a different "
                          "number of processors."));
  }

  NYPE = NPES / NXPE;

  output_progress.write(_("\tDomain split (NXPE={:d}, NYPE={:d}) into domains "
                          "(localNx={:d}, localNy={:d})\n"),
                        NXPE, NYPE, MX / NXPE, ny / NYPE);
}

void BoutMesh::setDerivedGridSizes() {
  // Check that nx is large enough
  if (nx <= 2 * MXG) {
    throw BoutException(_("Error: nx must be greater than 2 times MXG (2 * {:d})"), MXG);
  }

  GlobalNx = nx;
  GlobalNy = ny + 2 * MYG;
  GlobalNz = nz;

  // If we've got a second pair of diverator legs, we need an extra
  // pair of boundary regions
  if (jyseps1_2 != jyseps2_1) {
    GlobalNy += 2 * MYG;
  }

  // Set global grid sizes, excluding boundary points
  GlobalNxNoBoundaries = nx - 2*MXG;
  GlobalNyNoBoundaries = ny;
  GlobalNzNoBoundaries = nz;

  // Split MX points between NXPE processors
  // MXG at each end needed for edge boundary regions
  MX = nx - 2 * MXG;
  MXSUB = MX / NXPE;
  if ((MX % NXPE) != 0) {
    throw BoutException(_("Cannot split {:d} X points equally between {:d} processors\n"),
                        MX, NXPE);
  }

  // NOTE: No grid data reserved for Y boundary cells - copy from neighbours
  MY = ny;
  MYSUB = MY / NYPE;
  if ((MY % NYPE) != 0) {
    throw BoutException(
        _("\tERROR: Cannot split {:d} Y points equally between {:d} processors\n"), MY,
        NYPE);
  }

  MZ = nz;
  MZSUB = MZ / NZPE;
  if ((MZ % NZPE) != 0) {
    throw BoutException(
        _("\tERROR: Cannot split {:d} Z points equally between {:d} processors\n"), MZ,
        NZPE);
  }

  // Set global offsets
  OffsetX = PE_XIND * MXSUB;
  OffsetY = PE_YIND * MYSUB;
  OffsetZ = 0;

  // Number of grid cells on this processor is ng* = M*SUB + guard/boundary cells
  LocalNx = MXSUB + 2 * MXG;
  LocalNy = MYSUB + 2 * MYG;
  LocalNz = MZSUB + 2 * MZG;

  // Set local index ranges
  xstart = MXG;
  xend = MXG + MXSUB - 1;

  ystart = MYG;
  yend = MYG + MYSUB - 1;

  zstart = MZG;
  zend = MZG + MZSUB - 1;
}

int BoutMesh::load() {
  TRACE("BoutMesh::load()");

  output_progress << _("Loading mesh") << endl;

  // Use root level options
  auto& options = Options::root();

  //////////////
  // Number of processors

  MPI_Comm_size(BoutComm::get(), &NPES);
  MPI_Comm_rank(BoutComm::get(), &MYPE);

  //////////////
  // Grid sizes

  if (Mesh::get(nx, "nx") != 0) {
    throw BoutException(_("Mesh must contain nx"));
  }

  if (Mesh::get(ny, "ny") != 0) {
    throw BoutException(_("Mesh must contain ny"));
  }

  if (Mesh::get(nz, "nz") != 0) {
    // No "nz" variable in the grid file. Instead read MZ from options

    OPTION(options, MZ, 64);
    OPTION(options, nz, MZ);
    ASSERT0(nz == MZ);
    if (!is_pow2(nz)) {
      // Should be a power of 2 for efficient FFTs
      output_warn.write(
          _("WARNING: Number of toroidal points should be 2^n for efficient "
            "FFT performance -- consider changing MZ ({:d}) if using FFTs\n"),
          nz);
    }
  } else {
    MZ = nz;
    output_info.write(_("\tRead nz from input grid file\n"));
  }

  output_info << _("\tGrid size: ") << nx << " x " << ny << " x " << nz << endl;

  // Get guard cell sizes
  // Try to read from grid file first, then if not found
  // get from options
  if (Mesh::get(MXG, "MXG") != 0) {
    // Error code returned
    MXG = options["MXG"].doc("Number of guard cells on each side in X").withDefault(2);
  }
  ASSERT0(MXG >= 0);

  if (Mesh::get(MYG, "MYG") != 0) {
    MYG = options["MYG"].doc("Number of guard cells on each side in Y").withDefault(2);
  }
  ASSERT0(MYG >= 0);

  // For now only support no z-guard cells
  MZG = 0;
  ASSERT0(MZG >= 0);

  // For now don't parallelise z
  NZPE = 1;

  output_info << _("\tGuard cells (x,y,z): ") << MXG << ", " << MYG << ", " << MZG
              << std::endl;

  // separatrix location
  Mesh::get(ixseps1, "ixseps1", nx);
  Mesh::get(ixseps2, "ixseps2", nx);
  Mesh::get(jyseps1_1, "jyseps1_1", -1);
  Mesh::get(jyseps1_2, "jyseps1_2", ny / 2);
  Mesh::get(jyseps2_1, "jyseps2_1", jyseps1_2);
  Mesh::get(jyseps2_2, "jyseps2_2", ny - 1);
  Mesh::get(ny_inner, "ny_inner", jyseps2_1);

  // Check inputs
  setYDecompositionIndices(jyseps1_1, jyseps2_1, jyseps1_2, jyseps2_2, ny_inner);

  if (options.isSet("NXPE") or options.isSet("NYPE")) {
    chooseProcessorSplit(options);
  } else {
    findProcessorSplit();
  }

  // Get X and Y processor indices
  PE_YIND = MYPE / NXPE;
  PE_XIND = MYPE % NXPE;

  // Set the other grid sizes from nx, ny, nz
  setDerivedGridSizes();

  /// Get mesh options
  OPTION(options, IncIntShear, false);
  periodicX = options["periodicX"].doc("Make grid periodic in X?").withDefault(false);

  async_send = options["async_send"]
                   .doc("Whether to use asyncronous MPI sends")
                   .withDefault(false);

  if (options.isSet("zperiod")) {
    OPTION(options, zperiod, 1);
    ZMIN = 0.0;
    ZMAX = 1.0 / static_cast<BoutReal>(zperiod);
  } else {
    OPTION(options, ZMIN, 0.0);
    OPTION(options, ZMAX, 1.0);

    zperiod = ROUND(1.0 / (ZMAX - ZMIN));
  }

  ///////////////////// TOPOLOGY //////////////////////////
  /// Call topology to set layout of grid
  topology();

  TwistShift = options["TwistShift"]
                   .doc("Apply a Twist-Shift boundary using ShiftAngle?")
                   .withDefault(false);

  // Try to read the shift angle from the grid file
  // NOTE: All processors should know the twist-shift angle (for invert_parderiv)
  // NOTE: Always read ShiftAngle as Coordinates will use hasBranchCutLower and
  //       hasBranchCutUpper to set zShift for ShiftedMetric

  ShiftAngle.resize(LocalNx);

  if (!source->get(this, ShiftAngle, "ShiftAngle", LocalNx, getGlobalXIndex(0))) {
    ShiftAngle.clear();
  }

  if (TwistShift) {
    output_info.write("Applying Twist-Shift condition. Interpolation: FFT\n");
    if (ShiftAngle.empty()) {
      throw BoutException("ERROR: Twist-shift angle 'ShiftAngle' not found. "
          "Required when TwistShift==true.");
    }
  }

  //////////////////////////////////////////////////////
  /// Communicator

  createCommunicators();
  output_debug << "Got communicators" << endl;

  //////////////////////////////////////////////////////
  // Boundary regions
  createXBoundaries();
  createYBoundaries();

  auto possible_boundaries = getPossibleBoundaries();
  if (possible_boundaries.empty()) {
    output_info.write(_("No boundary regions; domain is periodic\n"));
  } else {
    output_info.write(_("Possible boundary regions are: "));
    for (const auto& boundary : possible_boundaries) {
      output_info.write("{}, ", boundary);
    }
  }

  if (!boundary.empty()) {
    output_info << _("Boundary regions in this processor: ");
    for (const auto &bndry : boundary) {
      output_info << bndry->label << ", ";
    }
    output_info << endl;
  } else {
    output_info << _("No boundary regions in this processor") << endl;
  }

  output_info << _("Constructing default regions") << endl;
  createDefaultRegions();

  // Add boundary regions
  addBoundaryRegions();

  // Initialize default coordinates
  getCoordinates();

  output_info.write(_("\tdone\n"));

  return 0;
}

void BoutMesh::createCommunicators() {
  MPI_Group group_world{};
  MPI_Comm_group(BoutComm::get(), &group_world); // Get the entire group

  //////////////////////////////////////////////////////
  /// Communicator in X

  MPI_Group group{};
  MPI_Comm comm_tmp{};

  int proc[3]; // Processor range

  for (int yp = 0; yp < NYPE; yp++) {
    proc[0] = PROC_NUM(0, yp);        // First
    proc[1] = PROC_NUM(NXPE - 1, yp); // Last
    proc[2] = 1;                      // stride

    output_debug << "XCOMM " << proc[0] << ", " << proc[1] << endl;

    if (MPI_Group_range_incl(group_world, 1, &proc, &group) != MPI_SUCCESS) {
      throw BoutException(
          "Could not create X communication group for yp={:d} (xind={:d},yind={:d})\n",
          yp, PE_XIND, PE_YIND);
    }
    if (MPI_Comm_create(BoutComm::get(), group, &comm_tmp) != MPI_SUCCESS) {
      throw BoutException(
          "Could not create X communicator for yp={:d} (xind={:d},yind={:d})\n", yp,
          PE_XIND, PE_YIND);
    }
    MPI_Group_free(&group);

    if (yp == PE_YIND) {
      // Should be in this group
      if (comm_tmp == MPI_COMM_NULL) {
        throw BoutException("X communicator null");
      }

      comm_x = comm_tmp;
    } else {
      if (comm_tmp != MPI_COMM_NULL) {
        throw BoutException("X communicator should be null");
      }
    }
  }

  //////////////////////////////////////////////////////
  /// Communicators for Y gather/scatter

  MPI_Group group_tmp1{};
  MPI_Group group_tmp2{};

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

      if (MPI_Group_range_incl(group_world, 1, &proc, &group) != MPI_SUCCESS) {
        throw BoutException("MPI_Group_range_incl failed for xp = {:d}", NXPE);
      }
      MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
      if (comm_tmp != MPI_COMM_NULL) {
        comm_outer = comm_tmp;
      }
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
      TRACE("Creating lower PF communicators for xp={:d}", i);

      output_debug << "Creating lower PF communicators for xp = " << i << endl;

      if (jyseps1_1 >= 0) {
        proc[0] = PROC_NUM(i, 0);
        proc[1] = PROC_NUM(i, YPROC(jyseps1_1));

        output_debug << "PF1 " << proc[0] << ", " << proc[1] << endl;

        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
      } else {
        group_tmp1 = MPI_GROUP_EMPTY;
      }

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
      TRACE("Creating upper PF communicators for xp={:d}", i);

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
      if (group_tmp1 != MPI_GROUP_EMPTY) {
        MPI_Group_free(&group_tmp1);
      }
      if (group_tmp2 != MPI_GROUP_EMPTY) {
        MPI_Group_free(&group_tmp2);
      }

      output_debug << "done upper PF\n";
    }

    // Core region
    TRACE("Creating core communicators");
    if (jyseps2_1 > jyseps1_1) {
      proc[0] = PROC_NUM(i, YPROC(jyseps1_1 + 1));
      proc[1] = PROC_NUM(i, YPROC(jyseps2_1));

      output_debug << "CORE1 " << proc[0] << ", " << proc[1] << endl;

      if ((proc[0] < 0) || (proc[1] < 0)) {
        throw BoutException("Invalid processor range for core processors");
      }
      MPI_Group_range_incl(group_world, 1, &proc, &group_tmp1);
    } else {
      // no core region between jyseps1_1 and jyseps2_1
      group_tmp1 = MPI_GROUP_EMPTY;
    }

    if (jyseps2_2 > jyseps1_2) {
      proc[0] = PROC_NUM(i, YPROC(jyseps1_2 + 1));
      proc[1] = PROC_NUM(i, YPROC(jyseps2_2));

      output_debug << "CORE2 " << proc[0] << ", " << proc[1] << endl;

      if ((proc[0] < 0) || (proc[1] < 0)) {
        group_tmp2 = MPI_GROUP_EMPTY;
      } else {
        MPI_Group_range_incl(group_world, 1, &proc, &group_tmp2);
      }
    } else {
      // no core region between jyseps1_2 and jyseps2_2
      group_tmp2 = MPI_GROUP_EMPTY;
    }

    MPI_Group_union(group_tmp1, group_tmp2, &group);
    MPI_Comm_create(BoutComm::get(), group, &comm_tmp);
    if (comm_tmp != MPI_COMM_NULL) {
      comm_inner = comm_tmp;

      if (ixseps_inner == ixseps_outer) {
        MPI_Comm_dup(comm_inner, &comm_middle);
      }
    }

    if (group_tmp1 != MPI_GROUP_EMPTY) {
      MPI_Group_free(&group_tmp1);
    }
    if (group_tmp2 != MPI_GROUP_EMPTY) {
      MPI_Group_free(&group_tmp2);
    }
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
        if (comm_tmp != MPI_COMM_NULL) {
          comm_middle = comm_tmp;
        }

        if (group_tmp1 != MPI_GROUP_EMPTY) {
          MPI_Group_free(&group_tmp1);
        }
        if (group_tmp2 != MPI_GROUP_EMPTY) {
          MPI_Group_free(&group_tmp2);
        }
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
        if (comm_tmp != MPI_COMM_NULL) {
          comm_middle = comm_tmp;
        }

        if (group_tmp1 != MPI_GROUP_EMPTY) {
          MPI_Group_free(&group_tmp1);
        }
        if (group_tmp2 != MPI_GROUP_EMPTY) {
          MPI_Group_free(&group_tmp2);
        }
        MPI_Group_free(&group);
      }
    }
  }
  MPI_Group_free(&group_world);
  // Now have communicators for all regions.
}

void BoutMesh::createXBoundaries() {
  // Need boundaries in X if not periodic and have X guard cells
  if (periodicX) {
    return;
  }
  if (MXG <= 0) {
    return;
  }

  if (PE_XIND == 0) {
    // Inner either core or PF

    // Get a global index in this processor
    const int yg = getGlobalYIndexNoBoundaries(MYG);

    if (((yg > jyseps1_1) and (yg <= jyseps2_1))
        or ((yg > jyseps1_2) and (yg <= jyseps2_2))) {
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

void BoutMesh::createYBoundaries() {
  if (MYG <= 0) {
    return;
  }
  // Need boundaries in Y

  // Alter x-limits so that y-boundary conditions set corner-boundary cells
  // i.e. if there is an x-boundary, include corner cells. If
  // include_corner_cells==false, this modification is disabled to match the behaviour
  // of BOUT++ up to v4.
  // Note that including the corner cells requires that the x-boundary conditions are
  // applied before the y-boundary conditions. This is ensured here in the
  // BOUT++-applied boundary conditions because the y-boundaries are added to the
  // 'boundary' vector after the x-boundaries, but beware **THE ORDER IS IMPORTANT**
  const int yboundary_xstart = (include_corner_cells and IDATA_DEST == -1) ? 0 : xstart;
  const int yboundary_xend =
      (include_corner_cells and ODATA_DEST == -1) ? LocalNx - 1 : xend;

  if ((UDATA_INDEST < 0) && (UDATA_XSPLIT > yboundary_xstart)) {
    boundary.push_back(
        new BoundaryRegionYUp("upper_target", yboundary_xstart, UDATA_XSPLIT - 1, this));
  }
  if ((UDATA_OUTDEST < 0) && (UDATA_XSPLIT <= yboundary_xend)) {
    boundary.push_back(
        new BoundaryRegionYUp("upper_target", UDATA_XSPLIT, yboundary_xend, this));
  }

  if ((DDATA_INDEST < 0) && (DDATA_XSPLIT > yboundary_xstart)) {
    boundary.push_back(new BoundaryRegionYDown("lower_target", yboundary_xstart,
                                               DDATA_XSPLIT - 1, this));
  }
  if ((DDATA_OUTDEST < 0) && (DDATA_XSPLIT <= yboundary_xend)) {
    boundary.push_back(
        new BoundaryRegionYDown("lower_target", DDATA_XSPLIT, yboundary_xend, this));
  }
}

std::set<std::string> BoutMesh::getPossibleBoundaries() const {
  // Result set: set so it automatically takes care of duplicates
  std::set<std::string> all_boundaries{};

  // Lambda that modifies `all_boundaries`
  const auto get_boundaries_on_different_rank = [mesh = this, &all_boundaries](
                                                    int x_rank, int y_rank) {
    // Don't try to check boundaries on unphysical processors
    if (x_rank < 0 or x_rank >= mesh->NXPE) {
      return;
    }
    if (y_rank < 0 or y_rank >= mesh->NYPE) {
      return;
    }

    // Make a copy of this mesh, EXCEPT we change the (X, Y) rank of the processor
    BoutMesh mesh_copy{mesh->GlobalNx, mesh->GlobalNyNoBoundaries, mesh->GlobalNz,
                       mesh->MXG, mesh->MYG, mesh->NXPE, mesh->NYPE, x_rank, y_rank,
                       mesh->symmetricGlobalX, mesh->symmetricGlobalY, mesh->periodicX,
                       mesh->ixseps1, mesh->ixseps2, mesh->jyseps1_1, mesh->jyseps2_1,
                       mesh->jyseps1_2, mesh->jyseps2_2, mesh->ny_inner};
    // We need to create the boundaries
    mesh_copy.createXBoundaries();
    mesh_copy.createYBoundaries();

    // Get the boundaries and shove their names into the set
    auto boundaries = mesh_copy.getBoundaries();
    std::transform(boundaries.begin(), boundaries.end(),
                   std::inserter(all_boundaries, all_boundaries.begin()),
                   [](BoundaryRegionBase* boundary) { return boundary->label; });
  };

  // This is sufficient to get the SOL boundary, if it exists
  get_boundaries_on_different_rank(NXPE - 1, 0);

  // Now we need to work out which Y rank at each of the possible
  // targets and branch cuts. This makes sure we get a processor on
  // each possible boundary. We only need to do this at PE_XIND==0, as
  // the other X boundary is covered by the case above
  for (const auto& y_index :
       {0, jyseps1_1, jyseps1_2, ny_inner - 1, ny_inner, jyseps2_1, jyseps2_2, ny - 1}) {
    get_boundaries_on_different_rank(0, YPROC(y_index));
  }

  return all_boundaries;
}

void BoutMesh::setShiftAngle(const std::vector<BoutReal>& shift_angle) {
  if (shift_angle.size() != static_cast<std::size_t>(LocalNx)) {
    throw BoutException("shift_angle vector wrong size: got {}, expected {}",
                        shift_angle.size(), LocalNx);
  }
  TwistShift = true;
  ShiftAngle = shift_angle;
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

void BoutMesh::post_receiveX(CommHandle &ch) {
  /// Post receive data from left (x-1)

  if (IDATA_DEST != -1) {
    mpi->MPI_Irecv(std::begin(ch.imsg_recvbuff),
                   msg_len(ch.var_list.get(), 0, MXG, 0,
                           ch.include_x_corners ? LocalNy : MYSUB),
                   PVEC_REAL_MPI_TYPE, IDATA_DEST, OUT_SENT_IN, BoutComm::get(),
                   &ch.request[4]);
  }

  // Post receive data from right (x+1)

  if (ODATA_DEST != -1) {
    mpi->MPI_Irecv(std::begin(ch.omsg_recvbuff),
                   msg_len(ch.var_list.get(), 0, MXG, 0,
                           ch.include_x_corners ? LocalNy : MYSUB),
                   PVEC_REAL_MPI_TYPE, ODATA_DEST, IN_SENT_OUT, BoutComm::get(),
                   &ch.request[5]);
  }
}

void BoutMesh::post_receiveY(CommHandle &ch) {
  BoutReal *inbuff;
  int len;

  /// Post receive data from above (y+1)

  len = 0;
  if (UDATA_INDEST != -1) {
    len = msg_len(ch.var_list.get(), 0, UDATA_XSPLIT, 0, MYG);
    mpi->MPI_Irecv(std::begin(ch.umsg_recvbuff), len, PVEC_REAL_MPI_TYPE, UDATA_INDEST,
                   IN_SENT_DOWN, BoutComm::get(), &ch.request[0]);
  }
  if (UDATA_OUTDEST != -1) {
    inbuff = &ch.umsg_recvbuff[len]; // pointer to second half of the buffer
    mpi->MPI_Irecv(inbuff, msg_len(ch.var_list.get(), UDATA_XSPLIT, LocalNx, 0, MYG),
                   PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, OUT_SENT_DOWN, BoutComm::get(),
                   &ch.request[1]);
  }

  /// Post receive data from below (y-1)

  len = 0;

  if (DDATA_INDEST != -1) { // If sending & recieving data from a processor
    len = msg_len(ch.var_list.get(), 0, DDATA_XSPLIT, 0, MYG);
    mpi->MPI_Irecv(std::begin(ch.dmsg_recvbuff), len, PVEC_REAL_MPI_TYPE, DDATA_INDEST,
                   IN_SENT_UP, BoutComm::get(), &ch.request[2]);
  }
  if (DDATA_OUTDEST != -1) {
    inbuff = &ch.dmsg_recvbuff[len];
    mpi->MPI_Irecv(inbuff, msg_len(ch.var_list.get(), DDATA_XSPLIT, LocalNx, 0, MYG),
                   PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, OUT_SENT_UP, BoutComm::get(),
                   &ch.request[3]);
  }
}

comm_handle BoutMesh::send(FieldGroup &g) {
  /// Start timer
  Timer timer("comms");

  if (include_corner_cells) {
    throw BoutException("Cannot use send() when include_corner_cells==true as it sends "
                        "in x- and y-directions simultaneously. Use sendX() and sendY() "
                        "instead.");
  }

  /// Work out length of buffer needed
  int xlen = msg_len(g.get(), 0, MXG, 0, MYSUB);
  int ylen = msg_len(g.get(), 0, LocalNx, 0, MYG);

  /// Get a communications handle of (at least) the needed size
  CommHandle *ch = get_handle(xlen, ylen);
  ch->var_list = g; // Group of fields to send

  sendX(g, ch, true);
  sendY(g, ch);

  return static_cast<void *>(ch);
}

comm_handle BoutMesh::sendX(FieldGroup &g, comm_handle handle, bool disable_corners) {
  /// Start timer
  Timer timer("comms");

  const bool with_corners = include_corner_cells and not disable_corners;

  CommHandle* ch = nullptr;
  if (handle == nullptr) {
    /// Work out length of buffer needed
    int xlen = msg_len(g.get(), 0, MXG, 0, with_corners ? LocalNy : MYSUB);

    /// Get a communications handle of (at least) the needed size
    ch = get_handle(xlen, 0);
    ch->var_list = g; // Group of fields to send
  } else {
    ch = static_cast<CommHandle*>(handle);
  }

  ch->include_x_corners = with_corners;

  /// Post receives
  post_receiveX(*ch);

  //////////////////////////////////////////////////

  /// Send to the left (x-1)

  if (IDATA_DEST != -1) {
    int len = pack_data(ch->var_list.get(), MXG, 2 * MXG, ch->include_x_corners ? 0 : MYG,
                        ch->include_x_corners ? LocalNy : MYG + MYSUB,
                        std::begin(ch->imsg_sendbuff));
    if (async_send) {
      mpi->MPI_Isend(std::begin(ch->imsg_sendbuff), len, PVEC_REAL_MPI_TYPE, IDATA_DEST,
                     IN_SENT_OUT, BoutComm::get(), &(ch->sendreq[4]));
    } else {
      mpi->MPI_Send(std::begin(ch->imsg_sendbuff), len, PVEC_REAL_MPI_TYPE, IDATA_DEST,
                    IN_SENT_OUT, BoutComm::get());
    }
  }

  /// Send to the right (x+1)

  if (ODATA_DEST != -1) {
    int len = pack_data(ch->var_list.get(), MXSUB, MXSUB + MXG, ch->include_x_corners ? 0 :
                        MYG, ch->include_x_corners ? LocalNy : MYG + MYSUB,
                        std::begin(ch->omsg_sendbuff));
    if (async_send) {
      mpi->MPI_Isend(std::begin(ch->omsg_sendbuff), len, PVEC_REAL_MPI_TYPE, ODATA_DEST,
                     OUT_SENT_IN, BoutComm::get(), &(ch->sendreq[5]));
    } else {
      mpi->MPI_Send(std::begin(ch->omsg_sendbuff), len, PVEC_REAL_MPI_TYPE, ODATA_DEST,
                    OUT_SENT_IN, BoutComm::get());
    }
  }

  /// Mark communication handle as in progress
  ch->in_progress = true;

  return static_cast<void *>(ch);
}

comm_handle BoutMesh::sendY(FieldGroup &g, comm_handle handle) {
  /// Start timer
  Timer timer("comms");

  CommHandle* ch;
  if (handle == nullptr) {
    /// Work out length of buffer needed
    int ylen = msg_len(g.get(), 0, LocalNx, 0, MYG);

    /// Get a communications handle of (at least) the needed size
    ch = get_handle(0, ylen);
    ch->var_list = g; // Group of fields to send
  } else {
    ch = static_cast<CommHandle*>(handle);
  }

  /// Post receives
  post_receiveY(*ch);

  //////////////////////////////////////////////////

  /// Send data going up (y+1)

  int len = 0;
  BoutReal* outbuff = nullptr;

  if (UDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(ch->var_list.get(), 0, UDATA_XSPLIT, MYSUB, MYSUB + MYG,
                    std::begin(ch->umsg_sendbuff));
    // Send the data to processor UDATA_INDEST

    if (async_send) {
      mpi->MPI_Isend(std::begin(ch->umsg_sendbuff), // Buffer to send
                     len,                           // Length of buffer in BoutReals
                     PVEC_REAL_MPI_TYPE,            // Real variable type
                     UDATA_INDEST,                  // Destination processor
                     IN_SENT_UP,                    // Label (tag) for the message
                     BoutComm::get(), &(ch->sendreq[0]));
    } else {
      mpi->MPI_Send(std::begin(ch->umsg_sendbuff), len, PVEC_REAL_MPI_TYPE, UDATA_INDEST,
                    IN_SENT_UP, BoutComm::get());
    }
  }
  if (UDATA_OUTDEST != -1) {             // if destination for outer x data
    outbuff = &(ch->umsg_sendbuff[len]); // A pointer to the start of the second part
                                         // of the buffer
    len =
        pack_data(ch->var_list.get(), UDATA_XSPLIT, LocalNx, MYSUB, MYSUB + MYG, outbuff);
    // Send the data to processor UDATA_OUTDEST
    if (async_send) {
      mpi->MPI_Isend(outbuff, len, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, OUT_SENT_UP,
                     BoutComm::get(), &(ch->sendreq[1]));
    } else {
      mpi->MPI_Send(outbuff, len, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, OUT_SENT_UP,
                    BoutComm::get());
    }
  }

  /// Send data going down (y-1)

  len = 0;
  if (DDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(ch->var_list.get(), 0, DDATA_XSPLIT, MYG, 2 * MYG,
                    std::begin(ch->dmsg_sendbuff));
    // Send the data to processor DDATA_INDEST
    if (async_send) {
      mpi->MPI_Isend(std::begin(ch->dmsg_sendbuff), len, PVEC_REAL_MPI_TYPE, DDATA_INDEST,
                     IN_SENT_DOWN, BoutComm::get(), &(ch->sendreq[2]));
    } else {
      mpi->MPI_Send(std::begin(ch->dmsg_sendbuff), len, PVEC_REAL_MPI_TYPE, DDATA_INDEST,
                    IN_SENT_DOWN, BoutComm::get());
    }
  }
  if (DDATA_OUTDEST != -1) {             // if destination for outer x data
    outbuff = &(ch->dmsg_sendbuff[len]); // A pointer to the start of the second part
                                         // of the buffer
    len = pack_data(ch->var_list.get(), DDATA_XSPLIT, LocalNx, MYG, 2 * MYG, outbuff);
    // Send the data to processor DDATA_OUTDEST

    if (async_send) {
      mpi->MPI_Isend(outbuff, len, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, OUT_SENT_DOWN,
                     BoutComm::get(), &(ch->sendreq[3]));
    } else {
      mpi->MPI_Send(outbuff, len, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, OUT_SENT_DOWN,
                    BoutComm::get());
    }
  }

  /// Mark communication handle as in progress
  ch->in_progress = true;

  /// Mark as y-communication
  ch->has_y_communication = true;

  return static_cast<void *>(ch);
}

int BoutMesh::wait(comm_handle handle) {
  TRACE("BoutMesh::wait(comm_handle)");

  if (handle == nullptr) {
    return 1;
  }

  auto* ch = static_cast<CommHandle*>(handle);

  if (!ch->in_progress) {
    return 2;
  }

  /// Start timer
  Timer timer("comms");

  ///////////// WAIT FOR DATA //////////////

  int ind, len;
  MPI_Status status;

  if (ch->var_list.empty()) {

    // Just waiting for a single MPI request
    mpi->MPI_Wait(ch->request, &status);
    free_handle(ch);

    return 0;
  }

  do {
    mpi->MPI_Waitany(6, ch->request, &ind, &status);
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
      unpack_data(ch->var_list.get(), 0, MXG, ch->include_x_corners ? 0 : MYG,
                  ch->include_x_corners ? LocalNy : MYG + MYSUB,
                  std::begin(ch->imsg_recvbuff));
      break;
    }
    case 5: { // outer
      unpack_data(ch->var_list.get(), MXSUB + MXG, MXSUB + 2 * MXG,
                  ch->include_x_corners ? 0 : MYG,
                  ch->include_x_corners ? LocalNy : MYG + MYSUB,
                  std::begin(ch->omsg_recvbuff));
      break;
    }
    }
    if (ind != MPI_UNDEFINED) {
      ch->request[ind] = MPI_REQUEST_NULL;
    }
  } while (ind != MPI_UNDEFINED);

  if (async_send) {
    /// Asyncronous sending: Need to check if sends have completed (frees MPI memory)
    MPI_Status async_status;

    if (UDATA_INDEST != -1) {
      mpi->MPI_Wait(ch->sendreq, &async_status);
    }
    if (UDATA_OUTDEST != -1) {
      mpi->MPI_Wait(ch->sendreq + 1, &async_status);
    }
    if (DDATA_INDEST != -1) {
      mpi->MPI_Wait(ch->sendreq + 2, &async_status);
    }
    if (DDATA_OUTDEST != -1) {
      mpi->MPI_Wait(ch->sendreq + 3, &async_status);
    }
    if (IDATA_DEST != -1) {
      mpi->MPI_Wait(ch->sendreq + 4, &async_status);
    }
    if (ODATA_DEST != -1) {
      mpi->MPI_Wait(ch->sendreq + 5, &async_status);
    }
  }

  if (ch->has_y_communication) {
    // TWIST-SHIFT CONDITION
    // Loop over 3D fields
    for (const auto& var : ch->var_list.field3d()) {
      if (var->requiresTwistShift(TwistShift)) {

        // Twist-shift only needed for field-aligned fields
        int jx = 0;
        int jy = 0;

        // Perform Twist-shift using shifting method
        if (var->getDirectionY() == YDirectionType::Aligned) {
          // Only variables in field-aligned coordinates need the twist-shift boundary
          // condition to be applied

          // Lower boundary
          if (TS_down_in && (DDATA_INDEST != -1)) {
            for (jx = 0; jx < DDATA_XSPLIT; jx++) {
              for (jy = 0; jy != MYG; jy++) {
                shiftZ(*var, jx, jy, ShiftAngle[jx]);
              }
            }
          }
          if (TS_down_out && (DDATA_OUTDEST != -1)) {
            for (jx = DDATA_XSPLIT; jx < LocalNx; jx++) {
              for (jy = 0; jy != MYG; jy++) {
                shiftZ(*var, jx, jy, ShiftAngle[jx]);
              }
            }
          }

          // Upper boundary
          if (TS_up_in && (UDATA_INDEST != -1)) {
            for (jx = 0; jx < UDATA_XSPLIT; jx++) {
              for (jy = LocalNy - MYG; jy != LocalNy; jy++) {
                shiftZ(*var, jx, jy, -ShiftAngle[jx]);
              }
            }
          }
          if (TS_up_out && (UDATA_OUTDEST != -1)) {
            for (jx = UDATA_XSPLIT; jx < LocalNx; jx++) {
              for (jy = LocalNy - MYG; jy != LocalNy; jy++) {
                shiftZ(*var, jx, jy, -ShiftAngle[jx]);
              }
            }
          }
        }
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

  MPI_Request request{};

  mpi->MPI_Isend(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(xproc, yproc), tag,
                 BoutComm::get(), &request);

  return request;
}

comm_handle BoutMesh::receiveFromProc(int xproc, int yproc, BoutReal *buffer, int size,
                                      int tag) {
  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(xproc, yproc), tag,
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

bool BoutMesh::firstX() const { return PE_XIND == 0; }

bool BoutMesh::lastX() const { return PE_XIND == NXPE - 1; }

int BoutMesh::sendXOut(BoutReal* buffer, int size, int tag) {
  if (PE_XIND == NXPE - 1) {
    return 1;
  }

  Timer timer("comms");

  mpi->MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND + 1, PE_YIND), tag,
                BoutComm::get());

  return 0;
}

int BoutMesh::sendXIn(BoutReal* buffer, int size, int tag) {
  if (PE_XIND == 0) {
    return 1;
  }

  Timer timer("comms");

  mpi->MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND - 1, PE_YIND), tag,
                BoutComm::get());

  return 0;
}

comm_handle BoutMesh::irecvXOut(BoutReal *buffer, int size, int tag) {
  if (PE_XIND == NXPE - 1)
    return nullptr;

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle *ch = get_handle(0, 0);

  mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND + 1, PE_YIND), tag,
                 BoutComm::get(), ch->request);

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvXIn(BoutReal* buffer, int size, int tag) {
  if (PE_XIND == 0) {
    return nullptr;
  }

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle* ch = get_handle(0, 0);

  mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, PROC_NUM(PE_XIND - 1, PE_YIND), tag,
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
  int xglobal = getGlobalXIndex(xpos);
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
  int xglobal = getGlobalXIndex(xpos);
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

int BoutMesh::sendYOutIndest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1) {
    return 1;
  }

  Timer timer("comms");

  if (UDATA_INDEST != -1) {
    mpi->MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_INDEST, tag, BoutComm::get());
  } else {
    throw BoutException("Expected UDATA_INDEST to exist, but it does not.");
  }
  return 0;
}

int BoutMesh::sendYOutOutdest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1) {
    return 1;
  }

  Timer timer("comms");

  if (UDATA_OUTDEST != -1) {
    mpi->MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, tag, BoutComm::get());
  } else {
    throw BoutException("Expected UDATA_OUTDEST to exist, but it does not.");
  }

  return 0;
}

int BoutMesh::sendYInIndest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == 0) {
    return 1;
  }

  Timer timer("comms");

  if (DDATA_INDEST != -1) {
    mpi->MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_INDEST, tag, BoutComm::get());
  } else {
    throw BoutException("Expected DDATA_INDEST to exist, but it does not.");
  }

  return 0;
}

int BoutMesh::sendYInOutdest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == 0) {
    return 1;
  }

  Timer timer("comms");

  if (DDATA_OUTDEST != -1) {
    mpi->MPI_Send(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, tag, BoutComm::get());
  } else {
    throw BoutException("Expected DDATA_OUTDEST to exist, but it does not.");
  }

  return 0;
}

comm_handle BoutMesh::irecvYOutIndest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1) {
    return nullptr;
  }

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle* ch = get_handle(0, 0);

  if (UDATA_INDEST != -1) {
    mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_INDEST, tag, BoutComm::get(),
                   ch->request);
  } else {
    throw BoutException("Expected UDATA_INDEST to exist, but it does not.");
  }

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvYOutOutdest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == NYPE - 1) {
    return nullptr;
  }

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle* ch = get_handle(0, 0);

  if (UDATA_OUTDEST != -1) {
    mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, UDATA_OUTDEST, tag, BoutComm::get(),
                   ch->request);
  } else {
    throw BoutException("Expected UDATA_OUTDEST to exist, but it does not.");
  }

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvYInIndest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == 0) {
    return nullptr;
  }

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle* ch = get_handle(0, 0);

  if (DDATA_INDEST != -1) {
    mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_INDEST, tag, BoutComm::get(),
                   ch->request);
  } else {
    throw BoutException("Expected DDATA_INDEST to exist, but it does not.");
  }

  ch->in_progress = true;

  return static_cast<comm_handle>(ch);
}

comm_handle BoutMesh::irecvYInOutdest(BoutReal* buffer, int size, int tag) {
  if (PE_YIND == 0) {
    return nullptr;
  }

  Timer timer("comms");

  // Get a communications handle. Not fussy about size of arrays
  CommHandle* ch = get_handle(0, 0);

  if (DDATA_OUTDEST != -1) {
    mpi->MPI_Irecv(buffer, size, PVEC_REAL_MPI_TYPE, DDATA_OUTDEST, tag, BoutComm::get(),
                   ch->request);
  } else {
    throw BoutException("Expected DDATA_OUTDEST to exist, but it does not.");
  }

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

int BoutMesh::PROC_NUM(int xind, int yind) const {
  if ((xind >= NXPE) || (xind < 0)) {
    return -1;
  }
  if ((yind >= NYPE) || (yind < 0)) {
    return -1;
  }

  return yind * NXPE + xind;
}

/// Returns the global X index given a local index
int BoutMesh::XGLOBAL(BoutReal xloc, BoutReal &xglo) const {
  xglo = xloc + PE_XIND * MXSUB;
  return static_cast<int>(xglo);
}

int BoutMesh::getGlobalXIndex(int xlocal) const { return xlocal + PE_XIND * MXSUB; }

int BoutMesh::getGlobalXIndexNoBoundaries(int xlocal) const {
  return xlocal + PE_XIND * MXSUB - MXG;
}

int BoutMesh::getLocalXIndex(int xglobal) const { return xglobal - PE_XIND * MXSUB; }

int BoutMesh::getLocalXIndexNoBoundaries(int xglobal) const {
  return xglobal - PE_XIND * MXSUB + MXG;
}

int BoutMesh::YGLOBAL(BoutReal yloc, BoutReal &yglo) const {
  yglo = yloc + PE_YIND * MYSUB - MYG;
  return static_cast<int>(yglo);
}

int BoutMesh::getGlobalYIndex(int ylocal) const {
  int yglobal =  ylocal + PE_YIND * MYSUB;
  if (jyseps1_2 > jyseps2_1 and PE_YIND*MYSUB + 2*MYG + 1 > ny_inner) {
    // Double null, and we are past the upper target
    yglobal += 2*MYG;
  }
  return yglobal;
}

int BoutMesh::getGlobalYIndexNoBoundaries(int ylocal) const {
  return ylocal + PE_YIND * MYSUB - MYG;
}

int BoutMesh::getLocalYIndex(int yglobal) const {
  int ylocal = yglobal - PE_YIND * MYSUB;
  if (jyseps1_2 > jyseps2_1 and PE_YIND * MYSUB + 2 * MYG + 1 > ny_inner) {
    // Double null, and we are past the upper target
    ylocal -= 2 * MYG;
  }
  return ylocal;
}

int BoutMesh::getLocalYIndexNoBoundaries(int yglobal) const {
  return yglobal - PE_YIND * MYSUB + MYG;
}

int BoutMesh::YGLOBAL(int yloc, int yproc) const { return yloc + yproc * MYSUB - MYG; }

int BoutMesh::YLOCAL(int yglo, int yproc) const { return yglo - yproc * MYSUB + MYG; }

int BoutMesh::getGlobalZIndex(int zlocal) const { return zlocal; }

int BoutMesh::getGlobalZIndexNoBoundaries(int zlocal) const { return zlocal; }

int BoutMesh::getLocalZIndex(int zglobal) const { return zglobal; }

int BoutMesh::getLocalZIndexNoBoundaries(int zglobal) const { return zglobal; }

int BoutMesh::YPROC(int yind) const {
  if ((yind < 0) || (yind >= ny)) {
    return -1;
  }
  return yind / MYSUB;
}

int BoutMesh::XPROC(int xind) const { return (xind >= MXG) ? (xind - MXG) / MXSUB : 0; }

/****************************************************************
 *                     TESTING UTILITIES
 ****************************************************************/

BoutMesh::BoutMesh(int input_nx, int input_ny, int input_nz, int mxg, int myg, int nxpe,
                   int nype, int pe_xind, int pe_yind, bool create_topology,
                   bool symmetric_X, bool symmetric_Y)
    : nx(input_nx), ny(input_ny), nz(input_nz), symmetricGlobalX(symmetric_X),
      symmetricGlobalY(symmetric_Y), MXG(mxg), MYG(myg), MZG(0) {
  maxregionblocksize = MAXREGIONBLOCKSIZE;

  PE_XIND = pe_xind;
  PE_YIND = pe_yind;
  NPES = nxpe * nype;
  MYPE = nxpe * pe_yind + pe_xind;

  ixseps1 = nx;
  ixseps2 = nx;
  jyseps1_1 = -1;
  jyseps1_2 = ny / 2;
  jyseps2_1 = jyseps1_2;
  jyseps2_2 = ny - 1;
  ny_inner = jyseps2_1;
  numberOfXPoints = 0;

  NXPE = nxpe;
  NYPE = nype;
  NZPE = 1;

  // Set the other grid sizes from nx, ny, nz
  setDerivedGridSizes();

  periodicX = false;

  ZMIN = 0.0;
  ZMAX = 1.0;
  zperiod = 1.0;

  if (not create_topology) {
    return;
  }

  // Call topology to set layout of grid
  topology();

  TwistShift = false;
  ShiftAngle.resize(LocalNx);

  ShiftAngle.clear();

  createDefaultRegions();
  addBoundaryRegions();
}

BoutMesh::BoutMesh(int input_nx, int input_ny, int input_nz, int mxg, int myg, int nxpe,
                   int nype, int pe_xind, int pe_yind, bool symmetric_X, bool symmetric_Y,
                   bool periodicX_, int ixseps1_, int ixseps2_, int jyseps1_1_,
                   int jyseps2_1_, int jyseps1_2_, int jyseps2_2_, int ny_inner_)
    : nx(input_nx), ny(input_ny), nz(input_nz), NPES(nxpe * nype),
      MYPE(nxpe * pe_yind + pe_xind), PE_YIND(pe_yind), NYPE(nype), NZPE(1),
      ixseps1(ixseps1_), ixseps2(ixseps2_), symmetricGlobalX(symmetric_X),
      symmetricGlobalY(symmetric_Y), MXG(mxg), MYG(myg), MZG(0) {
  NXPE = nxpe;
  PE_XIND = pe_xind;
  periodicX = periodicX_;
  setYDecompositionIndices(jyseps1_1_, jyseps2_1_, jyseps1_2_, jyseps2_2_, ny_inner_);
  setDerivedGridSizes();
  topology();
  createDefaultRegions();
  addBoundaryRegions();
}
/****************************************************************
 *                       CONNECTIONS
 ****************************************************************/

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
    if (PE_XIND == (NXPE - 1)) {
      ODATA_DEST = PROC_NUM(0, PE_YIND);
    }

    if (PE_XIND == 0) {
      IDATA_DEST = PROC_NUM(NXPE - 1, PE_YIND);
    }
  }
}

void BoutMesh::set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts) {
  if (xlt <= xge) {
    return;
  }

  if ((ypos1 < 0) || (ypos1 >= MY)) {
    output_warn.write("WARNING adding connection: poloidal index {:d} out of range\n",
                      ypos1);
    return;
  }
  if ((ypos2 < 0) || (ypos2 >= MY)) {
    output_warn.write("WARNING adding connection: poloidal index {:d} out of range\n",
                      ypos2);
    return;
  }

  const int ype1 = YPROC(ypos1);
  const int ype2 = YPROC(ypos2);

  /* y index within processors */
  const int yind1 = YLOCAL(ypos1, ype1);
  const int yind2 = YLOCAL(ypos2, ype2);

  /* Check which boundary the connection is on */
  int ypeup = 0;
  int ypedown = 0;
  if ((yind1 == MYG) && (yind2 == MYSUB + MYG - 1)) {
    ypeup = ype2;   /* processor sending data up (+ve y) */
    ypedown = ype1; /* processor sending data down (-ve y) */
  } else if ((yind2 == MYG) && (yind1 == MYSUB + MYG - 1)) {
    ypeup = ype1;
    ypedown = ype2;
  } else {
    throw BoutException(
        "ERROR adding connection: y index {:d} or {:d} not on processor boundary\n",
        ypos1, ypos2);
  }

  /* check the x ranges are possible */
  if ((xge != 0) && (xlt != MX)) {
    throw BoutException(
        "ERROR adding connection({:d},{:d},{:d},{:d}): can only divide X domain in 2\n",
        ypos1, ypos2, xge, xlt);
  }

  output_info.write(
      "Connection between top of Y processor {:d} and bottom of {:d} in range {:d} <= x < {:d}\n",
      ypeup, ypedown, xge, xlt);

  // Convert X coordinates into local indices

  xge = getLocalXIndex(xge);
  xlt = getLocalXIndex(xlt);

  if ((xge >= LocalNx) || (xlt <= 0)) {
    return; // Not in this x domain
  }

  if (xge < 0) {
    xge = 0;
  }
  if (xlt > LocalNx) {
    xlt = LocalNx;
  }

  if (MYPE == PROC_NUM(PE_XIND, ypeup)) { /* PROCESSOR SENDING +VE Y */
    /* Set the branch cut x position */
    if (xge <= MXG) {
      /* Connect on the inside */
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = PROC_NUM(PE_XIND, ypedown);
      if (UDATA_XSPLIT == LocalNx) {
        UDATA_OUTDEST = -1;
      }

      TS_up_in = ts; // Twist-shift

      output_info.write("=> This processor sending in up\n");
    } else {
      /* Connect on the outside */
      if (UDATA_XSPLIT <= 0) {
        UDATA_INDEST = UDATA_OUTDEST;
      }
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = PROC_NUM(PE_XIND, ypedown);
      if (UDATA_XSPLIT <= 0) {
        UDATA_INDEST = -1;
      }

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
      if (DDATA_XSPLIT == LocalNx) {
        DDATA_OUTDEST = -1;
      }

      TS_down_in = ts;

      output_info.write("=> This processor sending in down\n");
    } else {
      /* Connect on the outside */
      if (DDATA_XSPLIT <= 0) {
        DDATA_INDEST = DDATA_OUTDEST;
      }
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = PROC_NUM(PE_XIND, ypeup);
      if (DDATA_XSPLIT == 0) {
        DDATA_INDEST = -1;
      }

      TS_down_out = ts;

      output_info.write("=> This processor sending out down\n");
    }
  }
}

void BoutMesh::add_target(int ypos, int xge, int xlt) {
  if (xlt <= xge) {
    return;
  }

  if ((ypos < 0) || (ypos >= MY)) {
    output_warn.write("WARNING adding target: poloidal index {:d} out of range\n", ypos);
    return;
  }

  int ypeup = YPROC(ypos);
  int ypedown = YPROC(ypos + 1);
  if (ypeup == ypedown) {
    throw BoutException("Adding target at y={:d} in middle of processor {:d}\n", ypos,
                        ypeup);
  }

  output_info.write(
      "Target at top of Y processor {:d} and bottom of {:d} in range {:d} <= x < {:d}\n", ypeup,
      ypedown, xge, xlt);

  // Convert X coordinates into local indices
  xge = getLocalXIndex(xge);
  xlt = getLocalXIndex(xlt);
  if ((xge >= LocalNx) || (xlt <= 0)) {
    return; // Not in this x domain
  }

  if (MYPE == PROC_NUM(PE_XIND, ypeup)) {
    // Target on upper processor boundary
    if (xge <= MXG) {
      // Target on inside
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = -1;
      if (xlt >= LocalNx) {
        UDATA_OUTDEST = -1;
      }
      output_info.write("=> This processor has target upper inner\n");
    } else {
      // Target on outside
      if (UDATA_XSPLIT <= 0) {
        UDATA_INDEST = UDATA_OUTDEST;
      }
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = -1;
      if (xge <= 0) {
        UDATA_INDEST = -1;
      }
      output_info.write("=> This processor has target upper outer\n");
    }
  }
  if (MYPE == PROC_NUM(PE_XIND, ypedown)) {
    // Target on upper processor boundary
    if (xge <= MXG) {
      // Target on inside
      DDATA_XSPLIT = xlt;
      DDATA_INDEST = -1;
      if (xlt >= LocalNx) {
        DDATA_OUTDEST = -1;
      }
      output_info.write("=> This processor has target lower inner\n");
    } else {
      // Target on outside
      if (DDATA_XSPLIT <= 0) {
        DDATA_INDEST = DDATA_OUTDEST;
      }
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = -1;
      if (xge <= 0) {
        DDATA_INDEST = -1;
      }
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
    throw BoutException("\tTopology error: npes={:d} is not equal to NXPE*NYPE={:d}\n",
                        NPES, NXPE * NYPE);
  }
  if (MYSUB * NYPE != MY) {
    throw BoutException("\tTopology error: MYSUB[{:d}] * NYPE[{:d}] != MY[{:d}]\n", MYSUB,
                        NYPE, MY);
  }
  if (MXSUB * NXPE != MX) {
    throw BoutException("\tTopology error: MXSUB[{:d}] * NXPE[{:d}] != MX[{:d}]\n", MXSUB,
                        NXPE, MX);
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

  if ((ixseps_inner > 0)
      && (((PE_YIND * MYSUB > jyseps1_1) && (PE_YIND * MYSUB <= jyseps2_1))
          || ((PE_YIND * MYSUB > jyseps1_2) && (PE_YIND * MYSUB <= jyseps2_2)))) {
    MYPE_IN_CORE = true; /* processor is in the core */
  }

  if (DDATA_XSPLIT > LocalNx) {
    DDATA_XSPLIT = LocalNx;
  }
  if (UDATA_XSPLIT > LocalNx) {
    UDATA_XSPLIT = LocalNx;
  }

  // Print out settings
  output_info.write("\tMYPE_IN_CORE = {}\n", MYPE_IN_CORE);
  output_info.write("\tDXS = {:d}, DIN = {:d}. DOUT = {:d}\n", DDATA_XSPLIT, DDATA_INDEST,
                    DDATA_OUTDEST);
  output_info.write("\tUXS = {:d}, UIN = {:d}. UOUT = {:d}\n", UDATA_XSPLIT, UDATA_INDEST,
                    UDATA_OUTDEST);
  output_info.write("\tXIN = {:d}, XOUT = {:d}\n", IDATA_DEST, ODATA_DEST);

  output_info.write("\tTwist-shift: ");
  if (TS_down_in) {
    output_info.write("DI ");
  }
  if (TS_down_out) {
    output_info.write("DO ");
  }
  if (TS_up_in) {
    output_info.write("UI ");
  }
  if (TS_up_out) {
    output_info.write("UO ");
  }
  output_info.write("\n");
}

/****************************************************************
 *                     Communication handles
 ****************************************************************/

BoutMesh::CommHandle *BoutMesh::get_handle(int xlen, int ylen) {
  if (comm_list.empty()) {
    // Allocate a new CommHandle

    auto* ch = new CommHandle;
    for (auto& i : ch->request) {
      i = MPI_REQUEST_NULL;
    }

    if (ylen > 0) {
      ch->umsg_sendbuff.reallocate(ylen);
      ch->dmsg_sendbuff.reallocate(ylen);
      ch->umsg_recvbuff.reallocate(ylen);
      ch->dmsg_recvbuff.reallocate(ylen);
    }

    if (xlen > 0) {
      ch->imsg_sendbuff.reallocate(xlen);
      ch->omsg_sendbuff.reallocate(xlen);
      ch->imsg_recvbuff.reallocate(xlen);
      ch->omsg_recvbuff.reallocate(xlen);
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
    ch->umsg_sendbuff.reallocate(ylen);
    ch->dmsg_sendbuff.reallocate(ylen);
    ch->umsg_recvbuff.reallocate(ylen);
    ch->dmsg_recvbuff.reallocate(ylen);

    ch->ybufflen = ylen;
  }
  if (ch->xbufflen < xlen) {
    ch->imsg_sendbuff.reallocate(xlen);
    ch->omsg_sendbuff.reallocate(xlen);
    ch->imsg_recvbuff.reallocate(xlen);
    ch->omsg_recvbuff.reallocate(xlen);

    ch->xbufflen = xlen;
  }

  ch->in_progress = false;
  ch->include_x_corners = false;
  ch->has_y_communication = false;

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

/// For debugging purposes (when creating fake parallel meshes), make
/// the send and receive buffers share memory. This allows for
/// communications to be faked between meshes as though they were on
/// different processors.
void BoutMesh::overlapHandleMemory(BoutMesh* yup, BoutMesh* ydown, BoutMesh* xin,
                                   BoutMesh* xout) {
  const int xlen = LocalNy * LocalNz * MXG * 5, ylen = LocalNx * LocalNz * MYG * 5;

  CommHandle* ch = get_handle(xlen, ylen);
  if (yup != nullptr) {
    CommHandle* other = (yup == this) ? ch : yup->get_handle(xlen, ylen);
    if (other->dmsg_sendbuff.unique()) {
      ch->umsg_recvbuff = other->dmsg_sendbuff;
    }
    if (yup != this) {
      yup->free_handle(other);
    }
  }
  if (ydown != nullptr) {
    CommHandle* other = (ydown == this) ? ch : ydown->get_handle(xlen, ylen);
    if (other->umsg_sendbuff.unique()) {
      ch->dmsg_recvbuff = other->umsg_sendbuff;
    }
    if (ydown != this) {
      ydown->free_handle(other);
    }
  }
  if (xin != nullptr) {
    CommHandle* other = (xin == this) ? ch : xin->get_handle(xlen, ylen);
    if (other->omsg_sendbuff.unique()) {
      ch->imsg_recvbuff = other->omsg_sendbuff;
    }
    if (xin != this) {
      xin->free_handle(other);
    }
  }
  if (xout != nullptr) {
    CommHandle* other = (xout == this) ? ch : xout->get_handle(xlen, ylen);
    if (other->imsg_sendbuff.unique()) {
      ch->omsg_recvbuff = other->imsg_sendbuff;
    }
    if (xout != this) {
      xout->free_handle(other);
    }
  }
  free_handle(ch);
}

/****************************************************************
 *                   Communication utilities
 ****************************************************************/

int BoutMesh::pack_data(const std::vector<FieldData *> &var_list, int xge, int xlt, int yge,
                        int ylt, BoutReal *buffer) {

  int len = 0;

  /// Loop over variables
  for (const auto &var : var_list) {
    if (var->is3D()) {
      // 3D variable
      auto& var3d_ref = *dynamic_cast<Field3D*>(var);
      ASSERT2(var3d_ref.isAllocated());
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++) {
          for (int jz = 0; jz < LocalNz; jz++, len++) {
            buffer[len] = var3d_ref(jx, jy, jz);
          }
        }
      }
    } else {
      // 2D variable
      auto& var2d_ref = *dynamic_cast<Field2D*>(var);
      ASSERT2(var2d_ref.isAllocated());
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++, len++) {
          buffer[len] = var2d_ref(jx, jy);
        }
      }
    }
  }

  return (len);
}

int BoutMesh::unpack_data(const std::vector<FieldData *> &var_list, int xge, int xlt, int yge,
                          int ylt, BoutReal *buffer) {

  int len = 0;

  /// Loop over variables
  for (const auto &var : var_list) {
    if (var->is3D()) {
      // 3D variable
      auto& var3d_ref = *dynamic_cast<Field3D*>(var);
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++) {
          for (int jz = 0; jz < LocalNz; jz++, len++) {
            var3d_ref(jx, jy, jz) = buffer[len];
          }
        }
      }
    } else {
      // 2D variable
      auto& var2d_ref = *dynamic_cast<Field2D*>(var);
      for (int jx = xge; jx != xlt; jx++) {
        for (int jy = yge; jy < ylt; jy++, len++) {
          var2d_ref(jx, jy) = buffer[len];
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
  return MYPE_IN_CORE and (getGlobalXIndex(jx) < ixseps_inner);
}

bool BoutMesh::periodicY(int jx, BoutReal& ts) const {
  ts = 0.;
  if (periodicY(jx)) {
    if (TwistShift) {
      ts = ShiftAngle[jx];
    }
    return true;
  }
  return false;
}

int BoutMesh::numberOfYBoundaries() const {
  if (jyseps2_1 != jyseps1_2) {
    return 2;
  } else {
    return 1;
  }
}

std::pair<bool, BoutReal> BoutMesh::hasBranchCutLower(int jx) const {
  if ( (TS_down_in and DDATA_INDEST != -1 and jx < DDATA_XSPLIT)
      or (TS_down_out and DDATA_OUTDEST != -1 and jx >= DDATA_XSPLIT) ) {
    // this processor has branch cut at lower boundary for jx
    if (ShiftAngle.empty()) {
      // This function should only be called during initialization, so always check
      throw BoutException("BoutMesh failed to read ShiftAngle from the grid");
    }
    return std::make_pair(true, ShiftAngle[jx]);
  }

  return std::make_pair(false, 0.);
}


std::pair<bool, BoutReal> BoutMesh::hasBranchCutUpper(int jx) const {
  if ( (TS_up_in and UDATA_INDEST != -1 and jx < UDATA_XSPLIT)
      or (TS_up_out and UDATA_OUTDEST != -1 and jx >= UDATA_XSPLIT) ) {
    // this processor has branch cut at upper boundary for jx
    if (ShiftAngle.empty()) {
      // This function should only be called during initialization, so always check
      throw BoutException("BoutMesh failed to read ShiftAngle from the grid");
    }
    return std::make_pair(true, ShiftAngle[jx]);
  }

  return std::make_pair(false, 0.);
}

int BoutMesh::ySize(int xpos) const {
  int xglobal = getGlobalXIndex(xpos);
  int yglobal = getGlobalYIndexNoBoundaries(MYG);

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
  int xglobal = getGlobalXIndex(xpos);

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
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart)) {
      xs = DDATA_XSPLIT;
    }
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1)) {
      xe = DDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
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
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart)) {
      xs = DDATA_XSPLIT;
    }
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1)) {
      xe = DDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
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
  if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart)) {
    xs = DDATA_XSPLIT;
  }
  if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1)) {
    xe = DDATA_XSPLIT - 1;
  }

  if (xs < xstart) {
    xs = xstart;
  }
  if (xe > xend) {
    xe = xend;
  }

  if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
    // Include corner cells on x-boundary
    xs = 0;
  }
  if (include_corner_cells and lastX() and xe == xend and xs < xend) {
    // Include corner cells on x-boundary
    xe = LocalNx - 1;
  }

  addRegion3D("RGN_LOWER_Y", Region<Ind3D>(xs, xe, 0, ystart-1, 0, LocalNz-1,
                                           LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_LOWER_Y", Region<Ind2D>(xs, xe, 0, ystart-1, 0, 0,
                                           LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_LOWER_Y");
  
  // Upper Inner Y

  xs = 0;
  xe = LocalNx - 1;

  if (!lastY()) {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart)) {
      xs = UDATA_XSPLIT;
    }
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1)) {
      xe = UDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }
  } else {
    xs = -1;
    xe = -2;
  }

  if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
    // Include corner cells on x-boundary
    xs = 0;
  }
  if (include_corner_cells and lastX() and xe == xend and xs < xend) {
    // Include corner cells on x-boundary
    xe = LocalNx - 1;
  }
  
  addRegion3D("RGN_UPPER_INNER_Y", Region<Ind3D>(xs, xe, yend+1, LocalNy-1, 0, LocalNz-1,
                                                 LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_UPPER_INNER_Y", Region<Ind2D>(xs, xe, yend+1, LocalNy-1, 0, 0,
                                                 LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_UPPER_INNER_Y");

  // Upper Outer Y
  
  xs = 0;
  xe = LocalNx - 1;

  if (!lastY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart)) {
      xs = UDATA_XSPLIT;
    }
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1)) {
      xe = UDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
  }

  addRegion3D("RGN_UPPER_OUTER_Y", Region<Ind3D>(xs, xe, yend+1, LocalNy-1, 0, LocalNz-1,
                                                 LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_UPPER_OUTER_Y", Region<Ind2D>(xs, xe, yend+1, LocalNy-1, 0, 0,
                                                 LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_UPPER_OUTER_Y");

  // Upper Y

  xs = 0;
  xe = LocalNx - 1;
  if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart)) {
    xs = UDATA_XSPLIT;
  }
  if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1)) {
    xe = UDATA_XSPLIT - 1;
  }

  if (xs < xstart) {
    xs = xstart;
  }
  if (xe > xend) {
    xe = xend;
  }

  if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
    // Include corner cells on x-boundary
    xs = 0;
  }
  if (include_corner_cells and lastX() and xe == xend and xs < xend) {
    // Include corner cells on x-boundary
    xe = LocalNx - 1;
  }

  addRegion3D("RGN_UPPER_Y", Region<Ind3D>(xs, xe, yend+1, LocalNy-1, 0, LocalNz-1,
                                           LocalNy, LocalNz, maxregionblocksize));
  addRegion2D("RGN_UPPER_Y", Region<Ind2D>(xs, xe, yend+1, LocalNy-1, 0, 0,
                                           LocalNy, 1, maxregionblocksize));
  all_boundaries.emplace_back("RGN_UPPER_Y");
  
  // Inner X
  if(firstX() && !periodicX) {
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, xstart-1, ystart, yend, 0, LocalNz-1,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_INNER_X", Region<Ind2D>(0, xstart-1, ystart, yend, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
    addRegionPerp("RGN_INNER_X", Region<IndPerp>(0, xstart - 1, 0, 0, 0, LocalNz - 1, 1,
                                                 LocalNz, maxregionblocksize));
    all_boundaries.emplace_back("RGN_INNER_X");
    
    output_info.write("\tBoundary region inner X\n");
  } else {
    // Empty region
    addRegion3D("RGN_INNER_X", Region<Ind3D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_INNER_X", Region<Ind2D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
    addRegionPerp("RGN_INNER_X",
                  Region<IndPerp>(0, -1, 0, 0, 0, 0, 1, LocalNz, maxregionblocksize));
  }

  // Outer X
  if(lastX() && !periodicX) {
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(xend+1, LocalNx-1, ystart, yend, 0, LocalNz-1,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_OUTER_X", Region<Ind2D>(xend+1, LocalNx-1, ystart, yend, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
    addRegionPerp("RGN_OUTER_X",
                  Region<IndPerp>(xend + 1, LocalNx - 1, 0, 0, 0, LocalNz - 1, 1, LocalNz,
                                  maxregionblocksize));
    all_boundaries.emplace_back("RGN_OUTER_X");
    
    output_info.write("\tBoundary region outer X\n");
  } else {
    // Empty region
    addRegion3D("RGN_OUTER_X", Region<Ind3D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, LocalNz, maxregionblocksize));
    addRegion2D("RGN_OUTER_X", Region<Ind2D>(0, -1, 0, 0, 0, 0,
                                             LocalNy, 1, maxregionblocksize));
    addRegionPerp("RGN_OUTER_X",
                  Region<IndPerp>(0, -1, 0, 0, 0, 0, 1, LocalNz, maxregionblocksize));
  }

  // Join boundary regions together
  
  Region<Ind3D> bndry3d; // Empty
  for (const auto &region_name : all_boundaries) {
    bndry3d += getRegion3D(region_name);
  }
  bndry3d.unique(); // Ensure that the points are unique

  // Create a region which is all boundaries
  addRegion3D("RGN_BNDRY", bndry3d);

  // Create a region including all x-boundaries
  bndry3d = getRegion3D("RGN_NOBNDRY") + getRegion3D("RGN_INNER_X")
            + getRegion3D("RGN_OUTER_X");
  bndry3d.unique();
  addRegion3D("RGN_WITH_XBNDRIES", bndry3d);

  // Create a region including all y-boundaries
  bndry3d = getRegion3D("RGN_NOBNDRY") + getRegion3D("RGN_LOWER_Y")
            + getRegion3D("RGN_UPPER_Y");
  bndry3d.unique();
  addRegion3D("RGN_WITH_YBNDRIES", bndry3d);

  // Create a region including all boundaries
  bndry3d = getRegion3D("RGN_NOBNDRY") + getRegion3D("RGN_BNDRY");
  bndry3d.unique();
  addRegion3D("RGN_WITH_BNDRIES", bndry3d);

  Region<Ind2D> bndry2d; // Empty
  for (const auto &region_name : all_boundaries) {
    bndry2d += getRegion2D(region_name);
  }
  bndry2d.unique(); // Ensure that the points are unique

  // Create a region which is all boundaries
  addRegion2D("RGN_BNDRY", bndry2d);

  // Create a region including all x-boundaries
  bndry2d = getRegion2D("RGN_NOBNDRY") + getRegion2D("RGN_INNER_X")
            + getRegion2D("RGN_OUTER_X");
  bndry2d.unique();
  addRegion2D("RGN_WITH_XBNDRIES", bndry2d);

  // Create a region including all y-boundaries
  bndry2d = getRegion2D("RGN_NOBNDRY") + getRegion2D("RGN_LOWER_Y")
            + getRegion2D("RGN_UPPER_Y");
  bndry2d.unique();
  addRegion2D("RGN_WITH_YBNDRIES", bndry2d);

  // Create a region including all boundaries
  bndry2d = getRegion2D("RGN_NOBNDRY") + getRegion2D("RGN_BNDRY");
  bndry2d.unique();
  addRegion2D("RGN_WITH_BNDRIES", bndry2d);
}

RangeIterator BoutMesh::iterateBndryLowerInnerY() const {

  int xs = 0;
  int xe = LocalNx - 1;

  if (!firstY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart)) {
      xs = DDATA_XSPLIT;
    }
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1)) {
      xe = DDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
  }
  return RangeIterator(xs, xe);
}

RangeIterator BoutMesh::iterateBndryLowerOuterY() const {

  int xs = 0;
  int xe = LocalNx - 1;
  if (!firstY()) {
    if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart)) {
      xs = DDATA_XSPLIT;
    }
    if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1)) {
      xe = DDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
  } else {
    xs = -1;
    xe = -2;
  }
  return RangeIterator(xs, xe);
}

RangeIterator BoutMesh::iterateBndryLowerY() const {
  int xs = 0;
  int xe = LocalNx - 1;
  if ((DDATA_INDEST >= 0) && (DDATA_XSPLIT > xstart)) {
    xs = DDATA_XSPLIT;
  }
  if ((DDATA_OUTDEST >= 0) && (DDATA_XSPLIT < xend + 1)) {
    xe = DDATA_XSPLIT - 1;
  }

  if (xs < xstart) {
    xs = xstart;
  }
  if (xe > xend) {
    xe = xend;
  }

  if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
    // Include corner cells on x-boundary
    xs = 0;
  }
  if (include_corner_cells and lastX() and xe == xend and xs < xend) {
    // Include corner cells on x-boundary
    xe = LocalNx - 1;
  }

  return RangeIterator(xs, xe);
}

RangeIterator BoutMesh::iterateBndryUpperInnerY() const {
  int xs = 0;
  int xe = LocalNx - 1;

  if (!lastY()) {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart)) {
      xs = UDATA_XSPLIT;
    }
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1)) {
      xe = UDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
  } else {
    xs = -1;
    xe = -2;
  }
  return RangeIterator(xs, xe);
}

RangeIterator BoutMesh::iterateBndryUpperOuterY() const {
  int xs = 0;
  int xe = LocalNx - 1;

  if (!lastY()) {
    xs = -1;
    xe = -2;
  } else {
    if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart)) {
      xs = UDATA_XSPLIT;
    }
    if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1)) {
      xe = UDATA_XSPLIT - 1;
    }

    if (xs < xstart) {
      xs = xstart;
    }
    if (xe > xend) {
      xe = xend;
    }

    if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
      // Include corner cells on x-boundary
      xs = 0;
    }
    if (include_corner_cells and lastX() and xe == xend and xs < xend) {
      // Include corner cells on x-boundary
      xe = LocalNx - 1;
    }
  }
  return RangeIterator(xs, xe);
}

RangeIterator BoutMesh::iterateBndryUpperY() const {
  int xs = 0;
  int xe = LocalNx - 1;
  if ((UDATA_INDEST >= 0) && (UDATA_XSPLIT > xstart)) {
    xs = UDATA_XSPLIT;
  }
  if ((UDATA_OUTDEST >= 0) && (UDATA_XSPLIT < xend + 1)) {
    xe = UDATA_XSPLIT - 1;
  }

  if (xs < xstart) {
    xs = xstart;
  }
  if (xe > xend) {
    xe = xend;
  }

  if (include_corner_cells and firstX() and xs == xstart and xe > xstart) {
    // Include corner cells on x-boundary
    xs = 0;
  }
  if (include_corner_cells and lastX() and xe == xend and xs < xend) {
    // Include corner cells on x-boundary
    xe = LocalNx - 1;
  }

  return RangeIterator(xs, xe);
}

std::vector<BoundaryRegion *> BoutMesh::getBoundaries() { return boundary; }

std::vector<BoundaryRegionPar *> BoutMesh::getBoundariesPar() { return par_boundary; }

void BoutMesh::addBoundaryPar(BoundaryRegionPar *bndry) {
  output_info << "Adding new parallel boundary: " << bndry->label << endl;
  par_boundary.push_back(bndry);
}

Field3D BoutMesh::smoothSeparatrix(const Field3D& f) {
  Field3D result{emptyFrom(f)};
  if ((ixseps_inner > 0) && (ixseps_inner < nx - 1)) {
    if (XPROC(ixseps_inner) == PE_XIND) {
      int x = getLocalXIndex(ixseps_inner);
      for (int y = 0; y < LocalNy; y++) {
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x - 1, y, z));
        }
      }
    }
    if (XPROC(ixseps_inner - 1) == PE_XIND) {
      int x = getLocalXIndex(ixseps_inner - 1);
      for (int y = 0; y < LocalNy; y++) {
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x + 1, y, z));
        }
      }
    }
  }
  if ((ixseps_outer > 0) && (ixseps_outer < nx - 1) && (ixseps_outer != ixseps_inner)) {
    if (XPROC(ixseps_outer) == PE_XIND) {
      int x = getLocalXIndex(ixseps_outer);
      for (int y = 0; y < LocalNy; y++) {
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x - 1, y, z));
        }
      }
    }
    if (XPROC(ixseps_outer - 1) == PE_XIND) {
      int x = getLocalXIndex(ixseps_outer - 1);
      for (int y = 0; y < LocalNy; y++) {
        for (int z = 0; z < LocalNz; z++) {
          result(x, y, z) = 0.5 * (f(x, y, z) + f(x + 1, y, z));
        }
      }
    }
  }
  return result;
}

BoutReal BoutMesh::GlobalX(int jx) const {
  if (symmetricGlobalX) {
    // With this definition the boundary sits dx/2 away form the first/last inner points
    return (0.5 + getGlobalXIndex(jx) - (nx - MX) * 0.5) / static_cast<BoutReal>(MX);
  }
  return static_cast<BoutReal>(getGlobalXIndex(jx)) / static_cast<BoutReal>(MX);
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
    BoutReal yi = getGlobalYIndexNoBoundaries(jy);
    int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

    if (yi < ny_inner) {
      yi -= jyseps1_1 + 0.5;
    } else {
      // Result in core between 0.5 and 1.0
      yi -= jyseps1_1 + 0.5 + (jyseps1_2 - jyseps2_1);
    }
    return yi / nycore;
  }

  int ly = getGlobalYIndexNoBoundaries(jy); // global poloidal index across subdomains
  int nycore = (jyseps2_1 - jyseps1_1) + (jyseps2_2 - jyseps1_2);

  if (MYPE_IN_CORE) {
    // Turn ly into an index over the core cells only
    if (ly <= jyseps2_1) {
      ly -= jyseps1_1 + 1;
    } else {
      ly -= jyseps1_1 + 1 + (jyseps1_2 - jyseps2_1);
    }
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
  file.add(MZSUB, "MZSUB", false);
  file.add(PE_XIND, "PE_XIND", false);
  file.add(PE_YIND, "PE_YIND", false);
  file.add(MYPE, "MYPE", false);
  file.add(MXG, "MXG", false);
  file.add(MYG, "MYG", false);
  file.add(MZG, "MZG", false);
  file.add(nx, "nx", false);
  file.add(ny, "ny", false);
  file.add(nz, "nz", false);
  file.add(MZ, "MZ", false);
  file.add(NXPE, "NXPE", false);
  file.add(NYPE, "NYPE", false);
  file.add(NZPE, "NZPE", false);
  file.add(ZMAX, "ZMAX", false);
  file.add(ZMIN, "ZMIN", false);
  file.add(ixseps1, "ixseps1", false);
  file.add(ixseps2, "ixseps2", false);
  file.add(jyseps1_1, "jyseps1_1", false);
  file.add(jyseps1_2, "jyseps1_2", false);
  file.add(jyseps2_1, "jyseps2_1", false);
  file.add(jyseps2_2, "jyseps2_2", false);
  file.add(ny_inner, "ny_inner", false);

  getCoordinates()->outputVars(file);

  // Try and save some provenance tracking info that new enough versions of
  // hypnotoad provide in the grid file.
  // Note with current Datafile/DataFormat implementation, must not write an
  // empty string because it ends up as a null char* pointer, which causes a
  // segfault.
  if (this->get(grid_id, "grid_id") == 0 and not grid_id.empty()) {
    file.add(grid_id, "grid_id", false);
  }
  if (this->get(hypnotoad_version, "hypnotoad_version") == 0
      and not hypnotoad_version.empty()) {

    file.add(hypnotoad_version, "hypnotoad_version", false);
  }
  if (this->get(hypnotoad_git_hash, "hypnotoad_git_hash") == 0
      and not hypnotoad_git_hash.empty()) {

    file.add(hypnotoad_git_hash, "hypnotoad_git_hash", false);
  }
  if (this->get(hypnotoad_git_diff, "hypnotoad_git_diff") == 0
      and not hypnotoad_git_diff.empty()) {

    file.add(hypnotoad_git_diff, "hypnotoad_git_diff", false);
  }
  if (this->get(hypnotoad_geqdsk_filename, "hypnotoad_geqdsk_filename") == 0
      and not hypnotoad_geqdsk_filename.empty()) {

    file.add(hypnotoad_geqdsk_filename, "hypnotoad_geqdsk_filename", false);
  }
}
