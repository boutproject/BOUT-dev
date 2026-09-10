/// Integrated test for ordering equivalence between PetscCellMapping
/// and the SNES solver's globalIndex traversal.
///
/// For every locally owned evolving interior cell this program checks:
///
///   1. The local offset into mapOwnedInteriorCells equals the local SNES
///      index from globalIndex (ordering equivalence).
///
///   2. The shift (rowStart() - Istart) is the same constant for every cell
///      on this rank (scalar-shift property).
///
/// Results are written to the BOUT++ output so the runtest script can parse
/// them.  The program exits with a non-zero status if any check fails.

#include <bout/build_defines.hxx>

#if !BOUT_HAS_PETSC
// Without PETSc there is nothing to test.  Exit cleanly so CTest skips rather
// than fails.
#include <cstdlib>
int main() { return EXIT_SUCCESS; }
#else

#include <bout/bout.hxx>
#include <bout/boutcomm.hxx>
#include <bout/field3d.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/petsc_operators.hxx>
#include <bout/region.hxx>
#include <bout/utils.hxx>

#include <mpi.h>

#include <petscmat.h>
#include <petscsys.h>

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

int main(int argc, char** argv) {
  int total_mismatches = 0;
  int total_shift_fail = 0;

  // Initialise BOUT++ (reads BOUT.inp, sets up mesh, MPI, PETSc, etc.)
  BoutInitialise(argc, argv);
  {
    Mesh* mesh = bout::globals::mesh;

    // ── Build cell_number field ────────────────────────────────────────────────
    // Assign stored indices to interior cells exactly as the Python weights
    // module does: consecutive integers over the evolving region.
    Field3D cell_number{-1.0, mesh};
    {
      const int ngy = mesh->GlobalNyNoBoundaries;
      const int ngz = mesh->GlobalNzNoBoundaries;

      for (int x = mesh->xstart; x <= mesh->xend; ++x) {
        const int global_x =
            mesh->getGlobalXIndexNoBoundaries(x); // Exclude boundary cells
        for (int y = mesh->ystart; y <= mesh->yend; ++y) {
          const int global_y = mesh->getGlobalYIndexNoBoundaries(y);
          for (int z = 0; z < mesh->LocalNz; ++z) {
            const int global_z = mesh->getGlobalZIndexNoBoundaries(z);
            cell_number(x, y, z) = ((global_x * ngy) + global_y) * ngz + global_z;
          }
        }
      }
      // Communicate so guard cells are filled (needed by some mesh operations)
      mesh->communicate(cell_number);
    }

    const Field3D forward_cell_number{-1.0, mesh};
    const Field3D backward_cell_number{-1.0, mesh};

    // Count total evolving cells across all ranks
    int n_local_evolving = 0;
    for (const auto& i : mesh->getRegion3D("RGN_NOBNDRY")) {
      if (cell_number[i] >= 0) {
        ++n_local_evolving;
      }
    }
    int n_global_evolving = 0;
    MPI_Allreduce(&n_local_evolving, &n_global_evolving, 1, MPI_INT, MPI_SUM,
                  BoutComm::get());

    // ── Build PetscCellMapping ─────────────────────────────────────────────────
    auto mapping = std::make_shared<const PetscCellMapping>(
        cell_number, forward_cell_number, backward_cell_number, n_global_evolving);

    // ── Build globalIndex field (replicates FDJinitialise logic) ──────────────
    // globalIndex(0) gives 0-based local indices over RGN_NOBNDRY.
    // After shifting by Istart these become the SNES global indices.
    Field3D snes_index{-1.0, mesh};
    {
      int local_idx = 0;
      for (const auto& i : mesh->getRegion3D("RGN_NOBNDRY")) {
        if (cell_number[i] >= 0) {
          snes_index[i] = local_idx++;
        }
      }
    }

    // Determine SNES global offset for this rank by constructing a dummy PETSc
    // matrix of size n_global_evolving and reading its ownership range.
    Mat dummy{nullptr};
    MatCreate(BoutComm::get(), &dummy);
    MatSetSizes(dummy, n_local_evolving, n_local_evolving, PETSC_DETERMINE,
                PETSC_DETERMINE);
    MatSetType(dummy, MATMPIAIJ);
    MatSetUp(dummy);
    PetscInt Istart, Iend;
    MatGetOwnershipRange(dummy, &Istart, &Iend);
    MatDestroy(&dummy);

    // ── Compare traversal orderings ───────────────────────────────────────────
    // Walk mapOwnedInteriorCells and RGN_NOBNDRY in parallel, checking that
    // the local offset produced by each is identical for every cell.

    struct Mismatch {
      int x, y, z;
      PetscInt mapping_local; // local offset from mapOwnedInteriorCells
      int snes_local;         // local offset from globalIndex
    };
    std::vector<Mismatch> mismatches;

    // Collect (Ind3D -> mapping local offset) from mapOwnedInteriorCells
    // so we can look up by field index.
    // mapping_local[Ind3D linear] = local offset (0-based on this rank)
    const PetscInt rstart = mapping->rowStart();
    std::vector<int> mapping_local_by_ind(
        static_cast<std::size_t>(mesh->LocalNx * mesh->LocalNy * mesh->LocalNz), -1);
    {
      PetscInt local_offset = 0;
      mapping->mapOwnedInteriorCells(
          [&](PetscInt petsc_row, const Ind3D& i, int /*stored*/) {
            mapping_local_by_ind[i.ind] = static_cast<int>(petsc_row - rstart);
            ++local_offset;
          });
    }

    // Walk RGN_NOBNDRY and compare
    bool shift_initialised = false;
    PetscInt expected_shift = 0;
    bool shift_constant = true;

    int snes_local_idx = 0;
    for (const auto& i : mesh->getRegion3D("RGN_NOBNDRY")) {
      if (cell_number[i] < 0) {
        continue;
      }

      const int map_local = mapping_local_by_ind[i.ind];
      if (map_local < 0) {
        // Cell is in RGN_NOBNDRY but not visited by mapOwnedInteriorCells
        mismatches.push_back({i.x(), i.y(), i.z(), -1, snes_local_idx});
        ++snes_local_idx;
        continue;
      }

      // Check ordering equivalence: local offsets must match
      if (map_local != snes_local_idx) {
        mismatches.push_back(
            {i.x(), i.y(), i.z(), static_cast<PetscInt>(map_local), snes_local_idx});
      }

      // Check scalar-shift property: petsc_global - snes_global must be constant
      const PetscInt petsc_global = rstart + map_local;
      const PetscInt snes_global = Istart + snes_local_idx;
      const PetscInt shift = petsc_global - snes_global;
      if (!shift_initialised) {
        expected_shift = shift;
        shift_initialised = true;
      } else if (shift != expected_shift) {
        shift_constant = false;
      }

      ++snes_local_idx;
    }

    // ── Gather results across ranks ───────────────────────────────────────────
    const int rank_mismatches = static_cast<int>(mismatches.size());
    MPI_Allreduce(&rank_mismatches, &total_mismatches, 1, MPI_INT, MPI_SUM,
                  BoutComm::get());

    int rank_shift_fail = shift_constant ? 0 : 1;
    MPI_Allreduce(&rank_shift_fail, &total_shift_fail, 1, MPI_INT, MPI_SUM,
                  BoutComm::get());

    // ── Report ────────────────────────────────────────────────────────────────
    const int my_rank = BoutComm::rank();

    // Each rank writes its own diagnostics; rank 0 also writes the summary.
    for (const auto& m : mismatches) {
      output.write("MISMATCH rank={} cell=({},{},{}) mapping_local={} snes_local={}\n",
                   my_rank, m.x, m.y, m.z, m.mapping_local, m.snes_local);
    }
    if (!shift_constant) {
      output.write("SHIFT_NOT_CONSTANT rank={} expected_shift={}\n", my_rank,
                   expected_shift);
    }

    if (my_rank == 0) {
      output.write("total_mismatches={}\n", total_mismatches);
      output.write("total_shift_failures={}\n", total_shift_fail);
      if (total_mismatches == 0 && total_shift_fail == 0) {
        output.write("ordering_check=PASS\n");
      } else {
        output.write("ordering_check=FAIL\n");
      }
    }
  }
  BoutFinalise();

  return (total_mismatches == 0 && total_shift_fail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

#endif // BOUT_HAS_PETSC
