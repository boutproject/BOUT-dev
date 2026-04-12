#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/assert.hxx"
#include "bout/petsc_jacobian.hxx"
#include "bout/petsclib.hxx"

#include <petscmat.h>

void addOperatorSparsity(Mat Jfd, Mat sub, int out_var, int in_var) {
  // Infer nvars from global sizes
  PetscInt jfd_global{0}, sub_global{0};
  BOUT_DO_PETSC(MatGetSize(Jfd, &jfd_global, nullptr));
  BOUT_DO_PETSC(MatGetSize(sub, &sub_global, nullptr));

  ASSERT1(sub_global > 0);
  ASSERT1(jfd_global % sub_global == 0);

  const PetscInt nvars = jfd_global / sub_global;

  ASSERT1(out_var >= 0 && out_var < static_cast<int>(nvars));
  ASSERT1(in_var >= 0 && in_var < static_cast<int>(nvars));

  // Iterate over locally owned rows of sub and insert into Jfd
  PetscInt rstart{0}, rend{0};
  BOUT_DO_PETSC(MatGetOwnershipRange(sub, &rstart, &rend));

  const PetscScalar one = 1.0;

  for (PetscInt sub_row = rstart; sub_row < rend; ++sub_row) {
    PetscInt ncols{0};
    const PetscInt* sub_cols{nullptr};
    const PetscScalar* vals{nullptr};
    BOUT_DO_PETSC(MatGetRow(sub, sub_row, &ncols, &sub_cols, &vals));

    const PetscInt jfd_row = (sub_row * nvars) + out_var;

    for (PetscInt k = 0; k < ncols; ++k) {
      const PetscInt jfd_col = (sub_cols[k] * nvars) + in_var;
      BOUT_DO_PETSC(MatSetValues(Jfd, 1, &jfd_row, 1, &jfd_col, &one, INSERT_VALUES));
    }

    BOUT_DO_PETSC(MatRestoreRow(sub, sub_row, &ncols, &sub_cols, &vals));
  }
}

#endif // BOUT_HAS_PETSC
