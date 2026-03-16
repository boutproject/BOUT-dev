#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_operators.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"
#include <memory>
#include <string>
#include <utility>

Region<Ind3D> PetscMapping::create_region(const Field3D& cell_number) {
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, cell_number.getRegion("RGN_NOBNDRY")) {
    if (cell_number[i] > -1) {
      indices.push_back(i);
    }
  }
  return Region<Ind3D>(indices);
}

Region<Ind3D> PetscMapping::create_region_xin(const Field3D& cell_number) {
  Region<Ind3D>::RegionIndices xin_indices;
  const Mesh* mesh = cell_number.getMesh();
  if (mesh->firstX()) {
    for (int i = 0; i < mesh->xstart; ++i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = mesh->zstart; k <= mesh->zend; ++k) {
          if (cell_number(i, j, k) > -1) {
            xin_indices.push_back(cell_number.indexAt(i, j, k));
          }
        }
      }
    }
  }
  return Region<Ind3D>(xin_indices);
}

Region<Ind3D> PetscMapping::create_region_xout(const Field3D& cell_number) {
  Region<Ind3D>::RegionIndices xout_indices;
  const Mesh* mesh = cell_number.getMesh();
  if (mesh->lastX()) {
    for (int i = mesh->xend + 1; i < mesh->LocalNx; ++i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = mesh->zstart; k <= mesh->zend; ++k) {
          if (cell_number(i, j, k) > -1) {
            xout_indices.push_back(cell_number.indexAt(i, j, k));
          }
        }
      }
    }
  }
  return Region<Ind3D>(xout_indices);
}

PetscMapping::PetscMapping(const Field3D& cell_number, const Field3D& forward_cell_number,
                           const Field3D& backward_cell_number)
    : cell_number(cell_number), forward_cell_number(forward_cell_number),
      backward_cell_number(backward_cell_number),
      evolving_region(create_region(cell_number)),
      xin_region(create_region_xin(cell_number)),
      xout_region(create_region_xout(cell_number)),
      yup_region(create_region(forward_cell_number)),
      ydown_region(create_region(backward_cell_number)) {
  // Calculate size of each region
  const unsigned int nlocal = this->size();
  output.write("Local {} : evolving {} xin {} xout {} yup {} ydown {}", nlocal,
               evolving_region.size(), xin_region.size(), xout_region.size(),
               yup_region.size(), ydown_region.size());

  // Create a PETSc matrix
  // Note: Numbering is different from that used in PetscVector / PetscMatrix
  // because yup/down boundaries are included.

  // Renumbering matrix.
  // Maps PETSc row (or column) indices to the global indices used in
  // the mesh.  This is needed because the PETSc indices depend on the
  // number of processors.
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), &mat_mesh_to_petsc));
  BOUT_DO_PETSC(
      MatSetSizes(mat_mesh_to_petsc, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
  BOUT_DO_PETSC(MatSetType(mat_mesh_to_petsc, MATMPIAIJ));
  // Each row will have one non-zero entry, which could be in
  // either the "diagonal" or "off-diagonal" block.
  BOUT_DO_PETSC(MatMPIAIJSetPreallocation(mat_mesh_to_petsc, 1, nullptr, 1, nullptr));

  // Get the range of rows owned by this processor
  BOUT_DO_PETSC(MatGetOwnershipRange(mat_mesh_to_petsc, &row_start, &row_end));
  output.write("Local row range: {} -> {}\n", row_start, row_end);

  // Iterate through regions in this order
  this->map([&](PetscInt row, PetscInt col) {
    // `row` is the PETSc index; `col` is the Mesh index
    ASSERT1(row >= 0);
    ASSERT1(col >= 0);
    const PetscScalar ONE = 1.0;
    BOUT_DO_PETSC(MatSetValues(mat_mesh_to_petsc, 1, &row, 1, &col, &ONE, INSERT_VALUES));
  });
  BOUT_DO_PETSC(MatAssemblyBegin(mat_mesh_to_petsc, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(mat_mesh_to_petsc, MAT_FINAL_ASSEMBLY));

  // Transpose to map petsc indices to mesh indices
  MatTranspose(mat_mesh_to_petsc, MAT_INITIAL_MATRIX, &mat_petsc_to_mesh);
}

PetscMapping::~PetscMapping() {
  MatDestroy(&this->mat_mesh_to_petsc);
  MatDestroy(&this->mat_petsc_to_mesh);
}

PetscOperator::PetscOperator(PetscMappingPtr mapping, Array<int> rows, Array<int> cols,
                             Array<BoutReal> weights)
    : mapping(std::move(mapping)), rhs_vec(createVec(this->mapping->size())),
      result_vec(createVec(this->mapping->size())) {

  output.write("PetscOperator Array sizes {} {} {}\n", rows.size(), cols.size(),
               weights.size());

  Mat mat{nullptr};
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), &mat));
  BOUT_DO_PETSC(MatSetType(mat, MATMPIAIJ));
  const int nlocal = this->mapping->size();
  BOUT_DO_PETSC(MatSetSizes(mat, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));

  this->mapping->map([&](PetscInt row, PetscInt mesh_index) {
    if (mesh_index >= rows.size()) {
      return; // No weights -> skip
    }
    // Get the range of indices into columns and weights
    const int start_ind = rows[mesh_index];
    if (start_ind < 0) {
      return; // No entries
    }
    int end_ind = cols.size(); // End of the columns array
    for (int i = mesh_index + 1; i < rows.size(); ++i) {
      // rows[i] can be -1 if no weights
      if (rows[i] > -1) {
        // This is the next entry in the columns / weights array
        end_ind = rows[i];
        break;
      }
    }
    BOUT_DO_PETSC(MatSetValues(mat, 1, &row, end_ind - start_ind, &cols[start_ind],
                               &weights[start_ind], INSERT_VALUES));
  });
  BOUT_DO_PETSC(MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY));

  // Row indices are PETSc indices but columns are mesh indices
  // Multiply on the right by PetscToMesh.
  BOUT_DO_PETSC(MatMatMult(mat, this->mapping->getPetscToMesh(), MAT_INITIAL_MATRIX,
                           PETSC_DETERMINE, &this->mat_operator));
  // Destroy temporary matrix
  BOUT_DO_PETSC(MatDestroy(&mat));
}

void PetscOperator::copyToVec(PetscMappingPtr mapping, const Field3D& f, Vec vec) {
  ASSERT1(f.hasParallelSlices());
  ASSERT1(f.yup().isAllocated());
  ASSERT1(f.ydown().isAllocated());

  PetscScalar* x{nullptr};
  // Copy rows from Field3D cells
  BOUT_DO_PETSC(VecGetArray(vec, &x));
  mapping->map_local_field([&](PetscInt row, const Ind3D& i) { x[row] = f[i]; });

  // Copy Yup / Ydown values from boundaries
  const Field3D& yup = f.yup();
  mapping->map_local_yup([&](PetscInt row, const Ind3D& i) { x[row] = yup[i]; });

  const Field3D& ydown = f.ydown();
  mapping->map_local_ydown([&](PetscInt row, const Ind3D& i) { x[row] = ydown[i]; });

  BOUT_DO_PETSC(VecRestoreArray(vec, &x));
}

Vec PetscOperator::createVec(PetscInt local_size) {
  Vec vec{nullptr};
  BOUT_DO_PETSC(VecCreate(BoutComm::get(), &vec));
  BOUT_DO_PETSC(VecSetType(vec, VECMPI));
  BOUT_DO_PETSC(VecSetSizes(vec, local_size, PETSC_DETERMINE));
  return vec;
}

PetscOperator::~PetscOperator() {
  MatDestroy(&this->mat_operator);
  VecDestroy(&this->rhs_vec);
  VecDestroy(&this->result_vec);
}

PetscOperator PetscOperator::diagonal(PetscMappingPtr mapping, const Field3D& f) {
  Vec diag{createVec(mapping->size())};
  copyToVec(mapping, f, diag);

  Mat mat{nullptr};
  // Note: MatMatMult with diagonal and mpiaij not supported
  // -> Create MPIAIJ
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), &mat));
  BOUT_DO_PETSC(MatSetType(mat, MATMPIAIJ));
  const int nlocal = mapping->size();
  BOUT_DO_PETSC(MatSetSizes(mat, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
  BOUT_DO_PETSC(MatMPIAIJSetPreallocation(mat, 1, nullptr, 0, nullptr));
  BOUT_DO_PETSC(MatDiagonalSet(mat, diag, INSERT_VALUES));
  BOUT_DO_PETSC(MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY));

  BOUT_DO_PETSC(VecDestroy(&diag));

  return PetscOperator(std::move(mapping), mat);
}

/// Perform operation
Field3D PetscOperator::operator()(const Field3D& rhs) const {
  // Fill vec from rhs
  copyToVec(this->mapping, rhs, this->rhs_vec);

  // Perform Mat-Vec muliplication
  BOUT_DO_PETSC(MatMult(this->mat_operator, rhs_vec, result_vec));

  // Copy result_vec into a Field3D
  Field3D result{emptyFrom(rhs)}; // This allocates memory
  const PetscScalar* r{nullptr};
  BOUT_DO_PETSC(VecGetArrayRead(result_vec, &r));
  this->mapping->map_local_field(
      [&](PetscInt row, const Ind3D& i) { result[i] = r[row]; });
  BOUT_DO_PETSC(VecRestoreArrayRead(result_vec, &r));

  return result;
}

/// Operator composition
PetscOperator PetscOperator::operator*(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  Mat mat;
  BOUT_DO_PETSC(MatMatMult(this->mat_operator, rhs.mat_operator, MAT_INITIAL_MATRIX,
                           PETSC_DETERMINE, &mat));
  return PetscOperator(this->mapping, mat);
}

/// Operator addition
PetscOperator PetscOperator::operator+(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  Mat mat;
  BOUT_DO_PETSC(MatDuplicate(mat_operator, MAT_COPY_VALUES, &mat));
  BOUT_DO_PETSC(MatAXPY(mat, 1.0, rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
  return PetscOperator(this->mapping, mat);
}

/// Operator subtraction
PetscOperator PetscOperator::operator-(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  Mat mat;
  BOUT_DO_PETSC(MatDuplicate(mat_operator, MAT_COPY_VALUES, &mat));
  BOUT_DO_PETSC(MatAXPY(mat, -1.0, rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
  return PetscOperator(this->mapping, mat);
}

PetscOperator PetscOperator::transpose() const {
  Mat mat{nullptr};
  BOUT_DO_PETSC(MatTranspose(this->mat_operator, MAT_INITIAL_MATRIX, &mat));
  return PetscOperator(this->mapping, mat);
}

Field3D PetscOperators::meshGetField3D(Mesh* mesh, const std::string& name) const {
  Field3D result;
  if (mesh->get(result, name) != 0) {
    throw BoutException("PetscOperators requires field '{}'", name);
  }
  return result;
}

PetscOperators::PetscOperators(Mesh* mesh)
    : mesh(mesh), mapping(std::make_shared<const PetscMapping>(
                      meshGetField3D(mesh, "cell_number"),
                      meshGetField3D(mesh, "forward_cell_number"),
                      meshGetField3D(mesh, "backward_cell_number"))) {

  int mesh_total_cells;
  if (mesh->get(mesh_total_cells, "total_cells") == 0) {
    // Check total number of cells
    if (mesh_total_cells != mapping->size()) {
      throw BoutException("Total cells in mesh {} doesn't match mapping size {}",
                          mesh_total_cells, mapping->size());
    }
  }
}

PetscOperator PetscOperators::get(const std::string& name) const {
  return PetscOperator(this->mapping, this->meshGetArray<int>(this->mesh, name + "_rows"),
                       this->meshGetArray<int>(this->mesh, name + "_columns"),
                       this->meshGetArray<BoutReal>(this->mesh, name + "_weights"));
}

#endif // BOUT_HAS_PETSC
