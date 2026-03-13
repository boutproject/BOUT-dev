#include "bout/petsc_operators.hxx"
#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/region.hxx"
#include <memory>
#include <string>

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
          if (cell_number(i, j, k) > 0) {
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
          if (cell_number(i, j, k) > 0) {
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
    : cell_number(cell_number), evolving_region(create_region(cell_number)),
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
  MatCreate(BoutComm::get(), &mat_mesh_to_petsc);
  MatSetSizes(mat_mesh_to_petsc, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE);
  MatSetType(mat_mesh_to_petsc, MATMPIAIJ);
  // Each row will have one non-zero entry, which could be in
  // either the "diagonal" or "off-diagonal" block.
  MatMPIAIJSetPreallocation(mat_mesh_to_petsc, 1, nullptr, 1, nullptr);

  // Get the range of rows owned by this processor
  MatGetOwnershipRange(mat_mesh_to_petsc, &row_start, &row_end);
  output.write("Local row range: {} -> {}\n", row_start, row_end);

  // Iterate through regions in this order
  this->map([&](PetscInt row, PetscInt col) {
    // `row` is the PETSc index; `col` is the Mesh index
    const PetscScalar ONE = 1.0;
    MatSetValues(mat_mesh_to_petsc, 1, &row, 1, &col, &ONE, INSERT_VALUES);
  });
  MatAssemblyBegin(mat_mesh_to_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat_mesh_to_petsc, MAT_FINAL_ASSEMBLY);
}

PetscOperator::PetscOperator(PetscMappingPtr mapping, Array<int> rows, Array<int> cols,
                             Array<BoutReal> weights)
    : mapping(mapping) {

  MatCreate(BoutComm::get(), &mat_operator);
  int nlocal = mapping->size();
  MatSetSizes(mat_operator, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE);

  mapping->map_evolving([&](PetscInt row, Ind3D, PetscInt mesh_index) {
    if (mesh_index >= rows.size()) {
      return; // No weights -> skip
    }
    // Get the range of indices into columns and weights
    int start_ind = rows[mesh_index];
    int end_ind = (mesh_index + 1 < rows.size()) ? rows[mesh_index + 1] : cols.size();
    output.write("Mapping {} : {} -> {}\n", row, start_ind, end_ind);

    MatSetValues(mat_operator, 1, &row, end_ind - start_ind, &cols[start_ind],
                 &weights[start_ind], INSERT_VALUES);
  });
  MatAssemblyBegin(mat_operator, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat_operator, MAT_FINAL_ASSEMBLY);
}

/// Perform operation
Field3D PetscOperator::operator()(const Field3D& rhs) const {
  Vec rhs_vec;
  VecCreate(BoutComm::get(), &rhs_vec);
  VecSetSizes(rhs_vec, this->mapping->size(), PETSC_DETERMINE);

  // Fill vec from rhs
  PetscScalar* x;
  VecGetArray(rhs_vec, &x);
  throw BoutException("Not implemented");
  VecRestoreArray(rhs_vec, &x);

  // Perform Mat-Vec muliplication
  Vec result_vec;
  VecDuplicate(rhs_vec, &result_vec);
  MatMult(this->mat_operator, rhs_vec, result_vec);

  // Copy result_vec into a Field3D
  Field3D result;
  const PetscScalar* r;
  VecGetArrayRead(result_vec, &r);
  throw BoutException("Not implemented");
  VecRestoreArrayRead(result_vec, &r);

  return result;
}

/// Operator composition
PetscOperator PetscOperator::operator*(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  Mat mat;
  MatMatMult(this->mat_operator, rhs.mat_operator, MAT_INITIAL_MATRIX, PETSC_DETERMINE,
             &mat);
  return PetscOperator(this->mapping, mat);
}

/// Operator addition
PetscOperator PetscOperator::operator+(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  Mat mat;
  MatDuplicate(mat_operator, MAT_COPY_VALUES, &mat);
  MatAXPY(mat, 1.0, rhs.mat_operator, UNKNOWN_NONZERO_PATTERN);
  return PetscOperator(this->mapping, mat);
}

/// Operator subtraction
PetscOperator PetscOperator::operator-(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  Mat mat;
  MatDuplicate(mat_operator, MAT_COPY_VALUES, &mat);
  MatAXPY(mat, -1.0, rhs.mat_operator, UNKNOWN_NONZERO_PATTERN);
  return PetscOperator(this->mapping, mat);
}

Field3D PetscOperators::meshGetField3D(Mesh* mesh, const std::string& name) {
  Field3D result;
  if (mesh->get(result, name) != 0) {
    throw BoutException("PetscOperators requires field '{}'", name);
  }
  return result;
}

PetscOperators::PetscOperators(Mesh* mesh)
    : mapping(std::make_shared<const PetscMapping>(
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

  // Read forward operator
  PetscOperator forward(mapping, this->meshGetArray<int>(mesh, "forward_rows"),
                        this->meshGetArray<int>(mesh, "forward_columns"),
                        this->meshGetArray<BoutReal>(mesh, "forward_weights"));
}
