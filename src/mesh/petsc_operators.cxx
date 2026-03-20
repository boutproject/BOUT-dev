#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/array.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_interface.hxx"
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
  const unsigned int nlocal = this->localSize();
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
    : mapping(std::move(mapping)), mat_operator(new Mat{nullptr}),
      rhs_vec(createVec(this->mapping->localSize())),
      result_vec(createVec(this->mapping->localSize())) {

  output.write("PetscOperator Array sizes {} {} {}\n", rows.size(), cols.size(),
               weights.size());

  Mat mat{nullptr};
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), &mat));
  BOUT_DO_PETSC(MatSetType(mat, MATMPIAIJ));
  const int nlocal = this->mapping->localSize();
  BOUT_DO_PETSC(MatSetSizes(mat, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));

  // This reads CSR format but is defensive in handling negative indices and missing
  // final 'row' entry.
  this->mapping->map([&](PetscInt row, PetscInt mesh_index) {
    if (mesh_index >= rows.size()) {
      return; // No weights -> skip
    }
    // Get the range of indices into columns and weights
    const int start_ind = rows[mesh_index];
    if (start_ind < 0) {
      return; // No entries (non-standard CSR)
    }
    int end_ind =
        cols.size(); // End of the columns array (should be last element of rows)
    for (int i = mesh_index + 1; i < rows.size(); ++i) {
      // rows[i] can be -1 if no weights (non-standard CSR)
      if (rows[i] > -1) {
        // This is the next entry in the columns / weights array
        end_ind = rows[i];
        break;
      }
    }
    if (end_ind == start_ind) {
      // Empty row in CSR format
      return;
    }
    BOUT_DO_PETSC(MatSetValues(mat, 1, &row, end_ind - start_ind, &cols[start_ind],
                               &weights[start_ind], INSERT_VALUES));
  });
  BOUT_DO_PETSC(MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY));

  // Row indices are PETSc indices but columns are mesh indices
  // Multiply on the right by PetscToMesh.
  BOUT_DO_PETSC(MatMatMult(mat, this->mapping->getPetscToMesh(), MAT_INITIAL_MATRIX,
                           PETSC_DETERMINE, this->mat_operator.get()));
  // Destroy temporary matrix
  BOUT_DO_PETSC(MatDestroy(&mat));
}

void PetscOperator::copyToVec(const PetscMappingPtr& mapping, const Field3D& f, Vec vec) {
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

PetscOperator::UniqueVec PetscOperator::createVec(PetscInt local_size) {
  UniqueVec vec(new Vec{nullptr});
  BOUT_DO_PETSC(VecCreate(BoutComm::get(), vec.get()));
  BOUT_DO_PETSC(VecSetType(*vec, VECMPI));
  BOUT_DO_PETSC(VecSetSizes(*vec, local_size, PETSC_DETERMINE));
  return vec;
}

PetscOperator PetscOperator::diagonal(PetscMappingPtr mapping, const Field3D& f) {
  const UniqueVec diag{createVec(mapping->localSize())};
  copyToVec(mapping, f, *diag);

  UniqueMat mat(new Mat{nullptr});
  // Note: MatMatMult with diagonal and mpiaij not supported
  // -> Create MPIAIJ
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), mat.get()));
  BOUT_DO_PETSC(MatSetType(*mat, MATMPIAIJ));
  const PetscInt nlocal = mapping->localSize();
  BOUT_DO_PETSC(MatSetSizes(*mat, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE));
  BOUT_DO_PETSC(MatMPIAIJSetPreallocation(*mat, 1, nullptr, 0, nullptr));
  BOUT_DO_PETSC(MatDiagonalSet(*mat, *diag, INSERT_VALUES));
  BOUT_DO_PETSC(MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY));

  return PetscOperator(std::move(mapping), std::move(mat));
}

/// Perform operation
Field3D PetscOperator::operator()(const Field3D& rhs) const {
  // Fill vec from rhs
  copyToVec(this->mapping, rhs, *this->rhs_vec);

  // Perform Mat-Vec muliplication
  BOUT_DO_PETSC(MatMult(*this->mat_operator, *rhs_vec, *result_vec));

  // Copy result_vec into a Field3D
  Field3D result{emptyFrom(rhs)}; // This allocates memory
  const PetscScalar* r{nullptr};
  BOUT_DO_PETSC(VecGetArrayRead(*result_vec, &r));
  this->mapping->map_local_field(
      [&](PetscInt row, const Ind3D& i) { result[i] = r[row]; });
  BOUT_DO_PETSC(VecRestoreArrayRead(*result_vec, &r));

  return result;
}

/// Operator composition
PetscOperator PetscOperator::operator*(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  UniqueMat mat(new Mat{nullptr});
  BOUT_DO_PETSC(MatMatMult(*this->mat_operator, *rhs.mat_operator, MAT_INITIAL_MATRIX,
                           PETSC_DETERMINE, mat.get()));
  return PetscOperator(this->mapping, std::move(mat));
}

/// Operator addition
PetscOperator PetscOperator::operator+(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  UniqueMat mat(new Mat{nullptr});
  BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, mat.get()));
  BOUT_DO_PETSC(MatAXPY(*mat, 1.0, *rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
  return PetscOperator(this->mapping, std::move(mat));
}

/// Operator subtraction
PetscOperator PetscOperator::operator-(const PetscOperator& rhs) const {
  ASSERT0(this->mapping == rhs.mapping);
  UniqueMat mat(new Mat{nullptr});
  BOUT_DO_PETSC(MatDuplicate(*this->mat_operator, MAT_COPY_VALUES, mat.get()));
  BOUT_DO_PETSC(MatAXPY(*mat, -1.0, *rhs.mat_operator, UNKNOWN_NONZERO_PATTERN));
  return PetscOperator(this->mapping, std::move(mat));
}

PetscOperator PetscOperator::transpose() const {
  UniqueMat mat(new Mat{nullptr});
  BOUT_DO_PETSC(MatTranspose(*this->mat_operator, MAT_INITIAL_MATRIX, mat.get()));
  return PetscOperator(this->mapping, std::move(mat));
}

Field3D PetscOperators::meshGetField3D(Mesh* mesh, const std::string& name) {
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

  int mesh_total_cells{0};
  if (mesh->get(mesh_total_cells, "total_cells") == 0) {
    // Check total number of cells
    if (mesh_total_cells != mapping->globalSize()) {
      throw BoutException("Total cells in mesh {} doesn't match global mapping size {}",
                          mesh_total_cells, mapping->globalSize());
    }
  }
}

PetscOperator PetscOperators::get(const std::string& name) const {
  return PetscOperator(this->mapping, this->meshGetArray<int>(this->mesh, name + "_rows"),
                       this->meshGetArray<int>(this->mesh, name + "_columns"),
                       this->meshGetArray<BoutReal>(this->mesh, name + "_weights"));
}

PetscOperators::Parallel PetscOperators::getParallel() const {
  // Read maps from the mesh
  auto forward = this->get("forward");
  auto backward = this->get("backward");

  // ---- Construct Grad_par ----
  //
  // Create a parallel gradient operator by combining the parallel
  // length dl = dy * sqrt(g_22) with forward & backward operators
  auto* coords = this->mesh->getCoordinates();
  Field3D dl = coords->dy * sqrt(coords->g_22);
  dl.splitParallelSlices();
  dl.yup() = 0.0;
  dl.ydown() = 0.0;
  dl.applyParallelBoundary("parallel_neumann_o1");

  auto inv_2dl = 0.5 / dl;
  inv_2dl.splitParallelSlices();
  inv_2dl.yup() = 0.0;
  inv_2dl.ydown() = 0.0;
  inv_2dl.applyParallelBoundary("parallel_neumann_o1");

  auto inv_2dl_op = this->diagonal(inv_2dl);

  auto Grad_par = inv_2dl_op * (forward - backward);

  // ---- Construct Div_par ----
  //
  // Use the Support Operator Method (SOM) to calculate
  // Div_par from Grad_par and cell volumes.

  // Cell volume
  Field3D dV = coords->J * coords->dx * coords->dy * coords->dz;
  dV.splitParallelSlices();
  dV.yup() = 0.0;
  dV.ydown() = 0.0;
  dV.applyParallelBoundary("parallel_neumann_o1");

  // Fractional boundary cells
  const auto forward_boundary_fraction =
      meshGetField3D(this->mesh, "forward_boundary_fraction");
  const auto backward_boundary_fraction =
      meshGetField3D(this->mesh, "backward_boundary_fraction");

  BOUT_FOR(i, dV.getRegion("RGN_NOBNDRY")) {
    dV.yup()[i.yp()] *= forward_boundary_fraction[i];
    dV.ydown()[i.ym()] *= backward_boundary_fraction[i];
  }

  const auto dV_op = this->diagonal(dV);

  Field3D neg_inv_dV = -1. / dV;
  neg_inv_dV.splitParallelSlices();
  neg_inv_dV.yup() = 0.0;
  neg_inv_dV.ydown() = 0.0;
  neg_inv_dV.applyParallelBoundary("parallel_neumann_o1");
  auto neg_inv_dV_op = this->diagonal(neg_inv_dV);

  auto Div_par = neg_inv_dV_op * Grad_par.transpose() * dV_op;

  // ---- Construct Div_par_Grad_par ----
  //
  // Requires gradients between planes, at +1/2 and -1/2, and interpolation
  // operators to calculate quantities between cells.

  // Identity operator
  Field3D one{1.0};
  one.splitParallelSlices();
  one.yup() = 1.0;
  one.ydown() = 1.0;
  const auto identity = this->diagonal(one);

  // Interpolate at + 1/2
  auto interp_plus_op = (identity + forward) * 0.5;

  // dl averaged at +1/2
  const Field3D dl_plus = interp_plus_op(dl);
  Field3D inv_dl_plus = 1. / dl_plus;
  inv_dl_plus.splitParallelSlices();
  inv_dl_plus.yup() = 0.0;
  inv_dl_plus.ydown() = 0.0;
  inv_dl_plus.applyParallelBoundary("parallel_neumann_o1");
  const auto inv_dl_plus_op = this->diagonal(inv_dl_plus);

  // Gradient at + 1/2
  auto Grad_plus = inv_dl_plus_op * (forward - identity);

  // Divergence at -1/2
  auto Div_minus = neg_inv_dV_op * Grad_plus.transpose() * dV_op;

  // Interpolate at - 1/2
  auto interp_minus_op = (identity + backward) * 0.5;

  // dl averaged at -1/2
  const Field3D dl_minus = interp_minus_op(dl);
  Field3D inv_dl_minus = 1. / dl_minus;
  inv_dl_minus.splitParallelSlices();
  inv_dl_minus.yup() = 0.0;
  inv_dl_minus.ydown() = 0.0;
  inv_dl_minus.applyParallelBoundary("parallel_neumann_o1");
  const auto inv_dl_minus_op = this->diagonal(inv_dl_minus);

  // Gradient at - 1/2
  auto Grad_minus = inv_dl_minus_op * (identity - backward);

  // Divergence at +1/2
  auto Div_plus = neg_inv_dV_op * Grad_minus.transpose() * dV_op;

  // Div(Grad_par()) operator
  auto Div_par_Grad_par = ((Div_minus * Grad_plus) + (Div_plus * Grad_minus)) * 0.5;

  return Parallel{std::move(Grad_par),         std::move(Div_par),
                  std::move(Div_par_Grad_par), std::move(dV),
                  std::move(Grad_minus),       std::move(Grad_plus),
                  std::move(Div_minus),        std::move(Div_plus),
                  std::move(interp_minus_op),  std::move(interp_plus_op)};
}

Field3D PetscOperators::Parallel::Div_par_K_Grad_par(const Field3D& K,
                                                     const Field3D& f) const {

  // There are 4 matrix-vector products that could be performed in parallel.
  // The best way to optimize this is probably to pack f and K into one Vec,
  // and assemble a single sparse matrix that contains the four blocks:
  //
  //  (Grad_plus       0     )           (grad_plus  )
  //  (Grad_minus      0     ) (f)   ->  (grad_minus )
  //  (    0     Interp_plus ) (K)       (K_plus     )
  //  (    0     Interp_minus)           (K_minus    )
  //
  // For now we perform each operation in sequence.

  // Calculate gradients in + and - directions
  const Field3D grad_plus = this->Grad_plus(f);
  const Field3D grad_minus = this->Grad_minus(f);

  // Interpolate K to + and - locations
  const Field3D K_plus = this->Interp_plus(K);
  const Field3D K_minus = this->Interp_minus(K);

  // Calculate divergence
  return (this->Div_minus(K_plus * grad_plus) + this->Div_plus(K_minus * grad_minus))
         * 0.5;
}

#endif // BOUT_HAS_PETSC
