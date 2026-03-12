#include "bout/petsc_operators.hxx"
#include "bout/boutexception.hxx"
#include "bout/output_bout_types.hxx"

Region<Ind3D> PetscOperators::create_region(Field3D cell_number) {
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, cell_number.getRegion("RGN_NOBNDRY")) {
    if (cell_number[i] > -1) {
      indices.push_back(i);
    }
  }
  return Region<Ind3D>(indices);
}

Region<Ind3D> PetscOperators::create_region_xin(Field3D cell_number) {
  Region<Ind3D>::RegionIndices xin_indices;
  const Mesh* mesh = cell_number.getMesh();
  if (mesh->firstX()) {
    for (int i = 0; i < mesh->xstart; ++i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = mesh->zstart; k <= mesh->zend; ++k) {
          if (cell_number(i, j, k) > 0) {
            xin_indices.push_back(Ind3D{i, j, k});
          }
        }
      }
    }
  }
  return Region<Ind3D>(xin_indices);
}

Region<Ind3D> PetscOperators::create_region_xout(Field3D cell_number) {
  Region<Ind3D>::RegionIndices xout_indices;
  const Mesh* mesh = cell_number.getMesh();
  if (mesh->lastX()) {
    for (int i = mesh->xend + 1; i < mesh->LocalNx; ++i) {
      for (int j = mesh->ystart; j <= mesh->yend; ++j) {
        for (int k = mesh->zstart; k <= mesh->zend; ++k) {
          if (cell_number(i, j, k) > 0) {
            xout_indices.push_back(Ind3D{i, j, k});
          }
        }
      }
    }
  }
  return Region<Ind3D>(xout_indices);
}

PetscOperators::PetscOperators(Mesh* mesh) {

  // Get the cell center global numbering
  // Note: This is an integer but Field3D only stores BoutReal.
  Field3D cell_number;
  if (mesh->get(cell_number, "cell_number") != 0) {
    throw BoutException("PetscOperators requires field 'cell_number'");
  }

  // Forward Y boundary cell numbers.
  // >0 only if a forward boundary cell
  Field3D forward_cell_number;
  if (mesh->get(forward_cell_number, "forward_cell_number") != 0) {
    throw BoutException("PetscOperators requires field 'forward_cell_number'");
  }

  // Backward Y boundary cell numbers.
  // >0 only if a backward boundary cell
  Field3D backward_cell_number;
  if (mesh->get(backward_cell_number, "backward_cell_number") != 0) {
    throw BoutException("PetscOperators requires field 'backward_cell_number'");
  }

  // Calculate Regions for evolving variables and boundaries
  evolving_region = create_region(cell_number);
  xin_region = create_region_xin(cell_number);
  xout_region = create_region_xout(cell_number);
  yup_region = create_region(forward_cell_number);
  ydown_region = create_region(backward_cell_number);

  const std::size_t nlocal = evolving_region.size() + xin_region.size()
                             + xout_region.size() + yup_region.size()
                             + ydown_region.size();
  output.write("Local {} : evolving {} xin {} xout {} yup {} ydown {}", nlocal,
               evolving_region.size(), xin_region.size(), xout_region.size(),
               yup_region.size(), ydown_region.size());

  int mesh_total_cells;
  if (mesh->get(mesh_total_cells, "total_cells") == 0) {
    // Check total number of cells
    output.write("Total cells in mesh: {}\n", mesh_total_cells);
  }

  // Create a PETSc matrix
  // Note: Numbering is different from that used in PetscVector / PetscMatrix
  // because yup/down boundaries are included.

  MPI_Comm comm = BoutComm::get();

  // Renumbering matrix.
  // Maps PETSc row (or column) indices to the global indices used in
  // the mesh.  This is needed because the PETSc indices depend on the
  // number of processors.
  MatCreate(comm, &mat_mesh_to_petsc);
  MatSetSizes(mat_mesh_to_petsc, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE);
  MatSetType(mat_mesh_to_petsc, MATMPIAIJ);
  // Each row will have one non-zero entry, which could be in
  // either the "diagonal" or "off-diagonal" block.
  MatMPIAIJSetPreallocation(mat_mesh_to_petsc, 1, nullptr, 1, nullptr);

  // Get the range of rows owned by this processor
  PetscInt row_start, row_end;
  MatGetOwnershipRange(mat_mesh_to_petsc, &row_start, &row_end);
  output.write("Local row range: {} -> {}\n", row_start, row_end);

  PetscInt row = row_start;
  const PetscScalar ONE = 1.0;
  BOUT_FOR_SERIAL(i, evolving_region) {
    PetscInt col = ROUND(cell_number[i]);
    MatSetValues(mat_mesh_to_petsc, 1, &row, 1, &col, &ONE, INSERT_VALUES);
    ++row;
  }

  //MatCreate(comm, &mat_yforward);
  //MatSetSizes(mat_yforward, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE);
}
