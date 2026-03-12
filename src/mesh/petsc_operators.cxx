#include "bout/petsc_operators.hxx"
#include "bout/boutexception.hxx"
#include "bout/field3d.hxx"
#include "bout/output.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/region.hxx"

#include <functional>
#include <string>
#include <vector>

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
    : evolving_region(create_region(cell_number)),
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
  PetscInt row_start, row_end;
  MatGetOwnershipRange(mat_mesh_to_petsc, &row_start, &row_end);
  output.write("Local row range: {} -> {}\n", row_start, row_end);

  PetscInt row = row_start;
  const PetscScalar ONE = 1.0;
  // Iterate through regions in this order
  const std::vector<std::reference_wrapper<Region<Ind3D>>> regions = {
      evolving_region, xin_region, xout_region, yup_region, ydown_region};
  for (const auto& region : regions) {
    BOUT_FOR_SERIAL(i, region.get()) {
      const PetscInt col = ROUND(cell_number[i]);
      MatSetValues(mat_mesh_to_petsc, 1, &row, 1, &col, &ONE, INSERT_VALUES);
      ++row;
    }
  }
  MatAssemblyBegin(mat_mesh_to_petsc, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(mat_mesh_to_petsc, MAT_FINAL_ASSEMBLY);
}

Field3D PetscOperators::meshGetField3D(Mesh* mesh, const std::string& name) {
  Field3D result;
  if (mesh->get(result, name) != 0) {
    throw BoutException("PetscOperators requires field '{}'", name);
  }
  return result;
}

PetscOperators::PetscOperators(Mesh* mesh)
    : mapping(PetscMapping(meshGetField3D(mesh, "cell_number"),
                           meshGetField3D(mesh, "forward_cell_number"),
                           meshGetField3D(mesh, "backward_cell_number"))) {

  int mesh_total_cells;
  if (mesh->get(mesh_total_cells, "total_cells") == 0) {
    // Check total number of cells
    if (mesh_total_cells != mapping.size()) {
      throw BoutException("Total cells in mesh {} doesn't match mapping size {}",
                          mesh_total_cells, mapping.size());
    }
  }

  //MatCreate(comm, &mat_yforward);
  //MatSetSizes(mat_yforward, nlocal, nlocal, PETSC_DECIDE, PETSC_DECIDE);
}
