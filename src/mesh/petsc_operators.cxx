#include "bout/build_defines.hxx"

#if BOUT_HAS_PETSC

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutexception.hxx"
#include "bout/coordinates.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/output_bout_types.hxx"
#include "bout/petsc_interface.hxx"
#include "bout/petsc_operators.hxx"
#include "bout/petsclib.hxx"
#include "bout/region.hxx"
#include "bout/utils.hxx"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

PetscIndexMapping::~PetscIndexMapping() {
  if (mat_stored_to_petsc != nullptr) {
    MatDestroy(&mat_stored_to_petsc);
  }
  if (mat_petsc_to_stored != nullptr) {
    MatDestroy(&mat_petsc_to_stored);
  }
}

PetscInt PetscIndexMapping::storedToPetsc(int stored_index) const {
  auto it = local_stored_to_petsc.find(stored_index);
  if (it == local_stored_to_petsc.end()) {
    throw BoutException("Stored index {} is not owned locally", stored_index);
  }
  return it->second;
}

void PetscIndexMapping::buildPermutation(PetscInt nlocal, PetscInt nglobal,
                                         const std::vector<int>& stored_indices) {
  global_size = nglobal;
  local_stored_indices = stored_indices;

  // Renumbering matrix.
  // Maps PETSc row (or column) indices to the global indices stored in
  // the mesh.  This is needed because the PETSc indices depend on the
  // number of processors.
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), &mat_stored_to_petsc));
  BOUT_DO_PETSC(MatSetSizes(mat_stored_to_petsc, nlocal, nlocal, nglobal, nglobal));
  BOUT_DO_PETSC(MatSetType(mat_stored_to_petsc, MATMPIAIJ));

  // Each row will have one non-zero entry, which could be in
  // either the "diagonal" or "off-diagonal" block.
  BOUT_DO_PETSC(MatMPIAIJSetPreallocation(mat_stored_to_petsc, 1, nullptr, 1, nullptr));

  // Get the range of global rows owned by this processor
  BOUT_DO_PETSC(MatGetOwnershipRange(mat_stored_to_petsc, &row_start, &row_end));
  ASSERT1(row_end - row_start == nlocal);

  const PetscScalar one = 1.0;
  for (PetscInt i = 0; i < nlocal; ++i) {
    const PetscInt petsc_row = row_start + i;
    const PetscInt stored_col = stored_indices[i];
    ASSERT1(stored_col >= 0);
    local_stored_to_petsc[static_cast<int>(stored_col)] = petsc_row;
    BOUT_DO_PETSC(MatSetValues(mat_stored_to_petsc, 1, &petsc_row, 1, &stored_col, &one,
                               INSERT_VALUES));
  }
  BOUT_DO_PETSC(MatAssemblyBegin(mat_stored_to_petsc, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(mat_stored_to_petsc, MAT_FINAL_ASSEMBLY));
  // Transpose to map petsc indices to stored mesh indices
  BOUT_DO_PETSC(
      MatTranspose(mat_stored_to_petsc, MAT_INITIAL_MATRIX, &mat_petsc_to_stored));
}

Region<Ind3D> PetscCellMapping::create_region(const Field3D& cell_number) {
  Region<Ind3D>::RegionIndices indices;
  BOUT_FOR_SERIAL(i, cell_number.getRegion("RGN_NOBNDRY")) {
    if (cell_number[i] >= 0) {
      indices.push_back(i);
    }
  }
  return Region<Ind3D>(indices);
}

Region<Ind3D> PetscCellMapping::create_region_xin(const Field3D& cell_number) {
  Region<Ind3D>::RegionIndices indices;
  const auto& region = cell_number.getRegion("RGN_INNER_X");
  BOUT_FOR_SERIAL(i, region) {
    if (cell_number[i] >= 0) {
      indices.push_back(i);
    }
  }
  return Region<Ind3D>(indices);
}

Region<Ind3D> PetscCellMapping::create_region_xout(const Field3D& cell_number) {
  Region<Ind3D>::RegionIndices indices;
  const auto& region = cell_number.getRegion("RGN_OUTER_X");
  BOUT_FOR_SERIAL(i, region) {
    if (cell_number[i] >= 0) {
      indices.push_back(i);
    }
  }
  return Region<Ind3D>(indices);
}

PetscCellMapping::PetscCellMapping(const Field3D& cell_number,
                                   const Field3D& forward_cell_number,
                                   const Field3D& backward_cell_number, int total_cells)
    : cell_number(cell_number), forward_cell_number(forward_cell_number),
      backward_cell_number(backward_cell_number),
      evolving_region(create_region(cell_number)),
      xin_region(create_region_xin(cell_number)),
      xout_region(create_region_xout(cell_number)),
      yup_region(create_region(forward_cell_number)),
      ydown_region(create_region(backward_cell_number)) {
  std::vector<int> local_indices;
  const auto push_region = [&](const Region<Ind3D>& region, const Field3D& values) {
    BOUT_FOR_SERIAL(i, region) { local_indices.push_back(ROUND(values[i])); }
  };
  push_region(evolving_region, this->cell_number);
  push_region(xin_region, this->cell_number);
  push_region(xout_region, this->cell_number);
  push_region(yup_region, this->forward_cell_number);
  push_region(ydown_region, this->backward_cell_number);
  buildPermutation(static_cast<PetscInt>(local_indices.size()), total_cells,
                   local_indices);
}

PetscLegMapping::PetscLegMapping(int total_legs, std::vector<int> local_leg_indices) {
  std::sort(local_leg_indices.begin(), local_leg_indices.end());
  local_leg_indices.erase(std::unique(local_leg_indices.begin(), local_leg_indices.end()),
                          local_leg_indices.end());
  buildPermutation(static_cast<PetscInt>(local_leg_indices.size()), total_legs,
                   local_leg_indices);
}

PetscCellVector makePetscCellVector(const PetscCellMappingPtr& mapping,
                                    const Field3D& field) {
  ASSERT1(field.hasParallelSlices());
  ASSERT1(field.yup().isAllocated());
  ASSERT1(field.ydown().isAllocated());
  PetscCellVector out(mapping);
  PetscScalar* x{nullptr};
  BOUT_DO_PETSC(VecGetArray(out.raw(), &x));
  mapping->mapLocalField([&](PetscInt row, const Ind3D& i) { x[row] = field[i]; });
  const auto& yup = field.yup();
  mapping->mapLocalYup([&](PetscInt row, const Ind3D& i) { x[row] = yup[i]; });
  const auto& ydown = field.ydown();
  mapping->mapLocalYdown([&](PetscInt row, const Ind3D& i) { x[row] = ydown[i]; });
  BOUT_DO_PETSC(VecRestoreArray(out.raw(), &x));
  return out;
}

Field3D toField3D(const PetscCellVector& vec, const Field3D& prototype) {
  auto mapping = std::static_pointer_cast<const PetscCellMapping>(vec.getMapping());
  Field3D result{emptyFrom(prototype)};
  result.splitParallelSlices();
  result.yup().allocate();
  result.ydown().allocate();

  const PetscScalar* x{nullptr};
  BOUT_DO_PETSC(VecGetArrayRead(vec.raw(), &x));
  mapping->mapLocalField([&](PetscInt row, const Ind3D& i) { result[i] = x[row]; });
  mapping->mapLocalYup([&](PetscInt row, const Ind3D& i) { result.yup()[i] = x[row]; });
  mapping->mapLocalYdown(
      [&](PetscInt row, const Ind3D& i) { result.ydown()[i] = x[row]; });
  BOUT_DO_PETSC(VecRestoreArrayRead(vec.raw(), &x));

  return result;
}

Field3D PetscOperators::meshGetField3D(Mesh* mesh, const std::string& name) {
  Field3D result;
  if (mesh->get(result, name) != 0) {
    throw BoutException("PetscOperators requires field '{}'", name);
  }
  return result;
}

int PetscOperators::meshGetInt(const std::string& name) const {
  int result{0};
  if (mesh->get(result, name) != 0) {
    throw BoutException("PetscOperators requires int '{}'", name);
  }
  return result;
}

PetscOperators::PetscOperators(Mesh* mesh)
    : mesh(mesh),
      cell_mapping(std::make_shared<PetscCellMapping>(
          meshGetField3D(mesh, "cell_number"),
          meshGetField3D(mesh, "forward_cell_number"),
          meshGetField3D(mesh, "backward_cell_number"), meshGetInt("total_cells"))),
      forward_leg_interior_number(meshGetField3D(mesh, "forward_leg_interior_number")),
      forward_leg_boundary_number(meshGetField3D(mesh, "forward_leg_boundary_number")),
      backward_leg_interior_number(meshGetField3D(mesh, "backward_leg_interior_number")),
      backward_leg_boundary_number(meshGetField3D(mesh, "backward_leg_boundary_number")) {
  std::vector<int> local_forward_legs;
  std::vector<int> local_backward_legs;
  cell_mapping->mapOwnedInteriorCells([&](PetscInt, const Ind3D& i, int) {
    const int f_int = ROUND(forward_leg_interior_number[i]);
    const int f_bnd = ROUND(forward_leg_boundary_number[i]);
    const int b_int = ROUND(backward_leg_interior_number[i]);
    const int b_bnd = ROUND(backward_leg_boundary_number[i]);
    if (f_int >= 0) {
      local_forward_legs.push_back(f_int);
    }
    if (f_bnd >= 0) {
      local_forward_legs.push_back(f_bnd);
    }
    if (b_int >= 0) {
      local_backward_legs.push_back(b_int);
    }
    if (b_bnd >= 0) {
      local_backward_legs.push_back(b_bnd);
    }
  });

  forward_leg_mapping = std::make_shared<PetscLegMapping>(meshGetInt("n_forward_legs"),
                                                          std::move(local_forward_legs));
  backward_leg_mapping = std::make_shared<PetscLegMapping>(
      meshGetInt("n_backward_legs"), std::move(local_backward_legs));
}

PetscForwardOperator PetscOperators::forward() const {
  return PetscForwardOperator(
      forward_leg_mapping, cell_mapping, meshGetArray<int>("forward_rows"),
      meshGetArray<int>("forward_columns"), meshGetArray<BoutReal>("forward_weights"));
}

PetscBackwardOperator PetscOperators::backward() const {
  return PetscBackwardOperator(
      backward_leg_mapping, cell_mapping, meshGetArray<int>("backward_rows"),
      meshGetArray<int>("backward_columns"), meshGetArray<BoutReal>("backward_weights"));
}

template <typename LegTag>
PetscOperator<LegTag, CellSpaceTag>
PetscOperators::buildInjection(const Field3D& interior_leg_number,
                               const Field3D& boundary_leg_number,
                               std::shared_ptr<const PetscLegMapping> leg_mapping) const {
  using Op = PetscOperator<LegTag, CellSpaceTag>;
  bout::petsc::UniqueMat mat{new Mat{nullptr}};
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), mat.get()));
  BOUT_DO_PETSC(MatSetType(*mat, MATMPIAIJ));
  BOUT_DO_PETSC(MatSetSizes(*mat, leg_mapping->localSize(), cell_mapping->localSize(),
                            leg_mapping->globalSize(), cell_mapping->globalSize()));
  BOUT_DO_PETSC(MatMPIAIJSetPreallocation(*mat, 1, nullptr, 1, nullptr));
  const PetscScalar one = 1.0;
  cell_mapping->mapOwnedInteriorCells([&](PetscInt cell_row, const Ind3D& i, int) {
    const int interior = ROUND(interior_leg_number[i]);
    const int boundary = ROUND(boundary_leg_number[i]);
    if (interior >= 0) {
      const PetscInt leg_row = leg_mapping->storedToPetsc(interior);
      BOUT_DO_PETSC(MatSetValues(*mat, 1, &leg_row, 1, &cell_row, &one, INSERT_VALUES));
    }
    if (boundary >= 0) {
      const PetscInt leg_row = leg_mapping->storedToPetsc(boundary);
      BOUT_DO_PETSC(MatSetValues(*mat, 1, &leg_row, 1, &cell_row, &one, INSERT_VALUES));
    }
  });
  BOUT_DO_PETSC(MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY));
  return Op(leg_mapping, cell_mapping, std::move(mat));
}

template PetscForwardOperator PetscOperators::buildInjection<ForwardLegSpaceTag>(
    const Field3D&, const Field3D&, std::shared_ptr<const PetscLegMapping>) const;
template PetscBackwardOperator PetscOperators::buildInjection<BackwardLegSpaceTag>(
    const Field3D&, const Field3D&, std::shared_ptr<const PetscLegMapping>) const;

PetscForwardOperator PetscOperators::injectForward() const {
  return buildInjection<ForwardLegSpaceTag>(
      forward_leg_interior_number, forward_leg_boundary_number, forward_leg_mapping);
}

PetscBackwardOperator PetscOperators::injectBackward() const {
  return buildInjection<BackwardLegSpaceTag>(
      backward_leg_interior_number, backward_leg_boundary_number, backward_leg_mapping);
}

PetscForwardLegVector PetscOperators::forwardLegWeights() const {
  return loadLegWeights<PetscForwardLegVector>("forward_leg_weights",
                                               forward_leg_mapping);
}

PetscBackwardLegVector PetscOperators::backwardLegWeights() const {
  return loadLegWeights<PetscBackwardLegVector>("backward_leg_weights",
                                                backward_leg_mapping);
}

PetscCellOperator PetscOperators::diagonal(const Field3D& f) const {
  return PetscCellOperator::diagonal(cell_mapping, makePetscCellVector(cell_mapping, f));
}

PetscForwardLegOperator PetscOperators::diagonal(const PetscForwardLegVector& v) const {
  return PetscForwardLegOperator::diagonal(forward_leg_mapping, v);
}

PetscBackwardLegOperator PetscOperators::diagonal(const PetscBackwardLegVector& v) const {
  return PetscBackwardLegOperator::diagonal(backward_leg_mapping, v);
}

PetscOperators::Parallel PetscOperators::getParallel() const {
  // Set boundaries
  const PetscCellOperator& Neumann_fwd =
      makeNeumannOperator(cell_mapping, BoundaryDirection::Forward); // C <- C
  const PetscCellOperator& Neumann_bwd =
      makeNeumannOperator(cell_mapping, BoundaryDirection::Backward); // C <- C

  // Primitive operators from the mesh / metadata
  auto Forward = this->forward();             // L+ <- C
  auto Backward = this->backward();           // L- <- C
  auto Inject_plus = this->injectForward();   // L+ <- C
  auto Inject_minus = this->injectBackward(); // L- <- C

  // Leg weights
  const auto w_plus_vec = forwardLegWeights();
  const auto w_minus_vec = backwardLegWeights();
  const auto W_plus = diagonal(w_plus_vec);   // L+ <- L+
  const auto W_minus = diagonal(w_minus_vec); // L- <- L-

  // Weighted adjoints: C <- L+ and C <- L-
  // The Neumann_f/bwd operators set the boundary cells to which the legs are connected
  auto Restrict_minus = Neumann_fwd * Inject_plus.transpose() * W_plus;
  auto Restrict_plus = Neumann_bwd * Inject_minus.transpose() * W_minus;

  auto* coords = mesh->getCoordinates();

  // Parallel spacing in cell space
  Field3D dl = Coordinates::FieldMetric{coords->dy() * sqrt(coords->g_22())};
  dl.splitParallelSlices();
  dl.yup() = 0.0;
  dl.ydown() = 0.0;
  dl.applyParallelBoundary("parallel_neumann_o1");

  // Cell volume
  Field3D dV =
      Coordinates::FieldMetric{coords->J() * coords->dx() * coords->dy() * coords->dz()};
  dV.splitParallelSlices();
  dV.yup() = 0.0;
  dV.ydown() = 0.0;
  dV.applyParallelBoundary("parallel_neumann_o1");

  Field3D neg_inv_dV = -1.0 / dV;
  neg_inv_dV.splitParallelSlices();
  neg_inv_dV.yup() = 0.0;
  neg_inv_dV.ydown() = 0.0;
  neg_inv_dV.applyParallelBoundary("parallel_neumann_o1");

  const auto DV = diagonal(dV);
  const auto Neg_inv_dV = diagonal(neg_inv_dV);

  // Unweighted interpolation to +1/2 and -1/2 legs
  auto Interp_plus = (Forward + Inject_plus) * 0.5;    // L+ <- C
  auto Interp_minus = (Backward + Inject_minus) * 0.5; // L- <- C

  // Leg-centered dl
  const auto dl_plus = Interp_plus(dl);
  const auto dl_minus = Interp_minus(dl);

  const auto inv_dl_plus = PetscForwardLegVector::reciprocal(dl_plus);
  const auto inv_dl_minus = PetscBackwardLegVector::reciprocal(dl_minus);

  const auto Inv_dl_plus = diagonal(inv_dl_plus);   // L+ <- L+
  const auto Inv_dl_minus = diagonal(inv_dl_minus); // L- <- L-

  // Half-step gradients
  auto Grad_plus = Inv_dl_plus * (Forward - Inject_plus);     // L+ <- C
  auto Grad_minus = Inv_dl_minus * (Inject_minus - Backward); // L- <- C

  // Cell-centered central gradient
  // Apply Neumann so that the yup/ydown cells are set. This is needed
  // so that Div_par has the correct dependence on yup/down cells.
  auto Grad_par =
      ((Restrict_minus * Grad_plus) + (Restrict_plus * Grad_minus)) * 0.5; // C <- C

  // Leg-centered volumes from unweighted interpolation
  const auto dV_plus = Interp_plus(dV);   // L+
  const auto dV_minus = Interp_minus(dV); // L-

  // Physical leg volumes include the leg weights
  const auto dV_leg_plus = PetscForwardLegVector::pointwiseMultiply(dV_plus, w_plus_vec);
  const auto dV_leg_minus =
      PetscBackwardLegVector::pointwiseMultiply(dV_minus, w_minus_vec);

  const auto DV_leg_plus = diagonal(dV_leg_plus);   // L+ <- L+
  const auto DV_leg_minus = diagonal(dV_leg_minus); // L- <- L-

  // Support-operator divergence from half-step gradients
  auto Div_minus = Neg_inv_dV * Grad_plus.transpose() * DV_leg_plus;  // C <- L+
  auto Div_plus = Neg_inv_dV * Grad_minus.transpose() * DV_leg_minus; // C <- L-

  // Cell-centered divergence paired with the cell-centered gradient
  auto Div_par = Neg_inv_dV * Grad_par.transpose() * DV; // C <- C

  // Diffusion operator assembled from the half-step pieces
  auto Div_par_Grad_par =
      ((Div_minus * Grad_plus) + (Div_plus * Grad_minus)) * 0.5; // C <- C

  return Parallel{// Cell <- Cell operators
                  std::move(Grad_par), std::move(Div_par), std::move(Div_par_Grad_par),
                  std::move(dV),

                  // Leg <-> Cell operators
                  std::move(Grad_minus), std::move(Grad_plus), std::move(Div_minus),
                  std::move(Div_plus),

                  std::move(Backward), std::move(Forward), std::move(Inject_minus),
                  std::move(Inject_plus),

                  std::move(Interp_minus), std::move(Interp_plus),

                  std::move(Restrict_minus), std::move(Restrict_plus)};
}

Field3D PetscOperators::Parallel::Div_par_K_Grad_par(const Field3D& K,
                                                     const Field3D& f) const {
  const auto grad_plus = Grad_plus(f);
  const auto grad_minus = Grad_minus(f);
  const auto K_plus = Interp_plus(K);
  const auto K_minus = Interp_minus(K);

  const auto flux_plus = PetscForwardLegVector::pointwiseMultiply(K_plus, grad_plus);
  const auto flux_minus = PetscBackwardLegVector::pointwiseMultiply(K_minus, grad_minus);

  const auto div_minus = Div_minus(flux_plus);
  const auto div_plus = Div_plus(flux_minus);
  return (toField3D(div_minus, K) + toField3D(div_plus, K)) * 0.5;
}

PetscCellOperator makeNeumannOperator(const PetscCellMappingPtr& mapping,
                                      BoundaryDirection direction) {

  // The operator is built using stored indices for both rows and columns,
  // then post-multiplied by getPetscToStored() to convert columns to PETSc
  // ordering — exactly the same pattern as assembleFromCSR.

  const PetscInt n_global = mapping->globalSize();
  const PetscInt n_local = mapping->localSize();

  // Temporary matrix in stored-column ordering.
  // Each row has at most one nonzero (either diagonal or the interior
  // neighbour), so d_nz=1, o_nz=1 is the correct preallocation.
  bout::petsc::UniqueMat temp{new Mat{nullptr}};
  BOUT_DO_PETSC(MatCreate(BoutComm::get(), temp.get()));
  BOUT_DO_PETSC(MatSetType(*temp, MATMPIAIJ));
  BOUT_DO_PETSC(MatSetSizes(*temp, n_local, PETSC_DECIDE, n_global, n_global));
  BOUT_DO_PETSC(MatMPIAIJSetPreallocation(*temp, 1, nullptr, 1, nullptr));

  const PetscScalar one = 1.0;

  // Iterate over all locally owned cells using mapLocalField (interior +
  // X-boundary) and mapLocalYup / mapLocalYdown (parallel boundary virtual
  // cells).  For each cell we insert exactly one entry.
  //
  // mapLocalField visits cells whose stored index comes from cell_number.
  // mapLocalYup visits cells whose stored index comes from forward_cell_number.
  // mapLocalYdown visits cells whose stored index comes from backward_cell_number.
  //
  // For a forward boundary virtual cell at field index i.yp():
  //   row = stored index from forward_cell_number (the virtual cell)
  //   col = stored index from cell_number at the source cell i
  //   => virtual cell copies from the adjacent interior cell (Neumann)
  //
  // For a backward boundary virtual cell at field index i.ym():
  //   row = stored index from backward_cell_number (the virtual cell)
  //   col = stored index from cell_number at the source cell i
  //
  // For all other cells the row is the diagonal (identity).

  // Identity rows: all cells in the interior and X-boundary regions
  mapping->mapLocalField([&](PetscInt local_row, const Ind3D& i) {
    const PetscInt petsc_row = mapping->rowStart() + local_row;
    const int stored_col = ROUND(mapping->cell_number[i]);
    if (stored_col < 0) {
      throw BoutException("Check failed: stored_col, mapLocalField");
    }
    const PetscInt col = static_cast<PetscInt>(stored_col);
    BOUT_DO_PETSC(MatSetValues(*temp, 1, &petsc_row, 1, &col, &one, INSERT_VALUES));
  });

  if ((direction == BoundaryDirection::Forward)
      or (direction == BoundaryDirection::Both)) {
    // Set forward rows from interior columns
    mapping->mapLocalYup([&](PetscInt local_row, const Ind3D& i_yp) {
      const PetscInt petsc_row = mapping->rowStart() + local_row;
      const Ind3D i_src = i_yp.ym();
      const int stored_row_check = ROUND(mapping->forward_cell_number[i_src]);
      if (stored_row_check < 0) {
        throw BoutException("Check failed: stored_row, mapLocalYup");
      }
      const int stored_col = ROUND(mapping->cell_number[i_src]);
      if (stored_col < 0) {
        throw BoutException("Check failed: stored_col, mapLocalYup");
      }
      const PetscInt col = static_cast<PetscInt>(stored_col);
      BOUT_DO_PETSC(MatSetValues(*temp, 1, &petsc_row, 1, &col, &one, INSERT_VALUES));
    });
  } else {
    // Identity rows in forward boundary
    mapping->mapLocalYup([&](PetscInt local_row, const Ind3D& i_yp) {
      const PetscInt petsc_row = mapping->rowStart() + local_row;
      const Ind3D i_src = i_yp.ym();
      const int stored_col = ROUND(mapping->forward_cell_number[i_src]);
      if (stored_col < 0) {
        throw BoutException("Check failed: stored_col, mapLocalYup identity");
      }
      const PetscInt col = static_cast<PetscInt>(stored_col);
      BOUT_DO_PETSC(MatSetValues(*temp, 1, &petsc_row, 1, &col, &one, INSERT_VALUES));
    });
  }

  if ((direction == BoundaryDirection::Backward)
      or (direction == BoundaryDirection::Both)) {
    // Set backward rows from interior columns
    mapping->mapLocalYdown([&](PetscInt local_row, const Ind3D& i_ym) {
      const PetscInt petsc_row = mapping->rowStart() + local_row;
      const Ind3D i_src = i_ym.yp();
      const int stored_row_check = ROUND(mapping->backward_cell_number[i_src]);
      if (stored_row_check < 0) {
        throw BoutException("Check failed: stored_row, mapLocalYdown");
      }
      const int stored_col = ROUND(mapping->cell_number[i_src]);
      if (stored_col < 0) {
        throw BoutException("Check failed: stored_col, mapLocalYdown");
      }
      const PetscInt col = static_cast<PetscInt>(stored_col);
      BOUT_DO_PETSC(MatSetValues(*temp, 1, &petsc_row, 1, &col, &one, INSERT_VALUES));
    });
  } else {
    // Identity rows in backward boundary
    mapping->mapLocalYdown([&](PetscInt local_row, const Ind3D& i_ym) {
      const PetscInt petsc_row = mapping->rowStart() + local_row;
      const Ind3D i_src = i_ym.yp();
      const int stored_col = ROUND(mapping->backward_cell_number[i_src]);
      if (stored_col < 0) {
        throw BoutException("Check failed: stored_col, mapLocalYdown identity");
      }
      const PetscInt col = static_cast<PetscInt>(stored_col);
      BOUT_DO_PETSC(MatSetValues(*temp, 1, &petsc_row, 1, &col, &one, INSERT_VALUES));
    });
  }

  BOUT_DO_PETSC(MatAssemblyBegin(*temp, MAT_FINAL_ASSEMBLY));
  BOUT_DO_PETSC(MatAssemblyEnd(*temp, MAT_FINAL_ASSEMBLY));

  // Post-multiply by getPetscToStored() to convert stored column indices to
  // PETSc column indices — identical to the final step in assembleFromCSR.
  bout::petsc::UniqueMat out{new Mat{nullptr}};
  BOUT_DO_PETSC(MatMatMult(*temp, mapping->getPetscToStored(), MAT_INITIAL_MATRIX,
                           PETSC_DETERMINE, out.get()));

  return PetscCellOperator(mapping, mapping, std::move(out));
}

#endif
