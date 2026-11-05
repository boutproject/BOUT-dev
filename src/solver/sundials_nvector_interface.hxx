/**************************************************************************
 * Shared helpers for SUNDIALS N_Vector backends used by BOUT++ solvers.
 *
 * Copyright 2010-2026 BOUT++ contributors
 *
 * Contact: Ben Dudson, dudson2@llnl.gov
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

#ifndef BOUT_SUNDIALS_NVECTOR_INTERFACE_HXX
#define BOUT_SUNDIALS_NVECTOR_INTERFACE_HXX

#include "bout/bout_enum_class.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutcomm.hxx"
#include "bout/boutexception.hxx"
#include "bout/build_defines.hxx"
#include "bout/field2d.hxx"
#include "bout/globals.hxx"
#include "bout/region.hxx"
#include "bout/solver.hxx"
#include "bout/sundials_backports.hxx"
#include "bout/unused.hxx"

#if BOUT_HAS_SUNDIALS_MANYVECTOR
#include "nvector.hxx"
#endif

#include "sundials_cuda_packing.hxx"

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

/// Supported SUNDIALS ``N_Vector`` backends that can be selected at runtime.
BOUT_ENUM_CLASS(NVectorType, Sundials, ManyVector, Cuda);

/// Helper for moving solver state between BOUT++ fields and a chosen ``N_Vector``
/// backend without exposing backend-specific details in each solver wrapper.
class SundialsNVectorInterface {
public:
  SundialsNVectorInterface(Solver& solver_in, const sundials::Context& ctx_in,
                           NVectorType nvector_type_in,
                           std::shared_ptr<SundialsCudaPacking> cuda_packing_in = {})
      : solver(solver_in), ctx(ctx_in), nvector_type(nvector_type_in),
        cuda_packing(std::move(cuda_packing_in)) {}

  /// Return ``true`` when using the field-backed ManyVector implementation.
  bool use_manyvector() const { return nvector_type == NVectorType::ManyVector; }

  bool use_cuda() const { return nvector_type == NVectorType::Cuda; }

  /// Human-readable backend name for logs and error messages.
  const char* backend_name() const {
    if (use_cuda()) {
      return "SUNDIALS CUDA + MPIPlusX (GPU vector operations)";
    }
    return use_manyvector() ? "SUNDIALS ManyVector-backed custom" : "SUNDIALS parallel";
  }

  /// Reject backends unavailable in this build or solver wrapper.
  void ensure_manyvector_available() const {
    if (use_cuda()) {
#if BOUT_HAS_SUNDIALS_CUDA
      if (!cuda_packing) {
        throw BoutException("solver:nvector=cuda is currently supported only by CVODE");
      }
#else
      throw BoutException("solver:nvector=cuda requires BOUT++ CUDA/Umpire support and "
                          "SUNDIALS CUDA and MPIPlusX imported targets");
#endif
    }
    if (use_manyvector() and not BOUT_HAS_SUNDIALS_MANYVECTOR) {
      throw BoutException("SUNDIALS ManyVector backend requested, but BOUT++ was "
                          "built without SUNDIALS ManyVector support");
    }
  }

  /// Create the main solver state vector using the selected backend.
  N_Vector create_state_vector(int local_N, int neq) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      auto* vec = nvector_from_state();
      if (vec == nullptr) {
        throw BoutException("BOUT N_Vector failed\n");
      }
      return vec;
#endif
    }

#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      initialise_cuda_packing(local_N);
      auto* local = callWithSUNContext(N_VNew_Cuda, ctx, local_N);
      if (!local) {
        throw BoutException("SUNDIALS CUDA vector allocation failed");
      }
      try {
        cuda_packing->bind(local);
      } catch (...) {
        N_VDestroy(local);
        throw;
      }
      auto* vec = callWithSUNContext(N_VMake_MPIPlusX, ctx, BoutComm::get(), local);
      if (!vec) {
        N_VDestroy(local);
        throw BoutException("SUNDIALS MPIPlusX vector allocation failed");
      }
      // N_VMake_MPIPlusX borrows X; transfer ownership so N_VDestroy frees it.
      static_cast<N_VectorContent_MPIPlusX>(vec->content)->own_data = SUNTRUE;
      try {
        copy_state_to_vector(vec);
      } catch (...) {
        N_VDestroy(vec);
        throw;
      }
      return vec;
    }
#endif

    auto* vec = callWithSUNContext(N_VNew_Parallel, ctx, BoutComm::get(), local_N, neq);
    if (vec == nullptr) {
      throw BoutException("SUNDIALS memory allocation failed\n");
    }
    solver.save_vars(N_VGetArrayPointer(vec));
    return vec;
  }

  /// Clone a vector of the same backend/layout, with a descriptive error on failure.
  N_Vector clone_vector_like(N_Vector source, const char* description) const {
    auto* vec = N_VClone(source);
    if (vec == nullptr) {
      throw BoutException("{}\n", description);
    }
    return vec;
  }

  /// Copy the current solution state from ``u`` into the solver fields.
  void copy_state_from_vector(N_Vector u) const {
    ensure_manyvector_available();
#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      transfer_cuda(u, false, true);
      return;
    }
#endif

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_state(u);
      return;
#endif
    }

    solver.load_vars(N_VGetArrayPointer(u));
  }

  /// Copy the solver fields into the solution vector ``u``.
  void copy_state_to_vector(N_Vector u) const {
    ensure_manyvector_available();
#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      transfer_cuda(u, false, false);
      return;
    }
#endif

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_state(u);
      return;
#endif
    }

    solver.save_vars(N_VGetArrayPointer(u));
  }

  /// Copy time-derivative values from ``du`` into the solver derivative fields.
  void copy_deriv_from_vector(N_Vector du) const {
    ensure_manyvector_available();
#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      transfer_cuda(du, true, true);
      return;
    }
#endif

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_deriv(du);
      return;
#endif
    }

    solver.load_derivs(N_VGetArrayPointer(du));
  }

  /// Copy the solver derivative fields into ``du``.
  void copy_deriv_to_vector(N_Vector du) const {
    ensure_manyvector_available();
#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      transfer_cuda(du, true, false);
      return;
    }
#endif

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      swap_deriv(du);
      return;
#endif
    }

    solver.save_derivs(N_VGetArrayPointer(du));
  }

  /// Save callback outputs together. Prepare both sets of fields before any
  /// asynchronous reads, since vector conversion may replace component storage.
  void copy_state_and_deriv_to_vectors(N_Vector u, N_Vector du) const {
    ensure_manyvector_available();
#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda() && cuda_packing->async_transfers) {
      const bout::profiling::ProfileRange range(
          "SundialsNVectorInterface::save_cuda_pair");
      cuda_packing->wait_for_fields();
      // User RHS/preconditioner/Jacobian functions may use arbitrary streams.
      // Keep one producer fence for the pair rather than one per transfer.
      SundialsCudaPacking::check(cudaDeviceSynchronize());
      prepare_cuda_fields(false, false);
      prepare_cuda_fields(true, false);
      if (!solver.v2d.empty() || !solver.v3d.empty()) {
        SundialsCudaPacking::check(cudaDeviceSynchronize());
      }
      try {
        queue_cuda_transfer(u, false, false, true);
        queue_cuda_transfer(du, true, false, true);
        cuda_packing->record_field_reads();
      } catch (...) {
        cuda_packing->synchronize_noexcept();
        throw;
      }
      return;
    }
#endif
    copy_state_to_vector(u);
    copy_deriv_to_vector(du);
  }

  /// Finish pending field reads before monitors or other host field mutations.
  void wait_for_cuda_fields() const {
#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      cuda_packing->wait_for_fields();
    }
#endif
  }

  /// Fill every evolved entry in ``vec`` with per-field constants, respecting
  /// whether each field evolves boundary cells.
  void fill_vector_values(N_Vector vec, const std::vector<BoutReal>& f2d_values,
                          const std::vector<BoutReal>& f3d_values) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      std::size_t i = 0;
      for (std::size_t j = 0; j < f2d_values.size(); ++j, ++i) {
        fill_field(BoutNVector::get<Field2D>(vec, i), f2d_values[j],
                   solver.f2d[j].evolve_bndry);
      }
      for (std::size_t j = 0; j < f3d_values.size(); ++j, ++i) {
        fill_field(BoutNVector::get<Field3D>(vec, i), f3d_values[j],
                   solver.f3d[j].evolve_bndry);
      }
      return;
#endif
    }

#if BOUT_HAS_SUNDIALS_CUDA
    if (use_cuda()) {
      auto* local = N_VGetLocalVector_MPIPlusX(vec);
      auto* data = N_VGetHostArrayPointer_Cuda(local);
      int p = 0;
      for (const auto& i : bout::globals::mesh->getRegion2D("RGN_BNDRY")) {
        fill_vector_values_op(i, data, p, f2d_values, f3d_values, true);
      }
      for (const auto& i : bout::globals::mesh->getRegion2D("RGN_NOBNDRY")) {
        fill_vector_values_op(i, data, p, f2d_values, f3d_values, false);
      }
      cuda_packing->check_stream(local);
      SundialsCudaPacking::check(cudaMemcpyAsync(
          N_VGetDeviceArrayPointer_Cuda(local), data, p * sizeof(BoutReal),
          cudaMemcpyHostToDevice, cuda_packing->stream));
      // Setup vectors may be destroyed immediately by the caller.
      cuda_packing->synchronize();
      return;
    }
#endif

    BoutReal* option_data = N_VGetArrayPointer(vec);
    int p = 0;
    for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_BNDRY")) {
      fill_vector_values_op(i2d, option_data, p, f2d_values, f3d_values, true);
    }
    for (const auto& i2d : bout::globals::mesh->getRegion2D("RGN_NOBNDRY")) {
      fill_vector_values_op(i2d, option_data, p, f2d_values, f3d_values, false);
    }
  }

  /// Fill the IDA ``id`` vector with ``1`` for differential variables and ``0``
  /// for constrained algebraic variables.
  void fill_id_vector(N_Vector vec) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      std::size_t i = 0;
      for (const auto& field : solver.f2d) {
        fill_field(BoutNVector::get<Field2D>(vec, i), field.constraint ? 0.0 : 1.0,
                   field.evolve_bndry);
        ++i;
      }
      for (const auto& field : solver.f3d) {
        fill_field(BoutNVector::get<Field3D>(vec, i), field.constraint ? 0.0 : 1.0,
                   field.evolve_bndry);
        ++i;
      }
      return;
#endif
    }

    solver.set_id(N_VGetArrayPointer(vec));
  }

  /// Apply the IDA residual update ``residual -= du`` only on differential
  /// components selected by ``id``.
  void subtract_differential_term(N_Vector residual, N_Vector du, N_Vector id) const {
    ensure_manyvector_available();

    if (use_manyvector()) {
#if BOUT_HAS_SUNDIALS_MANYVECTOR
      std::size_t i = 0;
      for (const auto& field : solver.f2d) {
        subtract_field(BoutNVector::get<Field2D>(residual, i),
                       BoutNVector::get<Field2D>(du, i), BoutNVector::get<Field2D>(id, i),
                       field.evolve_bndry);
        ++i;
      }
      for (const auto& field : solver.f3d) {
        subtract_field(BoutNVector::get<Field3D>(residual, i),
                       BoutNVector::get<Field3D>(du, i), BoutNVector::get<Field3D>(id, i),
                       field.evolve_bndry);
        ++i;
      }
      return;
#endif
    }

    const auto length = N_VGetLocalLength_Parallel(id);
    const BoutReal* id_data = N_VGetArrayPointer(id);
    BoutReal* residual_data = N_VGetArrayPointer(residual);
    const BoutReal* du_data = N_VGetArrayPointer(du);
    for (int i = 0; i < length; ++i) {
      if (id_data[i] > 0.5) {
        residual_data[i] -= du_data[i];
      }
    }
  }

private:
  Solver& solver;
  const sundials::Context& ctx;
  NVectorType nvector_type;
  std::shared_ptr<SundialsCudaPacking> cuda_packing;

#if BOUT_HAS_SUNDIALS_CUDA
  void initialise_cuda_packing(int local_N) const {
    if (cuda_packing->ready) {
      return;
    }
    cuda_packing->fields.resize(solver.f2d.size() + solver.f3d.size());
    int slot = 0;
    const auto append = [&](const auto& region, bool boundary) {
      for (const auto& i : region) {
        for (std::size_t f = 0; f < solver.f2d.size(); ++f) {
          if (!boundary || solver.f2d[f].evolve_bndry) {
            cuda_packing->fields[f].append(i.ind, slot++);
          }
        }
        for (int z = 0; z < bout::globals::mesh->LocalNz; ++z) {
          for (std::size_t f = 0; f < solver.f3d.size(); ++f) {
            if (!boundary || solver.f3d[f].evolve_bndry) {
              const auto index = solver.f3d[f].var->getMesh()->ind2Dto3D(i, z);
              cuda_packing->fields[solver.f2d.size() + f].append(index.ind, slot++);
            }
          }
        }
      }
    };
    append(bout::globals::mesh->getRegion2D("RGN_BNDRY"), true);
    append(bout::globals::mesh->getRegion2D("RGN_NOBNDRY"), false);
    if (slot != local_N) {
      throw BoutException("SUNDIALS CUDA packing length mismatch");
    }
    for (auto& field : cuda_packing->fields) {
      field.upload();
    }
    cuda_packing->ready = true;
  }

  void prepare_cuda_fields(bool derivatives, bool load) const {
    // Preserve the same vector basis and field preparation as Solver::load/save_*.
    const auto prepare_vectors = [&](auto& vectors) {
      for (const auto& v : vectors) {
        auto* field = derivatives ? v.F_var : v.var;
        if (load) {
          field->covariant = v.covariant;
        } else if (v.covariant) {
          field->toCovariant();
        } else {
          field->toContravariant();
        }
      }
    };
    prepare_vectors(solver.v2d);
    prepare_vectors(solver.v3d);
  }

  void transfer_cuda(N_Vector vec, bool derivatives, bool load) const {
    const bout::profiling::ProfileRange range(
        load ? "SundialsNVectorInterface::load_cuda"
             : "SundialsNVectorInterface::save_cuda");
    cuda_packing->wait_for_fields();
    // Loads consume vectors on our shared stream. They do not need a device-wide
    // fence before preparation, which only sets vector metadata when loading.
    if (!cuda_packing->async_transfers || !load) {
      SundialsCudaPacking::check(cudaDeviceSynchronize());
    }
    prepare_cuda_fields(derivatives, load);
    if (!cuda_packing->async_transfers
        || (!load && (!solver.v2d.empty() || !solver.v3d.empty()))) {
      SundialsCudaPacking::check(cudaDeviceSynchronize());
    }
    try {
      queue_cuda_transfer(vec, derivatives, load);
      if (cuda_packing->async_transfers) {
        // Keep standalone saves synchronous too: callers may mutate fields as
        // soon as this method returns. Only paired callback saves defer the wait.
        cuda_packing->synchronize();
      } else {
        SundialsCudaPacking::check(cudaDeviceSynchronize());
      }
    } catch (...) {
      cuda_packing->synchronize_noexcept();
      throw;
    }
  }

  void queue_cuda_transfer(N_Vector vec, bool derivatives, bool load,
                           bool retain_fields = false) const {
    auto* local = N_VGetLocalVector_MPIPlusX(vec);
    cuda_packing->check_stream(local);
    auto* data = N_VGetDeviceArrayPointer_Cuda(local);
    std::size_t f = 0;
    for (auto& v : solver.f2d) {
      auto* field = derivatives ? v.F_var : v.var;
      if (load) {
        field->allocate();
      }
      if (!field->isAllocated()) {
        throw BoutException("Unallocated CUDA solver field");
      }
      if (retain_fields) {
        cuda_packing->retained_2d.push_back(*field);
      }
      cuda_packing->fields[f++].transfer(&(*field)(0, 0), data, load,
                                         cuda_packing->stream);
    }
    for (auto& v : solver.f3d) {
      auto* field = derivatives ? v.F_var : v.var;
      if (load) {
        field->allocate();
        field->setLocation(v.location);
      }
      if (!field->isAllocated()) {
        throw BoutException("Unallocated CUDA solver field");
      }
      if (derivatives && !load && field->getLocation() != v.var->getLocation()) {
        throw BoutException("CUDA solver derivative at wrong cell location");
      }
      if (retain_fields) {
        cuda_packing->retained_3d.push_back(*field);
      }
      cuda_packing->fields[f++].transfer(&(*field)(0, 0, 0), data, load,
                                         cuda_packing->stream);
    }
  }
#endif

  void fill_vector_values_op(Ind2D UNUSED(i2d), BoutReal* option_data, int& p,
                             const std::vector<BoutReal>& f2d_values,
                             const std::vector<BoutReal>& f3d_values, bool bndry) const {
    for (std::size_t i = 0; i < f2d_values.size(); ++i) {
      if (bndry && !solver.f2d[i].evolve_bndry) {
        continue;
      }
      option_data[p] = f2d_values[i];
      ++p;
    }

    for (int jz = 0; jz < bout::globals::mesh->LocalNz; ++jz) {
      for (std::size_t i = 0; i < f3d_values.size(); ++i) {
        if (bndry && !solver.f3d[i].evolve_bndry) {
          continue;
        }
        option_data[p] = f3d_values[i];
        ++p;
      }
    }
  }

#if BOUT_HAS_SUNDIALS_MANYVECTOR
  /// Build a field-backed ManyVector that aliases the solver state variables.
  N_Vector nvector_from_state() const {
    std::vector<N_Vector> subvectors;
    subvectors.reserve(solver.f2d.size() + solver.f3d.size());
    const auto inserter = std::back_inserter(subvectors);

    const auto field_to_nvector = [this](const auto& var_str) {
      return BoutNVector::create(ctx, *var_str.var, var_str.evolve_bndry);
    };
    std::transform(solver.f2d.cbegin(), solver.f2d.cend(), inserter, field_to_nvector);
    std::transform(solver.f3d.cbegin(), solver.f3d.cend(), inserter, field_to_nvector);
    return BoutNVector::create(ctx, subvectors);
  }

  /// Swap the active solver state fields with the subvectors contained in ``u``.
  void swap_state(const N_Vector u) const {
    std::size_t i = 0;
    for (auto& var_str : solver.f2d) {
      BoutNVector::swap(u, *var_str.var, i);
      ++i;
    }
    for (auto& var_str : solver.f3d) {
      BoutNVector::swap(u, *var_str.var, i);
      ++i;
    }
  }

  /// Swap the active solver derivative fields with the subvectors in ``du``.
  void swap_deriv(const N_Vector du) const {
    std::size_t i = 0;
    for (auto& var_str : solver.f2d) {
      BoutNVector::swap(du, *var_str.F_var, i);
      ++i;
    }
    for (auto& var_str : solver.f3d) {
      BoutNVector::swap(du, *var_str.F_var, i);
      ++i;
    }
  }

  template <typename FieldType>
  static void fill_field(FieldType& field, BoutReal value, bool evolve_bndry) {
    BOUT_FOR(i, field.getRegion(evolve_bndry ? RGN_ALL : RGN_NOBNDRY)) {
      field[i] = value;
    }
  }

  template <typename FieldType>
  static void subtract_field(FieldType& residual, const FieldType& du,
                             const FieldType& id, bool evolve_bndry) {
    BOUT_FOR(i, residual.getRegion(evolve_bndry ? RGN_ALL : RGN_NOBNDRY)) {
      if (id[i] > 0.5) {
        residual[i] -= du[i];
      }
    }
  }
#endif
};

#endif
