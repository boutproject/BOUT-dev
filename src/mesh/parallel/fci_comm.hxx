/**************************************************************************
 * Communication for Flux-coordinate Independent interpolation
 *
 **************************************************************************
 * Copyright 2025 BOUT++ contributors
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

#pragma once

#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/boutcomm.hxx"
#include "bout/field3d.hxx"
#include "bout/mesh.hxx"
#include "bout/region.hxx"
#include <map>
#include <memory>
#include <set>
#include <stddef.h>
#include <utility>
#include <vector>

/// GlobalField3DAccess is a class to set up the communication
/// patterns, to request abitrary data form the global
/// field. GlobalField3DAccessInstance is an instance. after
/// GlobalField3DAccess has communicated the pre-requested data. Only
/// data that has been pre-requested can be requested from
/// GlobalField3DAccessInstance after communication.
///
/// The usage looks a bit like this:
///
/// GlobalField3DAccess gfa;
/// // Request an abitrary number of global indices ``gi``:
/// gfa.request(gi)
/// // Communicate data
/// const auto data = gfa.communicate(f3d);
/// // Now data can be accesssed for all previously requested ``gi``s:
/// data[gi]
class GlobalField3DAccess;

namespace fci_comm {
struct ProcLocal {
  int proc;
  int index;
};

/// Class to convert global to local indices for 1D
/// given the global index, it returns the local index and the processor.
struct GlobalToLocal1D {
  GlobalToLocal1D(int mg, int npe, int localwith, bool periodic)
      : mg(mg), npe(npe), localwith(localwith), local(localwith - (2 * mg)),
        global(local * npe), globalwith(global + (2 * mg)), periodic(periodic){};
  ProcLocal convert(int id) const;
  int getLocalWith() const { return localwith; }
  int getGlobalWith() const { return globalwith; }
  int getNPE() const { return npe; }

private:
  int mg;
  int npe;
  int localwith;
  int local;
  int global;
  int globalwith;
  bool periodic;
};

/// Convert an x-y-z tupple to an Ind
template <class ind>
struct XYZ2Ind {
  XYZ2Ind(const int nx, const int ny, const int nz) : nx(nx), ny(ny), nz(nz) {}
  ind convert(const int x, const int y, const int z) const {
    return {z + ((y + x * ny) * nz), ny, nz};
  }
  ind operator()(const int x, const int y, const int z) const { return convert(x, y, z); }

private:
  int nx;
  int ny;
  int nz;
};
} // namespace fci_comm

class GlobalField3DAccessInstance {
public:
  const BoutReal& operator[](IndG3D ind) const;
  GlobalField3DAccessInstance(const GlobalField3DAccess* gfa,
                              std::vector<BoutReal>&& data)
      : gfa(gfa), data(std::move(data)){};

private:
  const GlobalField3DAccess* gfa;
  std::vector<BoutReal> data;
};

class GlobalField3DAccess {
public:
  friend class GlobalField3DAccessInstance;
  GlobalField3DAccess(Mesh* mesh)
      : mesh(mesh),
        global2local_x(mesh->xstart, mesh->getNXPE(), mesh->LocalNx, mesh->periodicX),
        global2local_y(mesh->ystart, mesh->getNYPE(), mesh->LocalNy, true),
        global2local_z(mesh->zstart, mesh->getNZPE(), mesh->LocalNz, true),
        xyzlocal(global2local_x.getLocalWith(), global2local_y.getLocalWith(),
                 global2local_z.getLocalWith()),
        xyzglobal(global2local_x.getGlobalWith(), global2local_y.getGlobalWith(),
                  global2local_z.getGlobalWith()),
        comm(BoutComm::get()) {
#ifdef _OPENMP
    openmp_ids.resize(omp_get_max_threads());
#endif
#if CHECK >= 2
    // We could also allow false, but then we would need to ensure it
    // is false everywhere.
    for (int x = 0; x < mesh->LocalNx; ++x) {
      ASSERT2(mesh->periodicY(x) == true);
    }
#endif
  };
  void request(IndG3D ind) {
    ASSERT2(is_setup == false);
#ifdef _OPENMP
    ASSERT2(openmp_ids.size() > static_cast<size_t>(omp_get_thread_num()));
    openmp_ids[omp_get_thread_num()].emplace(ind.ind);
#else
    ids.emplace(ind.ind);
#endif
  }

  GlobalField3DAccessInstance communicate(const Field3D& f) {
    return {this, communicate_data(f)};
  }
  std::unique_ptr<GlobalField3DAccessInstance> communicate_asPtr(const Field3D& f) {
    return std::make_unique<GlobalField3DAccessInstance>(this, communicate_data(f));
  }

private:
  void setup();
  void commCommLists();
  Mesh* mesh;
#ifdef _OPENMP
  // openmp thread-local variable
  std::vector<std::set<int>> openmp_ids;
#endif
  std::set<int> ids;
  std::map<int, int> mapping;
  bool is_setup{false};
  fci_comm::GlobalToLocal1D global2local_x;
  fci_comm::GlobalToLocal1D global2local_y;
  fci_comm::GlobalToLocal1D global2local_z;

public:
  fci_comm::XYZ2Ind<Ind3D> xyzlocal;
  fci_comm::XYZ2Ind<IndG3D> xyzglobal;

private:
  std::vector<std::vector<int>> toGet;
  std::vector<std::vector<int>> toSend;
  std::vector<int> getOffsets;
  int sendBufferSize{0};
  MPI_Comm comm;
  std::vector<BoutReal> communicate_data(const Field3D& f);
};

inline const BoutReal& GlobalField3DAccessInstance::operator[](IndG3D ind) const {
  auto it = gfa->mapping.find(ind.ind);
  ASSERT2(it != gfa->mapping.end());
  return data[it->second];
}
