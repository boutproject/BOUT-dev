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
#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <stddef.h>
#include <utility>
#include <vector>
class GlobalField3DAccess;

namespace fci_comm {
struct ProcLocal {
  int proc;
  int ind;
};
struct globalToLocal1D {
  const int mg;
  const int npe;
  const int localwith;
  const int local;
  const int global;
  const int globalwith;
  const bool periodic;
  globalToLocal1D(int mg, int npe, int localwith, bool periodic)
      : mg(mg), npe(npe), localwith(localwith), local(localwith - 2 * mg),
        global(local * npe), globalwith(global + 2 * mg), periodic(periodic) {};
  ProcLocal convert(int id) const {
    if (periodic) {
      while (id < mg) {
        id += global;
      }
      while (id >= global + mg) {
        id -= global;
      }
    }
    int idwo = id - mg;
    int proc = idwo / local;
    if (not periodic) {
      if (proc >= npe) {
        proc = npe - 1;
      }
    }
    int loc = id - local * proc;
#if CHECK > 1
    if ((loc < 0 or loc > localwith or proc < 0 or proc >= npe)
        or (periodic and (loc < mg or loc >= local + mg))) {
      printf("globalToLocal1D failure: %d %d, %d %d, %d %d %s\n", id, idwo, globalwith,
             npe, proc, loc, periodic ? "periodic" : "non-periodic");
      ASSERT0(false);
    }
#endif
    return {proc, loc};
  }
};
template <class ind>
struct XYZ2Ind {
  const int nx;
  const int ny;
  const int nz;
  ind convert(const int x, const int y, const int z) const {
    return {z + (y + x * ny) * nz, ny, nz};
  }
  ind operator()(const int x, const int y, const int z) const { return convert(x, y, z); }
  XYZ2Ind(const int nx, const int ny, const int nz) : nx(nx), ny(ny), nz(nz) {}
};
} // namespace fci_comm

class GlobalField3DAccessInstance {
public:
  const BoutReal& operator[](IndG3D ind) const;

  GlobalField3DAccessInstance(const GlobalField3DAccess* gfa,
                              const std::vector<BoutReal>&& data)
      : gfa(*gfa), data(std::move(data)) {};

private:
  const GlobalField3DAccess& gfa;
  const std::vector<BoutReal> data;
};

class GlobalField3DAccess {
public:
  friend class GlobalField3DAccessInstance;
  GlobalField3DAccess(Mesh* mesh)
      : mesh(mesh), g2lx(mesh->xstart, mesh->getNXPE(), mesh->LocalNx, false),
        g2ly(mesh->ystart, mesh->getNYPE(), mesh->LocalNy, true),
        g2lz(mesh->zstart, 1, mesh->LocalNz, true),
        xyzl(g2lx.localwith, g2ly.localwith, g2lz.localwith),
        xyzg(g2lx.globalwith, g2ly.globalwith, g2lz.globalwith), comm(BoutComm::get()) {
#ifdef _OPENMP
    o_ids.resize(omp_get_max_threads());
#endif
  };
  void get(IndG3D ind) {
    ASSERT2(is_setup == false);
#ifdef _OPENMP
    ASSERT2(o_ids.size() > static_cast<size_t>(omp_get_thread_num()));
    o_ids[omp_get_thread_num()].emplace(ind.ind);
#else
    ids.emplace(ind.ind);
#endif
  }

  void operator[](IndG3D ind) { return get(ind); }
  void setup() {
    ASSERT2(is_setup == false);
#ifdef _OPENMP
    for (auto& o_id : o_ids) {
      ids.merge(o_id);
    }
    o_ids.clear();
#endif
    toGet.resize(g2lx.npe * g2ly.npe * g2lz.npe);
    for (const auto id : ids) {
      IndG3D gind{id, g2ly.globalwith, g2lz.globalwith};
      const auto pix = g2lx.convert(gind.x());
      const auto piy = g2ly.convert(gind.y());
      const auto piz = g2lz.convert(gind.z());
      ASSERT3(piz.proc == 0);
      toGet[piy.proc * g2lx.npe + pix.proc].push_back(
          xyzl.convert(pix.ind, piy.ind, piz.ind).ind);
    }
    for (auto& v : toGet) {
      std::sort(v.begin(), v.end());
    }
    commCommLists();
    {
      int offset = 0;
      for (auto get : toGet) {
        getOffsets.push_back(offset);
        offset += get.size();
      }
      getOffsets.push_back(offset);
    }
    for (const auto id : ids) {
      IndG3D gind{id, g2ly.globalwith, g2lz.globalwith};
      const auto pix = g2lx.convert(gind.x());
      const auto piy = g2ly.convert(gind.y());
      const auto piz = g2lz.convert(gind.z());
      ASSERT3(piz.proc == 0);
      const auto proc = piy.proc * g2lx.npe + pix.proc;
      const auto& vec = toGet[proc];
      const auto tofind = xyzl.convert(pix.ind, piy.ind, piz.ind).ind;
      auto it = std::lower_bound(vec.begin(), vec.end(), tofind);
      ASSERT3(it != vec.end());
      ASSERT3(*it == tofind);
      mapping[id] = std::distance(vec.begin(), it) + getOffsets[proc];
    }
    is_setup = true;
  }
  GlobalField3DAccessInstance communicate(const Field3D& f) {
    return {this, communicate_data(f)};
  }
  std::unique_ptr<GlobalField3DAccessInstance> communicate_asPtr(const Field3D& f) {
    return std::make_unique<GlobalField3DAccessInstance>(this, communicate_data(f));
  }

private:
  void commCommLists() {
    toSend.resize(toGet.size());
    std::vector<int> toGetSizes(toGet.size(), -1);
    std::vector<int> toSendSizes(toSend.size(), -1);
#if CHECK > 3
    {
      int thisproc;
      MPI_Comm_rank(comm, &thisproc);
      ASSERT0(thisproc == mesh->getYProcIndex() * g2lx.npe + mesh->getXProcIndex());
    }
#endif
    std::vector<MPI_Request> reqs(toSend.size());
    for (size_t proc = 0; proc < toGet.size(); ++proc) {
      auto ret = MPI_Irecv(static_cast<void*>(&toSendSizes[proc]), 1, MPI_INT, proc, 666,
                           comm, &reqs[proc]);
      ASSERT0(ret == MPI_SUCCESS);
    }
    for (size_t proc = 0; proc < toGet.size(); ++proc) {
      toGetSizes[proc] = toGet[proc].size();
      auto ret =
          MPI_Send(static_cast<void*>(&toGetSizes[proc]), 1, MPI_INT, proc, 666, comm);
      ASSERT0(ret == MPI_SUCCESS);
    }
    std::vector<MPI_Request> reqs2(toSend.size());
    int cnt = 0;
    for ([[maybe_unused]] auto dummy : reqs) {
      int ind{0};
      auto ret = MPI_Waitany(reqs.size(), &reqs[0], &ind, MPI_STATUS_IGNORE);
      ASSERT0(ret == MPI_SUCCESS);
      ASSERT3(ind != MPI_UNDEFINED);
      ASSERT2(static_cast<size_t>(ind) < toSend.size());
      ASSERT3(toSendSizes[ind] >= 0);
      if (toSendSizes[ind] == 0) {
        continue;
      }
      sendBufferSize += toSendSizes[ind];
      toSend[ind].resize(toSendSizes[ind], -1);

      ret = MPI_Irecv(static_cast<void*>(toSend[ind].data()), toSend[ind].size(), MPI_INT,
                      ind, 666 * 666, comm, reqs2.data() + cnt++);
      ASSERT0(ret == MPI_SUCCESS);
    }
    for (size_t proc = 0; proc < toGet.size(); ++proc) {
      if (toGet[proc].size() != 0) {
        const auto ret = MPI_Send(static_cast<void*>(toGet[proc].data()),
                                  toGet[proc].size(), MPI_INT, proc, 666 * 666, comm);
        ASSERT0(ret == MPI_SUCCESS);
      }
    }
    for (int c = 0; c < cnt; c++) {
      int ind{0};
      const auto ret = MPI_Waitany(cnt, reqs2.data(), &ind, MPI_STATUS_IGNORE);
      ASSERT0(ret == MPI_SUCCESS);
      ASSERT3(ind != MPI_UNDEFINED);
    }
  }
  Mesh* mesh;
#ifdef _OPENMP
  std::vector<std::set<int>> o_ids;
#endif
  std::set<int> ids;
  std::map<int, int> mapping;
  bool is_setup{false};
  const fci_comm::globalToLocal1D g2lx;
  const fci_comm::globalToLocal1D g2ly;
  const fci_comm::globalToLocal1D g2lz;

public:
  const fci_comm::XYZ2Ind<Ind3D> xyzl;
  const fci_comm::XYZ2Ind<IndG3D> xyzg;

private:
  std::vector<std::vector<int>> toGet;
  std::vector<std::vector<int>> toSend;
  std::vector<int> getOffsets;
  int sendBufferSize{0};
  MPI_Comm comm;
  std::vector<BoutReal> communicate_data(const Field3D& f) {
    ASSERT2(is_setup);
    ASSERT2(f.getMesh() == mesh);
    std::vector<BoutReal> data(getOffsets.back());
    std::vector<BoutReal> sendBuffer(sendBufferSize);
    std::vector<MPI_Request> reqs(toSend.size());
    int cnt1 = 0;
    for (size_t proc = 0; proc < toGet.size(); ++proc) {
      if (toGet[proc].size() == 0) {
        continue;
      }
      auto ret =
          MPI_Irecv(static_cast<void*>(data.data() + getOffsets[proc]),
                    toGet[proc].size(), MPI_DOUBLE, proc, 666, comm, reqs.data() + cnt1);
      ASSERT0(ret == MPI_SUCCESS);
      cnt1++;
    }
    int cnt = 0;
    for (size_t proc = 0; proc < toGet.size(); ++proc) {
      if (toSend[proc].size() == 0) {
        continue;
      }
      const void* start = static_cast<void*>(sendBuffer.data() + cnt);
      for (auto i : toSend[proc]) {
        sendBuffer[cnt++] = f[Ind3D(i)];
      }
      auto ret = MPI_Send(start, toSend[proc].size(), MPI_DOUBLE, proc, 666, comm);
      ASSERT0(ret == MPI_SUCCESS);
    }
    for (int j = 0; j < cnt1; ++j) {
      int ind{0};
      auto ret = MPI_Waitany(cnt1, reqs.data(), &ind, MPI_STATUS_IGNORE);
      ASSERT0(ret == MPI_SUCCESS);
      ASSERT3(ind != MPI_UNDEFINED);
      ASSERT3(ind >= 0);
      ASSERT3(ind < cnt1);
    }
    return data;
  }
};
