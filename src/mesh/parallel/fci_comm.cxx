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

#include "fci_comm.hxx"
#include "bout/assert.hxx"
#include "bout/bout_types.hxx"
#include "bout/field3d.hxx"
#include "bout/region.hxx"

#include <algorithm>
#include <cstddef>
#include <mpi.h>
#include <vector>

fci_comm::ProcLocal fci_comm::GlobalToLocal1D::convert(int id) const {
  if (periodic) {
    while (id < mg) {
      id += global;
    }
    while (id >= global + mg) {
      id -= global;
    }
  }
  const int idwo = id - mg;
  int proc = idwo / local;
  if (not periodic) {
    if (proc >= npe) {
      proc = npe - 1;
    }
  }
  const int loc = id - (local * proc);
#if CHECK > 1
  ASSERT1(loc >= 0);
  ASSERT1(loc <= localwith);
  ASSERT1(proc >= 0);
  ASSERT1(proc < npe);
  if (periodic and (loc < mg or loc >= local + mg)) {
    throw BoutException(
        "GlobalToLocal1D failure - expected {} < {} < {} because we are periodic\n", mg,
        loc, local + mg);
  }
#endif
  return {proc, loc};
}

void GlobalField3DAccess::setup() {
  // We need to send a list of data to every processor of which data
  // we want We also get a list of all the data that every other
  // processor wants Some of theses lists may be empty. We still need
  // to send them, later we can skip that.  We also compute where the
  // requested data will be stored later on. This is currently
  // implemented as a map, to memory efficient, as the data is
  // sparse. We could also store the mapping as a dense array, if it
  // turns out this lookup is not fast enough, but that may limit
  // scaling at some point.
  ASSERT2(is_setup == false);
#ifdef _OPENMP
  for (auto& o_id : o_ids) {
    ids.merge(o_id);
  }
  o_ids.clear();
#endif
  toGet.resize(static_cast<size_t>(global2local_x.getNPE() * global2local_y.getNPE()
                                   * global2local_z.getNPE()));
  for (const auto id : ids) {
    const IndG3D gind{id, global2local_y.getGlobalWith(), global2local_z.getGlobalWith()};
    const auto pix = global2local_x.convert(gind.x());
    const auto piy = global2local_y.convert(gind.y());
    const auto piz = global2local_z.convert(gind.z());
    ASSERT3(piz.proc == 0);
    toGet[mesh->getProcIndex(pix.proc, piy.index, piz.index)].push_back(
        xyzlocal.convert(pix.ind, piy.ind, piz.ind).ind);
  }
  for (auto& v : toGet) {
    std::sort(v.begin(), v.end());
  }
  commCommLists();
  {
    int offset = 0;
    for (const auto& get : toGet) {
      getOffsets.push_back(offset);
      offset += get.size();
    }
    getOffsets.push_back(offset);
  }
  for (const auto id : ids) {
    const IndG3D gind{id, global2local_y.getGlobalWith(), global2local_z.getGlobalWith()};
    const auto pix = global2local_x.convert(gind.x());
    const auto piy = global2local_y.convert(gind.y());
    const auto piz = global2local_z.convert(gind.z());
    ASSERT3(piz.proc == 0);
    const auto proc = mesh->getProcIndex(pix.proc, piy.index, piz.index);
    const auto& vec = toGet[proc];
    const auto tofind = xyzlocal.convert(pix.ind, piy.ind, piz.ind).ind;
    auto it = std::lower_bound(vec.begin(), vec.end(), tofind);
    ASSERT3(it != vec.end());
    ASSERT3(*it == tofind);
    mapping[id] = std::distance(vec.begin(), it) + getOffsets[proc];
  }
  is_setup = true;
}

void GlobalField3DAccess::commCommLists() {
  toSend.resize(toGet.size());
  std::vector<int> toGetSizes(toGet.size(), -1);
  std::vector<int> toSendSizes(toSend.size(), -1);
#if CHECK > 3
  {
    int thisproc;
    MPI_Comm_rank(comm, &thisproc);
    ASSERT0(thisproc == mesh.getProcIndex(pix.proc, piy.index, piz.index));
  }
#endif
  std::vector<MPI_Request> reqs(toSend.size());
  for (size_t proc = 0; proc < toGet.size(); ++proc) {
    auto ret = MPI_Irecv(&toSendSizes[proc], 1, MPI_INT, proc, 666, comm, &reqs[proc]);
    ASSERT0(ret == MPI_SUCCESS);
  }
  for (size_t proc = 0; proc < toGet.size(); ++proc) {
    toGetSizes[proc] = toGet[proc].size();
    auto ret = MPI_Send(&toGetSizes[proc], 1, MPI_INT, proc, 666, comm);
    ASSERT0(ret == MPI_SUCCESS);
  }
  std::vector<MPI_Request> reqs2(toSend.size());
  int cnt = 0;
  for ([[maybe_unused]] auto dummy : reqs) {
    int ind{0};
    auto ret = MPI_Waitany(reqs.size(), reqs.data(), &ind, MPI_STATUS_IGNORE);
    ASSERT0(ret == MPI_SUCCESS);
    ASSERT3(ind != MPI_UNDEFINED);
    ASSERT2(static_cast<size_t>(ind) < toSend.size());
    ASSERT3(toSendSizes[ind] >= 0);
    if (toSendSizes[ind] == 0) {
      continue;
    }
    sendBufferSize += toSendSizes[ind];
    toSend[ind].resize(toSendSizes[ind], -1);

    ret = MPI_Irecv(toSend[ind].data(), toSend[ind].size(), MPI_INT, ind, 666 * 666, comm,
                    reqs2.data() + cnt++);
    ASSERT0(ret == MPI_SUCCESS);
  }
  for (size_t proc = 0; proc < toGet.size(); ++proc) {
    if (!toGet[proc].empty()) {
      const auto ret = MPI_Send(toGet[proc].data(), toGet[proc].size(), MPI_INT, proc,
                                666 * 666, comm);
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

std::vector<BoutReal> GlobalField3DAccess::communicate_data(const Field3D& f) {
  // Ensure setup is called, to setup communication pattern
  if (not is_setup) {
    setup();
  }
  // We now send the previosly requested data to every processor, that wanted some.
  // We also get the data we requested.
  ASSERT2(f.getMesh() == mesh);
  std::vector<BoutReal> data(getOffsets.back());
  std::vector<BoutReal> sendBuffer(sendBufferSize);
  std::vector<MPI_Request> reqs(toSend.size());
  int cnt1 = 0;
  for (size_t proc = 0; proc < toGet.size(); ++proc) {
    if (toGet[proc].empty()) {
      continue;
    }
    auto ret = MPI_Irecv(data.data() + getOffsets[proc], toGet[proc].size(), MPI_DOUBLE,
                         proc, 666, comm, reqs.data() + cnt1);
    ASSERT0(ret == MPI_SUCCESS);
    cnt1++;
  }
  int cnt = 0;
  for (size_t proc = 0; proc < toGet.size(); ++proc) {
    if (toSend[proc].empty()) {
      continue;
    }
    const void* start = sendBuffer.data() + cnt;
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
