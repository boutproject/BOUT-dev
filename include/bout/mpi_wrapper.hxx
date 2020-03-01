/*!************************************************************************
 * Wrappers around MPI routines, which can be overridden for testing
 * purposes
 *
 *
 **************************************************************************
 * Copyright 2019 C. MacMackin
 *
 * Contact: Chris MacMackin, chris.macmackin@ukaea.uk
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

class MpiWrapper;

#ifndef __MPIWRAPPER_H__
#define __MPIWRAPPER_H__

#include <mpi.h>

/// Provides wrappers around MPI functions, taking the same names. These can
/// then be overloaded for testing purposes.
class MpiWrapper {
public:
  virtual ~MpiWrapper(){};

  virtual int MPI_Abort(MPI_Comm comm, int errorcode) {
    return ::MPI_Abort(comm, errorcode);
  }

  virtual int MPI_Allreduce(const void* sendbuf, void* recvbuf, int count,
                            MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    return ::MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm);
  }

  virtual int MPI_Barrier(MPI_Comm comm) { return ::MPI_Barrier(comm); }

  virtual int MPI_Comm_create(MPI_Comm comm, MPI_Group group, MPI_Comm* newcomm) {
    return ::MPI_Comm_create(comm, group, newcomm);
  }

  virtual int MPI_Comm_dup(MPI_Comm comm, MPI_Comm* newcomm) {
    return ::MPI_Comm_dup(comm, newcomm);
  }

  virtual int MPI_Comm_free(MPI_Comm* comm) { return ::MPI_Comm_free(comm); }

  virtual int MPI_Comm_group(MPI_Comm comm, MPI_Group* group) {
    return ::MPI_Comm_group(comm, group);
  }

  virtual int MPI_Group_range_incl(MPI_Group group, int n, int ranges[][3],
                                   MPI_Group* newgroup) {
    return ::MPI_Group_range_incl(group, n, ranges, newgroup);
  }

  virtual int MPI_Comm_rank(MPI_Comm comm, int* rank) {
    return ::MPI_Comm_rank(comm, rank);
  }

  virtual int MPI_Comm_size(MPI_Comm comm, int* size) {
    return ::MPI_Comm_size(comm, size);
  }

  virtual int MPI_Comm_split(MPI_Comm comm, int color, int key, MPI_Comm* newcomm) {
    return ::MPI_Comm_split(comm, color, key, newcomm);
  }

  virtual int MPI_Gatherv(const void* sendbuf, int sendcount, MPI_Datatype sendtype,
                          void* recvbuf, const int* recvcounts, const int* displs,
                          MPI_Datatype recvtype, int root, MPI_Comm comm) {
    return ::MPI_Gatherv(sendbuf, sendcount, sendtype, recvbuf, recvcounts, displs,
                         recvtype, root, comm);
  }

  virtual int MPI_Group_union(MPI_Group group1, MPI_Group group2, MPI_Group* newgroup) {
    return ::MPI_Group_union(group1, group2, newgroup);
  }

  virtual int MPI_Group_free(MPI_Group* group) { return ::MPI_Group_free(group); }

  virtual int MPI_Irecv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
                        MPI_Comm comm, MPI_Request* request) {
    return ::MPI_Irecv(buf, count, datatype, source, tag, comm, request);
  }

  virtual int MPI_Isend(const void* buf, int count, MPI_Datatype datatype, int dest,
                        int tag, MPI_Comm comm, MPI_Request* request) {
    return ::MPI_Isend(buf, count, datatype, dest, tag, comm, request);
  }

  virtual int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag,
                       MPI_Comm comm, MPI_Status* status) {
    return ::MPI_Recv(buf, count, datatype, source, tag, comm, status);
  }

  virtual int MPI_Scan(const void* sendbuf, void* recvbuf, int count,
                       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    return ::MPI_Scan(sendbuf, recvbuf, count, datatype, op, comm);
  }

  virtual int MPI_Scatterv(const void* sendbuf, const int* sendcounts, const int* displs,
                           MPI_Datatype sendtype, void* recvbuf, int recvcount,
                           MPI_Datatype recvtype, int root, MPI_Comm comm) {
    return ::MPI_Scatterv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount,
                          recvtype, root, comm);
  }

  virtual int MPI_Send(const void* buf, int count, MPI_Datatype datatype, int dest,
                       int tag, MPI_Comm comm) {
    return ::MPI_Send(buf, count, datatype, dest, tag, comm);
  }

  virtual int MPI_Type_commit(MPI_Datatype* datatype) {
    return ::MPI_Type_commit(datatype);
  }

  virtual int MPI_Type_free(MPI_Datatype* datatype) { return ::MPI_Type_free(datatype); }

  virtual int MPI_Type_vector(int count, int blocklength, int stride,
                              MPI_Datatype oldtype, MPI_Datatype* newtype) {
    return ::MPI_Type_vector(count, blocklength, stride, oldtype, newtype);
  }

  virtual int MPI_Wait(MPI_Request* request, MPI_Status* status) {
    return ::MPI_Wait(request, status);
  }

  virtual int MPI_Waitall(int count, MPI_Request array_of_requests[],
                          MPI_Status array_of_statuses[]) {
    return ::MPI_Waitall(count, array_of_requests, array_of_statuses);
  }

  virtual int MPI_Waitany(int count, MPI_Request array_of_requests[], int* indx,
                          MPI_Status* status) {
    return ::MPI_Waitany(count, array_of_requests, indx, status);
  }

  virtual double MPI_Wtime() { return ::MPI_Wtime(); }
};

#endif // __MPIWRAPPER_H__
