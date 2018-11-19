/***********************************************************************
 * Non-blocking collective operations (gather, scatter)
 * 
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact Ben Dudson, bd512@york.ac.uk
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
 ***********************************************************************/

#ifndef __COMM_GROUP_H__
#define __COMM_GROUP_H__

#include "mpi.h"

namespace comm_group {
  /// Communication handle
  typedef struct {
    // Communication info
    int root, myrank, nprocs;
    
    // Handles
    int nreq;     ///< Number of MPI requests
    MPI_Request *request;
    bool current; ///< True if currently in progress
  }Comm_handle_t;
  
  /// Begin a gather operation
  bool Comm_gather_start(void *local, int nlocal, MPI_Datatype type,
			 void *data, 
			 int root, MPI_Comm comm,
			 Comm_handle_t *handle);
  
  bool Comm_scatter_start( void *sendbuf, int sendcnt, MPI_Datatype type, 
			   void *recvbuf, 
			   int root, MPI_Comm comm,
			   Comm_handle_t *handle);
  
  /// Wait for a single operation to finish
  bool Comm_wait(Comm_handle_t *handle);
  
  /// Wait for all the communications to finish
  bool Comm_wait_all(int n, Comm_handle_t *handles);
}

using comm_group::Comm_handle_t;
using comm_group::Comm_gather_start;
using comm_group::Comm_scatter_start;
using comm_group::Comm_wait;
using comm_group::Comm_wait_all;

#endif // __COMM_GROUP_H__
