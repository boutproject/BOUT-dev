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

#include <comm_group.hxx>
#include <string.h>
#include <stdlib.h>

#include <output.hxx>
#include <options.hxx>
#include <msg_stack.hxx>

namespace comm_group {

  const int COMM_GROUP_TAG = 31415;
  
  static bool initialised = false;
  static bool nonblock;
  
  void Comm_initialise()
  {
    if(initialised)
      return;
    
    output.write("Initialising group comms\n");
    Options *options = Options::getRoot();
    options = options->getSection("comms");
    options->get("group_nonblock", nonblock, true);

    initialised = true;
  }

  bool Comm_gather_start(void *local, int nlocal, MPI_Datatype type,
			 void *data, 
			 int root, MPI_Comm comm,
			 Comm_handle_t *handle)
  {
    TRACE("Comm_gather_start(%d -> %d)", nlocal, root);
   
    Comm_initialise();
 
    // Put communication info into handle
    MPI_Comm_size(comm, &handle->nprocs);
    MPI_Comm_rank(comm, &handle->myrank);
    handle->root = root;

    if(!nonblock) {
      // Blocking comms using the MPI call
      // May be faster for different architectures, and useful
      // for debugging
      
      MPI_Gather(local, nlocal, type,
		 data, nlocal, type,
		 root, comm);
      
      handle->current = true;
      
      return (root == handle->myrank);
    }
    
    int type_size; // Get size of the datatype
    MPI_Type_size(type, &type_size);
    
    if(root == handle->myrank) {
      // All arriving on this processor. Post receives

      handle->request = (MPI_Request *)malloc(sizeof(MPI_Request) * (handle->nprocs - 1));
      handle->nreq = handle->nprocs-1;
      
      MPI_Request *r = handle->request;
      for(int p=0;p<handle->nprocs; p++) {
	if(p != root) {
	  MPI_Irecv((void*) ( ((char*) data) + p*nlocal*type_size),
		    nlocal,
		    type,
		    p,
		    COMM_GROUP_TAG,
		    comm,
		    r);
	  r++; // Next request
	}
      }
      
      // Copy local data into output
      memcpy((void*) (((char*) data) + handle->myrank*nlocal*type_size),
	     local, 
	     nlocal*type_size);
    }else {
      // Sending to root processor

      handle->request = (MPI_Request *)malloc(sizeof(MPI_Request));
      handle->nreq = 1;

      MPI_Isend(local,
		nlocal,
		type,
		root,
		COMM_GROUP_TAG,
		comm,
		handle->request);
    }
    
    handle->current = true; // Mark as in progress
    
    return (root == handle->myrank);
  }
  
  /// start a scatter operation
  /*!
   * @param[in]  sendbuf   Buffer of data to be send (only used on root process)
   * @param[in]  sendcnt   Number of elements to be sent to each process (NOT TOTAL!)
   * @param[in]  type      Data type
   * @param[out] recvbuf   Where to store the data
   * @param[in]  root      Process rank where data originates
   * @param[in]  comm      Communicator
   * @param[out] handle    Used for completion later
   */ 
  bool Comm_scatter_start( void *sendbuf, int sendcnt, MPI_Datatype type, 
			   void *recvbuf, 
			   int root, MPI_Comm comm,
			   Comm_handle_t *handle)
  {
    TRACE("Comm_scatter_start(%d -> %d)", root, sendcnt);

    Comm_initialise();

    MPI_Comm_size(comm, &handle->nprocs);
    MPI_Comm_rank(comm, &handle->myrank);
    handle->root = root;
    
    if(!nonblock) {
      MPI_Scatter ( sendbuf, sendcnt, type, recvbuf, sendcnt, type, root, comm );
      handle->current = true;
      
      return (root == handle->myrank);
    }
    
    int type_size; // Get size of the datatype
    MPI_Type_size(type, &type_size);
    
    if(root == handle->myrank) {
      // Sending from this processor
      
      handle->request = (MPI_Request*) malloc(sizeof(MPI_Request)*(handle->nprocs-1));
      handle->nreq = handle->nprocs-1;
      
      MPI_Request *r = handle->request;
      for(int p=0;p<handle->nprocs; p++) {
	if(p != root) {
	  
	  MPI_Isend((void*) ( ((char*) sendbuf) + p*sendcnt*type_size),
		    sendcnt,
		    type,
		    p,
		    COMM_GROUP_TAG,
		    comm,
		    r);
	  r++; // Next request
	}
      }
      
      // Copy local data
      memcpy(recvbuf,
	     (void*) (((char*) sendbuf) + root*sendcnt*type_size),
	     sendcnt*type_size);
    }else {
      // Receiving one message from root
      
      handle->request = (MPI_Request*) malloc(sizeof(MPI_Request));
      handle->nreq = 1;

      MPI_Irecv(recvbuf,
		sendcnt,
		type,
		root,
		COMM_GROUP_TAG,
		comm,
		handle->request);
    }

    handle->current = true;
    
    return (root == handle->myrank);
  }
  
  
  bool Comm_wait(Comm_handle_t *handle)
  {
    TRACE("Comm_gather_wait(%d)", handle->root);

    if(handle == NULL)
      return false;
    if(!handle->current)
      return false;

    if(!nonblock) { // Already done communication
      handle->current = false;
      return (handle->root == handle->myrank);
    }

    MPI_Status *s = (MPI_Status*) malloc(sizeof(MPI_Status)*(handle->nreq));
    MPI_Waitall(handle->nreq, 
		handle->request,
		s);

    free(s); 
    free(handle->request);
    handle->current = false;
    
    return (handle->root == handle->myrank);
  }
  
  bool Comm_wait_all(int n, Comm_handle_t *handles)
  {
    bool newdata = false;
    
    for(int i=0; i<n;i++) {
      newdata |= Comm_wait(handles + i);
    }
    return newdata;
  }
}
