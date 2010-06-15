/**************************************************************************
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
 **************************************************************************/

class Communicator;

#ifndef __COMMUNICATOR_H__
#define __COMMUNICATOR_H__

#include "mpi.h"

#include "globals.h"
#include "field_data.h"

#include <vector>

/// Class for parallel communication
/*!
 * Fields and vectors can be added, and then guard cells transmitted between processors
 * \author B.Dudson, University of York
 * \date Sept 2007
 * \warning Once an object has been added, it should not be deleted until communicator
 * is finished (stores a pointer to the object internally).
 *
 * \note July 2008: Modified to communicate in X and Y. Generalised to use the FieldData
 * interface. Changed to use MPI_Isend instead of MPI_Send, and MPI_Waitany instead of MPI_Wait
 * for (hopefully) faster communications.
 */
class Communicator {
 public:
  /// Constructor
  Communicator();
  /// Copy constructor
  Communicator(Communicator &copy);
  /// Destructor
  ~Communicator();
  
  /// Add a data field to be communicated
  void add(FieldData &f);
  
  /// Remove all fields
  void clear();

  /// Tells MPI to expect a message, and sends data to processors
  void send();
  /// Waits for data from processors. This should always be called AFTER send().
  void receive();

  /// Perform communications. Same as send() then receive();
  void run();

  /// Elapsed wall-time. Used to keep track of time spent communicating
  static real wtime;
 private:
  
  void post_receive();
  
  std::vector<FieldData*> var_list; ///< Array of fields to communicate

  int xbufflen, ybufflen;  ///< Length of the buffers used to send/receive (in reals)
  real *umsg_sendbuff, *dmsg_sendbuff, *imsg_sendbuff, *omsg_sendbuff;
  real *umsg_recvbuff, *dmsg_recvbuff, *imsg_recvbuff, *omsg_recvbuff;

  /// Take data from objects and put into a buffer
  int pack_data(int xge, int xlt, int yge, int ylt, real *buffer);
  /// Copy data from a buffer back into the fields
  int unpack_data(int xge, int xlt, int yge, int ylt, real *buffer);
  /// Calculates the size of a message for a given x and y range
  int msg_len(int xge, int xlt, int yge, int ylt);
  
  /// Array of request handles for MPI
  MPI_Request request[6];
  /// Record whether sends have been performed previously
  bool send_cur;

  // Communication options, from BOUT.inp
  static bool options_set; ///< Prevents options being read each time
  static bool async_send; ///< Switch to asyncronous sends (ISend, not Send)
  static bool pre_post; ///< Post receives early. May speed up comms.

  /// When using pre_post, need to make an exception for first time
  bool first_time;

  /// Need to make sure sends complete before buffers are overwritten
  bool send_pending;
  /// Array of send requests (non-blocking send)
  MPI_Request sendreq[6];
};


#endif // __COMMUNICATOR_H__
