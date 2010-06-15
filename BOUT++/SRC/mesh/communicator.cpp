/**************************************************************************
 * Class for parallel communication of scalar and vector fields
 * 
 * 
 * Changelog: 
 * 
 * 2008-07 Ben Dudson <bd512@york.ac.uk>
 *   * Extended for X domain decomposition
 *
 * 2007-09 Ben Dudson <bd512@york.ac.uk>
 *   * Initial version
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
 **************************************************************************/

#include "communicator.h"

#include <stdlib.h>
#include <string.h>

#include "mpi.h"

// This was defined in nvector.h
#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

// Print detailed timing
//#define PRINT_TIME

/**************************************************************************
 * TAG NAMES - label MPI messages 
 **************************************************************************/

const int IN_SENT_UP    = 0;
const int OUT_SENT_UP   = 1;
const int IN_SENT_DOWN  = 2;
const int OUT_SENT_DOWN = 3;
// X communication signals
const int IN_SENT_OUT = 4;
const int OUT_SENT_IN  = 5;

real Communicator::wtime = 0.0;

/**************************************************************************
 * Global communicator options
 *
 * These are set in BOUT.inp, and are common to all Communicator objects
 **************************************************************************/

bool Communicator::options_set = false;
bool Communicator::async_send = false;
bool Communicator::pre_post = false;

/**************************************************************************
 * Constructor / Destructor
 **************************************************************************/

Communicator::Communicator()
{
  var_list.clear();
  
  send_cur = send_pending = false;

  first_time = true;

  for(int i=0;i<6;i++) {
    request[i] = MPI_REQUEST_NULL;
    sendreq[i] = MPI_REQUEST_NULL;
  }

  ybufflen = xbufflen = 0;
}

Communicator::Communicator(Communicator &copy)
{
  /// Copy the arrays of pointers to objects
  
  var_list = copy.var_list;

  ybufflen = copy.ybufflen;
  xbufflen = copy.xbufflen;

  if(ybufflen > 0) {
    umsg_sendbuff = new real[ybufflen];
    dmsg_sendbuff = new real[ybufflen];
    umsg_recvbuff = new real[ybufflen];
    dmsg_recvbuff = new real[ybufflen];
  }
  
  if(xbufflen > 0) {
    imsg_sendbuff = new real[xbufflen];
    omsg_sendbuff = new real[xbufflen];
    imsg_recvbuff = new real[xbufflen];
    omsg_recvbuff = new real[xbufflen];
  }

  send_cur = send_pending = false;

  first_time = true;

  for(int i=0;i<6;i++) {
    request[i] = MPI_REQUEST_NULL;
    sendreq[i] = MPI_REQUEST_NULL;
  }
}

Communicator::~Communicator()
{
  /// Free pointers to fields

  var_list.clear();

  /// Free message buffers
  
  if(ybufflen > 0) {
    delete[] umsg_sendbuff;
    delete[] umsg_recvbuff;
    delete[] dmsg_sendbuff;
    delete[] dmsg_recvbuff;
  }
  
  if(xbufflen > 0) {
    delete[] imsg_sendbuff;
    delete[] imsg_recvbuff;
    delete[] omsg_sendbuff;
    delete[] omsg_recvbuff;
  }
}

/************************************************************************//**
 * Add fields
 **************************************************************************/

void Communicator::add(FieldData &f)
{
  if(send_cur) { // Currently mid-send
    bout_error("ERROR: Variable added to Communicator between send and receive\n");
  }

  if(!options_set) {
    output.write("Initialising comms\n");
    options.setSection("comms");
    options.get("async", async_send, false);
    options.get("pre_post", pre_post, false); 
    options_set = true;
  }
  
  // Add to the list
  var_list.push_back(&f);
}

void Communicator::clear()
{
  /// Free pointers to fields

  var_list.clear();
  
  /// Reset message flags
  
  send_cur = send_pending = false;

  for(int i=0;i<6;i++) {
    request[i] = MPI_REQUEST_NULL;
    sendreq[i] = MPI_REQUEST_NULL;
  }
}

/**************************************************************************
 * Post receives. Usually this is done just before sending, but
 * setting pre_post=true moves this to just after receiving the previous
 * messages (may be faster).
 **************************************************************************/

void Communicator::post_receive()
{
  real *inbuff;
  int len;
  
  /// Post receive data from above (y+1)

  len = 0;
  if(UDATA_INDEST != -1) {
    len = msg_len(0, UDATA_XSPLIT, 0, MYG);
    MPI_Irecv(umsg_recvbuff,
	      len,
	      PVEC_REAL_MPI_TYPE,
	      UDATA_INDEST,
	      IN_SENT_DOWN,
	      MPI_COMM_WORLD,
	      &request[0]);
  }
  if(UDATA_OUTDEST != -1) {
    inbuff = &umsg_recvbuff[len]; // pointer to second half of the buffer
    MPI_Irecv(inbuff,
	      msg_len(UDATA_XSPLIT, ngx, 0, MYG),
	      PVEC_REAL_MPI_TYPE,
	      UDATA_OUTDEST,
	      OUT_SENT_DOWN,
	      MPI_COMM_WORLD,
	      &request[1]);
  }
  
  /// Post receive data from below (y-1)

  len = 0;

  if(DDATA_INDEST != -1) { // If sending & recieving data from a processor
    len = msg_len(0, DDATA_XSPLIT, 0, MYG);
    MPI_Irecv(dmsg_recvbuff, 
	      len,
	      PVEC_REAL_MPI_TYPE,
	      DDATA_INDEST,
	      IN_SENT_UP,
	      MPI_COMM_WORLD,
	      &request[2]);
  }
  if(DDATA_OUTDEST != -1) {
    inbuff = &dmsg_recvbuff[len];
    MPI_Irecv(inbuff,
	      msg_len(DDATA_XSPLIT, ngx, 0, MYG),
	      PVEC_REAL_MPI_TYPE,
	      DDATA_OUTDEST,
	      OUT_SENT_UP,
	      MPI_COMM_WORLD,
	      &request[3]);
  }

  /// Post receive data from left (x-1)
  
  if(IDATA_DEST != -1) {
    MPI_Irecv(imsg_recvbuff,
	      msg_len(0, MXG, 0, MYSUB),
	      PVEC_REAL_MPI_TYPE,
	      IDATA_DEST,
	      OUT_SENT_IN,
	      MPI_COMM_WORLD,
	      &request[4]);
  }

  // Post receive data from right (x+1)

  if(ODATA_DEST != -1) {
    MPI_Irecv(omsg_recvbuff,
	      msg_len(0, MXG, 0, MYSUB),
	      PVEC_REAL_MPI_TYPE,
	      ODATA_DEST,
	      IN_SENT_OUT,
	      MPI_COMM_WORLD,
	      &request[5]);
  }
}

/**************************************************************************
 * Main communication routines
 **************************************************************************/

void Communicator::send()
{
  real *outbuff;
  int len;
  real t;
  
#ifdef PRINT_TIME
  real t2, t3;
  output.write("SEND: ");
#endif

  if(send_cur) {
    // Sending twice in a row!
    bout_error("Error: Communicator send called twice without receive\n");
  }

  /// Record starting wall-time
  t = MPI_Wtime();

  if(async_send && send_pending) {
    /// Asyncronous sending: Need to check if previous sends have completed

    MPI_Status status;

    if(UDATA_INDEST != -1)
      MPI_Wait(sendreq, &status);
    if(UDATA_OUTDEST != -1)
      MPI_Wait(sendreq+1, &status);
    if(DDATA_INDEST != -1)
      MPI_Wait(sendreq+2, &status);
    if(DDATA_OUTDEST != -1)
      MPI_Wait(sendreq+3, &status);
    if(IDATA_DEST != -1)
      MPI_Wait(sendreq+4, &status);
    if(ODATA_DEST != -1)
      MPI_Wait(sendreq+5, &status);
  }

#ifdef PRINT_TIME
  output.write("Wait: %e ", (t2 = MPI_Wtime()) - t);
#endif  

  ////////////////// ALLOCATE BUFFERS //////////////////

  /// work out how many reals need to be sent
  len = msg_len(0, ngx, 0, MYG);
  /// Make sure buffers are the correct size
  if(ybufflen < len) {
    if(ybufflen != 0) {
      delete[] umsg_sendbuff;
      delete[] umsg_recvbuff;
      delete[] dmsg_sendbuff;
      delete[] dmsg_recvbuff;
    }
    
    umsg_sendbuff = new real[len];
    umsg_recvbuff = new real[len];
    dmsg_sendbuff = new real[len];
    dmsg_recvbuff = new real[len];
    ybufflen = len;
  }

  len = msg_len(0, MXG, 0, MYSUB);

  if(xbufflen != len) {
    if(xbufflen != 0) {
      delete[] imsg_sendbuff;
      delete[] imsg_recvbuff;
      delete[] omsg_sendbuff;
      delete[] omsg_recvbuff;
    }
    
    imsg_sendbuff = new real[len];
    imsg_recvbuff = new real[len];
    omsg_sendbuff = new real[len];
    omsg_recvbuff = new real[len];
    xbufflen = len;
  }

  ////////////// POST RECEIVES ////////////////

  if((!pre_post) || first_time)
    post_receive();
  first_time = false; // after the first time, receives are posted the previous time around

#ifdef PRINT_TIME
  output.write("Post: %e ", (t3 = MPI_Wtime()) - t2);
#endif

  ////////////////// SEND DATA ////////////////////

  /// Send data going up (y+1)
  
  len = 0;
  if(UDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(0, UDATA_XSPLIT, MYSUB, MYSUB+MYG, umsg_sendbuff);
    // Send the data to processor UDATA_INDEST

    if(async_send) {
      MPI_Isend(umsg_sendbuff,   // Buffer to send
		len,             // Length of buffer in reals
		PVEC_REAL_MPI_TYPE,  // Real variable type
		UDATA_INDEST,        // Destination processor
		IN_SENT_UP,          // Label (tag) for the message
		MPI_COMM_WORLD,
		&sendreq[0]);
    }else
      MPI_Send(umsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       UDATA_INDEST,
	       IN_SENT_UP,
	       MPI_COMM_WORLD);
  }
  if(UDATA_OUTDEST != -1) { // if destination for outer x data
    outbuff = &umsg_sendbuff[len]; // A pointer to the start of the second part
                                   // of the buffer 
    len = pack_data(UDATA_XSPLIT, ngx, MYSUB, MYSUB+MYG, outbuff);
    // Send the data to processor UDATA_OUTDEST
    if(async_send) {
      MPI_Isend(outbuff, 
		len, 
		PVEC_REAL_MPI_TYPE,
		UDATA_OUTDEST,
		OUT_SENT_UP,
		MPI_COMM_WORLD,
		&sendreq[1]);
    }else
      MPI_Send(outbuff, 
	       len, 
	       PVEC_REAL_MPI_TYPE,
	       UDATA_OUTDEST,
	       OUT_SENT_UP,
	       MPI_COMM_WORLD);
  }
    
  /// Send data going down (y-1)

  len = 0;
  if(DDATA_INDEST != -1) { // If there is a destination for inner x data
    len = pack_data(0, DDATA_XSPLIT, MYG, 2*MYG, dmsg_sendbuff);    
    // Send the data to processor DDATA_INDEST
    if(async_send) {
      MPI_Isend(dmsg_sendbuff, 
		len,
		PVEC_REAL_MPI_TYPE,
		DDATA_INDEST,
		IN_SENT_DOWN,
		MPI_COMM_WORLD,
		&sendreq[2]);
    }else
      MPI_Send(dmsg_sendbuff, 
	       len,
	       PVEC_REAL_MPI_TYPE,
	       DDATA_INDEST,
	       IN_SENT_DOWN,
	       MPI_COMM_WORLD);
  }
  if(DDATA_OUTDEST != -1) { // if destination for outer x data
    outbuff = &dmsg_sendbuff[len]; // A pointer to the start of the second part
			           // of the buffer
    len = pack_data(DDATA_XSPLIT, ngx, MYG, 2*MYG, outbuff);
    // Send the data to processor DDATA_OUTDEST

    if(async_send) {
      MPI_Isend(outbuff,
		len,
		PVEC_REAL_MPI_TYPE,
		DDATA_OUTDEST,
		OUT_SENT_DOWN,
		MPI_COMM_WORLD,
		&sendreq[3]);
    }else
      MPI_Send(outbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       DDATA_OUTDEST,
	       OUT_SENT_DOWN,
	       MPI_COMM_WORLD);
  }

  /// Send to the left (x-1)
  
  if(IDATA_DEST != -1) {
    len = pack_data(MXG, 2*MXG, MYG, MYG+MYSUB, imsg_sendbuff);
    if(async_send) {
      MPI_Isend(imsg_sendbuff,
		len,
		PVEC_REAL_MPI_TYPE,
		IDATA_DEST,
		IN_SENT_OUT,
		MPI_COMM_WORLD,
		&sendreq[4]);
    }else
      MPI_Send(imsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       IDATA_DEST,
	       IN_SENT_OUT,
	       MPI_COMM_WORLD);
  }

  /// Send to the right (x+1)

  if(ODATA_DEST != -1) {
    len = pack_data(MXSUB, MXSUB+MXG, MYG, MYG+MYSUB, omsg_sendbuff);
    if(async_send) {
      MPI_Isend(omsg_sendbuff,
		len,
		PVEC_REAL_MPI_TYPE,
		ODATA_DEST,
		OUT_SENT_IN,
		MPI_COMM_WORLD,
		&sendreq[5]);
    }else
      MPI_Send(omsg_sendbuff,
	       len,
	       PVEC_REAL_MPI_TYPE,
	       ODATA_DEST,
	       OUT_SENT_IN,
	       MPI_COMM_WORLD);
  }
  
  send_cur = true; // Mark a send in process
  send_pending = true;

#ifdef PRINT_TIME
  output.write("Send: %e\n", (t2 = MPI_Wtime()) - t3);
#endif

  /// Add the elapsed wall-time to wtime
  wtime += MPI_Wtime() - t;
}

void Communicator::receive()
{
  MPI_Status status;
  int len;
  real t;

  if(!send_cur) {
    output.write("Error: Communicator receive called without call to send first\n");
    exit(1);
  }

  t = MPI_Wtime();

  ///////////// WAIT FOR DATA //////////////

#ifdef PRINT_TIME
  output.write("RECEIVING: ");
#endif  
  
  int ind;
  do {
    MPI_Waitany(6, request, &ind, &status);
    switch(ind) {
    case 0: { // Up, inner
      unpack_data(0, UDATA_XSPLIT, MYSUB+MYG, MYSUB+2*MYG, umsg_recvbuff);
#ifdef PRINT_TIME
      output.write(" UI");
#endif
      break;
    }
    case 1: { // Up, outer
      len = msg_len(0, UDATA_XSPLIT, 0, MYG);
      unpack_data(UDATA_XSPLIT, ngx, MYSUB+MYG, MYSUB+2*MYG, &umsg_recvbuff[len]);
#ifdef PRINT_TIME
      output.write(" UO");
#endif
      break;
    }
    case 2: { // Down, inner
      unpack_data(0, DDATA_XSPLIT, 0, MYG, dmsg_recvbuff);
#ifdef PRINT_TIME
      output.write(" DI");
#endif
      break;
    }
    case 3: { // Down, outer
      len = msg_len(0, DDATA_XSPLIT, 0, MYG);
      unpack_data(DDATA_XSPLIT, ngx, 0, MYG, &dmsg_recvbuff[len]);
#ifdef PRINT_TIME
      output.write(" DO");
#endif
      break;
    }
    case 4: { // inner
      unpack_data(0, MXG, MYG, MYG+MYSUB, imsg_recvbuff);
#ifdef PRINT_TIME
      output.write(" I");
#endif
      break;
    }
    case 5: { // outer
      unpack_data(MXSUB+MXG, MXSUB+2*MXG, MYG, MYG+MYSUB, omsg_recvbuff);
#ifdef PRINT_TIME
      output.write(" O");
#endif
      break;
    }
    }
    if(ind != MPI_UNDEFINED)
      request[ind] = MPI_REQUEST_NULL;
    
  }while(ind != MPI_UNDEFINED);
  
  send_cur = false; // Finished this send-receive pair

#ifdef PRINT_TIME
  output.write("RECV: %e\n", MPI_Wtime() - t);
#endif

  // If using pre_post, post receives for next time
  if(pre_post)
    post_receive();

  // TWIST-SHIFT CONDITION
  if(TwistShift && (TwistOrder == 0)) {
    int jx, jy;
    
    // Perform Twist-shift using shifting method (rather than in SetStencil)
    for(std::vector<FieldData*>::iterator it = var_list.begin(); it != var_list.end(); it++)
      if((*it)->is3D()) {
	
	// Lower boundary

	if(TS_down_in && (DDATA_INDEST  != -1)) {
	  for(jx=0;jx<DDATA_XSPLIT;jx++)
	    for(jy=0;jy != MYG; jy++)
	      (*it)->ShiftZ(jx, jy, ShiftAngle[jx]);
      
	}
	if(TS_down_out && (DDATA_OUTDEST  != -1)) {
	  for(jx=DDATA_XSPLIT;jx<ngx; jx++)
	    for(jy=0;jy != MYG; jy++)
	      (*it)->ShiftZ(jx, jy, ShiftAngle[jx]);
	  
	}
	
	// Upper boundary
	
	if(TS_up_in && (UDATA_INDEST  != -1)) {
	  for(jx=0;jx<UDATA_XSPLIT; jx++)
	    for(jy=ngy-MYG;jy != ngy; jy++)
	      (*it)->ShiftZ(jx, jy, -ShiftAngle[jx]);
	  
	}
	if(TS_up_out && (UDATA_OUTDEST  != -1)) {
	  for(jx=UDATA_XSPLIT;jx<ngx; jx++)
	    for(jy=ngy-MYG;jy != ngy; jy++)
	      (*it)->ShiftZ(jx, jy, -ShiftAngle[jx]);
	  
	}
      }
  }

  wtime += MPI_Wtime() - t;

#ifdef CHECK
  // Keeping track of whether communications have been done
  for(std::vector<FieldData*>::iterator it = var_list.begin(); it != var_list.end(); it++)
    (*it)->doneComms();
#endif
}

void Communicator::run()
{
  send();
  receive();
}


/************************************************************************//**
 * Pack and unpack data from buffers
 **************************************************************************/

int Communicator::pack_data(int xge, int xlt, int yge, int ylt, real *buffer)
{
  int jx, jy, jz;
  int len = 0;
  std::vector<FieldData*>::iterator it;
  
  for(jx=xge; jx != xlt; jx++) {
    
    /// Loop over variables
    for(it = var_list.begin(); it != var_list.end(); it++) {
      if((*it)->is3D()) {
	// 3D variable
	
	for(jy=yge;jy != ylt;jy++)
	  for(jz=0;jz != ncz;jz++)
	    len += (*it)->getData(jx,jy,jz,buffer+len);
	
      }else {
	// 2D variable
	for(jy=yge;jy < ylt;jy++)
	  len += (*it)->getData(jx,jy,0,buffer+len);
      }
    }
    
  }
  
  return(len);
}

int Communicator::unpack_data(int xge, int xlt, int yge, int ylt, real *buffer)
{
  int jx, jy, jz;
  int len = 0;
  std::vector<FieldData*>::iterator it;

  //output.write("Unpacking for %d <= x < %d\n", xge, xlt);

  for(jx=xge; jx != xlt; jx++) {

    /// Loop over variables
    for(it = var_list.begin(); it != var_list.end(); it++) {
      if((*it)->is3D()) {
	// 3D variable
   
	for(jy=yge;jy != ylt;jy++)
	  for(jz=0;jz != ncz;jz++) {
	    //output.write("Setting %d,%d,%d\n", jx, jy, jz);
	    len += (*it)->setData(jx,jy,jz,buffer+len);
	  }
	
      }else {
	// 2D variable
	for(jy=yge;jy < ylt;jy++)
	  len += (*it)->setData(jx,jy,0,buffer+len);
      }
    }
    
  }
  
  return(len);
}

int Communicator::msg_len(int xge, int xlt, int yge, int ylt)
{
  int len = 0;

  /// Loop over variables
  for(std::vector<FieldData*>::iterator it = var_list.begin(); it != var_list.end(); it++) {
    if((*it)->is3D()) {
      len += (xlt - xge) * (ylt - yge) * ncz * (*it)->realSize();
    }else
      len += (xlt - xge) * (ylt - yge) * (*it)->realSize();
  }
  
  return len;
}
