/*!************************************************************************
* MPI Communicator for BOUT++ representation
*
* 
*
**************************************************************************
* Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
*
* Contact: Ben Dudson, bd512@york.ac.uk
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

class BoutComm;

#ifndef __BOUTCOMM_H__
#define __BOUTCOMM_H__

#include <mpi.h>

#include "bout/sys/variant.hxx"
#include <vector>

/// Class representing a communication method
class Communicator {
public:
  Communicator() = default;
  virtual ~Communicator() = default;
  
  /// Handles for communication requests
  /// Note: monostate used as null reference (default constructor)
  using Request = bout::utils::variant<bout::utils::monostate, int, MPI_Request>;

  const Request REQUEST_NULL {}; ///< Default initialised to monostate (null)
  const int UNDEFINED = -1; ///< Index returned by waitAny
  
  /// Return the rank of this processor
  /// Note: Can't use "rank" for now because BoutComm has static function of that name
  virtual int commRank() = 0;

  /// Number of processors in this communicator
  virtual int commSize() = 0;

  /// Non-blocking receive
  /// If the type is unknown, it is communicated as bytes
  template<typename T>
  Request irecv(T* storage, int length, int destination, int identifier) {
    return irecvBytes(static_cast<void*>(storage), length * sizeof(T), destination, identifier);
  }

  /// Non-blocking receive
  virtual Request irecvBytes(void *data, int length, int destination, int identifier) = 0;
  
  /// Non-blocking send
  virtual void sendBytes(void *data, int length, int destination, int identifier) = 0;

  /// Wait for any of a vector of requests. Returns the index or -1 (UNDEFINED) if no valid requests.
  /// An integer is returned rather than an iterator, because the user may associate the request index
  /// with an array index.
  virtual int waitAny(const std::vector<Request> &requests) = 0; 
};

class SingleProcessorCommunicator : public Communicator {
public:
  int commRank() override {return 1;}
  int commSize() override {return 1;}
  Request irecvBytes(void *data, int length, int destination, int identifier) override;
  void sendBytes(void *data, int length, int destination, int identifier) override;
  int waitAny(const std::vector<Request> &requests) override;
private:
  struct RequestInfo {
    void *data;
    int length;
    bool arrived;
  };
  int next_request_id {0}; // A counter used to make unique requests
  std::map<int, std::list<RequestInfo>> posted_requests; // Currently waiting requests
};

/// Class to represent the 'global' communicator
class BoutComm : public Communicator {
public:
  ~BoutComm();
  /// Get a pointer to the only instance
  static BoutComm* getInstance();

  /// Shortcut method
  static MPI_Comm get();

  static void setArgs(int &c, char** &v);
  
  static void cleanup();

  static int rank(); ///< Rank: my processor number
  static int size(); ///< Size: number of processors

  // Setting options
  void setComm(MPI_Comm c);

  // Getters
  MPI_Comm getComm();
  bool isSet();
  
  // Communicator interface
  int commRank() override {return rank();}
  int commSize() override {return size();}
  Request irecvBytes(void *data, int length, int destination, int identifier) override;
  void sendBytes(void *data, int length, int destination, int identifier) override;
  int waitAny(const std::vector<Request> &requests) override;
  
 private:
  BoutComm();
  
  int *pargc; char ***pargv; ///< Command-line arguments. These can be modified by MPI init, so pointers are used
  bool hasBeenSet;
  MPI_Comm comm;
  
  static BoutComm* instance; ///< The only instance of this class (Singleton)

};



#endif // __BOUTCOMM_H__
