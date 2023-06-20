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

/// Class to represent the 'global' communicator
class BoutComm {
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

 private:
  BoutComm();

  int* pargc{nullptr};
  char*** pargv{nullptr}; ///< Command-line arguments. These can be modified by MPI init,
                          ///< so pointers are used
  bool hasBeenSet{false};
  MPI_Comm comm;
  
  static BoutComm* instance; ///< The only instance of this class (Singleton)

};

#endif // __BOUTCOMM_H__
