/**************************************************************************
 * Maintains a list of values which are printed every time-step
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

#include "mpi.h"

#include <diagnos.hxx>
#include <utils.hxx>
#include <globals.hxx>
#include <output.hxx>
#include <msg_stack.hxx>

#define PVEC_REAL_MPI_TYPE MPI_DOUBLE

/// Global variable initialisation
bool Diagnos::init = false;
bool Diagnos::global_vals = false;

Diagnos::Diagnos()
{
  if(!init) {
    output.write("Initialising diagnostics\n");
    Options *options = Options::getRoot();
    options = options->getSection("diagnos");
    options->get("global", global_vals, true);
    
    init = true;
  }
}

Diagnos::~Diagnos()
{

}

void Diagnos::add(FieldData &f, DIAGNOS_FUNC func, int x, int y, int z, int component, const char* label)
{
  diag_item i;
  
  i.func = func;
  
  i.label = copy_string(label);

  i.var = &f;
  i.x = x;
  i.y = y;
  i.z = z;

//  int nr = i.var->BoutRealSize(); // Number of BoutReals per point

  i.component = component;

  item.push_back(i);
}

/// Calculate the values and return in an array
const vector< BoutReal > Diagnos::run()
{
#ifdef CHECK
  msg_stack.push("Diagnos::run\n");
#endif

  vector< BoutReal > result;
  
  for(std::vector< diag_item >::iterator it = item.begin(); it != item.end(); it++) {
    result.push_back(run(*it));
  }

#ifdef CHECK
  msg_stack.pop();
#endif

  return result;
}

BoutReal Diagnos::run(const diag_item &i)
{
  /*
  switch(i.func) {
  case DIAG_INDX: { ///< Index 
    if(global_vals) {
      int nr = i.var->BoutRealSize(); // Number of BoutReals per point
      if(nr <= 0)
	return 0.0;
      
      static BoutReal *rptr;
      static int rlen = 0;
      
      if(rlen < nr) {
	if(rlen > 0)
	  delete[] rptr;
	rptr = new BoutReal[nr];
	rlen = nr;
      }

      // Use a global index
      int np = PROC_NUM(i.x, i.y);
      if((np < 0) || (np >= NPES))
	return 0.0;
      
      BoutReal val = 0.0;
      if(MYPE == np) {
	i.var->getData(XLOCAL(i.x), YLOCAL(i.y), i.z, rptr);
	int c = i.component;
	if((i.component < 0) || (i.component >= nr))
	  c = 0;
	val = rptr[c];
      }
      
      MPI_Bcast ( &val, 1, PVEC_REAL_MPI_TYPE, np, BoutComm::get()); 
      
      return val;
    }else {
      // Local index
      
      if((i.x < 0) || (i.x >= mesh->ngx) ||
	 (i.y < 0) || (i.y >= mesh->ngy) ||
	 (i.z < 0) || (i.z >= mesh->ngz)) {
	return 0.0;
      }
    }
    break;
  }
  case DIAG_MAX: {
    
    
    
    if(global_vals) {
      
    }else {
      
    }
    
    break;
  }
  default:
    break;
  };

  */
  return 0.0;
}
