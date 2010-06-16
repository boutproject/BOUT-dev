/*!************************************************************************
 * Topology - detect topology of grid
 * Set communication values for Y communications:
 * UDATA_XSPLIT, UDATA_INDEST, UDATA_OUTDEST
 * DDATA_XSPLIT, DDATA_INDEST, DDATA_OUTDEST
 *
 * and X communications:
 * IDATA_DEST, ODATA_DEST
 * 
 * Set whether processor is in core: MYPE_IN_CORE
 * 
 * Set separatrix locations:
 * ixseps_inner, ixseps_outer, ixseps_upper, ixseps_lower
 * based on uedge.grd values:
 * ixseps1, ixseps2, jyseps*_*
 * 
 * NOTE: THIS CODE NEEDS TO BE MADE MORE GENERAL
 *       SORT OUT!!
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
 *************************************************************************/

#define TOPGLOBORIGIN
#include "mesh_topology.h"
#include "globals.h"

#include <stdlib.h>

/****************************************************************
 *                 GRID INDEX ROUTINES
 * 
 * These routines translate between local and global coordinates
 * This so that more flexible indexing (e.g. not equal on processors)
 * can be implemented later (maybe)
 ****************************************************************/

/// Returns the processor number, given X and Y processor indices.
/*!
 * If out of range returns -1 (no processor)
 */
int PROC_NUM(int xind, int yind)
{
  if((xind >= NXPE) || (xind < 0))
    return -1;
  if((yind >= NYPE) || (yind < 0))
    return -1;
  
  return yind * NXPE + xind;
}

/// Returns true if the given grid-point coordinates are in this processor
bool IS_MYPROC(int xind, int yind)
{
  return ((xind / MXSUB) == PE_XIND) && ((yind / MYSUB) == PE_YIND);
}

/// Returns the global X index given a local index
int XGLOBAL(int xloc)
{
  return xloc + PE_XIND * MXSUB;
}

/// Returns a local X index given a global index
int XLOCAL(int xglo)
{
  return xglo - PE_XIND * MXSUB;
}

/// Returns the global Y index given a local index
int YGLOBAL(int yloc)
{
  return yloc + PE_YIND*MYSUB - MYG;
}

/// Global Y index given local index and processor
int YGLOBAL(int yloc, int yproc)
{
  return yloc + yproc*MYSUB - MYG;
}

/// Returns a local Y index given a global index
int YLOCAL(int yglo)
{
  return yglo - PE_YIND*MYSUB + MYG;
}

int YLOCAL(int yglo, int yproc)
{
  return yglo - yproc*MYSUB + MYG;
}

/// Return the Y processor number given a global Y index
int YPROC(int yind)
{
  return yind / MYSUB;
}

/// Return the X processor number given a global X index
int XPROC(int xind)
{
  return (xind >= MXG) ? (xind - MXG) / MXSUB : 0;
}

/*
bool PE_IN_UPPER_PF(int yproc)
{
  return (PE_YIND > YPROC(jyseps1_2)) && (PE_YIND < YPROC(jyseps2_1+1));
}

bool PE_IN_LOWER_PF(int yproc)
{
  return (yproc <= YPROC(jyseps1_1)) || (yproc > YPROC(jyseps2_2));
					 }

bool PE_IN_CORE(int yproc)
{
  return (!PE_IN_LOWER_PF(yproc)) && (!PE_IN_UPPER_PF(yproc));
}

bool PE_IN_MY_REGION(int yproc)
{
  if(PE_IN_UPPER_PF(MYPE)) {
    return PE_IN_UPPER_PF(yproc);    
  }else if(PE_IN_LOWER_PF(MYPE)) {
    return PE_IN_LOWER_PF(yproc);
  }

  return PE_IN_CORE(yproc);
}
*/

