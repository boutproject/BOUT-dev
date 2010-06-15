/**************************************************************************
 * Topology - detect topology of grid
 * Set communication values:
 * UDATA_XSPLIT, UDATA_INDEST, UDATA_OUTDEST
 * DDATA_XSPLIT, DDATA_INDEST, DDATA_OUTDEST
 *
 * Set whether processor is in core: MYPE_IN_CORE
 * 
 * based on uedge.grd values:
 * ixseps1, ixseps2, jyseps*_*
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

#ifndef __TOPOLOGY_H__
#define __TOPOLOGY_H__

#include "globals.h"

int PROC_NUM(int xind, int yind); // (PE_XIND, PE_YIND) -> MYPE
int XGLOBAL(int xloc);
int XLOCAL(int xglo);
int YGLOBAL(int yloc);
int YGLOBAL(int yloc, int yproc);
int YLOCAL(int yglo);
int YLOCAL(int yglo, int yproc);
int YPROC(int yind);
int XPROC(int xind);

void add_boundary(int direction, int ypos, int xge, int xlt);
void topology();

// BOUNDARY POSITIONS

typedef struct { // Definition of a boundary
  int direction; // +1 -> boundary +ve y (e.g. MYPE=NYPE-1) , -1 -> -ve y (e.g. MYPE = 0)
  int ypos;      // y index of boundary
  int xge;       // minimum x of boundary
  int xlt;       // x less than this value   (xge <= x < xlt)
}TBoundary;

#ifdef TOPGLOBORIGIN
#define GLOBAL
#else
#define GLOBAL extern
#endif

GLOBAL int N_BOUNDARY;
GLOBAL TBoundary **BOUNDARY; // Array of pointers to TBoundary structures

#undef GLOBAL

#endif // __TOPOLOGY_H__
