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

#endif // __TOPOLOGY_H__
