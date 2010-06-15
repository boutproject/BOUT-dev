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


/****************************************************************
 *                       BOUNDARIES
 ****************************************************************/

void add_boundary_local(int direction, int ypos, int xge, int xlt);

/// Add a boundary in global indices
void add_boundary(int direction, int ypos, int xge, int xlt)
{
  // Calculate Y processor this boundary appears on
  if( (((float) ypos) / ((float) MYSUB)) == PE_YIND) {

    // In this Y range. Check X indices
    
    // Shift into local indices
    int loc_xge = XLOCAL(xge);
    int loc_xlt = XLOCAL(xlt);

    if((loc_xge >= mesh->ngx) || (loc_xlt <= 0)) {
      return; // Not in this x domain
    }
    
    if(loc_xge < 0)   loc_xge = 0;
    if(loc_xlt > mesh->ngx) loc_xlt = mesh->ngx;
  
    /// Send to routine in local indices
    add_boundary_local(direction, YLOCAL(ypos), loc_xge, loc_xlt);
  }
}

/// Add a boundary in the local processor's indices
void add_boundary_local(int direction, int ypos, int xge, int xlt)
{
  if(N_BOUNDARY <= 0) {
    BOUNDARY = (TBoundary**) malloc(sizeof(TBoundary*));
    N_BOUNDARY = 0;
  }else {
    BOUNDARY = (TBoundary**) realloc(BOUNDARY, sizeof(TBoundary*)*(N_BOUNDARY+1));
  }

  BOUNDARY[N_BOUNDARY] = (TBoundary*) malloc(sizeof(TBoundary));
  BOUNDARY[N_BOUNDARY]->direction = direction;
  BOUNDARY[N_BOUNDARY]->ypos = ypos;
  BOUNDARY[N_BOUNDARY]->xge = xge;
  BOUNDARY[N_BOUNDARY]->xlt = xlt;

  N_BOUNDARY++;
}

/****************************************************************
 *                       CONNECTIONS
 ****************************************************************/

/// Connection initialisation: Set processors in a simple 2D grid
void default_connections()
{
  DDATA_XSPLIT = UDATA_XSPLIT = 0;  // everything by default outside (arb. choice)
  DDATA_INDEST = UDATA_INDEST = -1; // since nothing inside

  DDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND - 1);
  UDATA_OUTDEST = PROC_NUM(PE_XIND, PE_YIND + 1);
  
  IDATA_DEST = PROC_NUM(PE_XIND - 1, PE_YIND);
  ODATA_DEST = PROC_NUM(PE_XIND + 1, PE_YIND);
  
  TS_up_in = TS_up_out = TS_down_in = TS_down_out = false; // No twist-shifts

  /// Check if X is periodic
  bool xperiodic;
  options.setSection(NULL);
  options.get("xperiodic", xperiodic, false);
  if(xperiodic) {
    if(PE_XIND == (NXPE-1))
      ODATA_DEST = PROC_NUM(0, PE_YIND);
    
    if(PE_XIND == 0)
      IDATA_DEST = PROC_NUM(NXPE-1, PE_YIND);
  }
}

/// Add a topology connection
/*!
 * Set ypos1 and ypos2 to be neighbours in the range xge <= x < xlt.
 * Optional argument ts sets whether to use twist-shift condition
 */
void set_connection(int ypos1, int ypos2, int xge, int xlt, bool ts = false)
{
  int ype1, ype2; // the two Y processor indices
  int ypeup, ypedown;
  int yind1, yind2;

  if(xlt <= xge)
    return;

  if((ypos1 < 0) || (ypos1 >= MY)) {
    output.write("WARNING adding connection: poloidal index %d out of range\n", ypos1);
    return;
  }
  if((ypos2 < 0) || (ypos2 >= MY)) {
    output.write("WARNING adding connection: poloidal index %d out of range\n", ypos2);
    return;
  }
  
  ype1 = YPROC(ypos1);
  ype2 = YPROC(ypos2);

  /* y index within processors */
  yind1 = YLOCAL(ypos1, ype1);
  yind2 = YLOCAL(ypos2, ype2);

  /* Check which boundary the connection is on */
  if((yind1 == MYG) && (yind2 == MYSUB+MYG-1)) {
    ypeup = ype2; /* processor sending data up (+ve y) */
    ypedown = ype1; /* processor sending data down (-ve y) */
  }else if((yind2 == MYG) && (yind1 == MYSUB+MYG-1)) {
    ypeup = ype1;
    ypedown = ype2;
  }else {
    output.write("ERROR adding connection: y index %d or %d not on processor boundary\n", ypos1, ypos2);
    exit(1);
  }

  /* check the x ranges are possible */
  if((xge != 0) && (xlt != MX)) {
    output.write("ERROR adding connection(%d,%d,%d,%d): can only divide X domain in 2\n",
		 ypos1, ypos2, xge, xlt);
    exit(1);
  }

  output.write("Connection between top of Y processor %d and bottom of %d in range %d <= x < %d\n",
	       ypeup, ypedown, xge, xlt);

  // Convert X coordinates into local indices

  xge = XLOCAL(xge);
  xlt = XLOCAL(xlt);
  
  if(( xge >= mesh->ngx ) || ( xlt <= 0 )) {
    return; // Not in this x domain
  }
  
  if(xge < 0)   xge = 0;
  if(xlt > mesh->ngx) xlt = mesh->ngx;

  if(MYPE == PROC_NUM(PE_XIND, ypeup)) { /* PROCESSOR SENDING +VE Y */
    /* Set the branch cut x position */
    if(xge <= MXG) {
      /* Connect on the inside */
      UDATA_XSPLIT = xlt;
      UDATA_INDEST = PROC_NUM(PE_XIND, ypedown);
      if(UDATA_XSPLIT == mesh->ngx)
	UDATA_OUTDEST = -1;

      TS_up_in = ts; // Twist-shift
      
      output.write("=> This processor sending in up\n");
    }else {
      /* Connect on the outside */
      UDATA_XSPLIT = xge;
      UDATA_OUTDEST = PROC_NUM(PE_XIND, ypedown);
      if(UDATA_XSPLIT == 0)
	UDATA_INDEST = -1;

      TS_up_out = ts;
      output.write("=> This processor sending out up\n");
    }
  }
  
  if(MYPE == PROC_NUM(PE_XIND, ypedown)) { /* PROCESSOR SENDING -VE Y */
    /* Set the branch cut x position */
    if(xge <= MXG) {
      /* Connect on the inside */
      DDATA_XSPLIT = xlt;
      DDATA_INDEST = PROC_NUM(PE_XIND, ypeup);
      if(DDATA_XSPLIT == mesh->ngx)
	DDATA_OUTDEST = -1;

      TS_down_in = ts;

      output.write("=> This processor sending in down\n");
    }else {
      /* Connect on the outside */
      DDATA_XSPLIT = xge;
      DDATA_OUTDEST = PROC_NUM(PE_XIND, ypeup);
      if(DDATA_XSPLIT == 0)
	DDATA_INDEST = -1;

      TS_down_out = ts;

      output.write("=> This processor sending out down\n");
    }
  }
}

/****************************************************************
 *                MAIN TOPOLOGY FUNCTION
 ****************************************************************/

void topology()
{
  // Default settings
  N_BOUNDARY = 0;   // No boundaries 

  // Perform checks common to all topologies

  if (NPES != NXPE*NYPE) {
    output.write("\tTopology error: npes=%d is not equal to NXPE*NYPE=%d\n",
		 NPES,NXPE*NYPE);
    exit(1);
  }
  if(MYSUB * NYPE != MY) {
    output.write("\tTopology error: MYSUB[%d] * NYPE[%d] != MY[%d]\n",MYSUB,NYPE,MY);
    exit(1);
  }
  if(MXSUB * NXPE != MX) {
    output.write("\tTopology error: MXSUB[%d] * NXPE[%d] != MX[%d]\n",MXSUB,NXPE,MX);
    exit(1);
  }

  if((NXPE > 1) && (MXSUB < MXG)) {
    output.write("\tERROR: Grid X size must be >= guard cell size\n");
    exit(1);
  }
  if(MYSUB < MYG) {
    output.write("\tERROR: Grid Y size must be >= guard cell size\n");
    exit(1);
  }
  
  if(jyseps2_1 == jyseps1_2) {
    /********* SINGLE NULL OPERATION *************/
    output.write("\tEQUILIBRIUM IS SINGLE NULL (SND) \n");

    /* Set separatrices - x location all the same */
    ixseps_inner = ixseps_outer = ixseps_upper = ixseps_lower = ixseps1;

    default_connections();
    set_connection(jyseps1_1+1, jyseps2_2, 0, ixseps1, true); // Twist-shift this connection
    set_connection(jyseps1_1, jyseps2_2+1, 0, ixseps1); // No twist-shift in PF region

    /* Set boundaries */

    /* format: direction, y position, minimum x, maximum x+1 */
    add_boundary(-1, 0, 0, MX); /* MYPE = 0 */
    add_boundary(1, ny-1, 0, MX); /* MYPE = NYPE-1 */
    
  }else {
    /*************** DOUBLE NULL OPERATION *******************/
    /* UPPER LEGS: Do not have to be the same length as each
       other or lower legs, but do have to have an integer number
       of processors */
    if((ny_inner-jyseps2_1-1) % MYSUB != 0) {
      output.write("\tTopology error: Upper inner leg does not have integer number of processors\n");
      exit(1);
    }
    if((jyseps1_2-ny_inner+1) % MYSUB != 0) {
      output.write("\tTopology error: Upper outer leg does not have integer number of processors\n");
    }

    if(ixseps1 == ixseps2) {
      /*************** CONNECTED (balanced) DOUBLE NULL ******************/
      output.write("\tEQUILIBRIUM IS CONNECTED DOUBLE NULL (CDND)\n");
      /* all separatrix indices the same */
      ixseps_inner = ixseps_outer = ixseps_upper = ixseps_lower = ixseps1;
        
    }else if(ixseps1 < ixseps2){
      /*************** LOWER DOUBLE NULL **********************/
      output.write("\tEQUILIBRIUM IS LOWER DOUBLE NULL (LDND)\n");
      ixseps_inner = ixseps_lower = ixseps1;
      ixseps_outer = ixseps_upper = ixseps2;
    }else {
      /*************** UPPER DOUBLE NULL **********************/
      output.write("\tEQUILIBRIUM IS UPPER DOUBLE NULL (UDND)\n");
      ixseps_inner = ixseps_upper = ixseps2;
      ixseps_outer = ixseps_lower = ixseps1;
    }
    
    /* Following code works for any Double Null */

    /********* DND CONNECTIONS **********/
    default_connections();
    /* Lower x-point */
    set_connection(jyseps1_1+1, jyseps2_2  , 0, ixseps_lower, ixseps1 <= ixseps2); /* Core */
    set_connection(jyseps1_1  , jyseps2_2+1, 0, ixseps_lower); /* PF   */
    /* Upper x-point */
    set_connection(jyseps2_1  , jyseps1_2+1, 0, ixseps_upper, ixseps1 > ixseps2); /* Core */
    set_connection(jyseps2_1+1, jyseps1_2  , 0, ixseps_upper); /* PF   */
    
    /********* DND STRIKE POINTS ********/
    /* Lower strike-points */
    add_boundary(-1, 0         , 0, MX); /* inner */
    add_boundary( 1, ny-1      , 0, MX); /* outer */
    /* Upper strike-points */
    add_boundary( 1, ny_inner-1, 0, MX); /* inner */
    add_boundary(-1, ny_inner  , 0, MX); /* outer */
  }

  MYPE_IN_CORE = 0; // processor not in core
  if( (ixseps_inner > 0) && ( ((PE_YIND*MYSUB > jyseps1_1) && (PE_YIND*MYSUB <= jyseps2_1)) || ((PE_YIND*MYSUB > jyseps1_2) && (PE_YIND*MYSUB <= jyseps2_2)) ) ) {
    MYPE_IN_CORE = 1; /* processor is in the core */
  }
  
  ///////////////////////////////////////////////
  // Create MPI groups
  
  MPI_Group group_world, comm_x_group, comm_y_group;
  MPI_Comm_group(MPI_COMM_WORLD, &group_world); // Group of all processors
  
  int *ranks = new int[NPES];
  // X communication group.
  
  for(int i=0;i<NXPE;i++)
    ranks[i] = PROC_NUM(i, PE_YIND); // All ranks at this Y
  MPI_Group_incl(group_world, NXPE, ranks, &comm_x_group);
  
  for(int i=0;i<NYPE;i++)
    ranks[i] = PROC_NUM(PE_XIND, i); // All ranks at this X
  MPI_Group_incl(group_world, NYPE, ranks, &comm_y_group);
  
  
  /*
  // Y communication group. Only defined in core; others more complicated
  if((PE_YIND <= YPROC(jyseps1_1)) || (PE_YIND > YPROC(jyseps2_2))) {
    // In lower PF region
    
    int npf1 = (jyseps1_1+1)/MYSUB;
    int npf2 = (MY - jyseps2_2 - 1)/MYSUB;
    
    ranks = new int[npf1 + npf2];
    for(int i=0;i<npf1;i++)
      ranks[i] = PROC_NUM(PE_XIND, i);
    for(int i=0;i<npf2;i++)
      ranks[npf1 + i] = PROC_NUM(PE_XIND, YPROC(jyseps2_2+1)+i);
    
    MPI_Group_incl(group_world, n, ranks, &comm_y_group);
    
  }else if((PE_YIND > YPROC(jyseps1_2)) && (PE_YIND < YPROC(jyseps2_1+1))) {
    // Upper PF region
    
  }else {
    // In core
    
  }
  */
  
  // Create communicators
  MPI_Comm_create(MPI_COMM_WORLD, comm_x_group, &comm_x);
  MPI_Comm_create(MPI_COMM_WORLD, comm_y_group, &comm_y);

  delete[] ranks;

  
  /****************** LIMITERS **********************
   * To add a limiter between y position 40 and 41  *
   * and radial position 25 <= x < MX use:          *
add_boundary( 1, 40, 25, MX);
add_boundary(-1, 41, 25, MX);
   * Determining which processor this               *
   * applies to is done automatically. The first    *
   * argument, direction = -1 if going towards      *
   * the boundary is in the -ve y direction.        *
   * If 40 is outboard, The first command           *
   * adds the geometric upper edge of the limiter,  *
   * whilst the second adds the lower edge          *
   **************************************************/

  // Print out settings
  output.write("\tMYPE_IN_CORE = %d\n", MYPE_IN_CORE);
  output.write("\tDXS = %d, DIN = %d. DOUT = %d\n", 
	       DDATA_XSPLIT, DDATA_INDEST, DDATA_OUTDEST);
  output.write("\tUXS = %d, UIN = %d. UOUT = %d\n", 
	       UDATA_XSPLIT, UDATA_INDEST, UDATA_OUTDEST);
  output.write("\tXIN = %d, XOUT = %d\n",
	       IDATA_DEST, ODATA_DEST);

  output.write("\tTwist-shift: ");
  if(TS_down_in)
    output.write("DI ");
  if(TS_down_out)
    output.write("DO ");
  if(TS_up_in)
    output.write("UI ");
  if(TS_up_out)
    output.write("UO ");
  output.write("\n");
}
