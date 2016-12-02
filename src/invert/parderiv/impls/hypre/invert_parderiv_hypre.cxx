/************************************************************************
 * Inversion of parallel derivatives
 * 
 * Inverts a matrix of the form 
 *
 * A + B * Grad2_par2
 * 
 * - Each flux surface can be solved independently
 * - By taking FFT in Z direction, toroidal modes can be
 *   solved independently
 * - Result is a set of complex band-diagonal matrices to invert
 * 
 * Author: Ben Dudson, University of York, Oct 2011
 * 
 * Known issues:
 * ------------
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
 ************************************************************************/

#include "mpi.h"

#include <invert_parderiv_hypre.hxx>
#include <utils.hxx>
#include <fft.hxx>

InvertParHypre::InvertParHypre(Mesh *mi) {
  // Store pointer to the mesh
  m = mi;
  
  // Allocate complex arrays for input rhs and result
  crhs = cmatrix(m->ngy, (m->LocalNz)/2 + 1);
  cresult = cmatrix(m->ngy, (m->LocalNz)/2 + 1);
  
  // Create Hypre structured grids
  grid = new HYPRE_StructGrid[m->ngx];
  
  // Define discretisation stencil
  HYPRE_StructStencilCreate(1, 10, &stencil);
  
  // Define the geometry of the stencil
  int offsets[5][1] = {{-2},{-1},{0},{1},{2}};
  for(int entry=0;entry<5;entry++)
    HYPRE_StructStencilSetElement(stencil, entry, offsets[entry]);

  // Allocate matrices
  mat = new HYPRE_StructMatrix*[m->ngx];
  mat[0] = new HYPRE_StructMatrix[m->ngx * ( (m->LocalNz)/2 + 1 )];
  for(int i=1;i<m->ngx;i++)
    mat[i] = mat[i-1] + (m->LocalNz)/2 + 1;

  // Loop over surfaces
  SurfaceIter *surf = m->iterateSurfaces();
  for(surf->first(); !surf->isDone(); surf->next()) {
    // Create empty 1D grid
    
    int np;
    MPI_Comm_size(surf->communicator(), &np);
    int myp;
    MPI_Comm_rank(surf->communicator(), &myp);
    
    // Create an empty grid object
    HYPRE_StructGridCreate(surf->communicator(), 1, grid+surf->xpos);
    
    // Set local extents
    int ys = m->ystart;
    int ye = m->yend;
    if(!surf->closed()) {
      if(myp == 0) {
        // On lower boundary
        ys = 0;
      }else if(myp == (np-1)) {
        // On Upper boundary
        ye = m->ngy - 1;
      }
    }
    int nyloc = ye - ys + 1; // Number of y points on this processor
    
    int ilower[1], iupper[1]; // Turn into global indices
    ilower[0] = surf->yGlobal(ys);
    iupper[0] = surf->yGlobal(ye);
    
    HYPRE_StructGridSetExtents(grid[surf->xpos], ilower, iupper);
    
    // Set periodicity
    if(surf->closed()) {
      int period = surf->ySize();
      HYPRE_StructGridSetPeriodic(grid[surf->xpos], &period);
    }
    
    // Collective call to finalise the grid assembly
    HYPRE_StructGridAssemble(grid[surf->xpos]);
    
    // Set up the matrices, one for each k
    int nk = (m->LocalNz)/2 + 1;
    
    int stencil_indices[5] = {0,1,2,3,4}; // Labels for stencil entries
    int *values = new int[nyloc * 5];
    for(int k=0;k<nk;k++) {
      HYPRE_StructMatrixCreate(surf->communicator(), grid[surf->xpos], stencil, &mat[surf->xpos][k]);
      HYPRE_StructMatrixInitialize(mat[surf->xpos][k]);
      
      for(int y=0;y<nyloc; y++) {
        int yc = y + ys;

        BoutReal dy = m->dy[surf->xpos][yc];

        values[5*y] = A[surf->xpos][yc] + B[surf->xpos][yc] * 1./12. ;
        
        
      }
    }
    
  }
}

InvertParHypre::~InvertParHypre() {
  free_cmatrix(crhs);
  free_cmatrix(cresult);
}

const Field3D InvertParHypre::solve(const Field3D &f) {
  return solve(f, f / A);
}

const Field3D InvertParHypre::solve(const Field3D &f, const Field3D &start) {
  Field3D result;
  result.allocate();
  
  // Decide on x range. Solve in boundary conditions
  int xs = (m->firstX()) ? 0 : 2;
  int xe = (m->lastX()) ? m->ngx-1 : (m->ngx-3);

  int nxsolve = xe - xs + 1; // Number of X points to solve    
  int nylocal = m->yend - m->ystart + 1;
  
  // Loop over surfaces
  SurfaceIter *surf = m->iterateSurfaces();
  for(surf->first(); !surf->isDone(); surf->next()) {
    // Take FFT, shifting in Z to get into orthogonal coordinates
    for(int y=0;y<m->ngy;y++)
      ZFFT(f[surf->xpos][y], m->zShift[surf->xpos][y], crhs[y]);
    
    // Perform inversion
    
    
    // FFT back
    for(int y=0;y<m->ngy;y++)
      ZFFT_rev(cresult[y], m->zShift[surf->xpos][y], result[surf->xpos][y]);
  }
  
}

