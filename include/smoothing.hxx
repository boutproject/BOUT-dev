/**************************************************************
 * \file smoothing.hxx
 * 
 * Smoothing operators
 *
 **************************************************************
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
 **************************************************************/

#ifndef __SMOOTHING_H__
#define __SMOOTHING_H__

#include "field3d.hxx"

/// Smooth in X using simple 1-2-1 filter
const Field3D smooth_x(const Field3D &f);

/// Smooth in Y using 1-2-1 filter
const Field3D smooth_y(const Field3D &f);

/// Smooth using a stencil in X and Y
const Field3D smoothXY(const Field3D &f);

/*! 
 * Average over Y
 *
 * Issues
 * ======
 * 
 * Important: Only works if there are no branch cuts
 *
 * Assumes every processor has the same domain shape
 * 
 */
const Field2D averageY(const Field2D &f);

/*!
 * Average in Y
 * 
 * Issues
 * ======
 * 
 * Important: Only works if there are no branch cuts
 *
 * Creates static arrays
 * 
 * Not thread safe
 *
 * Assumes every processor has the same domain shape
 * 
 */
const Field3D averageY(const Field3D &f);

/// Average over X
const Field2D averageX(const Field2D &f);
const Field3D averageX(const Field3D &f);

/*!
  Volume integral of Field2D variable
  Developed by T. Rhee and S. S. Kim
  
  Issues
  ======
  
  Assumes every processor has the same domain shape
  
  Will only work if X communicator is constant in Y
  so no processor/branch cuts in X
 */
BoutReal Average_XY(const Field2D &var);

/// Volume integral of Field2D variable
/// which uses Average_XY
BoutReal Vol_Integral(const Field2D &var);

/// Nonlinear filtering to remove grid-scale noise in X
/*!
  From a paper:

  W.Shyy et. al. JCP 102 (1) September 1992 page 49

  "On the Suppression of Numerical Oscillations Using a Non-Linear Filter"
  
 */
const Field3D nl_filter_x(const Field3D &f, BoutReal w=1.0);

/// Nonlinear filtering to remove grid-scale noise in Y
/*!
  From a paper:

  W.Shyy et. al. JCP 102 (1) September 1992 page 49

  "On the Suppression of Numerical Oscillations Using a Non-Linear Filter"
  
 */
const Field3D nl_filter_y(const Field3D &f, BoutReal w=1.0);

/// Nonlinear filtering to remove grid-scale noise in Z
/*!
  From a paper:

  W.Shyy et. al. JCP 102 (1) September 1992 page 49

  "On the Suppression of Numerical Oscillations Using a Non-Linear Filter"
  
 */
const Field3D nl_filter_z(const Field3D &f, BoutReal w=1.0);

/// Nonlinear filtering to remove grid-scale noise in X,Y and Z
/*!
  From a paper:

  W.Shyy et. al. JCP 102 (1) September 1992 page 49

  "On the Suppression of Numerical Oscillations Using a Non-Linear Filter"
  
 */
const Field3D nl_filter(const Field3D &f, BoutReal w=1.0);

#endif // __SMOOTHING_H__
