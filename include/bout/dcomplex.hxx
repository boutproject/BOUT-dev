/*!
 * \file dcomplex.hxx
 * \brief Complex number class definition
 *
 * Changelog
 * ---------
 *
 * 2015-03-09 Ben Dudson <bd512@york.ac.uk>
 *    o Removed, redefined in terms of std::complex
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
 */
#ifndef __DCOMPLEX_H__
#define __DCOMPLEX_H__

#include <complex>
#include "bout_types.hxx"

using dcomplex = std::complex<BoutReal>;

const dcomplex Im(0,1); // 1i 

/// Complex type for passing data to/from FORTRAN
struct fcmplx {
  BoutReal r, i;
};

#endif // __DCOMPLEX_H__
