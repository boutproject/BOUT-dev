/*******************************************************************************
 * FFT routines
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
 *******************************************************************************/

#ifndef __FFT_H__
#define __FFT_H__

#include "dcomplex.h"

void cfft(dcomplex *cv, int length, int isign);
void ZFFT(dcomplex *cv, real zoffset, int isign, bool shift = true);

// more optimised code (for reals)

void rfft(real *in, int length, dcomplex *out);
void irfft(dcomplex *in, int length, real *out);

void ZFFT(real *in, real zoffset, dcomplex *cv, bool shift = true);
void ZFFT_rev(dcomplex *cv, real zoffset, real *out, bool shift = true);

#endif // __FFT_H__
