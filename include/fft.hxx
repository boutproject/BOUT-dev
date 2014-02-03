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

#include "dcomplex.hxx"

void cfft(dcomplex *cv, int length, int isign);
void ZFFT(dcomplex *cv, BoutReal zoffset, int isign, bool shift = true);

// more optimised code (for BoutReals)

void rfft(BoutReal *in, int length, dcomplex *out);
void irfft(dcomplex *in, int length, BoutReal *out);

void ZFFT(BoutReal *in, BoutReal zoffset, dcomplex *cv, bool shift = true);
void ZFFT_rev(dcomplex *cv, BoutReal zoffset, BoutReal *out, bool shift = true);
// Discrete Sine Transform
void DST(BoutReal *in, int length, dcomplex *out);
void DST_rev(dcomplex *in, int length, BoutReal *out);

#endif // __FFT_H__
