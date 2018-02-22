/*!******************************************************************************
 * \file fft.hxx
 *
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

/*!
 * Returns the fft of a real signal using fftw_forward
 *
 * The fftw_forward returns
 * out_k = sum_{j=0}^(length-1) in_j*exp(-2*pi*j*k*sqrt(-1)/length)
 *
 * Thus, out_k must be divided by 'length' in order for DFT[IDFT[in]] = in
 * where IDFT is the inverse fourier transform.
 * See the the fftw user manual for details.
 *
 * \param[in] in  Pointer to the 1D array to take the fourier transform of
 * \param[in] length Number of points in the input array
 *
 * \param[out] Pointer to the complex 1D array which is the FFT of in
 */
void rfft(const BoutReal *in, int length, dcomplex *out);

/*!
 * Take the inverse fft of signal where the outputs are only reals.
 *
 * This is done through a call to fftw_plan_dft_c2r_1d
 * which is calling fftw_backwards.
 *
 * That is
 * out_k = sum_{j=0}^(length-1) in_j*exp(2*pi*j*k*sqrt(-1)/length)
 *
 * See the the fftw user manual for details.
 *
 * \param[in] in Pointer to the 1D array to take the inverse fourier transform
 *               of
 * \param[in] length Number of points in the input array
 *
 * \param[out] Pointer to the complex 1D array which is IFFTed
 */
void irfft(const dcomplex *in, int length, BoutReal *out);

/*!
 * Discrete Sine Transform
 *
 * in and out arrays must both be of the same length
 */
void DST(const BoutReal *in, int length, dcomplex *out);

/*!
 * Inverse Discrete Sine Transform
 *
 * in and out arrays must both be of the same length
 */
void DST_rev(dcomplex *in, int length, BoutReal *out);

#endif // __FFT_H__
