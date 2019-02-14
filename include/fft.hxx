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
#include <bout/array.hxx>

class Options;

namespace bout {
namespace fft {

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
 * \param[in] in     Pointer to the 1D array to take the fourier transform of
 * \param[in] length Number of points in the input array
 * \param[out] out   Pointer to the complex 1D array which is the FFT of in
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
 * \param[in] in     Pointer to the 1D array to take the inverse fourier
 *                   transform of
 * \param[in] length Number of points in the input array
 * \param[out] out   Pointer to the complex 1D array which is IFFTed
 */
void irfft(const dcomplex *in, int length, BoutReal *out);

/*!
 * Discrete Sine Transform
 *
 * \p in and \p out arrays must both be of the same \p length
 */
void DST(const BoutReal *in, int length, dcomplex *out);

/*!
 * Inverse Discrete Sine Transform
 *
 * \p in and \p out arrays must both be of the same \p length
 */
void DST_rev(dcomplex *in, int length, BoutReal *out);

/// Should the FFT functions find and use an optimised plan?
void fft_init(bool fft_measure);
/// Should the FFT functions find and use an optimised plan?
///
/// If \p options is not nullptr, it should contain a bool called
/// "fftw_measure". If it is nullptr, use the global `Options` root
void fft_init(Options* options = nullptr);

/// Returns the fft of a real signal \p in using fftw_forward
Array<dcomplex> rfft(const Array<BoutReal>& in);

/// Take the inverse fft of signal \p in where the outputs are only reals.
/// Requires the \p length of the original real signal
///
/// \p length is required because input signals to the forward
/// transform of length `n` and `n + 1` both produce ffts of length
/// `(n / 2) + 1` -- i.e. it's not possible to recover the length of
/// the original signal from the fft alone.
///
/// Expects that `in.size() == (length / 2) + 1`
Array<BoutReal> irfft(const Array<dcomplex>& in, int length);

} // namespace fft
} // namespace bout

// Legacy non-namespaced versions

inline void rfft(const BoutReal *in, int length, dcomplex *out) {
  return bout::fft::rfft(in, length, out);
}

inline void irfft(const dcomplex *in, int length, BoutReal *out) {
  return bout::fft::irfft(in, length, out);
}

inline void DST(const BoutReal *in, int length, dcomplex *out) {
  return bout::fft::DST(in, length, out);
}

inline void DST_rev(dcomplex *in, int length, BoutReal *out) {
  return bout::fft::DST_rev(in, length, out);
}

#endif // __FFT_H__
