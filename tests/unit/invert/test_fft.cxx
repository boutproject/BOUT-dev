#include "gtest/gtest.h"

#include "bout/array.hxx"
#include "bout/constants.hxx"
#include "dcomplex.hxx"
#include "fft.hxx"
#include "test_extras.hxx"

#include <algorithm>
#include <iostream>
#include <numeric>

TEST(FFTTest, rfft) {

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int size{8};

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> indices(size);
  std::iota(indices.begin(), indices.end(), 0);

  Array<BoutReal> input(size);
  std::transform(indices.begin(), indices.end(), input.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  Array<dcomplex> output(size);

  // Compute forward real FFT
  rfft(input.begin(), size, output.begin());

  // Real part should have exactly one frequency
  Array<BoutReal> expected_real(size);
  std::fill(expected_real.begin(), expected_real.end(), 0.);
  expected_real[2] = 0.5;

  // Imaginary part should have exactly one frequency
  Array<BoutReal> expected_imag(size);
  std::fill(expected_imag.begin(), expected_imag.end(), 0.);
  expected_imag[1] = -0.5;

  for (int i = 0; i < size; ++i) {
    EXPECT_NEAR(real(output[i]), expected_real[i], BoutRealTolerance);
    EXPECT_NEAR(imag(output[i]), expected_imag[i], BoutRealTolerance);
  }
}

TEST(FFTTest, irfft) {

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int size{8};

  // Make a spectrum with two frequencies
  Array<dcomplex> input(size);
  std::fill(input.begin(), input.end(), dcomplex{0., 0.});
  input[1] = dcomplex{0., -0.5};
  input[2] = dcomplex{0.5, 0.};

  Array<BoutReal> output(size);

  // Compute inverse real FFT
  irfft(input.begin(), size, output.begin());

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> indices(size);
  std::iota(indices.begin(), indices.end(), 0);

  Array<BoutReal> expected(size);
  std::transform(indices.begin(), indices.end(), expected.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  for (int i = 0; i < size; ++i) {
    EXPECT_NEAR(output[i], expected[i], BoutRealTolerance);
  }
}
