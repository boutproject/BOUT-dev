#include "gtest/gtest.h"

#include "dcomplex.hxx"
#include "fft.hxx"
#include "test_extras.hxx"
#include "bout/array.hxx"
#include "bout/constants.hxx"

#include <algorithm>
#include <iostream>
#include <numeric>

TEST(FFTTest, rfft) {

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int size{8};
  constexpr int nmodes{(size / 2) + 1};

  // Make grid indices from [0, size - 1]
  Array<BoutReal> indices{size};
  std::iota(indices.begin(), indices.end(), 0.0);

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> input{size};
  std::transform(indices.begin(), indices.end(), input.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  Array<dcomplex> output{nmodes};

  // Compute forward real FFT
  rfft(input.begin(), size, output.begin());

  // Real part should have exactly one frequency
  Array<BoutReal> expected_real{nmodes};
  std::fill(expected_real.begin(), expected_real.end(), 0.);
  expected_real[2] = 0.5;

  // Imaginary part should have exactly one frequency
  Array<BoutReal> expected_imag{nmodes};
  std::fill(expected_imag.begin(), expected_imag.end(), 0.);
  expected_imag[1] = -0.5;

  EXPECT_EQ(output.size(), nmodes);

  for (int i = 0; i < nmodes; ++i) {
    EXPECT_NEAR(real(output[i]), expected_real[i], BoutRealTolerance);
    EXPECT_NEAR(imag(output[i]), expected_imag[i], BoutRealTolerance);
  }
}

TEST(FFTTest, irfft) {

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int nmodes{5};
  constexpr int size{(nmodes - 1) * 2};

  // Make a spectrum with two frequencies
  Array<dcomplex> input{nmodes};
  std::fill(input.begin(), input.end(), dcomplex{0., 0.});
  input[1] = dcomplex{0., -0.5};
  input[2] = dcomplex{0.5, 0.};

  Array<BoutReal> output{size};

  // Compute inverse real FFT
  irfft(input.begin(), size, output.begin());

  // Make grid indices from [0, size - 1]
  Array<BoutReal> indices{size};
  std::iota(indices.begin(), indices.end(), 0.0);

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> expected{size};
  std::transform(indices.begin(), indices.end(), expected.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  EXPECT_EQ(output.size(), size);

  for (int i = 0; i < size; ++i) {
    EXPECT_NEAR(output[i], expected[i], BoutRealTolerance);
  }
}

TEST(FFTTest, rfftWithArray) {

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int size{8};
  constexpr int nmodes{(size / 2) + 1};

  // Make grid indices from [0, size - 1]
  Array<BoutReal> indices{size};
  std::iota(indices.begin(), indices.end(), 0.0);

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> input{size};
  std::transform(indices.begin(), indices.end(), input.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  // Compute forward real FFT
  Array<dcomplex> output = bout::fft::rfft(input);

  // Real part should have exactly one frequency
  Array<BoutReal> expected_real{nmodes};
  std::fill(expected_real.begin(), expected_real.end(), 0.);
  expected_real[2] = 0.5;

  // Imaginary part should have exactly one frequency
  Array<BoutReal> expected_imag{nmodes};
  std::fill(expected_imag.begin(), expected_imag.end(), 0.);
  expected_imag[1] = -0.5;

  EXPECT_EQ(output.size(), nmodes);

  for (int i = 0; i < nmodes; ++i) {
    EXPECT_NEAR(real(output[i]), expected_real[i], BoutRealTolerance);
    EXPECT_NEAR(imag(output[i]), expected_imag[i], BoutRealTolerance);
  }
}

TEST(FFTTest, irfftWithArray) {

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int nmodes{5};
  constexpr int size{(nmodes - 1) * 2};

  // Make a spectrum with two frequencies
  Array<dcomplex> input{nmodes};
  std::fill(input.begin(), input.end(), dcomplex{0., 0.});
  input[1] = dcomplex{0., -0.5};
  input[2] = dcomplex{0.5, 0.};

  // Compute inverse real FFT
  Array<BoutReal> output = bout::fft::irfft(input);

  // Make grid indices from [0, size - 1]
  Array<BoutReal> indices{size};
  std::iota(indices.begin(), indices.end(), 0.0);

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> expected{size};
  std::transform(indices.begin(), indices.end(), expected.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  EXPECT_EQ(output.size(), size);

  for (int i = 0; i < size; ++i) {
    EXPECT_NEAR(output[i], expected[i], BoutRealTolerance);
  }
}

TEST(FFTTest, RoundTrip) {
  // Checks normalisation is == 1

  // Make sure fft functions are quiet by setting fft_measure to false
  bout::fft::fft_init(false);

  constexpr int size{8};

  // Make grid indices from [0, size - 1]
  Array<BoutReal> indices{size};
  std::iota(indices.begin(), indices.end(), 0.0);

  // Calculate sin(x) + cos(2x) on [0, 2pi]
  Array<BoutReal> input{size};
  std::transform(indices.begin(), indices.end(), input.begin(),
                 [](BoutReal i) -> BoutReal {
                   return std::sin(i * TWOPI / size) + std::cos(2. * i * TWOPI / size);
                 });

  Array<BoutReal> output{bout::fft::irfft(bout::fft::rfft(input))};

  EXPECT_EQ(output.size(), input.size());

  for (int i = 0; i < size; ++i) {
    EXPECT_NEAR(output[i], input[i], BoutRealTolerance);
  }
}
