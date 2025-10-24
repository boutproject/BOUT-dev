/*
 * Implements the shifted metric method for parallel derivatives
 *
 * By default fields are stored so that X-Z are orthogonal,
 * and so not aligned in Y.
 *
 */

#include "bout/parallel_boundary_region.hxx"
#include "bout/paralleltransform.hxx"
#include <bout/boundary_region.hxx>
#include <bout/constants.hxx>
#include <bout/fft.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
#include <bout/sys/timer.hxx>

#include <cmath>

#if BOUT_HAS_CUDA
#include <bout/twiddle.hxx>
#include <cuda_runtime.h>
#endif

ShiftedMetric::ShiftedMetric(Mesh& m, CELL_LOC location_in, Field2D zShift_,
                             BoutReal zlength_in, Options* opt)
    : ParallelTransform(m, opt), location(location_in), zShift(std::move(zShift_)),
      zlength(zlength_in) {
  ASSERT1(zShift.getLocation() == location);
  // check the coordinate system used for the grid data source
  ShiftedMetric::checkInputGrid();

  cachePhases();
}

void ShiftedMetric::checkInputGrid() {
  std::string parallel_transform;
  if (mesh.isDataSourceGridFile()
      and !mesh.get(parallel_transform, "parallel_transform")) {
    if (parallel_transform != "shiftedmetric") {
      throw BoutException("Incorrect parallel transform type '" + parallel_transform
                          + "' used to generate metric components for ShiftedMetric. "
                            "Should be 'shiftedmetric'.");
    }
  } // else: parallel_transform variable not found in grid input, indicates older input
  //       file or grid from options so must rely on the user having ensured the type is
  //       correct
}

void ShiftedMetric::outputVars(Options& output_options) {
  Timer time("io");
  const std::string loc_string =
      (location == CELL_CENTRE) ? "" : "_" + toString(location);

  output_options["zShift" + loc_string].force(zShift, "ShiftedMetric");
}

void ShiftedMetric::cachePhases() {
  // If we wanted to be efficient we could move the following cached phase setup
  // into the relevant shifting routines (with static bool first protection)
  // so that we only calculate the phase if we actually call a relevant shift
  // routine -- however as we're only going to do this initialisation once I
  // think it's cleaner to put it in the constructor here.

  // As we're attached to a mesh we can expect the z direction to
  // not change once we've been created so precalculate the complex
  // phases used in transformations
  nmodes = mesh.LocalNz / 2 + 1;

  // Allocate storage for our 3d phase information.
  fromAlignedPhs = Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);
  toAlignedPhs = Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);

  // To/From field aligned phases
  BOUT_FOR(i, mesh.getRegion2D("RGN_ALL")) {
    int ix = i.x();
    int iy = i.y();
    for (int jz = 0; jz < nmodes; jz++) {
      BoutReal kwave = jz * 2.0 * PI / zlength; // wave number is 1/[rad]
      fromAlignedPhs(ix, iy, jz) =
          dcomplex(cos(kwave * zShift[i]), -sin(kwave * zShift[i]));
      toAlignedPhs(ix, iy, jz) = dcomplex(cos(kwave * zShift[i]), sin(kwave * zShift[i]));
    }
  }

  // Allocate space for parallel slice caches: y-guard cells in each
  // direction
  parallel_slice_phases.resize(mesh.ystart * 2);

  // Careful with the indices/offsets! Offsets are 1-indexed (as 0
  // would be the original slice), and Mesh::ystart is the number of
  // guard cells. The parallel slice vector stores the offsets as
  //    {+1, ..., +n, -1, ..., -n}
  // Once parallel_slice_phases is initialised though, each element
  // stores its phase and offset, so we don't need to faff about after
  // this
  for (int i = 0; i < mesh.ystart; ++i) {
    // NOTE: std::vector constructor here takes a **copy** of the
    // Array! We *must* call `Array::ensureUnique` on each element
    // before using it!
    parallel_slice_phases[i].phase_shift =
        Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);
    parallel_slice_phases[i].y_offset = i + 1;

    // Backwards parallel slices
    parallel_slice_phases[mesh.ystart + i].phase_shift =
        Tensor<dcomplex>(mesh.LocalNx, mesh.LocalNy, nmodes);
    parallel_slice_phases[mesh.ystart + i].y_offset = -(i + 1);
  }

  // Parallel slice phases -- note we don't shift in the boundaries/guards
  for (auto& slice : parallel_slice_phases) {
    BOUT_FOR(i, mesh.getRegion2D("RGN_NOY")) {

      int ix = i.x();
      int iy = i.y();
      BoutReal slice_shift = zShift[i] - zShift[i.yp(slice.y_offset)];

      for (int jz = 0; jz < nmodes; jz++) {
        // wave number is 1/[rad]
        BoutReal kwave = jz * 2.0 * PI / zlength;

        slice.phase_shift(ix, iy, jz) =
            dcomplex(cos(kwave * slice_shift), -sin(kwave * slice_shift));
      }
    }
  }
}

/*!
 * Shift the field so that X-Z is not orthogonal,
 * and Y is then field aligned.
 */
Field3D ShiftedMetric::toFieldAligned(const Field3D& f, const std::string& region) {
  ASSERT2(f.getDirectionY() == YDirectionType::Standard);
  return shiftZ(f, toAlignedPhs, YDirectionType::Aligned, region);
}
FieldPerp ShiftedMetric::toFieldAligned(const FieldPerp& f, const std::string& region) {
  ASSERT2(f.getDirectionY() == YDirectionType::Standard);
  // In principle, other regions are possible, but not yet implemented
  ASSERT2(region == "RGN_NOX");
  return shiftZ(f, toAlignedPhs, YDirectionType::Aligned, region);
}

/*!
 * Shift back, so that X-Z is orthogonal,
 * but Y is not field aligned.
 */
Field3D ShiftedMetric::fromFieldAligned(const Field3D& f, const std::string& region) {
  ASSERT2(f.getDirectionY() == YDirectionType::Aligned);
  return shiftZ(f, fromAlignedPhs, YDirectionType::Standard, region);
}
FieldPerp ShiftedMetric::fromFieldAligned(const FieldPerp& f, const std::string& region) {
  ASSERT2(f.getDirectionY() == YDirectionType::Aligned);
  // In principle, other regions are possible, but not yet implemented
  ASSERT2(region == "RGN_NOX");
  return shiftZ(f, fromAlignedPhs, YDirectionType::Standard, region);
}

Field3D ShiftedMetric::shiftZ(const Field3D& f, const Tensor<dcomplex>& phs,
                              const YDirectionType y_direction_out,
                              const std::string& region) const {
  ASSERT1(f.getMesh() == &mesh);
  ASSERT1(f.getLocation() == location);

  if (mesh.LocalNz == 1) {
    // Shifting does not change the array values
    Field3D result = copy(f).setDirectionY(y_direction_out);
    return result;
  }

  Field3D result{emptyFrom(f).setDirectionY(y_direction_out)};

  BOUT_FOR(i, mesh.getRegion2D(toString(region))) {
    shiftZ(&f(i, 0), &phs(i.x(), i.y(), 0), &result(i, 0));
  }

  return result;
}

FieldPerp ShiftedMetric::shiftZ(const FieldPerp& f, const Tensor<dcomplex>& phs,
                                const YDirectionType y_direction_out,
                                const std::string& UNUSED(region)) const {
  ASSERT1(f.getMesh() == &mesh);
  ASSERT1(f.getLocation() == location);

  if (mesh.LocalNz == 1) {
    // Shifting does not change the array values
    FieldPerp result = copy(f).setDirectionY(y_direction_out);
    return result;
  }

  FieldPerp result{emptyFrom(f).setDirectionY(y_direction_out)};

  int y = f.getIndex();
  // Note that this loop is essentially hardcoded to be RGN_NOX
  for (int i = mesh.xstart; i <= mesh.xend; ++i) {
    shiftZ(&f(i, 0), &phs(i, y, 0), &result(i, 0));
  }

  return result;
}

void ShiftedMetric::shiftZ(const BoutReal* in, const dcomplex* phs, BoutReal* out) const {
#if BOUT_HAS_UMPIRE
  // TODO: This static keyword is a hotfix and should be removed in
  //      future iterations. It is here because otherwise many allocations
  //      lead to very poor performance
  static Array<dcomplex> cmplx(nmodes);
#warning static hotfix used in ShiftedMetric::shiftZ. Not thread-safe.
#else
  Array<dcomplex> cmplx(nmodes);
#endif

  // Take forward FFT
  rfft(in, mesh.LocalNz, &cmplx[0]);

  // Following is an algorithm approach to write a = a*b where a and b are
  // vectors of dcomplex.
  //  std::transform(cmplxOneOff.begin(),cmplxOneOff.end(), ptr.begin(),
  //		 cmplxOneOff.begin(), std::multiplies<dcomplex>());

  for (int jz = 1; jz < nmodes; jz++) {
    cmplx[jz] *= phs[jz];
  }

  irfft(&cmplx[0], mesh.LocalNz, out); // Reverse FFT
}

#if BOUT_HAS_CUDA
// Bit-reversal
__device__ inline unsigned int bit_reverse(unsigned int x, unsigned int log2n) {
  unsigned int result = 0;
#pragma unroll
  for (unsigned int i = 0; i < log2n; i++) {
    result = (result << 1) | (x & 1);
    x >>= 1;
  }
  return result;
}

// Block-level cooperative FFT
// Multiple threads cooperate on each FFT using shared memory
template <int NZ, int FFTS_PER_BLOCK = 4>
__global__ void fft_block_cooperative(const BoutReal** __restrict__ in,
                                      BoutReal** __restrict__ out,
                                      const double2** __restrict__ blocks_phs,
                                      const int nbatches, const int nblocks) {

  constexpr int LOG2_NZ = __builtin_ctz(NZ);
  constexpr double INV_NZ = 1.0 / (double)NZ;
  constexpr int NMODES = (NZ / 2) + 1;

  // Shared memory for FFTS_PER_BLOCK FFTs
  // Each FFT needs NZ complex values
  __shared__ double2 shared_fft[FFTS_PER_BLOCK][NZ];

  // Select twiddles based on size
  const double2* twiddles;
  if constexpr (NZ == 16) {
    twiddles = c_twiddle_16;
  } else if constexpr (NZ == 64) {
    twiddles = c_twiddle_64;
  } else if constexpr (NZ == 128) {
    twiddles = c_twiddle_128;
  } else if constexpr (NZ == 256) {
    twiddles = c_twiddle_256;
  } else if constexpr (NZ == 512) {
    twiddles = c_twiddle_512;
  } else {
    static_assert(NZ == 16 || NZ == 64 || NZ == 128 || NZ == 256 || NZ == 512,
                  "Unsupported NZ");
  }

  // Each block processes FFTS_PER_BLOCK FFTs
  const int fft_id_in_block =
      threadIdx.y; // Which FFT this thread works on (0 to FFTS_PER_BLOCK-1)
  const int global_fft_id = blockIdx.x * FFTS_PER_BLOCK + fft_id_in_block;

  if (global_fft_id >= nblocks * nbatches)
    return;

  const int block = global_fft_id / nbatches;
  const int batch = global_fft_id % nbatches;

  const double* __restrict__ in_line = in[block] + batch * NZ;
  double* __restrict__ out_line = out[block] + batch * NZ;
  const double2* __restrict__ phs = blocks_phs[block];

  // Thread ID within the FFT computation
  const int tid = threadIdx.x;
  const int threads_per_fft = blockDim.x; // All threads in x-dimension work on same FFT

  // ===== LOAD INPUT WITH BIT-REVERSAL =====
  // Each thread loads some elements (strided)
  for (int i = tid; i < NZ; i += threads_per_fft) {
    const unsigned int rev_i = bit_reverse(i, LOG2_NZ);
    shared_fft[fft_id_in_block][rev_i].x = in_line[i];
    shared_fft[fft_id_in_block][rev_i].y = 0.0;
  }
  __syncthreads();

  // ===== FORWARD FFT: Cooley-Tukey DIT in Shared Memory =====
  for (int stage = 0; stage < LOG2_NZ; ++stage) {
    const int m = 1 << (stage + 1);
    const int m_half = m >> 1;

    // Each thread processes multiple butterflies
    for (int k = tid; k < NZ / 2; k += threads_per_fft) {
      const int butterfly_group = k / m_half;
      const int j = k % m_half;
      const int idx_top = butterfly_group * m + j;
      const int idx_bot = idx_top + m_half;

      // Twiddle factor
      const int twiddle_k = (j * NZ) / m;
      const double wr = twiddles[twiddle_k].x;
      const double wi = twiddles[twiddle_k].y;

      // Load from shared memory
      const double top_r = shared_fft[fft_id_in_block][idx_top].x;
      const double top_i = shared_fft[fft_id_in_block][idx_top].y;
      const double bot_r = shared_fft[fft_id_in_block][idx_bot].x;
      const double bot_i = shared_fft[fft_id_in_block][idx_bot].y;

      // Butterfly: t = W * bottom
      const double t_r = wr * bot_r - wi * bot_i;
      const double t_i = wr * bot_i + wi * bot_r;

      // Write back
      shared_fft[fft_id_in_block][idx_top].x = top_r + t_r;
      shared_fft[fft_id_in_block][idx_top].y = top_i + t_i;
      shared_fft[fft_id_in_block][idx_bot].x = top_r - t_r;
      shared_fft[fft_id_in_block][idx_bot].y = top_i - t_i;
    }
    __syncthreads();
  }

  // ===== APPLY PHASE SHIFT =====
  for (int k = tid; k < NMODES; k += threads_per_fft) {
    const double2 ph = phs[batch * NMODES + k];
    const double real = shared_fft[fft_id_in_block][k].x;
    const double imag = shared_fft[fft_id_in_block][k].y;
    shared_fft[fft_id_in_block][k].x = real * ph.x - imag * ph.y;
    shared_fft[fft_id_in_block][k].y = real * ph.y + imag * ph.x;
  }

  for (int k = tid + NMODES; k < NZ; k += threads_per_fft) {
    if (k >= NMODES) {
      const int kk = NZ - k;
      const double2 tmp = phs[batch * NMODES + kk];
      const double real = shared_fft[fft_id_in_block][k].x;
      const double imag = shared_fft[fft_id_in_block][k].y;
      shared_fft[fft_id_in_block][k].x = real * tmp.x + imag * tmp.y;
      shared_fft[fft_id_in_block][k].y = -real * tmp.y + imag * tmp.x;
    }
  }
  __syncthreads();

  // ===== INVERSE FFT: Conjugate, FFT, Conjugate =====
  // Conjugate input
  for (int i = tid; i < NZ; i += threads_per_fft) {
    shared_fft[fft_id_in_block][i].y = -shared_fft[fft_id_in_block][i].y;
  }
  __syncthreads();

  // Bit-reverse with standard swap to avoid temp array
  // This is tricky but saves memory
  for (int i = tid; i < NZ / 2; i += threads_per_fft) {
    const unsigned int rev_i = bit_reverse(i, LOG2_NZ);
    if (i < rev_i) { // Only swap once per pair
      double2 temp = shared_fft[fft_id_in_block][i];
      shared_fft[fft_id_in_block][i] = shared_fft[fft_id_in_block][rev_i];
      shared_fft[fft_id_in_block][rev_i] = temp;
    }
  }
  __syncthreads();

  // Forward FFT again (for inverse)
  for (int stage = 0; stage < LOG2_NZ; ++stage) {
    const int m = 1 << (stage + 1);
    const int m_half = m >> 1;

    for (int k = tid; k < NZ / 2; k += threads_per_fft) {
      const int butterfly_group = k / m_half;
      const int j = k % m_half;
      const int idx_top = butterfly_group * m + j;
      const int idx_bot = idx_top + m_half;

      const int twiddle_k = (j * NZ) / m;
      const double wr = twiddles[twiddle_k].x;
      const double wi = twiddles[twiddle_k].y;

      const double top_r = shared_fft[fft_id_in_block][idx_top].x;
      const double top_i = shared_fft[fft_id_in_block][idx_top].y;
      const double bot_r = shared_fft[fft_id_in_block][idx_bot].x;
      const double bot_i = shared_fft[fft_id_in_block][idx_bot].y;

      const double t_r = wr * bot_r - wi * bot_i;
      const double t_i = wr * bot_i + wi * bot_r;

      shared_fft[fft_id_in_block][idx_top].x = top_r + t_r;
      shared_fft[fft_id_in_block][idx_top].y = top_i + t_i;
      shared_fft[fft_id_in_block][idx_bot].x = top_r - t_r;
      shared_fft[fft_id_in_block][idx_bot].y = top_i - t_i;
    }
    __syncthreads();
  }

  // Store output (conjugate and normalize)
  for (int i = tid; i < NZ; i += threads_per_fft) {
    out_line[i] = shared_fft[fft_id_in_block][i].x * INV_NZ;
  }
}

// Launcher for block-level cooperative FFT
static void shiftZ_block_fft(const int Nz, const BoutReal** in, BoutReal** out,
                             const double2** phs, int nblocks, int nbatches,
                             cudaStream_t stream = 0) {
  if ((Nz & (Nz - 1)) != 0) {
    fprintf(stderr, "Error: Nz=%d must be power of 2\n", Nz);
    return;
  }

  const int total_ffts = nblocks * nbatches;

  if (Nz == 16) {
    constexpr int FFTS_PER_BLOCK = 16;
    constexpr int THREADS_PER_FFT = 16;

    dim3 block(THREADS_PER_FFT, FFTS_PER_BLOCK);
    dim3 grid((total_ffts + FFTS_PER_BLOCK - 1) / FFTS_PER_BLOCK);

    fft_block_cooperative<16, FFTS_PER_BLOCK>
        <<<grid, block, 0, stream>>>(in, out, phs, nbatches, nblocks);
  } else if (Nz == 64) {
    constexpr int FFTS_PER_BLOCK = 4;
    constexpr int THREADS_PER_FFT = 64;

    dim3 block(THREADS_PER_FFT, FFTS_PER_BLOCK);
    dim3 grid((total_ffts + FFTS_PER_BLOCK - 1) / FFTS_PER_BLOCK);

    fft_block_cooperative<64, FFTS_PER_BLOCK>
        <<<grid, block, 0, stream>>>(in, out, phs, nbatches, nblocks);

  } else if (Nz == 128) {
    constexpr int FFTS_PER_BLOCK = 2;
    constexpr int THREADS_PER_FFT = 128;

    dim3 block(THREADS_PER_FFT, FFTS_PER_BLOCK);
    dim3 grid((total_ffts + FFTS_PER_BLOCK - 1) / FFTS_PER_BLOCK);

    fft_block_cooperative<128, FFTS_PER_BLOCK>
        <<<grid, block, 0, stream>>>(in, out, phs, nbatches, nblocks);

  } else if (Nz == 256) {
    constexpr int FFTS_PER_BLOCK = 1;
    constexpr int THREADS_PER_FFT = 256;

    dim3 block(THREADS_PER_FFT, FFTS_PER_BLOCK);
    dim3 grid(total_ffts);

    fft_block_cooperative<256, FFTS_PER_BLOCK>
        <<<grid, block, 0, stream>>>(in, out, phs, nbatches, nblocks);

  } else if (Nz == 512) {
    constexpr int FFTS_PER_BLOCK = 1;
    constexpr int THREADS_PER_FFT = 512;

    dim3 block(THREADS_PER_FFT, FFTS_PER_BLOCK);
    dim3 grid(total_ffts);

    fft_block_cooperative<512, FFTS_PER_BLOCK>
        <<<grid, block, 0, stream>>>(in, out, phs, nbatches, nblocks);
  } else {
    throw std::runtime_error("Unsupported Nz " + std::to_string(Nz) + " for block FFT");
  }

  cudaError_t err = cudaGetLastError();
  if (err != cudaSuccess) {
    throw std::runtime_error(std::string("Block FFT failed: ") + cudaGetErrorString(err));
  }
}
#endif

void ShiftedMetric::calcParallelSlices(Field3D& f) {
  if (f.getDirectionY() == YDirectionType::Aligned) {
    // Cannot calculate parallel slices for field-aligned fields, so return without
    // setting yup or ydown
    return;
  }

  f.splitParallelSlices();

#if BOUT_HAS_CUDA
  auto& region = mesh.getRegion2D("RGN_NOY");
  static size_t nblocks = region.getBlocks().size();
  if (nblocks != region.getBlocks().size()) {
    throw BoutException("Number of blocks changed in ShiftedMetric::calcParallelSlices");
  }
  static Array<const BoutReal*> blocks_in(nblocks);
  static Array<BoutReal*> blocks_out(nblocks);
  static Array<const double2*> phs_in(nblocks);

  for (const auto& phase : parallel_slice_phases) {
    auto& f_slice = f.ynext(phase.y_offset);
    f_slice.allocate();

    static struct StreamRAII {
      cudaStream_t stream = 0;
      StreamRAII() {
        if (cudaStreamCreate(&stream) != cudaSuccess) {
          throw BoutException("Failed to create CUDA stream");
        }
      }

      cudaStream_t get() const { return stream; }

      void synchronize() const { cudaStreamSynchronize(stream); }

      ~StreamRAII() { cudaStreamDestroy(stream); }
    } stream;
    size_t block_idx = 0;
    int nbatches =
        region.getBlocks().cbegin()->second.ind - region.getBlocks().cbegin()->first.ind;

    for (auto block = region.getBlocks().cbegin(), end = region.getBlocks().cend();
         block < end; ++block) {
      auto idx_s = block->first;
      auto idx_e = block->second;
      int inner_nbatches = idx_e.ind - idx_s.ind;
      if (inner_nbatches != nbatches) {
        throw BoutException(
            "Non-uniform number of batches in ShiftedMetric::calcParallelSlices");
      }
      const int ix = idx_s.x();
      const int iy = idx_s.y();
      const int iy_offset = iy + phase.y_offset;

      blocks_in[block_idx] = &f(ix, iy_offset, 0);
      blocks_out[block_idx] = &f_slice(ix, iy_offset, 0);
      phs_in[block_idx] = reinterpret_cast<const double2*>(&phase.phase_shift(ix, iy, 0));

      block_idx++;
    }

    shiftZ_block_fft(mesh.LocalNz, &blocks_in[0], &blocks_out[0], &phs_in[0], nblocks,
                     nbatches, stream.get());

    stream.synchronize();
  }
#else
  for (const auto& phase : parallel_slice_phases) {
    auto& f_slice = f.ynext(phase.y_offset);
    f_slice.allocate();

    BOUT_FOR(i, mesh.getRegion2D("RGN_NOY")) {
      const int ix = i.x();
      const int iy = i.y();
      const int iy_offset = iy + phase.y_offset;
      shiftZ(&(f(ix, iy_offset, 0)), &(phase.phase_shift(ix, iy, 0)),
             &(f_slice(ix, iy_offset, 0)));
    }
  }
#endif
}

std::vector<Field3D>
ShiftedMetric::shiftZ(const Field3D& f,
                      const std::vector<ParallelSlicePhase>& phases) const {
  ASSERT1(f.getMesh() == &mesh);
  ASSERT1(f.getLocation() == location);
  ASSERT1(f.getDirectionY() == YDirectionType::Standard);

  const int nmodes = mesh.LocalNz / 2 + 1;

  // FFT in Z of input field at each (x, y) point
  Matrix<Array<dcomplex>> f_fft(mesh.LocalNx, mesh.LocalNy);
  f_fft = Array<dcomplex>(nmodes);

  BOUT_FOR(i, mesh.getRegion2D("RGN_ALL")) {
    int ix = i.x();
    int iy = i.y();
    f_fft(ix, iy).ensureUnique();
    rfft(&f(i, 0), mesh.LocalNz, f_fft(ix, iy).begin());
  }

  std::vector<Field3D> results{};

  for (const auto& phase : phases) {
    auto& current_result = results.emplace_back(&mesh);
    current_result.allocate();
    current_result.setLocation(f.getLocation());

    BOUT_FOR(i, mesh.getRegion2D("RGN_NOY")) {
      // Deep copy the FFT'd field
      int ix = i.x();
      int iy = i.y();

      Array<dcomplex> shifted_temp(f_fft(ix, iy + phase.y_offset));
      shifted_temp.ensureUnique();

      for (int jz = 1; jz < nmodes; ++jz) {
        shifted_temp[jz] *= phase.phase_shift(ix, iy, jz);
      }

      irfft(shifted_temp.begin(), mesh.LocalNz, &current_result(i.yp(phase.y_offset), 0));
    }
  }

  return results;
}
