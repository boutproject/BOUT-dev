/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_ExecutionPolicy
#define included_ExecutionPolicy

#include "RAJA/RAJA.hpp"

namespace bout {

namespace policy {
struct base {};
struct sequential : base {};
struct parallel : base {};
}

namespace detail {

template <typename pol>
struct policy_traits {};

template <>
struct policy_traits<policy::sequential> {
   using Policy = RAJA::seq_exec;

   using Policy1d = RAJA::KernelPolicy<
      RAJA::statement::For<0, RAJA::seq_exec,
         RAJA::statement::Lambda<0>
      >
   >;

   using Policy2d = RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::seq_exec,
         RAJA::statement::For<0, RAJA::seq_exec,
            RAJA::statement::Lambda<0>
         >
      >
   >;

   using Policy3d = RAJA::KernelPolicy<
      RAJA::statement::For<2, RAJA::seq_exec,
         RAJA::statement::For<1, RAJA::seq_exec,
            RAJA::statement::For<0, RAJA::seq_exec,
               RAJA::statement::Lambda<0>
            >
         >
      >
   >;

   using ReductionPolicy = RAJA::seq_reduce;
};

#if defined(HAVE_CUDA)

template <>
struct policy_traits<policy::parallel> {
   using Policy = RAJA::cuda_exec_async<128>;

   using Policy1d = RAJA::KernelPolicy<
      RAJA::statement::CudaKernelAsync<
         RAJA::statement::Tile<0, RAJA::statement::tile_fixed<128>, RAJA::cuda_block_x_loop,
            RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
               RAJA::statement::Lambda<0>
            >
         >
      >
   >;

   using Policy2d = RAJA::KernelPolicy<
      RAJA::statement::CudaKernelAsync<
         RAJA::statement::Tile<1, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_y_loop,
            RAJA::statement::Tile<0, RAJA::statement::tile_fixed<16>, RAJA::cuda_block_x_loop,
               RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                  RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                     RAJA::statement::Lambda<0>
                  >
               >
            >
         >
      >
   >;

   using Policy3d = RAJA::KernelPolicy<
      RAJA::statement::CudaKernelAsync<
         RAJA::statement::Tile<2, RAJA::statement::tile_fixed<8>, RAJA::cuda_block_z_loop,
            RAJA::statement::Tile<1, RAJA::statement::tile_fixed<8>, RAJA::cuda_block_y_loop,
               RAJA::statement::Tile<0, RAJA::statement::tile_fixed<8>, RAJA::cuda_block_x_loop,
                  RAJA::statement::For<2, RAJA::cuda_thread_z_loop,
                     RAJA::statement::For<1, RAJA::cuda_thread_y_loop,
                        RAJA::statement::For<0, RAJA::cuda_thread_x_loop,
                           RAJA::statement::Lambda<0>
                        >
                     >
                  >
               >
            >
         >
      >
   >;

   using ReductionPolicy = RAJA::cuda_reduce;
};

#else

// TODO: Make this an OpenMP policy if that is defined
template <>
struct policy_traits<policy::parallel> {
   using Policy = RAJA::loop_exec;

   using Policy1d = RAJA::KernelPolicy<
      RAJA::statement::For<0, RAJA::loop_exec,
         RAJA::statement::Lambda<0>
      >
   >;

   using Policy2d = RAJA::KernelPolicy<
      RAJA::statement::For<1, RAJA::loop_exec,
         RAJA::statement::For<0, RAJA::loop_exec,
            RAJA::statement::Lambda<0>
         >
      >
   >;

   using Policy3d = RAJA::KernelPolicy<
      RAJA::statement::For<2, RAJA::loop_exec,
         RAJA::statement::For<1, RAJA::loop_exec,
            RAJA::statement::For<0, RAJA::loop_exec,
               RAJA::statement::Lambda<0>
            >
         >
      >
   >;

   using ReductionPolicy = RAJA::seq_reduce;
};

#endif // HAVE_CUDA

} // namespace detail

} // namespace bout

#endif
