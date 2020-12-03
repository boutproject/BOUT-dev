/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_Collectives
#define included_Collectives

#if defined(BOUT_HAS_RAJA)

#include "RAJA/RAJA.hpp"

#include "bout/sys/ExecutionPolicy.h"

namespace bout {

template <typename policy>
inline void
synchronize() {}

#if defined(BOUT_USE_CUDA)
template<>
inline void
synchronize<policy::parallel>()
{
   RAJA::synchronize<RAJA::cuda_synchronize>();
}
#endif

inline void
parallel_synchronize() { synchronize<policy::parallel>(); }


// Reductions see https://raja.readthedocs.io/en/master/feature/reduction.html for these options
enum class Reduction {
   Sum,
   Min,
   Max,
   MinLoc,
   MaxLoc
};

template<typename Policy, Reduction R, typename TYPE = double>
struct reduction_variable;

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Sum, TYPE> {
   using type = RAJA::ReduceSum<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Min, TYPE> {
   using type = RAJA::ReduceMin<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::Max, TYPE> {
   using type = RAJA::ReduceMax<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::MinLoc, TYPE> {
   using type = RAJA::ReduceMinLoc<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<typename Policy, typename TYPE>
struct reduction_variable<Policy, Reduction::MaxLoc, TYPE> {
   using type = RAJA::ReduceMaxLoc<typename detail::policy_traits<Policy>::ReductionPolicy, TYPE>;
};

template<Reduction R, typename TYPE = double>
using parallel_reduction_variable = reduction_variable<tbox::policy::parallel, R, TYPE>;

template<typename Policy, Reduction R, typename TYPE = double>
using reduction_variable_t = typename reduction_variable<Policy, R, TYPE>::type;

template<Reduction R, typename TYPE = double>
using parallel_reduction_variable_t = typename parallel_reduction_variable<R, TYPE>::type;

} // namespace bout
#endif
#endif
