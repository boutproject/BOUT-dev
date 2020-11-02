/*************************************************************************
 *
 * This file is part of the SAMRAI distribution.  For full copyright
 * information, see COPYRIGHT and LICENSE.
 *
 * Copyright:     (c) 1997-2018 Lawrence Livermore National Security, LLC
 * Description:   Class to record statistics during program execution.
 *
 ************************************************************************/

#ifndef included_ForAll
#define included_ForAll

#include "BOUT_config.h"

#if defined(HAVE_RAJA)
#include "RAJA/RAJA.hpp"

#include "bout/sys/ExecutionPolicy.h"

#include <type_traits>
#include <tuple>
#include <cstdlib>  // for std::size_t

using Index_type = std::ptrdiff_t;

namespace bout
{

/*
 * Macros defined for host-device compilation.
 */
#define BOUT_INLINE inline

#if defined(HAVE_CUDA) && defined(__CUDACC__)
#define BOUT_HOST_DEVICE __host__ __device__
#else
#define BOUT_HOST_DEVICE
#endif

/*
parallel_for_all()  version that picks parallel policy (GPU if ENABLE_CUDA=ON)
for_all<policy>()      version that takes a custom RAJA policy (could be either host or device)
*/

namespace detail
{

template <typename T>
struct function_traits : function_traits<decltype(&T::operator())> {
};

// free function
template <typename R, typename... Args>
struct function_traits<R(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// pointer to function
template <typename R, typename... Args>
struct function_traits<R (*)(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// member function
template <typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...)> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

// const member function
template <typename T, typename R, typename... Args>
struct function_traits<R (T::*)(Args...) const> {
   using result_type = R;
   using argument_types = std::tuple<Args...>;
   enum { argument_count = sizeof...(Args) };
};

inline RAJA::RangeSegment make_range(const Index_type& ifirst, const Index_type& ilast, std::size_t index)
{
   return RAJA::RangeSegment(ifirst, ilast + 1);
}

template <int ArgumentCount>
struct for_all {
};

template <>
struct for_all<1> {
   template <typename Policy, typename LoopBody,
             typename std::enable_if<std::is_base_of<bout::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const Index_type& ifirst, const Index_type& ilast, LoopBody body)
   {
      RAJA::kernel<typename bout::detail::policy_traits<Policy>::Policy1d>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0)),
          body);
   }

   template <typename Policy, typename LoopBody,
             typename std::enable_if<!std::is_base_of<bout::policy::base, Policy>::value, int>::type = 0>
   inline static void eval(const Index_type& ifirst, const Index_type& ilast, LoopBody body)
   {
      RAJA::kernel<Policy>(
          RAJA::make_tuple(make_range(ifirst, ilast, 0)),
          body);
   }
};

// todo setup 2 and 3d variants 

}  // namespace detail

// does NOT include end
template <typename Policy, typename LoopBody,
          typename std::enable_if<std::is_base_of<bout::policy::base, Policy>::value, int>::type = 0>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<typename bout::detail::policy_traits<Policy>::Policy>(RAJA::RangeSegment(begin, end), body);
}

template <typename Policy, typename LoopBody,
          typename std::enable_if<!std::is_base_of<bout::policy::base, Policy>::value, int>::type = 0>
inline void for_all(int begin, int end, LoopBody body)
{
   RAJA::forall<Policy>(RAJA::RangeSegment(begin, end), body);
}

// does NOT include end
template <typename LoopBody>
inline void parallel_for_all(int begin, int end, LoopBody body)
{
   for_all<bout::policy::parallel>(begin, end, body);
}


}  // namespace bout

#endif

#endif  // included_ForAll
