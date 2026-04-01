#pragma once

#include "bout/build_defines.hxx"

#if BOUT_HAS_CALIPER
#include <caliper/cali-manager.h>
#include <caliper/cali.h>

#else

#define CALI_CXX_MARK_FUNCTION
#define CALI_MARK_BEGIN(...)
#define CALI_MARK_END(...)

#endif
