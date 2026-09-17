#include "bout/build_config.hxx"
#include "gtest/gtest.h"

#if BOUT_HAS_RAJA
#include "bout/rajalib.hxx"
#include <type_traits>

#if BOUT_HAS_HIP
static_assert(std::is_same_v<EXEC_POL, RAJA::hip_exec<256>>);
#endif

TEST(RajaExecution, HostDeviceRoundTrip) {
  Array<int> values(17);
  for (int i = 0; i < values.size(); ++i) {
    values[i] = i;
  }
  int* data = &values[0];
  // EXEC_POL is synchronous: results must be visible to the host on return.
  RAJA::forall<EXEC_POL>(RAJA::RangeSegment(0, values.size()),
                         [=] RAJA_DEVICE(int i) { data[i] = 3 * data[i] + 1; });
  for (int i = 0; i < values.size(); ++i) {
    EXPECT_EQ(values[i], 3 * i + 1);
  }
}
#endif
