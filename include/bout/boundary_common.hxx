#pragma once

#include <bout/bout_types.hxx>

#include <cstdint>

namespace bout {
namespace boundary {
/// Types of free boundary condition
/// ================== ===================================================
/// Name               Description
/// ================== ===================================================
/// ``limited``        use exponential if decreasing, otherwise Neumanm
/// ``exponential``    use exponential extrapolation
/// ``linear``         use linear extrapolation
/// ================== ===================================================
enum class BoundaryFreeExtrapolation : std::uint8_t { limited, exponential, linear };

// Limited free gradient of log of a quantity
// This ensures that the guard cell values remain positive
// while also ensuring that the quantity never increases
//
//  fm  fc | fp
//         ^ boundary
//
// exp( 2*log(fc) - log(fm) )
inline BoutReal limitFreeScale(BoutReal fm, BoutReal fc,
                               bout::boundary::BoundaryFreeExtrapolation mode) {
  if ((fm < fc) && (mode == bout::boundary::BoundaryFreeExtrapolation::limited)) {
    return fc; // Neumann rather than increasing into boundary
  }
  if (fm < 1e-10) {
    return fc; // Low / no density condition
  }

  BoutReal fp = 0;
  switch (mode) {
  case bout::boundary::BoundaryFreeExtrapolation::limited:
  case bout::boundary::BoundaryFreeExtrapolation::exponential:
    fp = fc * fc / fm; // Exponential
    break;
  case bout::boundary::BoundaryFreeExtrapolation::linear:
    fp = (2.0 * fc) - fm; // Linear
    break;
  }

#if CHECKLEVEL >= 2
  if (!std::isfinite(fp)) {
    throw BoutException("SheathBoundary limitFree: {}, {} -> {}", fm, fc, fp);
  }
#endif

  return fp;
}

} // namespace boundary
} // namespace bout
