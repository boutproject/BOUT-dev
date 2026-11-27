#pragma once

namespace bout::details {
/// Modulo that ensures the result is always positive, even when x is negative
constexpr int mod(int x, int y) {
  const int rem = x % y;
  return rem < 0 ? rem + y : rem;
}
} // namespace bout::details
