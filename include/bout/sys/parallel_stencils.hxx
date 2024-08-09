#pragma once

namespace parallel_stencil {
// generated by src/mesh/parallel_boundary_stencil.cxx.py
inline BoutReal pow(BoutReal val, int exp) {
  // constexpr int expval = exp;
  // static_assert(expval == 2 or expval == 3, "This pow is only for exponent 2 or 3");
  if (exp == 2) {
    return val * val;
  }
  ASSERT3(exp == 3);
  return val * val * val;
}
inline BoutReal dirichlet_o1(BoutReal UNUSED(spacing0), BoutReal value0) {
  return value0;
}
inline BoutReal dirichlet_o2(BoutReal spacing0, BoutReal value0, BoutReal spacing1,
                             BoutReal value1) {
  return (spacing0 * value1 - spacing1 * value0) / (spacing0 - spacing1);
}
inline BoutReal neumann_o2(BoutReal UNUSED(spacing0), BoutReal value0, BoutReal spacing1,
                           BoutReal value1) {
  return -spacing1 * value0 + value1;
}
inline BoutReal dirichlet_o3(BoutReal spacing0, BoutReal value0, BoutReal spacing1,
                             BoutReal value1, BoutReal spacing2, BoutReal value2) {
  return (pow(spacing0, 2) * spacing1 * value2 - pow(spacing0, 2) * spacing2 * value1
          - spacing0 * pow(spacing1, 2) * value2 + spacing0 * pow(spacing2, 2) * value1
          + pow(spacing1, 2) * spacing2 * value0 - spacing1 * pow(spacing2, 2) * value0)
         / ((spacing0 - spacing1) * (spacing0 - spacing2) * (spacing1 - spacing2));
}
inline BoutReal neumann_o3(BoutReal spacing0, BoutReal value0, BoutReal spacing1,
                           BoutReal value1, BoutReal spacing2, BoutReal value2) {
  return (2 * spacing0 * spacing1 * value2 - 2 * spacing0 * spacing2 * value1
          + pow(spacing1, 2) * spacing2 * value0 - pow(spacing1, 2) * value2
          - spacing1 * pow(spacing2, 2) * value0 + pow(spacing2, 2) * value1)
         / ((spacing1 - spacing2) * (2 * spacing0 - spacing1 - spacing2));
}
} // namespace parallel_stencil