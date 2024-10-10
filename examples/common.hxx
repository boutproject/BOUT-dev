
#ifndef BOUT_COMMON_HXX
#define BOUT_COMMON_HXX

#include "bout/bout.hxx"

inline Coordinates* tokamak_coordinates(Coordinates* coord, const FieldMetric& Rxy,
                                        const FieldMetric& Bpxy, const FieldMetric& hthe,
                                        const FieldMetric& I, const FieldMetric& B,
                                        const FieldMetric& Btxy) {

  const auto g11 = SQ(Rxy * Bpxy);
  const auto g22 = 1.0 / SQ(hthe);
  const auto g33 = SQ(I) * g11 + SQ(B) / g11;
  const auto g12 = 0.0;
  const auto g13 = -I * g11;
  const auto g23 = -sbp * Btxy / (hthe * Bpxy * Rxy);

  const auto g_11 = 1.0 / g11 + SQ(I * Rxy);
  const auto g_22 = SQ(B * hthe / Bpxy);
  const auto g_33 = Rxy * Rxy;
  const auto g_12 = sbp * Btxy * hthe * I * Rxy / Bpxy;
  const auto g_13 = I * Rxy * Rxy;
  const auto g_23 = Btxy * hthe * Rxy / Bpxy;

  coord->setMetricTensor(ContravariantMetricTensor(g11, g22, g33, g12, g13, g23),
                         CovariantMetricTensor(g_11, g_22, g_33, g_12, g_13, g_23));

  coord->setJ(hthe / Bpxy);
  coord->setBxy(B);

  return coord;
}


#endif //BOUT_COMMON_HXX
