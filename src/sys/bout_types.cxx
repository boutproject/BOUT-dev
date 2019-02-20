#include <bout_types.hxx>
#include <bout/assert.hxx>
#include <map>

template <typename T>
const std::string& saveat(const std::map<T, std::string>& mymap, T t) {
  auto found = mymap.find(t);
  if (found == mymap.end()) {
    throw BoutException("Did not find enum %d", static_cast<int>(t));
  }
  return found->second;
}

const std::string& CELL_LOC_STRING(CELL_LOC location) {
  const static std::map<CELL_LOC, std::string> CELL_LOCtoString = {
      ENUMSTR(CELL_DEFAULT), ENUMSTR(CELL_CENTRE), ENUMSTR(CELL_XLOW),
      ENUMSTR(CELL_YLOW),    ENUMSTR(CELL_ZLOW),   ENUMSTR(CELL_VSHIFT)};

  return saveat(CELL_LOCtoString, location);
}

const std::string& DIFF_METHOD_STRING(DIFF_METHOD location) {
  const static std::map<DIFF_METHOD, std::string> DIFF_METHODtoString = {
      {DIFF_DEFAULT, "DEFAULT"}, {DIFF_U1, "U1"},   {DIFF_U2, "U2"},      {DIFF_U3, "U3"},
      {DIFF_C2, "C2"},           {DIFF_C4, "C4"},   {DIFF_S2, "S2"},      {DIFF_W2, "W2"},
      {DIFF_W3, "W3"},           {DIFF_FFT, "FFT"}, {DIFF_SPLIT, "SPLIT"}};

  return saveat(DIFF_METHODtoString, location);
}

const std::string& REGION_STRING(REGION region) {
  ASSERT2(region >= 0 && region <= 4);
  const static std::map<REGION, std::string> REGIONtoString = {
      ENUMSTR(RGN_ALL), ENUMSTR(RGN_NOBNDRY), ENUMSTR(RGN_NOX), ENUMSTR(RGN_NOY),
      ENUMSTR(RGN_NOZ)};
  return saveat(REGIONtoString, region);
}

const std::string& DIRECTION_STRING(DIRECTION direction) {
  const static std::map<DIRECTION, std::string> DIRECTIONtoString = {
      {DIRECTION::X, "X"},
      {DIRECTION::Y, "Y"},
      {DIRECTION::Z, "Z"},
      {DIRECTION::YAligned, "Y - field aligned"},
      {DIRECTION::YOrthogonal, "Y - orthogonal"}};

  return saveat(DIRECTIONtoString, direction);
}

const std::string& STAGGER_STRING(STAGGER stagger) {
  const static std::map<STAGGER, std::string> STAGGERtoString = {
      {STAGGER::None, "No staggering"},
      {STAGGER::C2L, "Centre to Low"},
      {STAGGER::L2C, "Low to Centre"}};

  return saveat(STAGGERtoString, stagger);
}

const std::string& DERIV_STRING(DERIV deriv) {
  const static std::map<DERIV, std::string> DERIVtoString = {
      {DERIV::Standard, "Standard"},
      {DERIV::StandardSecond, "Standard -- second order"},
      {DERIV::StandardFourth, "Standard -- fourth order"},
      {DERIV::Upwind, "Upwind"},
      {DERIV::Flux, "Flux"}};

  return saveat(DERIVtoString, deriv);
}
