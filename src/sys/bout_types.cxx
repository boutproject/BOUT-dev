#include <boutexception.hxx>
#include <msg_stack.hxx>
#include <bout_types.hxx>
#include <boutexception.hxx>
#include <map>

template <typename T>
const std::string& safeAt(const std::map<T, std::string>& mymap, T t) {
  AUTO_TRACE();
  auto found = mymap.find(t);
  if (found == mymap.end()) {
    throw BoutException("Did not find enum %d", static_cast<int>(t));
  }
  return found->second;
}

const std::string& CELL_LOC_STRING(CELL_LOC location) {
  AUTO_TRACE();
  const static std::map<CELL_LOC, std::string> CELL_LOCtoString = {
      ENUMSTR(CELL_DEFAULT), ENUMSTR(CELL_CENTRE), ENUMSTR(CELL_XLOW),
      ENUMSTR(CELL_YLOW),    ENUMSTR(CELL_ZLOW),   ENUMSTR(CELL_VSHIFT)};

  return safeAt(CELL_LOCtoString, location);
}

const std::string& DIFF_METHOD_STRING(DIFF_METHOD location) {
  AUTO_TRACE();
  const static std::map<DIFF_METHOD, std::string> DIFF_METHODtoString = {
      {DIFF_DEFAULT, "DEFAULT"}, {DIFF_U1, "U1"},   {DIFF_U2, "U2"},      {DIFF_U3, "U3"},
      {DIFF_C2, "C2"},           {DIFF_C4, "C4"},   {DIFF_S2, "S2"},      {DIFF_W2, "W2"},
      {DIFF_W3, "W3"},           {DIFF_FFT, "FFT"}, {DIFF_SPLIT, "SPLIT"}};

  return safeAt(DIFF_METHODtoString, location);
}

const std::string& REGION_STRING(REGION region) {
  AUTO_TRACE();
  const static std::map<REGION, std::string> REGIONtoString = {
      ENUMSTR(RGN_ALL), ENUMSTR(RGN_NOBNDRY), ENUMSTR(RGN_NOX), ENUMSTR(RGN_NOY),
      ENUMSTR(RGN_NOZ)};
  return safeAt(REGIONtoString, region);
}

const std::string& DIRECTION_STRING(DIRECTION direction) {
  AUTO_TRACE();
  const static std::map<DIRECTION, std::string> DIRECTIONtoString = {
      {DIRECTION::X, "X"},
      {DIRECTION::Y, "Y"},
      {DIRECTION::Z, "Z"},
      {DIRECTION::YAligned, "Y - field aligned"},
      {DIRECTION::YOrthogonal, "Y - orthogonal"}};

  return safeAt(DIRECTIONtoString, direction);
}

void swap(DirectionTypes& first, DirectionTypes& second) {
  DirectionTypes temp = first;
  first = second;
  second = temp;
}

bool compatibleDirections(const DirectionTypes& d1, const DirectionTypes& d2) {
  if (d1.y == d2.y && d1.z == d2.z) {
    // direction types are the same, most common case, return immediately
    return true;
  }

  if (d2.z == ZDirectionType::Average && d2.y == YDirectionType::Standard
      && (d1.y == YDirectionType::Standard || d1.y == YDirectionType::Aligned)) {
    // If a field has ZDirectionType::Average, then it's compatible with
    // YDirectionType::Aligned as well as YDirectionType::Standard
    return true;
  }

  if (d1.z == ZDirectionType::Average && d1.y == YDirectionType::Standard
      && (d2.y == YDirectionType::Standard || d2.y == YDirectionType::Aligned)) {
    // If a field has ZDirectionType::Average, then it's compatible with
    // YDirectionType::Aligned as well as YDirectionType::Standard
    return true;
  }

  // No compatible cases found
  return false;
}

const std::string& STAGGER_STRING(STAGGER stagger) {
  AUTO_TRACE();
  const static std::map<STAGGER, std::string> STAGGERtoString = {
      {STAGGER::None, "No staggering"},
      {STAGGER::C2L, "Centre to Low"},
      {STAGGER::L2C, "Low to Centre"}};

  return safeAt(STAGGERtoString, stagger);
}

const std::string& DERIV_STRING(DERIV deriv) {
  AUTO_TRACE();
  const static std::map<DERIV, std::string> DERIVtoString = {
      {DERIV::Standard, "Standard"},
      {DERIV::StandardSecond, "Standard -- second order"},
      {DERIV::StandardFourth, "Standard -- fourth order"},
      {DERIV::Upwind, "Upwind"},
      {DERIV::Flux, "Flux"}};

  return safeAt(DERIVtoString, deriv);
}
