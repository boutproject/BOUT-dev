#include <bout_types.hxx>
#include <boutexception.hxx>
#include <map>
#include <msg_stack.hxx>

namespace {
template <typename T>
const std::string& safeAt(const std::map<T, std::string>& mymap, T t) {
  AUTO_TRACE();
  auto found = mymap.find(t);
  if (found == mymap.end()) {
    throw BoutException("Did not find enum %d", static_cast<int>(t));
  }
  return found->second;
}

template <typename T>
const T& safeAt(const std::map<std::string, T>& mymap, const std::string& s) {
  AUTO_TRACE();
  auto found = mymap.find(s);
  if (found == mymap.end()) {
    throw BoutException("Did not find enum %s", s.c_str());
  }
  return found->second;
}
}

std::string toString(CELL_LOC location) {
  AUTO_TRACE();
  const static std::map<CELL_LOC, std::string> CELL_LOCtoString = {
      ENUMSTR(CELL_DEFAULT), ENUMSTR(CELL_CENTRE), ENUMSTR(CELL_XLOW),
      ENUMSTR(CELL_YLOW),    ENUMSTR(CELL_ZLOW),   ENUMSTR(CELL_VSHIFT)};

  return safeAt(CELL_LOCtoString, location);
}

CELL_LOC CELL_LOCFromString(const std::string& location_string) {
  AUTO_TRACE();
  const static std::map<std::string, CELL_LOC> stringtoCELL_LOC = {
    STRENUM(CELL_DEFAULT), STRENUM(CELL_CENTRE), STRENUM(CELL_XLOW),
    STRENUM(CELL_YLOW),    STRENUM(CELL_ZLOW),   STRENUM(CELL_VSHIFT)};

  return safeAt(stringtoCELL_LOC, location_string);
}

std::string toString(DIFF_METHOD location) {
  AUTO_TRACE();
  const static std::map<DIFF_METHOD, std::string> DIFF_METHODtoString = {
      {DIFF_DEFAULT, "DEFAULT"}, {DIFF_U1, "U1"},   {DIFF_U2, "U2"},      {DIFF_U3, "U3"},
      {DIFF_C2, "C2"},           {DIFF_C4, "C4"},   {DIFF_S2, "S2"},      {DIFF_W2, "W2"},
      {DIFF_W3, "W3"},           {DIFF_FFT, "FFT"}, {DIFF_SPLIT, "SPLIT"}};

  return safeAt(DIFF_METHODtoString, location);
}

std::string toString(REGION region) {
  AUTO_TRACE();
  const static std::map<REGION, std::string> REGIONtoString = {
      ENUMSTR(RGN_ALL), ENUMSTR(RGN_NOBNDRY), ENUMSTR(RGN_NOX), ENUMSTR(RGN_NOY),
      ENUMSTR(RGN_NOZ)};
  return safeAt(REGIONtoString, region);
}

std::string toString(DIRECTION direction) {
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

bool areDirectionsCompatible(const DirectionTypes& d1, const DirectionTypes& d2) {
  if (d1.y == d2.y && d1.z == d2.z) {
    // direction types are the same, most common case, return immediately
    return true;
  }

  if (d2.z == ZDirectionType::Average && d2.y == YDirectionType::Standard
      && (d1.y == YDirectionType::Standard || d1.y == YDirectionType::Aligned)
      && d1.z == ZDirectionType::Standard) {
    // If d2 has ZDirectionType::Average, then it's compatible with d1 having
    // YDirectionType::Aligned as well as YDirectionType::Standard.  If d1 has
    // YDirectionType::Aligned, it should always have ZDirectionType::Standard,
    // and if d1 has ZDirectionType::Average it must have
    // YDirectionType::Standard and have been caught in the first condition
    // where d1 and d2 are identical, so only allow
    // 'd1.z == ZDirectionType::Standard' here.
    return true;
  }

  if (d1.z == ZDirectionType::Average && d1.y == YDirectionType::Standard
      && (d2.y == YDirectionType::Standard || d2.y == YDirectionType::Aligned)
      && d2.z == ZDirectionType::Standard) {
    // If d1 has ZDirectionType::Average, then it's compatible with d2 having
    // YDirectionType::Aligned as well as YDirectionType::Standard.  If d2 has
    // YDirectionType::Aligned, it should always have ZDirectionType::Standard,
    // and if d2 has ZDirectionType::Average it must have
    // YDirectionType::Standard and have been caught in the first condition
    // where d1 and d2 are identical, so only allow
    // 'd2.z == ZDirectionType::Standard' here.
    return true;
  }

  // No compatible cases found
  return false;
}

std::string toString(STAGGER stagger) {
  AUTO_TRACE();
  const static std::map<STAGGER, std::string> STAGGERtoString = {
      {STAGGER::None, "No staggering"},
      {STAGGER::C2L, "Centre to Low"},
      {STAGGER::L2C, "Low to Centre"}};

  return safeAt(STAGGERtoString, stagger);
}

std::string toString(DERIV deriv) {
  AUTO_TRACE();
  const static std::map<DERIV, std::string> DERIVtoString = {
      {DERIV::Standard, "Standard"},
      {DERIV::StandardSecond, "Standard -- second order"},
      {DERIV::StandardFourth, "Standard -- fourth order"},
      {DERIV::Upwind, "Upwind"},
      {DERIV::Flux, "Flux"}};

  return safeAt(DERIVtoString, deriv);
}

std::string toString(YDirectionType d) {
  AUTO_TRACE();
  const static std::map<YDirectionType, std::string> YDirectionTypeToString = {
      {YDirectionType::Standard, "Standard"},
      {YDirectionType::Aligned, "Aligned"}};

  return safeAt(YDirectionTypeToString, d);
}

YDirectionType YDirectionTypeFromString(const std::string& y_direction_string) {
  AUTO_TRACE();
  const static std::map<std::string, YDirectionType> stringToYDirectionType = {
    {"Standard", YDirectionType::Standard},
    {"Aligned", YDirectionType::Aligned}};

  return safeAt(stringToYDirectionType, y_direction_string);
}

std::string toString(ZDirectionType d) {
  AUTO_TRACE();
  const static std::map<ZDirectionType, std::string> ZDirectionTypeToString = {
      {ZDirectionType::Standard, "Standard"},
      {ZDirectionType::Average, "Average"}};

  return safeAt(ZDirectionTypeToString, d);
}

ZDirectionType ZDirectionTypeFromString(const std::string& z_direction_string) {
  AUTO_TRACE();
  const static std::map<std::string, ZDirectionType> stringToZDirectionType = {
    {"Standard", ZDirectionType::Standard},
    {"Average", ZDirectionType::Average}};

  return safeAt(stringToZDirectionType, z_direction_string);
}
