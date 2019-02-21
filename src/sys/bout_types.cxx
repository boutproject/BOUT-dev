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
      {DIRECTION::YOrthogonal, "Y - orthogonal"},
      {DIRECTION::Special, "Special"}};

  return safeAt(DIRECTIONtoString, direction);
}

void swap(DIRECTION& first, DIRECTION& second) {
  DIRECTION temp = first;
  first = second;
  second = temp;
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

bool isXDirectionType(DIRECTION x) {
  switch (x) {
  case (DIRECTION::X):
    return true;
  default:
    return false;
  }
}

bool isYDirectionType(DIRECTION y) {
  switch (y) {
  case (DIRECTION::Y):
  case (DIRECTION::YAligned):
  case (DIRECTION::YOrthogonal):
    return true;
  default:
    return false;
  }
}

bool isZDirectionType(DIRECTION z) {
  switch (z) {
  case (DIRECTION::Z):
    return true;
  default:
    return false;
  }
}

bool compatibleDirections(DIRECTION d1, DIRECTION d2) {
  switch (d1) {
  case (DIRECTION::X): {
    switch (d2) {
    case (DIRECTION::X):
    case (DIRECTION::Special):
      return true;
    default:
      return false;
    }
  }
  case (DIRECTION::Y):
    switch (d2) {
    case (DIRECTION::Y):
    case (DIRECTION::YAligned):
    case (DIRECTION::YOrthogonal):
    case (DIRECTION::Special):
      return true;
    default:
      return false;
    }
  case (DIRECTION::YAligned):
    switch (d2) {
    case (DIRECTION::Y):
    case (DIRECTION::YAligned):
    case (DIRECTION::Special):
      return true;
    default:
      return false;
    }
  case (DIRECTION::YOrthogonal):
    switch (d2) {
    case (DIRECTION::Y):
    case (DIRECTION::YOrthogonal):
    case (DIRECTION::Special):
      return true;
    default:
      return false;
    }
  case (DIRECTION::Z):
    switch (d2) {
    case (DIRECTION::Z):
    case (DIRECTION::Special):
      return true;
    default:
      return false;
    }
  case (DIRECTION::Special):
    return true;
  default:
    // Shouldn't reach this due to checks at start but in case
    // of future changes good to handle here.
    throw BoutException("Invalid y direction value");
  }
}
