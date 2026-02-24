///
/// Implement {fmt} formatters for BOUT++ types
///

#ifndef OUTPUT_BOUT_TYPES_H
#define OUTPUT_BOUT_TYPES_H

#include "bout/bout_types.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx" // IWYU pragma: keep
#include "bout/region.hxx"
#include "bout/traits.hxx"

#include "fmt/base.h"
#include <fmt/format.h>

#include <type_traits>

template <IND_TYPE N>
struct fmt::formatter<SpecificInd<N>> {
  // Presentation format: 'c' - components, 'i' - index.
  char presentation = 'c';

  // Parses format specifications of the form ['c' | 'i'].
  constexpr auto parse(format_parse_context& ctx) {
    const auto* it = ctx.begin();
    const auto* end = ctx.end();
    if (it != end && (*it == 'c' || *it == 'i')) {
      presentation = *it++;
    }

    // Check if reached the end of the range:
    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    // Return an iterator past the end of the parsed range:
    return it;
  }

  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  auto format(const SpecificInd<N>& ind, FormatContext& ctx) const {
    // ctx.out() is an output iterator to write to.
    if (presentation == 'c') {
      switch (N) {
      case IND_TYPE::IND_2D:
        return format_to(ctx.out(), "({}, {})", ind.x(), ind.y());
      case IND_TYPE::IND_3D:
        return format_to(ctx.out(), "({}, {}, {})", ind.x(), ind.y(), ind.z());
      case IND_TYPE::IND_PERP:
        return format_to(ctx.out(), "({}, {})", ind.x(), ind.z());
      }
    }
    return format_to(ctx.out(), "({})", ind.ind);
  }
};

/// Formatter for Fields
template <class T>
struct fmt::formatter<T, std::enable_if_t<bout::utils::is_Field_v<T>, char>> : fmt::formatter<BoutReal> {
  auto format(const T& f, format_context& ctx) const -> format_context::iterator {
    const auto* mesh = f.getMesh();
    int previous_x = 0;
    int previous_y = 0;
    int previous_z = 0;

    BOUT_FOR(i, f.getRegion("RGN_ALL")) {
      const auto ix = i.x();
      const auto iy = i.y();
      const auto iz = i.z();

      if (iz > previous_z) {
        format_to(ctx.out(), " ");
      }
      if (iy > previous_y) {
        format_to(ctx.out(), "\n");
      }
      if (ix > previous_x) {
        format_to(ctx.out(), "\n\n");
      }

      format_to(ctx.out(), "{:c}: ", i);
      formatter<BoutReal>::format(f[i], ctx);
      format_to(ctx.out(), ";");
      previous_x = ix;
      previous_y = iy;
      previous_z = iz;
    }
    return format_to(ctx.out(), "");
  }
};

#endif // OUTPUT_BOUT_TYPES_H
