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

#include <fmt/base.h>
#include <fmt/format.h>

#include <cstddef>
#include <string>
#include <string_view>
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
///
/// Format specification:
///
/// - ``n``: Don't show indices
/// - ``r'<region name>'``: Use given region (default: ``RGN_ALL``)
template <class T>
struct fmt::formatter<T, std::enable_if_t<bout::utils::is_Field_v<T>, char>> {
private:
  fmt::formatter<BoutReal> underlying;

  static constexpr auto default_region = "RGN_ALL";
  std::string_view region = default_region;

  bool show_indices = true;

public:
  constexpr auto parse(format_parse_context& ctx) {
    const auto* it = ctx.begin();
    const auto* end = ctx.end();

    if (it == end) {
      return underlying.parse(ctx);
    }

    while (it != end and *it != ':' and *it != '}') {
      // Other cases handled explicitly below
      // NOLINTNEXTLINE(bugprone-switch-missing-default-case)
      switch (*it) {
      case 'r':
        ++it;
        if (*it != '\'') {
          throw fmt::format_error("invalid format for Field");
        }
        {
          const auto* rgn_start = ++it;
          std::size_t size = 0;
          while (*it != '\'') {
            ++size;
            ++it;
          }
          region = std::string_view(rgn_start, size);
        }
        ++it;
        break;
      case 'n':
        show_indices = false;
        ++it;
        break;
      }
    }

    if (it != end and *it != '}') {
      if (*it != ':') {
        throw fmt::format_error("invalid format specifier");
      }
      ++it;
    }

    ctx.advance_to(it);
    return underlying.parse(ctx);
  }

  auto format(const T& f, format_context& ctx) const -> format_context::iterator {
    const auto* mesh = f.getMesh();

    const auto rgn = f.getRegion(std::string(region));
    const auto i = rgn.begin();
    int previous_x = i->x();
    int previous_y = i->y();
    int previous_z = i->z();

    BOUT_FOR(i, rgn) {
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

      if (show_indices) {
        format_to(ctx.out(), "{:c}: ", i);
      }
      underlying.format(f[i], ctx);
      format_to(ctx.out(), ";");
      previous_x = ix;
      previous_y = iy;
      previous_z = iz;
    }
    return format_to(ctx.out(), "");
  }
};

#endif // OUTPUT_BOUT_TYPES_H
