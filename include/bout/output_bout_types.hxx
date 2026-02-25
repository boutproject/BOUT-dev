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
#include <fmt/color.h>
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

namespace bout {
namespace details {
template <class T>
/// Transpose a region so that it iterates in Z first, then Y, then X
///
/// Caution: this is the most inefficient memory order!
auto region_transpose(const Region<T>& region) -> Region<T>;

auto colour(BoutReal value, BoutReal min, BoutReal max) -> fmt::text_style;
} // namespace details
} // namespace bout

/// Formatter for Fields
///
/// Format specification:
///
/// - ``n``: Don't show indices
/// - ``r'<region name>'``: Use given region (default: ``RGN_ALL``)
/// - ``T``: Transpose field so X is first dimension
/// - ``#``: Plot slices as 2D heatmap
template <class T>
struct fmt::formatter<T, std::enable_if_t<bout::utils::is_Field_v<T>, char>> {
private:
  fmt::formatter<BoutReal> underlying;

  static constexpr auto default_region = "RGN_ALL";
  std::string_view region = default_region;

  bool show_indices = true;
  bool transpose = false;
  bool plot = false;

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
      case 'T':
        transpose = true;
        ++it;
        break;
      case '#':
        plot = true;
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
    using namespace bout::details;

    const auto* mesh = f.getMesh();

    const auto rgn_str = std::string{region};
    const auto rgn_ = f.getRegion(rgn_str);
    const auto rgn = transpose ? region_transpose(rgn_) : rgn_;

    const auto i = rgn.begin();
    int previous_x = i->x();
    int previous_y = i->y();
    int previous_z = i->z();

    // Range of the data for plotting
    BoutReal plot_min = 0.0;
    BoutReal plot_max = 0.0;
    if (plot) {
      plot_min = min(f, false, rgn_str);
      plot_max = max(f, false, rgn_str);
    }

    // Separators
    const auto* const block_sep = "\n\n";
    const auto* const item_sep = plot ? "" : " ";

    BOUT_FOR(i, rgn) {
      const auto ix = i.x();
      const auto iy = i.y();
      const auto iz = i.z();

      if (iz > previous_z) {
        format_to(ctx.out(), transpose ? block_sep : item_sep);
      }
      if (iy > previous_y) {
        format_to(ctx.out(), "\n");
      }
      if (ix > previous_x) {
        format_to(ctx.out(), transpose ? item_sep : block_sep);
      }

      if (show_indices) {
        format_to(ctx.out(), "{:c}: ", i);
      }
      if (plot) {
        format_to(ctx.out(), "{}", styled("█", colour(f[i], plot_min, plot_max)));
      } else {
        underlying.format(f[i], ctx);
        format_to(ctx.out(), ";");
      }
      previous_x = ix;
      previous_y = iy;
      previous_z = iz;
    }
    return format_to(ctx.out(), "");
  }
};

#endif // OUTPUT_BOUT_TYPES_H
