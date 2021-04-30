///
/// Implement {fmt} formatters for BOUT++ types
///

#ifndef OUTPUT_BOUT_TYPES_H
#define OUTPUT_BOUT_TYPES_H

#include <fmt/format.h>

#include "output.hxx"
#include "bout/region.hxx"

template <IND_TYPE N>
struct fmt::formatter<SpecificInd<N>> {
  // Presentation format: 'c' - components, 'i' - index.
  char presentation = 'c';

  // Parses format specifications of the form ['c' | 'i'].
  constexpr auto parse(format_parse_context& ctx) {
    auto it = ctx.begin(), end = ctx.end();
    if (it != end && (*it == 'c' || *it == 'i')) presentation = *it++;

    // Check if reached the end of the range:
    if (it != end && *it != '}')
      throw format_error("invalid format");

    // Return an iterator past the end of the parsed range:
    return it;
  }

  // Formats the point p using the parsed format specification (presentation)
  // stored in this formatter.
  template <typename FormatContext>
  auto format(const Ind3D& ind, FormatContext& ctx) {
    // ctx.out() is an output iterator to write to.
    if (presentation == 'c') {
      return format_to(ctx.out(),
                       "({}, {}, {})",
                       ind.x(), ind.y(), ind.z());
    }
    return format_to(ctx.out(),
                     "({})",
                     ind.ind);
  }
};

#endif // OUTPUT_BOUT_TYPES_H
