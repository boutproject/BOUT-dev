///
/// Implement {fmt} formatters for BOUT++ types
///

#ifndef OUTPUT_BOUT_TYPES_H
#define OUTPUT_BOUT_TYPES_H

#include "fmt/base.h"
#include <fmt/format.h>

#include "bout/bout_types.hxx"
#include "bout/field2d.hxx"
#include "bout/mesh.hxx"
#include "bout/output.hxx" // IWYU pragma: keep
#include "bout/region.hxx"

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

class Field2D;
class Field3D;
class FieldPerp;

template <>
struct fmt::formatter<Field2D> : fmt::formatter<BoutReal> {
  auto format(const Field2D& f, format_context& ctx) const -> format_context::iterator {
    const auto* mesh = f.getMesh();
    for (int ix = 0; ix < mesh->LocalNx; ++ix) {
      for (int jy = 0; jy < mesh->LocalNy; ++jy) {
        format_to(ctx.out(), "({}, {}): ", ix, jy);
        formatter<BoutReal>::format(f(ix, jy), ctx);
        format_to(ctx.out(), (jy < mesh->LocalNy - 1) ? "; " : ";");
      }
      format_to(ctx.out(), "\n");
    }
    return format_to(ctx.out(), "\n");
  }
};

template <>
struct fmt::formatter<Field3D> : fmt::formatter<BoutReal> {
  auto format(const Field3D& f, format_context& ctx) const -> format_context::iterator {
    const auto* mesh = f.getMesh();
    for (int ix = 0; ix < mesh->LocalNx; ++ix) {
      for (int jy = 0; jy < mesh->LocalNy; ++jy) {
        for (int kz = 0; kz < mesh->LocalNz; ++kz) {
          format_to(ctx.out(), "({}, {}, {}): ", ix, jy, kz);
          formatter<BoutReal>::format(f(ix, jy, kz), ctx);
          format_to(ctx.out(), (kz < mesh->LocalNz - 1) ? "; " : ";");
        }
        format_to(ctx.out(), "\n");
      }
      format_to(ctx.out(), "\n");
    }
    return format_to(ctx.out(), "\n");
  }
};

template <>
struct fmt::formatter<FieldPerp> : fmt::formatter<BoutReal> {
  auto format(const FieldPerp& f, format_context& ctx) const -> format_context::iterator {
    const auto* mesh = f.getMesh();
    for (int ix = 0; ix < mesh->LocalNx; ++ix) {
      for (int kz = 0; kz < mesh->LocalNz; ++kz) {
        format_to(ctx.out(), "({}, {}): ", ix, kz);
        formatter<BoutReal>::format(f(ix, kz), ctx);
        format_to(ctx.out(), (kz < mesh->LocalNz - 1) ? "; " : ";");
      }
      format_to(ctx.out(), "\n");
    }
    return format_to(ctx.out(), "\n");
  }
};

#endif // OUTPUT_BOUT_TYPES_H
