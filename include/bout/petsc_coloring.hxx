#include "bout/build_config.hxx"

#ifndef PETSC_COLORING_H
#define PETSC_COLORING_H

#if BOUT_HAS_PETSC

#include <cmath>
#include <functional>
#include <vector>

#include <petsc.h>

class ColoringStencil {
private:
  bool static isInSquare(int const i, int const j, int const n_square) {
    return std::abs(i) <= n_square && std::abs(j) <= n_square;
  }
  bool static isInCross(int const i, int const j, int const n_cross) {
    if (i == 0) {
      return std::abs(j) <= n_cross;
    }
    if (j == 0) {
      return std::abs(i) <= n_cross;
    }
    return false;
  }
  bool static isInTaxi(int const i, int const j, int const n_taxi) {
    return std::abs(i) + std::abs(j) <= n_taxi;
  }

public:
  auto static getOffsets(int n_square, int n_taxi, int n_cross) {
    ASSERT2(n_square >= 0 && n_cross >= 0 && n_taxi >= 0
            && n_square + n_cross + n_taxi > 0);
    auto inside = [&](int i, int j) {
      return isInSquare(i, j, n_square) || isInTaxi(i, j, n_taxi)
             || isInCross(i, j, n_cross);
    };
    std::vector<std::pair<int, int>> xy_offsets;
    auto loop_bound = std::max({n_square, n_taxi, n_cross});
    for (int i = -loop_bound; i <= loop_bound; ++i) {
      for (int j = -loop_bound; j <= loop_bound; ++j) {
        if (inside(i, j)) {
          xy_offsets.emplace_back(i, j);
        }
      }
    }
    return xy_offsets;
  }
};

struct PetscPreconditioner {
  /// Represents a RHS function f(x): (Vec x, Vec f) -> PetscErrorCode
  using FormFunction = std::function<PetscErrorCode(Vec, Vec)>;

  /// @param[in] local_index  Local indices, starting from 0 on each processor
  /// @param[in] form_function   The operator to precondition.
  ///                   This can be a linearised and simplified form.
  PetscPreconditioner(Options& options, int nlocal, Field3D local_index,
                      FormFunction form_function);

  /// This is public because it is called from C callback functions
  FormFunction form_function;

  /// Return the finite difference Jacobian
  Mat jacobian() { return Jfd; }

private:
  Mesh* mesh;

  Mat Jfd; ///< Finite Difference Jacobian

  bool use_coloring;
  MatFDColoring fdcoloring{nullptr}; ///< Matrix coloring context
                                     ///< Jacobian evaluation

  void updateColoring(); ///< Updates the coloring using Jfd
};

#else // BOUT_HAS_PETSC

#endif // BOUT_HAS_PETSC

#endif // PETSC_COLORING_H
