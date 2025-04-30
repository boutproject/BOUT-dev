#include "bout/build_config.hxx"

#if BOUT_HAS_PETSC

#include <cmath>
#include <vector>


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



#endif // BOUT_HAS_PETSC
