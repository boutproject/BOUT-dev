#include <bout/yboundary_regions.hxx>

class yboundary_example_legacy {
public:
  yboundary_example_legacy(Options* opt, const Field3D& N, const Field3D& V)
      : N(N), V(V) {
    Options& options = *opt;
    lower_y = options["lower_y"].doc("Boundary on lower y?").withDefault<bool>(lower_y);
    upper_y = options["upper_y"].doc("Boundary on upper y?").withDefault<bool>(upper_y);
  }

  void rhs() {
    BoutReal totalFlux = 0;
    if (lower_y) {
      for (RangeIterator r = mesh->iterateBndryLowerY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Calculate flux through surface [normalised m^-2 s^-1],
          // should be positive since V < 0.0
          BoutReal flux =
              -0.5 * (N(r.ind, mesh->ystart, jz) + N(r.ind, mesh->ystart - 1, jz)) * 0.5
              * (V(r.ind, mesh->ystart, jz) + V(r.ind, mesh->ystart - 1, jz));
          totalFlux += flux;
        }
      }
    }
    if (upper_y) {
      for (RangeIterator r = mesh->iterateBndryUpperY(); !r.isDone(); r++) {
        for (int jz = 0; jz < mesh->LocalNz; jz++) {
          // Calculate flux through surface [normalised m^-2 s^-1],
          // should be positive since V < 0.0
          BoutReal flux = -0.5 * (N(r.ind, mesh->yend, jz) + N(r.ind, mesh->yend + 1, jz))
                          * 0.5
                          * (V(r.ind, mesh->yend, jz) + V(r.ind, mesh->yend + 1, jz));
          totalFlux += flux;
        }
      }
    }
  }

private:
  bool lower_y{true};
  bool upper_y{true};
  const Field3D& N;
  const Field3D& V;
}

class yboundary_example {
public:
  yboundary_example(Options* opt, const Field3D& N, const Field3D& V) : N(N), V(V) {
    // Set what kind of yboundaries you want to include
    yboundary.init(opt);
  }

  void rhs() {
    BoutReal totalFlux = 0;
    yboundary.iter_pnts([&](auto& pnt) {
      BoutReal flux = pnt.interpolate_sheath_o1(N) * pnt.interpolate_sheath_o1(V);
    });
  }

private:
  YBoundary ybounday;
  const Field3D& N;
  const Field3D& V;
};
