#include <bout/physicsmodel.hxx>
#include <cmath>
#include <type_traits>

// Mixed 2D/3D fields with different evolved boundaries exercise the complete
// interleaved packing map. Nonconstant data detect swapped points or fields.
class TestCvodeCuda : public PhysicsModel {
  Field2D a, b, initial_a, initial_b;
  Field3D c, d, initial_c, initial_d;

  template <class Field>
  static void initialise(Field& field, BoutReal offset) {
    field.allocate();
    BOUT_FOR_SERIAL(i, field.getRegion(RGN_ALL)) { field[i] = offset + 0.001 * i.ind; }
  }

  template <class Field>
  static void verify(const Field& field, const Field& initial, BoutReal rate,
                     BoutReal time, bool boundary) {
    const auto check = [&](const auto& i) {
      const auto expected = initial[i] * std::exp(-rate * time);
      if (!std::isfinite(field[i]) || std::abs(field[i] - expected) > 2e-7 * expected) {
        throw BoutException("CUDA CVODE analytic mismatch at index {}: {} vs {}", i.ind,
                            field[i], expected);
      }
    };
    BOUT_FOR_SERIAL(i, field.getRegion(RGN_NOBNDRY)) { check(i); }
    if (boundary) {
      // Only physical evolved boundaries belong to the solver vector. MPI
      // halos are not evolved, even when evolve_bndry=true.
      for (const auto& i : field.getMesh()->getRegion2D("RGN_BNDRY")) {
        if constexpr (std::is_same_v<Field, Field2D>) {
          check(i);
        } else {
          for (int z = 0; z < field.getMesh()->LocalNz; ++z) {
            check(field.getMesh()->ind2Dto3D(i, z));
          }
        }
      }
    }
  }

  int init(bool UNUSED(restarting)) override {
    solver->add(a, "a");
    solver->add(b, "b");
    solver->add(c, "c");
    solver->add(d, "d");
    initialise(a, 1.0);
    initialise(b, 2.0);
    initialise(c, 3.0);
    initialise(d, 4.0);
    initial_a = a;
    initial_b = b;
    initial_c = c;
    initial_d = d;
    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(a) = -0.1 * a;
    ddt(b) = -0.2 * b;
    ddt(c) = -0.3 * c;
    ddt(d) = -0.4 * d;
    return 0;
  }

  int outputMonitor(BoutReal time, int UNUSED(iter), int UNUSED(nout)) override {
    verify(a, initial_a, 0.1, time, false);
    verify(b, initial_b, 0.2, time, true);
    verify(c, initial_c, 0.3, time, false);
    verify(d, initial_d, 0.4, time, true);
    return 0;
  }
};

BOUTMAIN(TestCvodeCuda);
