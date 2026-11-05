#include "../../../src/solver/sundials_cuda_packing.hxx"
#include <bout/boutcomm.hxx>
#include <bout/physicsmodel.hxx>
#include <bout/sundials_backports.hxx>
#include <bout/vector2d.hxx>
#include <bout/vector3d.hxx>
#include <cmath>
#include <memory>
#include <type_traits>

// Mixed 2D/3D fields with different evolved boundaries exercise the complete
// interleaved packing map. Nonconstant data detect swapped points or fields.
class TestCvodeCuda : public PhysicsModel {
  Field2D a, b, initial_a, initial_b;
  Field3D c, d, initial_c, initial_d;
  Vector2D e, initial_e;
  Vector3D f, initial_f;
  bool vector_fields = false;
  bool user_callbacks = false;
  int preconditioner_calls = 0;
  int jacobian_calls = 0;

  // Exercise real CUDA stream ordering independently of CVODE's tolerance and
  // convergence decisions. No implicit default-stream ordering is available.
  static void check_transfer_ordering() {
    auto comm = BoutComm::get();
    auto context = createSUNContext(comm);
    SundialsCudaPacking packing;
    using VectorPtr =
        std::unique_ptr<std::remove_pointer_t<N_Vector>, decltype(&N_VDestroy)>;
    VectorPtr u(callWithSUNContext(N_VNew_Cuda, context, 2), N_VDestroy);
    if (!u) {
      throw BoutException("CUDA ordering test vector allocation failed");
    }
    packing.bind(u.get());
    VectorPtr v(N_VClone(u.get()), N_VDestroy);
    if (!v) {
      throw BoutException("CUDA ordering test vector clone failed");
    }
    packing.check_stream(v.get());

    Field2D field;
    field.allocate();
    SundialsCudaPacking::FieldMap map;
    map.append(0, 1);
    map.append(1, 0);
    map.upload();
    auto* field_data = &field(0, 0);
    N_VConst(3.0, u.get());
    N_VScale(2.0, u.get(), v.get());
    map.transfer(field_data, N_VGetDeviceArrayPointer_Cuda(v.get()), true,
                 packing.stream);
    packing.synchronize();
    if (field_data[0] != 6.0 || field_data[1] != 6.0) {
      throw BoutException("CUDA packing did not follow vector streaming operations");
    }

    field_data[0] = 5.0;
    field_data[1] = 7.0;
    packing.retained_2d.push_back(field);
    map.transfer(field_data, N_VGetDeviceArrayPointer_Cuda(v.get()), false,
                 packing.stream);
    packing.record_field_reads();
    N_VScale(2.0, v.get(), u.get());
    // This wait covers packing reads, not the following vector scale.
    packing.wait_for_fields();
    field_data[0] = 100.0;
    field_data[1] = 200.0;
    map.transfer(field_data, N_VGetDeviceArrayPointer_Cuda(u.get()), true,
                 packing.stream);
    packing.synchronize();
    if (field_data[0] != 10.0 || field_data[1] != 14.0
        || N_VDotProd(u.get(), u.get()) != 296.0 || !packing.retained_2d.empty()) {
      throw BoutException("CUDA packing event/reduction ordering failed");
    }

    VectorPtr other(callWithSUNContext(N_VNew_Cuda, context, 2), N_VDestroy);
    if (!other) {
      throw BoutException("CUDA ordering test vector allocation failed");
    }
    bool rejected = false;
    try {
      packing.check_stream(other.get());
    } catch (const BoutException&) {
      rejected = true;
    }
    if (!rejected) {
      throw BoutException("CUDA packing accepted an unordered vector stream");
    }
    SundialsCudaPacking::FieldMap empty;
    empty.transfer(nullptr, nullptr, false, packing.stream);
  }

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
    check_transfer_ordering();
    vector_fields = Options::root()["test"]["vector_fields"].withDefault(false);
    user_callbacks = Options::root()["test"]["user_callbacks"].withDefault(false);
    if (user_callbacks) {
      setPrecon(&TestCvodeCuda::preconditioner);
      setJacobian(&TestCvodeCuda::jacobian);
    }
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
    if (vector_fields) {
      solver->add(e, "e");
      solver->add(f, "f");
      initialise(e.x, 5.0);
      initialise(e.y, 6.0);
      initialise(e.z, 7.0);
      initialise(f.x, 8.0);
      initialise(f.y, 9.0);
      initialise(f.z, 10.0);
      initial_e = e;
      initial_f = f;
    }
    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(a) = -0.1 * a;
    ddt(b) = -0.2 * b;
    ddt(c) = -0.3 * c;
    ddt(d) = -0.4 * d;
    if (vector_fields) {
      // Force preparation to convert the state back to its registered basis.
      e.toContravariant();
      f.toContravariant();
      ddt(e) = -0.5 * e;
      ddt(f) = -0.6 * f;
    }
    return 0;
  }

  template <class Factor>
  void scale_derivatives(Factor factor) {
    ddt(a) *= factor(0.1);
    ddt(b) *= factor(0.2);
    ddt(c) *= factor(0.3);
    ddt(d) *= factor(0.4);
    if (vector_fields) {
      ddt(e) = factor(0.5) * ddt(e);
      ddt(f) = factor(0.6) * ddt(f);
    }
  }

  int preconditioner(BoutReal UNUSED(time), BoutReal gamma, BoutReal UNUSED(delta)) {
    ++preconditioner_calls;
    scale_derivatives([gamma](BoutReal rate) { return 1.0 / (1.0 + gamma * rate); });
    return 0;
  }

  int jacobian(BoutReal UNUSED(time)) {
    ++jacobian_calls;
    scale_derivatives([](BoutReal rate) { return -rate; });
    return 0;
  }

  int timestepMonitor(BoutReal UNUSED(time), BoutReal UNUSED(dt)) override {
    // Reuse derivative storage on the CPU immediately after CVODE returns.
    // Queued callback saves must have stopped reading it before this point.
    const auto overwrite = [](auto& field) {
      BOUT_FOR_SERIAL(i, field.getRegion(RGN_ALL)) { field[i] = -123.0; }
    };
    overwrite(ddt(a));
    overwrite(ddt(b));
    overwrite(ddt(c));
    overwrite(ddt(d));
    if (vector_fields) {
      overwrite(ddt(e).x);
      overwrite(ddt(e).y);
      overwrite(ddt(e).z);
      overwrite(ddt(f).x);
      overwrite(ddt(f).y);
      overwrite(ddt(f).z);
    }
    return 0;
  }

  int outputMonitor(BoutReal time, int UNUSED(iter), int UNUSED(nout)) override {
    if (user_callbacks && time > 0.0
        && (preconditioner_calls == 0 || jacobian_calls == 0)) {
      throw BoutException("CUDA CVODE user callbacks were not exercised");
    }
    verify(a, initial_a, 0.1, time, false);
    verify(b, initial_b, 0.2, time, true);
    verify(c, initial_c, 0.3, time, false);
    verify(d, initial_d, 0.4, time, true);
    if (vector_fields) {
      e.toCovariant();
      f.toCovariant();
      verify(e.x, initial_e.x, 0.5, time, false);
      verify(e.y, initial_e.y, 0.5, time, false);
      verify(e.z, initial_e.z, 0.5, time, false);
      verify(f.x, initial_f.x, 0.6, time, false);
      verify(f.y, initial_f.y, 0.6, time, false);
      verify(f.z, initial_f.z, 0.6, time, false);
    }
    return 0;
  }
};

BOUTMAIN(TestCvodeCuda);
