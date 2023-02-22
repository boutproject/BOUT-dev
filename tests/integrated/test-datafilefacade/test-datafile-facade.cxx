#include "bout/physicsmodel.hxx"
#include "bout/sys/timer.hxx"
#include "bout/vector2d.hxx"
#include "bout/vector3d.hxx"

class TestDataFileFacade : public PhysicsModel {
  int init([[maybe_unused]] bool restarting) override {
    solver->add(f3d, "f3d");
    solver->add(f2d, "f2d");

    f2d = 0.1;
    f3d = 0.2;

    fperp.setIndexFromGlobal(2);
    fperp = 1.1;
    SAVE_ONCE(fperp);

    v2d_contravariant.covariant = false;
    v2d_contravariant.x = 2.2;
    v2d_contravariant.y = 3.3;
    v2d_contravariant.z = 4.4;

    v3d_contravariant.covariant = false;
    v3d_contravariant.x = 5.5;
    v3d_contravariant.y = 6.6;
    v3d_contravariant.z = 7.7;
    SAVE_ONCE(v2d_contravariant, v3d_contravariant);

    v2d_covariant.x = 12.12;
    v2d_covariant.y = 13.13;
    v2d_covariant.z = 14.14;

    v3d_covariant.x = 15.15;
    v3d_covariant.y = 16.16;
    v3d_covariant.z = 17.17;
    SAVE_ONCE(v2d_covariant, v3d_covariant);

    SAVE_ONCE(integer, boolean, real, name);

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(f3d) = 0.;
    ddt(f2d) = 0.;
    return 0;
  }

  Field3D f3d;
  Field2D f2d;
  FieldPerp fperp;
  Vector2D v2d_contravariant;
  Vector3D v3d_contravariant;
  Vector2D v2d_covariant;
  Vector3D v3d_covariant;
  int integer = 42;
  bool boolean = true;
  BoutReal real = 3.14;
  std::string name = "test string";
};

BOUTMAIN(TestDataFileFacade);
