#include "bout/physicsmodel.hxx"

class TestRestartIO : public PhysicsModel {
  int init(bool UNUSED(restarting)) override {
    solver->add(f3d, "f3d");
    solver->add(f2d, "f2d");
    dump.addRepeat(fperp_lower, "fperp_lower");
    dump.addRepeat(fperp_upper, "fperp_upper");
    restart.addOnce(fperp_lower, "fperp_lower");
    restart.addOnce(fperp_upper, "fperp_upper");

    dump.addOnce(f3d, "f3d_once");
    dump.addOnce(f2d, "f2d_once");
    dump.addOnce(fperp_lower, "fperp_lower_once");
    dump.addOnce(fperp_upper, "fperp_upper_once");

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(f3d) = 0.;
    ddt(f2d) = 0.;
    return 0;
  }

  Field3D f3d;
  Field2D f2d;
  // fperp_lower is at yindex_global=0.
  // fperp_upper is at yindex_global=16, it is included to make sure the test does not
  // pass only for the special case of the FieldPerp being present on prcossor number 0.
  FieldPerp fperp_lower, fperp_upper;
};

BOUTMAIN(TestRestartIO);
