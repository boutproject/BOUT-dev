#include "bout/physicsmodel.hxx"

class TestRestartIO : public PhysicsModel {
  int init(bool restarting) override {
    solver->add(f3d, "f3d");
    solver->add(f2d, "f2d");

    if (restarting) {
      fperp_lower = readFromRestartFile("fperp_lower").as<FieldPerp>();
      fperp_upper = readFromRestartFile("fperp_upper").as<FieldPerp>();
    }

    return 0;
  }

  int rhs(BoutReal UNUSED(time)) override {
    ddt(f3d) = 0.;
    ddt(f2d) = 0.;
    return 0;
  }

  void outputVars(Options& options) override {
    options["fperp_lower"].assignRepeat(fperp_lower);
    options["fperp_upper"].assignRepeat(fperp_upper);
    options["f3d_once"] = f3d;
    options["f2d_once"] = f2d;
    options["fperp_lower_once"] = fperp_lower;
    options["fperp_upper_once"] = fperp_upper;
  }

  void restartVars(Options& restart) override {
    restart["fperp_lower"] = fperp_lower;
    restart["fperp_upper"] = fperp_upper;
  }

  Field3D f3d;
  Field2D f2d;
  // fperp_lower is at yindex_global=0.
  // fperp_upper is at yindex_global=16, it is included to make sure the test does not
  // pass only for the special case of the FieldPerp being present on prcossor number 0.
  FieldPerp fperp_lower, fperp_upper;
};

BOUTMAIN(TestRestartIO);
