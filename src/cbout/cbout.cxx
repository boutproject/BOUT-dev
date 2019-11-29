#include <cbout/bout.h>

#include <bout.hxx>
#include <field3d.hxx>

//////////////////////////////////////////////////////////////
// Initialisation

extern "C" int bout_initialise(int argc, char** argv) {
  return BoutInitialise(argc, argv);
}

extern "C" int bout_finalise() {
  return BoutFinalise(true);
}

//////////////////////////////////////////////////////////////
// Fields

extern "C" Field3D* Field3D_new() {
  return new Field3D();
}

extern "C" Field3D* Field3D_new_zerofrom(Field3D* field) {
  return new Field3D(zeroFrom(*field));
}

extern "C" Field3D* Field3D_new_emptyfrom(Field3D* field) {
  return new Field3D(emptyFrom(*field));
}

// Destroying

extern "C" void Field3D_delete(Field3D* field) {
  delete field;
}

// manipulating

extern "C" void Field3D_axpy(Field3D* result, BoutReal a, Field3D* x, Field3D* y) {
  (*result) = a * (*x) + (*y);
}

extern "C" void Field3D_scale(Field3D* x, BoutReal alpha) {
  (*x) = (*x) * alpha;
}

extern "C" int Field3D_getarray(Field3D* x, BoutReal **data) {
  try {
    *data = (*x)(0,0); // Pointer to start of internal data
    return BOUT_SUCCESS;
  } catch(BoutException &e) {
    return 1;
  }
}

//////////////////////////////////////////////////////////////
// Physics Models

namespace {
  class CAPIModel : public PhysicsModel {
  protected:
    /// Initialise
    int init(bool restarting) override {
      
    }
    /// Calculate time derivatives
    int rhs(BoutReal t) override {
      
    }
  };
}

