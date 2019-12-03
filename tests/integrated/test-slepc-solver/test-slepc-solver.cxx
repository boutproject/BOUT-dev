// A simple linear problem with one eigenvalue

#include <bout/physicsmodel.hxx>

class TestEigenSolver : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override {
    solver->add(field, "f");
    return 0;
  }
  int rhs(BoutReal time) override {
    mesh->communicate(field);
    ddt(field) = field * exp(time);
    return 0;
  }
private:
  Field3D field;
};

BOUTMAIN(TestEigenSolver);
