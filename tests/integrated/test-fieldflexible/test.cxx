#include <bout/physicsmodel.hxx>
#include <flexible.hxx>

class Test : public PhysicsModel {
private:
  Field3D n;

protected:
  int init(bool restart) override {
    n = 1;
    SOLVE_FOR(n);
    // test functions
    Field3D a = n;
    Flexible<Field3D> af = a;
    Field3D b = af.get(CELL_CENTRE);
    b = af.get(CELL_YLOW);
    output << strLocation(b.getLocation());
    n.setLocation(CELL_ZLOW);
    n *= af;
    // Test of cast
    Field3D ai = af;
    return 0;
  }

  int rhs(BoutReal time) override {
    ddt(n) = 0;
    return 0;
  }
};

BOUTMAIN(Test);
