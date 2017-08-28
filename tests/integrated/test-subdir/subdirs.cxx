#include <bout/physicsmodel.hxx>

// functions I need in subprojects
#include "fuu/fuu.hxx"
void bar();

class SubDirs : public PhysicsModel {
private:
  Field3D n;
  
protected:
  int init(bool restart) override {
    n=1;
    SOLVE_FOR(n);
    // test functions
    fuu();
    bar();
    return 0;
  }
  
  int rhs(BoutReal time) override {
    ddt(n)=0;
    return 0;
  }
};

BOUTMAIN(SubDirs);

