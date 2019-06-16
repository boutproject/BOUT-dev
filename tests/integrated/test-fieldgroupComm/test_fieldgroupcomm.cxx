#include <bout/physicsmodel.hxx>


class TestFieldGroupComm : public PhysicsModel{
protected:
  int init(bool UNUSED(restarting)) {
    //Create identical fields
    solver->add(fld1,"fld1");
    solver->add(fld2,"fld2");
    solver->add(fld3,"fld3");

    //Create different communicators
    comm1.add(fld1);
    comm2.add(fld2,fld2);
    comm3.add(fld3);

    return 0;
  }

  int rhs(BoutReal UNUSED(t)) {
    //Communicate -- ideally would all produce the same result
    //1. As it should be
    mesh->communicate(comm1);
    //2. Duplicate entry
    mesh->communicate(comm2);
    //3. Twice with single entry
    mesh->communicate(comm3);
    mesh->communicate(comm3);

    ddt(fld1) = Grad_par(fld1);
    ddt(fld2) = Grad_par(fld2);
    ddt(fld3) = Grad_par(fld3);
    return 0;
  }

private:
  Field3D fld1, fld2, fld3;
  FieldGroup comm1, comm2, comm3;
};

BOUTMAIN(TestFieldGroupComm);
