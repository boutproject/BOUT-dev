#include <bout/physicsmodel.hxx>
#include <bout/sys/timer.hxx>

class CommsSpeed : public PhysicsModel {
public:
  Field3D f;

  int init(bool) {

    SOLVE_FOR(f);

    return 0;
  }

  int rhs(BoutReal) {

    auto mpi = mesh->getMpi();

    mpi.MPI_Barrier(BoutComm::get());
    {
      Timer timer("pre-synchronised communications");
      mesh->communicate(f);
    }

    ddt(f) = 1.;

    return 0;
  }

  int outputMonitor(BoutReal, int, int) {

    output<<endl<<"just comms "<<Timer::getTime("pre-synchronised communications")<<endl;

    return 0;
  }
};

BOUTMAIN(CommsSpeed);
