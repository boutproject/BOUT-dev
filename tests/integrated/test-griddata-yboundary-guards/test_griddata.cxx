#include <bout.hxx>
#include <mpi.h>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Datafile df(Options::getRoot()->getSection("output"));

  Field2D test;
  mesh->get(test, "test");

  dump.add(test, "test");

  dump.write();

  MPI_Barrier(BoutComm::get());
  
  BoutFinalise();
  return 0;
}
