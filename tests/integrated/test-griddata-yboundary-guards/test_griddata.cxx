#include <bout.hxx>
#include <mpi.h>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field2D test;
  bout::globals::mesh->get(test, "test");

  bout::globals::dump.add(test, "test");

  bout::globals::dump.write();

  MPI_Barrier(BoutComm::get());
  
  BoutFinalise();
  return 0;
}
