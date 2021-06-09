#include <bout.hxx>
#include <mpi.h>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);

  Field2D test;
  bout::globals::mesh->get(test, "test");

  Options dump;
  dump["test"] = test;
  bout::writeDefaultOutputFile(dump);

  MPI_Barrier(BoutComm::get());
  
  BoutFinalise();
  return 0;
}
