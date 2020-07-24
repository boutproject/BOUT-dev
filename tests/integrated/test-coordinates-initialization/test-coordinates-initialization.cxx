/*
 * Test of initialization of Coordinates objects in Mesh::coords_map
 */

#include "bout.hxx"
#include "optionsreader.hxx"

int main() {

  // Initialize options, needed to load mesh from BOUT.inp
  Options *options = Options::getRoot();
  OptionsReader *reader = OptionsReader::getInstance();
  reader->read(options, "data/BOUT.inp");

  bout::globals::mpi = new MpiWrapper();

  // Initialize a mesh
  auto mesh = Mesh::create();
  mesh->load();

  // Test CELL_CENTRE
  Field3D f(0., mesh);
  f.applyBoundary("neumann");

  // Synchronise all processors
  MPI_Barrier(BoutComm::get());

  // Test CELL_YLOW
  Field3D f_ylow(0., mesh);
  f_ylow.setLocation(CELL_YLOW);
  f_ylow.applyBoundary("neumann");

  // Synchronise all processors
  MPI_Barrier(BoutComm::get());

  MPI_Finalize();

  return 0;
}
