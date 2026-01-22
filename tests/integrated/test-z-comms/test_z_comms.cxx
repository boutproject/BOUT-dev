#include "bout/assert.hxx"
#include "bout/bout.hxx"
#include "bout/field2d.hxx"
#include "bout/field3d.hxx"
#include "bout/globals.hxx"
#include "bout/mesh.hxx"
#include "bout/options.hxx"
#include "bout/output.hxx"

#include <memory>

int main(int argc, char** argv) {
  BoutComm::setArgs(argc, argv);

  const auto rank = BoutComm::rank();
  const auto size = BoutComm::size();

  if (rank == 0) {
    Output::getInstance()->enable();
  } else {
    Output::getInstance()->disable();
  }

  auto& root = Options::root();
  root["mesh"] = {
      {"MXG", 1}, {"MYG", 1}, {"MZG", 1}, {"nx", 3}, {"ny", 1}, {"nz", size},
  };

  root["NZPE"] = size;
  root["output"]["enabled"] = false;
  root["restart_files"]["enabled"] = false;
  root["datadir"] = "data";

  auto mpi_wrapper = std::make_unique<MpiWrapper>();
  bout::globals::mpi = mpi_wrapper.get();

  bout::globals::mesh = Mesh::create();
  bout::globals::mesh->load();

  Field3D field = rank;

  BOUT_FOR(i, field.getRegion("RGN_ZGUARDS")) { field[i] = -1; }

  bout::globals::mesh->communicate(field);

  if (rank == 0) {
    ASSERT0(field(1, 1, 0) == size - 1);
  }
  if (rank == size) {
    ASSERT0(field(1, 1, 2) == 0);
  }

  BoutComm::cleanup();
}
