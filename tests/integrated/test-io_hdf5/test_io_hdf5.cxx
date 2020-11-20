/*
 * Input/Output regression test
 *
 * Read from and write to data files to check that
 * the I/O routines are working.
 *
 * Test evolving and non-evolving variables
 */

#include <bout.hxx>

using bout::globals::dump;
using bout::globals::mesh;

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  // Variables to be read and written
  int ivar, ivar_evol;
  std::vector<int> ivar_vec = {1, 2, 3};
  std::vector<int> ivar_vec_evol = {4, 5, 6};
  std::vector<char> cvar_vec = {'a', 'b'};
  std::vector<char> cvar_vec_evol = {'c', 'd'};
  std::string svar = "ijkl";
  std::string svar_evol = "mnopq";
  BoutReal rvar, rvar_evol;
  bool bvar, bvar_evol;
  Field2D f2d;
  Field3D f3d;
  // fperp is at yindex_global=0.
  // fperp2 is at yindex_global=11, it is included to make sure the test does not pass
  // only for the special case of the FieldPerp being present on processor number 0.
  FieldPerp fperp, fperp2, fperp2_evol;
  Vector2D v2d;
  Vector3D v3d;
  bool check_double_add = Options::root()["check_double_add"].withDefault(false);
  auto check_incorrect_add = Options::root()["check_incorrect_add"].withDefault("none");

  f2d = 0.0;
  f3d = 0.0;
  fperp = 0.0;
  fperp2 = 0.0;

  // Read data from grid file
  mesh->get(ivar, "ivar");
  mesh->get(rvar, "rvar");
  mesh->get(bvar, "bvar");
  mesh->get(f2d, "f2d");
  mesh->get(f3d, "f3d");
  mesh->get(fperp, "fperp");
  mesh->get(fperp2, "fperp2");

  // Non-evolving variables
  dump.add(ivar, "ivar", false);
  dump.add(ivar_vec, "ivar_vec");
  dump.add(cvar_vec, "cvar_vec");
  dump.add(svar, "svar");
  dump.add(rvar, "rvar", false);
  dump.add(bvar, "bvar", false);
  dump.add(f2d, "f2d", false);
  dump.add(f3d, "f3d", false);
  dump.add(fperp, "fperp", false);
  dump.add(fperp2, "fperp2", false);

  // Evolving variables
  dump.add(ivar_evol, "ivar_evol", true);
  dump.add(ivar_vec_evol, "ivar_vec_evol", true);
  dump.add(cvar_vec_evol, "cvar_vec_evol", true);
  dump.add(svar_evol, "svar_evol", true);
  dump.add(rvar_evol, "rvar_evol", true);
  dump.add(bvar_evol, "bvar_evol", true);
  dump.add(v2d, "v2d_evol", true);
  dump.add(v3d, "v3d_evol", true);
  dump.add(fperp2_evol, "fperp2_evol", true);

  if (check_double_add) {
    // Add all variables twice to check this does not cause an error
    dump.add(ivar, "ivar", false);
    dump.add(ivar_vec, "ivar_vec");
    dump.add(cvar_vec, "cvar_vec");
    dump.add(svar, "svar");
    dump.add(rvar, "rvar", false);
    dump.add(bvar, "bvar", false);
    dump.add(f2d, "f2d", false);
    dump.add(f3d, "f3d", false);
    dump.add(fperp, "fperp", false);
    dump.add(fperp2, "fperp2", false);
    dump.add(ivar_evol, "ivar_evol", true);
    dump.add(ivar_vec_evol, "ivar_vec_evol");
    dump.add(cvar_vec_evol, "cvar_vec_evol");
    dump.add(svar_evol, "svar_evol");
    dump.add(rvar_evol, "rvar_evol", true);
    dump.add(bvar_evol, "bvar_evol", true);
    dump.add(v2d, "v2d_evol", true);
    dump.add(v3d, "v3d_evol", true);
    dump.add(fperp2_evol, "fperp2_evol", true);
  }

  // Cases to check expected fails
  if (check_incorrect_add == "ivar") {
    int dummy = 0;
    dump.add(dummy, "ivar", false);
  } else if (check_incorrect_add == "ivar_vec") {
    std::vector<int> dummy = {-1};
    dump.add(dummy, "ivar_vec", false);
  } else if (check_incorrect_add == "cvar_vec") {
    std::vector<char> dummy = {'z'};
    dump.add(dummy, "cvar_vec", false);
  } else if (check_incorrect_add == "svar") {
    std::string dummy = "y";
    dump.add(dummy, "svar", false);
  } else if (check_incorrect_add == "rvar") {
    BoutReal dummy = 0.0;
    dump.add(dummy, "rvar", false);
  } else if (check_incorrect_add == "bvar") {
    bool dummy = false;
    dump.add(dummy, "bvar", false);
  } else if (check_incorrect_add == "f2d") {
    Field2D dummy = 0.0;
    dump.add(dummy, "f2d", false);
  } else if (check_incorrect_add == "f3d") {
    Field3D dummy = 0.0;
    dump.add(dummy, "f3d", false);
  } else if (check_incorrect_add == "fperp") {
    FieldPerp dummy = 0.0;
    dump.add(dummy, "fperp", false);
  } else if (check_incorrect_add == "ivar_evol") {
    int dummy = 0;
    dump.add(dummy, "ivar_evol", true);
  } else if (check_incorrect_add == "ivar_vec_evol") {
    std::vector<int> dummy = {-1};
    dump.add(dummy, "ivar_vec_evol", false);
  } else if (check_incorrect_add == "cvar_vec_evol") {
    std::vector<char> dummy = {'z'};
    dump.add(dummy, "cvar_vec_evol", false);
  } else if (check_incorrect_add == "svar_evol") {
    std::string dummy = "y";
    dump.add(dummy, "svar_evol", false);
  } else if (check_incorrect_add == "rvar_evol") {
    BoutReal dummy = 0.0;
    dump.add(dummy, "rvar_evol", true);
  } else if (check_incorrect_add == "bvar_evol") {
    bool dummy = false;
    dump.add(dummy, "bvar_evol", true);
  } else if (check_incorrect_add == "v2d_evol") {
    Vector2D dummy;
    dump.add(dummy, "v2d_evol", true);
  } else if (check_incorrect_add == "v3d_evol") {
    Vector3D dummy;
    dump.add(dummy, "v3d_evol", true);
  }

  int MYPE;
  MPI_Comm_rank(BoutComm::get(), &MYPE);

  bvar_evol = bvar;
  for(int i=0;i<3;i++) {
    ivar_evol = ivar + i;
    ivar_vec_evol[0] += i; ivar_vec_evol[1] += i; ivar_vec_evol[2] += i;
    cvar_vec_evol[0] += i; cvar_vec_evol[1] += i;
    rvar_evol = rvar + 0.5 * i;
    bvar_evol = !bvar_evol;
    v2d.x = v2d.y = v2d.z = f2d;
    v3d.x = v3d.y = v3d.z = f3d;
    fperp2_evol = fperp2;

    dump.write();
  }

  dump.close(); // Ensure data is written

  // Need to wait for all processes to finish writing
  MPI_Barrier(BoutComm::get());

  /// Finished, tidy up and free memory
  BoutFinalise();

  return 0;
}
