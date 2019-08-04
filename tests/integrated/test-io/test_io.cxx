/*
 * Input/Output regression test
 *
 * Read from and write to data files to check that
 * the I/O routines are working.
 *
 * Test evolving and non-evolving variables
 */

#include <bout.hxx>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  // Variables to be read and written
  int ivar, ivar_evol;
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
  dump.add(rvar, "rvar", false);
  dump.add(bvar, "bvar", false);
  dump.add(f2d, "f2d", false);
  dump.add(f3d, "f3d", false);
  dump.add(fperp, "fperp", false);
  dump.add(fperp2, "fperp2", false);

  // Evolving variables
  dump.add(ivar_evol, "ivar_evol", true);
  dump.add(rvar_evol, "rvar_evol", true);
  dump.add(bvar_evol, "bvar_evol", true);
  dump.add(v2d, "v2d_evol", true);
  dump.add(v3d, "v3d_evol", true);
  dump.add(fperp2_evol, "fperp2_evol", true);

  int MYPE;
  MPI_Comm_rank(BoutComm::get(), &MYPE);

  for(int i=0;i<3;i++) {
    ivar_evol = ivar + i;
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
