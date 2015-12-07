/*
 * Input/Output regression test
 * 
 * Read from and write to data files to check that 
 * the I/O routines are working. 
 *
 * Test evolving and non-evolving variables
 */

#include <bout.hxx>

int main(int argc, char**argv) {

  /// Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);
  
  // Variables to be read and written
  int ivar, ivar_evol;
  BoutReal rvar, rvar_evol; 
  Field2D f2d;
  Field3D f3d;
  Vector2D v2d;
  Vector3D v3d;
  
  f2d = 0.0;
  f3d = 0.0;

  // Read data from grid file
  mesh->get(ivar, "ivar");
  mesh->get(rvar, "rvar");
  mesh->get(f2d, "f2d");
  mesh->get(f3d, "f3d");
  
  // Non-evolving variables
  dump.add(ivar, "ivar", 0);
  dump.add(rvar, "rvar", 0);
  dump.add(f2d, "f2d", 0);
  dump.add(f3d, "f3d", 0);
  
  // Evolving variables
  dump.add(ivar_evol, "ivar_evol", 1);
  dump.add(rvar_evol, "rvar_evol", 1);
  dump.add(v2d, "v2d_evol", 1);
  dump.add(v3d, "v3d_evol", 1);

  int MYPE;
  MPI_Comm_rank(BoutComm::get(), &MYPE);
  
  for(int i=0;i<3;i++) {
    ivar_evol = ivar + i;
    rvar_evol = rvar + 0.5 * i;
    v2d.x = v2d.y = v2d.z = f2d;
    v3d.x = v3d.y = v3d.z = f3d;
    
    dump.write();
  }
  
  dump.close(); // Ensure data is written

  // Need to wait for all processes to finish writing
  MPI_Barrier(BoutComm::get());

  /// Finished, tidy up and free memory
  BoutFinalise();
  
  return 0;
}
