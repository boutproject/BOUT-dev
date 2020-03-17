/*
 * Perpendicular Laplacian test
 * 
 * Tests the accuracy of the Delp2 operator on both uniform and
 * non-uniform grids. The same coefficients are used in both 
 * Delp2 and Laplacian inversion.
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>

int physics_init(bool UNUSED(restarting)) {
  Field3D input, reference, result;
  
  GRID_LOAD(input);                  // Read input from file
  GRID_LOAD(reference);              // Reference solution
  
  result = Delp2(input);             // Calculate perpendicular Laplacian
  result.applyBoundary("dirichlet"); // Apply Dirichlet boundary conditions
  
  SAVE_ONCE(input);
  SAVE_ONCE(result);                 // Queue the output variable to be written
  SAVE_ONCE(reference);              // Write the reference solution too
  

  dump.write();
  dump.close();
  
  MPI_Barrier(BoutComm::get());
  
  output << "\nFinished running test. Triggering error to quit\n\n";
  
  return 1;
}

int physics_run(BoutReal UNUSED(t)) {
  // Doesn't do anything
  return 1;
}
