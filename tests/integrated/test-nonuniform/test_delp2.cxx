/*
 * Perpendicular Laplacian test
 * 
 * Tests the accuracy of the Delp2 operator on both uniform and
 * non-uniform grids. The same coefficients are used in both 
 * Delp2 and Laplacian inversion.
 * 
 */

#include <bout.hxx>
#include <bout/physicsmodel.hxx>

class Test_delp2 : public PhysicsModel {
protected:
  int init(bool UNUSED(restarting)) override;
  int rhs(BoutReal UNUSED(t)) override;
};


int Test_delp2::init(bool UNUSED(restarting)) {
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

int Test_delp2::rhs(BoutReal UNUSED(t)) {
  // Doesn't do anything
  return 1;
}


BOUTMAIN(Test_delp2)
