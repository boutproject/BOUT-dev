/*
 * Laplacian inversion
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>

int physics_init(bool restarting) {
  FieldFactory f(mesh);
  
  Options *options = Options::getRoot();

  Field3D input = f.create3D("(1-gauss(x-0.5,0.2))*gauss(z-pi)");
  Field2D a = f.create2D("gauss(x)");
  Field2D c = f.create2D("sin(x) * gauss(x-0.5)");
  SAVE_ONCE3(input, a, c);
  
  // Create two solvers, using different options
  class Laplacian *solver1 = Laplacian::create(options->getSection("solver1"));
  class Laplacian *solver2 = Laplacian::create(options->getSection("solver2"));
  
  Field3D result1 = solver1->solve(input);
  Field3D result2 = solver2->solve(input);

  SAVE_ONCE2(result1, result2);
  
  dump.write();
  dump.close();
  
  output << "\nFinished running test. Triggering error to quit\n\n";
  
  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data
  
  return 1;
}

int physics_run(BoutReal t) {
  // Doesn't do anything
  return 1;
}
