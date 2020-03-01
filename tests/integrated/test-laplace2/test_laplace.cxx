/*
 * Laplacian inversion
 * 
 */

#include <bout.hxx>
#include <boutmain.hxx>
#include <invert_laplace.hxx>
#include <field_factory.hxx>

int physics_init(bool UNUSED(restarting)) {
  FieldFactory f(mesh);
  
  Options *options = Options::getRoot();

  // Read strings containing coefficients
  std::string in, acoef, ccoef;
  OPTION(options, in, "(1-gauss(x-0.5,0.2))*gauss(z-pi)");
  OPTION(options, acoef, "gauss(x)");
  OPTION(options, ccoef, "sin(x) * gauss(x-0.5)");  

  // Create the coefficients
  Field3D input = f.create3D(in);
  Field2D a = f.create2D(acoef);
  Field3D c = f.create3D(ccoef);
  SAVE_ONCE3(input, a, c);
  
  // Create two solvers, using different options
  auto solver1 = Laplacian::create(options->getSection("solver1"));
  auto solver2 = Laplacian::create(options->getSection("solver2"));
  
  solver1->setCoefA(a);
  solver1->setCoefC(c);

  solver2->setCoefA(a);
  solver2->setCoefC(c);

  Field3D result1 = solver1->solve(input);
  Field3D result2 = solver2->solve(input, result1);

  SAVE_ONCE2(result1, result2);
 
  Field3D check1 = a*result1 + Delp2(result1);
  check1.applyBoundary("dirichlet");
  Field3D check2 = a*result2 + Delp2(result2);
  check2.applyBoundary("dirichlet");

  SAVE_ONCE2(check1, check2);
 
  dump.write();
  dump.close();
  
  output << "\nFinished running test. Triggering error to quit\n\n";
  
  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data
  
  return 1;
}

int physics_run(BoutReal UNUSED(t)) {
  // Doesn't do anything
  return 1;
}
