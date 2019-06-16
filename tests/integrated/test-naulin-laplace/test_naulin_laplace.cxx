/**************************************************************************
 * Testing Perpendicular Laplacian inversion using LaplaceNaulin solver
 *
 **************************************************************************
 * Copyright 2013 J.T. Omotani
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 *
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <bout.hxx>
#include <bout/constants.hxx>
#include <field_factory.hxx>
#include <boutexception.hxx>
#include <options.hxx>
#include <invert_laplace.hxx>
#include <cmath>
#include <derivs.hxx>
#include "../../../src/invert/laplace/impls/naulin/naulin_laplace.hxx"

BoutReal max_error_at_ystart(const Field3D &error);
Field3D this_Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g);

int main(int argc, char** argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  Options* options = Options::getRoot()->getSection("laplace");
  LaplaceNaulin invert(options);

  // Solving equations of the form d*Grad_perp2(f) + 1/c*Grad_perp(c).Grad_perp(f) + a*f = b for various boundary conditions
  Field3D f1,a1,b1,c1,d1,sol1,bcheck1;
  Field3D absolute_error1;
  BoutReal max_error1; //Output of test

  dump.add(mesh->getCoordinates()->G1,"G1");
  dump.add(mesh->getCoordinates()->G3,"G3");

  ////////////////////////////////////////////////////////////////////////////////
  // Test 1: zero-value Dirichlet boundaries
  f1 = FieldFactory::get()->create3D("f1:function", Options::getRoot(), mesh);
  d1 = FieldFactory::get()->create3D("d1:function", Options::getRoot(), mesh);
  c1 = FieldFactory::get()->create3D("c1:function", Options::getRoot(), mesh);
  a1 = FieldFactory::get()->create3D("a1:function", Options::getRoot(), mesh);

  b1 = d1*Delp2(f1) + this_Grad_perp_dot_Grad_perp(c1,f1)/c1 + a1*f1;
  sol1 = 0.;

  invert.setInnerBoundaryFlags(0);
  invert.setOuterBoundaryFlags(0);
  invert.setCoefA(a1);
  invert.setCoefC(c1);
  invert.setCoefD(d1);

  try {
    sol1 = invert.solve(b1);
    mesh->communicate(sol1);
    checkData(sol1);
    bcheck1 = d1*Delp2(sol1) + this_Grad_perp_dot_Grad_perp(c1,sol1)/c1 + a1*sol1;
    absolute_error1 = f1-sol1;
    max_error1 = max_error_at_ystart(abs(absolute_error1, "RGN_NOBNDRY"));
  } catch (BoutException &err) {
    output << "BoutException occured in invert->solve(b1): " << err.what() << endl
           << "Laplacian inversion failed to converge (probably)" << endl;
    max_error1 = -1;
    sol1 = -1.;
    bcheck1 = -1.;
    absolute_error1 = -1.;
  }

  output<<endl<<"Test 1: zero Dirichlet"<<endl;
  output<<"Magnitude of maximum absolute error is "<<max_error1<<endl;
  output<<"Solver took "<<invert.getMeanIterations()<<" iterations to converge"<<endl;

  dump.add(a1,"a1");
  dump.add(b1,"b1");
  dump.add(c1,"c1");
  dump.add(d1,"d1");
  dump.add(f1,"f1");
  dump.add(sol1,"sol1");
  dump.add(bcheck1,"bcheck1");
  dump.add(absolute_error1,"absolute_error1");
  dump.add(max_error1,"max_error1");

  ////////////////////////////////////////////////////////////////////////////////

  invert.resetMeanIterations();
  Field3D f2,a2,b2,c2,d2,sol2,bcheck2;
  Field3D absolute_error2;
  BoutReal max_error2; //Output of test
  // Test 2: zero-value Neumann boundaries
  f2 = FieldFactory::get()->create3D("f2:function", Options::getRoot(), mesh);
  d2 = FieldFactory::get()->create3D("d2:function", Options::getRoot(), mesh);
  c2 = FieldFactory::get()->create3D("c2:function", Options::getRoot(), mesh);
  a2 = FieldFactory::get()->create3D("a2:function", Options::getRoot(), mesh);

  b2 = d2*Delp2(f2) + this_Grad_perp_dot_Grad_perp(c2,f2)/c2 + a2*f2;
  sol2 = 0.;

  invert.setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  invert.setOuterBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD);
  invert.setCoefA(a2);
  invert.setCoefC(c2);
  invert.setCoefD(d2);

  try {
    sol2 = invert.solve(b2);
    mesh->communicate(sol2);
    bcheck2 = d2*Delp2(sol2) + this_Grad_perp_dot_Grad_perp(c2,sol2)/c2 + a2*sol2;
    absolute_error2 = f2-sol2;
    max_error2 = max_error_at_ystart(abs(absolute_error2, "RGN_NOBNDRY"));
  } catch (BoutException &err) {
    output << "BoutException occured in invert->solve(b2): " << err.what() << endl
           << "Laplacian inversion failed to converge (probably)" << endl;
    max_error2 = -1;
    sol2 = -1.;
    bcheck2 = -1.;
    absolute_error2 = -1.;
  }

  output<<endl<<"Test 2: zero Neumann"<<endl;
  output<<"Magnitude of maximum absolute error is "<<max_error2<<endl;
  output<<"Solver took "<<invert.getMeanIterations()<<" iterations to converge"<<endl;

  dump.add(a2,"a2");
  dump.add(b2,"b2");
  dump.add(c2,"c2");
  dump.add(d2,"d2");
  dump.add(f2,"f2");
  dump.add(sol2,"sol2");
  dump.add(bcheck2,"bcheck2");
  dump.add(absolute_error2,"absolute_error2");
  dump.add(max_error2,"max_error2");

  ////////////////////////////////////////////////////////////////////////////////

  invert.resetMeanIterations();
  Field3D f3,a3,b3,c3,d3,sol3,bcheck3;
  Field3D absolute_error3;
  BoutReal max_error3; //Output of test
  // Test 3: set-value Dirichlet boundaries
  f3 = FieldFactory::get()->create3D("f3:function", Options::getRoot(), mesh);
  d3 = FieldFactory::get()->create3D("d3:function", Options::getRoot(), mesh);
  c3 = FieldFactory::get()->create3D("c3:function", Options::getRoot(), mesh);
  a3 = FieldFactory::get()->create3D("a3:function", Options::getRoot(), mesh);

  b3 = d3*Delp2(f3) + this_Grad_perp_dot_Grad_perp(c3,f3)/c3 + a3*f3;
  sol3 = 0.;

  invert.setInnerBoundaryFlags(INVERT_SET);
  invert.setOuterBoundaryFlags(INVERT_SET);
  invert.setCoefA(a3);
  invert.setCoefC(c3);
  invert.setCoefD(d3);

  // make field to pass in boundary conditions
  Field3D x0 = 0.;
  if (mesh->firstX())
    for (int k=0;k<mesh->LocalNz;k++)
      x0(mesh->xstart-1,mesh->ystart,k) = 0.5*(f3(mesh->xstart-1,mesh->ystart,k)+f3(mesh->xstart,mesh->ystart,k));
  if (mesh->lastX())
    for (int k=0;k<mesh->LocalNz;k++)
      x0(mesh->xend+1,mesh->ystart,k) = 0.5*(f3(mesh->xend+1,mesh->ystart,k)+f3(mesh->xend,mesh->ystart,k));

  try {
    sol3 = invert.solve(b3, x0);
    mesh->communicate(sol3);
    bcheck3 = d3*Delp2(sol3) + this_Grad_perp_dot_Grad_perp(c3,f3)/c3 + a3*sol3;
    absolute_error3 = f3-sol3;
    max_error3 = max_error_at_ystart(abs(absolute_error3, "RGN_NOBNDRY"));
  } catch (BoutException &err) {
    output << "BoutException occured in invert->solve(b3): " << err.what() << endl
           << "Laplacian inversion failed to converge (probably)" << endl;
    max_error3 = -1;
    sol3 = -1.;
    bcheck3 = -1.;
    absolute_error3 = -1.;
  }

  output<<endl<<"Test 3: set Dirichlet"<<endl;
  output<<"Magnitude of maximum absolute error is "<<max_error3<<endl;
  output<<"Solver took "<<invert.getMeanIterations()<<" iterations to converge"<<endl;

  dump.add(a3,"a3");
  dump.add(b3,"b3");
  dump.add(c3,"c3");
  dump.add(d3,"d3");
  dump.add(f3,"f3");
  dump.add(sol3,"sol3");
  dump.add(bcheck3,"bcheck3");
  dump.add(absolute_error3,"absolute_error3");
  dump.add(max_error3,"max_error3");

  ////////////////////////////////////////////////////////////////////////////////

  invert.resetMeanIterations();
  Field3D f4,a4,b4,c4,d4,sol4,bcheck4;
  Field3D absolute_error4;
  BoutReal max_error4; //Output of test
  // Test 4: set-value Neumann boundaries
  f4 = FieldFactory::get()->create3D("f4:function", Options::getRoot(), mesh);
  d4 = FieldFactory::get()->create3D("d4:function", Options::getRoot(), mesh);
  c4 = FieldFactory::get()->create3D("c4:function", Options::getRoot(), mesh);
  a4 = FieldFactory::get()->create3D("a4:function", Options::getRoot(), mesh);

  b4 = d4*Delp2(f4) + this_Grad_perp_dot_Grad_perp(c4,f4)/c4 + a4*f4;
  sol4 = 0.;

  invert.setInnerBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD + INVERT_SET);
  invert.setOuterBoundaryFlags(INVERT_DC_GRAD + INVERT_AC_GRAD + INVERT_SET);
  invert.setCoefA(a4);
  invert.setCoefC(c4);
  invert.setCoefD(d4);

  // make field to pass in boundary conditions
  x0 = 0.;
  if (mesh->firstX())
    for (int k=0;k<mesh->LocalNz;k++)
      x0(mesh->xstart-1,mesh->ystart,k) = (f4(mesh->xstart,mesh->ystart,k)-f4(mesh->xstart-1,mesh->ystart,k))
                                        /mesh->getCoordinates()->dx(mesh->xstart,mesh->ystart)
                                        /sqrt(mesh->getCoordinates()->g_11(mesh->xstart,mesh->ystart));
  if (mesh->lastX())
    for (int k=0;k<mesh->LocalNz;k++)
      x0(mesh->xend+1,mesh->ystart,k) = (f4(mesh->xend+1,mesh->ystart,k)-f4(mesh->xend,mesh->ystart,k))
                                        /mesh->getCoordinates()->dx(mesh->xend,mesh->ystart)
                                        /sqrt(mesh->getCoordinates()->g_11(mesh->xend,mesh->ystart));

  try {
    sol4 = invert.solve(b4, x0);
    mesh->communicate(sol4);
    bcheck4 = d4*Delp2(sol4) + this_Grad_perp_dot_Grad_perp(c4,sol4)/c4 + a4*sol4;
    absolute_error4 = f4-sol4;
    max_error4 = max_error_at_ystart(abs(absolute_error4, "RGN_NOBNDRY"));
  } catch (BoutException &err) {
    output << "BoutException occured in invert->solve(b4): " << err.what() << endl
           << "Laplacian inversion failed to converge (probably)" << endl;
    max_error4 = -1;
    sol4 = -1.;
    bcheck4 = -1.;
    absolute_error4 = -1.;
  }

  output<<endl<<"Test 4: set Neumann"<<endl;
  output<<"Magnitude of maximum absolute error is "<<max_error4<<endl;
  output<<"Solver took "<<invert.getMeanIterations()<<" iterations to converge"<<endl;

  dump.add(a4,"a4");
  dump.add(b4,"b4");
  dump.add(c4,"c4");
  dump.add(d4,"d4");
  dump.add(f4,"f4");
  dump.add(sol4,"sol4");
  dump.add(bcheck4,"bcheck4");
  dump.add(absolute_error4,"absolute_error4");
  dump.add(max_error4,"max_error4");

  ////////////////////////////////////////////////////////////////////////////////

  output << "\nFinished running test.\n\n";

  dump.write();

  MPI_Barrier(BoutComm::get()); // Wait for all processors to write data

  BoutFinalise();
  return 0;

}

Field3D this_Grad_perp_dot_Grad_perp(const Field3D &f, const Field3D &g) {
  Field3D result = mesh->getCoordinates()->g11 * ::DDX(f) * ::DDX(g) + mesh->getCoordinates()->g33 * ::DDZ(f) * ::DDZ(g)
                   + mesh->getCoordinates()->g13 * (DDX(f)*DDZ(g) + DDZ(f)*DDX(g));
  
  return result;
}

BoutReal max_error_at_ystart(const Field3D &error) {

  BoutReal local_max_error = error(mesh->xstart, mesh->ystart, 0);

  for (int jx=mesh->xstart; jx<=mesh->xend; jx++)
    for (int jz=0; jz<mesh->LocalNz; jz++)
      if (local_max_error<error(jx, mesh->ystart, jz)) local_max_error=error(jx, mesh->ystart, jz);
  
  BoutReal max_error;

  MPI_Allreduce(&local_max_error, &max_error, 1, MPI_DOUBLE, MPI_MAX, BoutComm::get());

  return max_error;
}
