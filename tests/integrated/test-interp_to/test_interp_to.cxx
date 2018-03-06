/*
 * FieldFactory regression test
 * 
 * Test the FieldFactory class
 *
 */

#include <bout.hxx>
#include <field_factory.hxx>
#include <interpolation.hxx>

int main(int argc, char **argv) {

  // Initialise BOUT++, setting up mesh
  BoutInitialise(argc, argv);

  FieldFactory f(mesh);
  
  Field2D f2d = f.create2D("sin(3*y)*cos(2*pi*x)");
  Field2D f2d_xlow = interp_to(f2d, CELL_XLOW);
  Field2D f2d_xlow_nobndry = interp_to(f2d, CELL_XLOW, RGN_NOBNDRY);
  Field2D f2d_ylow = interp_to(f2d, CELL_YLOW);
  Field2D f2d_ylow_nobndry = interp_to(f2d, CELL_YLOW, RGN_NOBNDRY);
  SAVE_ONCE5(f2d, f2d_xlow, f2d_xlow_nobndry, f2d_ylow, f2d_ylow_nobndry);

  Field3D f3d = f.create3D("gauss(x-0.5,0.2)*gauss(y-pi)*sin(z)");
  SAVE_ONCE(f3d);
  Field3D f3d_xlow = interp_to(f3d, CELL_XLOW);
  Field3D f3d_xlow_nobndry = interp_to(f3d, CELL_XLOW, RGN_NOBNDRY);
  Field3D f3d_ylow = interp_to(f3d, CELL_YLOW);
  Field3D f3d_ylow_nobndry = interp_to(f3d, CELL_YLOW, RGN_NOBNDRY);
  Field3D f3d_zlow = interp_to(f3d, CELL_ZLOW);
  Field3D f3d_zlow_nobndry = interp_to(f3d, CELL_ZLOW, RGN_NOBNDRY);
  SAVE_ONCE6(f3d_xlow, f3d_xlow_nobndry, f3d_ylow, f3d_ylow_nobndry, f3d_zlow, f3d_zlow_nobndry);

  // Write data to file
  dump.write();
  dump.close();
  
  // Need to wait for all processes to finish writing
  MPI_Barrier(BoutComm::get());

  /// Finished, tidy up and free memory
  BoutFinalise();

  return 0;
}
