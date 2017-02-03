#include <bout.hxx>

#include <bout/assert.hxx>

int main(int argc, char** argv) {
  BoutInitialise(argc, argv);
  
  Field2D f2d = 2;
  
  // Test 2D field
  f2d += 2;
  ASSERT0(fabs(f2d(2,2) - 4) < 1e-10);
  f2d -= 1;
  ASSERT0((f2d(2,1) - 3) < 1e-10);
  f2d *= f2d;
  ASSERT0((f2d(1,2) - 9) < 1e-10);
  f2d /= 2;
  ASSERT0((f2d(1,1) - 4.5) < 1e-10);
  f2d += f2d*3;
  ASSERT0(fabs(f2d(1,1) - 18) < 1e-10);
  
  // Test 3D field
  
  Field3D f3d = 3;
  f2d = 2;
  
  f3d += 2;
  ASSERT0(fabs(f3d(2,2,0) - 5) < 1e-10);
  f3d *= 1.5;
  ASSERT0(fabs(f3d(2,2,0) - 7.5) < 1e-10);
  f3d -= 3;
  ASSERT0(fabs(f3d(2,2,0) - 4.5) < 1e-10);
  f3d /= 2;
  ASSERT0(fabs(f3d(2,2,0) - 2.25) < 1e-10);
  f3d += f3d;
  ASSERT0(fabs(f3d(2,2,0) - 4.5) < 1e-10);
  f3d /= f3d;
  ASSERT0(fabs(f3d(2,2,0) - 1.0) < 1e-10);
  f3d -= f2d;
  ASSERT0(fabs(f3d(2,2,0) + 1.0) < 1e-10);
  f3d *= f2d;
  ASSERT0(fabs(f3d(2,2,0) + 2.0) < 1e-10);
  f3d = pow(f3d, f2d);
  ASSERT0(fabs(f3d(2,2,0) - 4.0) < 1e-10);
 
  // Some iterator tests
 
  //for(auto i : IndexRange{0,10,0,1,0,1})
  //  output.write("%d,%d,%d\n", i.x, i.y, i.z);

  Field3D result;
  result.allocate();
  
  for(auto i : result.region(RGN_NOY))
    result[i] = f3d[i.yp()] - 2*f3d[i] + f3d[i.ym()];
  
  BoutFinalise();
  return 0;
}
