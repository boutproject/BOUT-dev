#include <boutmain.hxx>
#include <derivs.hxx>
#include <field_factory.hxx>

#include <assert.h>

Field3D n;
Field2D n2;

Field3D xlow3;
Field2D xlow2;


int physics_init(bool restart) {
  xlow3.setLocation(CELL_XLOW);
  //xlow2.setLocation(CELL_XLOW);
  SOLVE_FOR(n);
  SOLVE_FOR(n2);
  //SOLVE_FOR(xlow2);
  SOLVE_FOR(xlow3);
  return 0;
}


int physics_run(BoutReal time) {
  ddt(n)=0;
  ddt(n2)=0;
  //ddt(xlow2)=0;
  ddt(xlow3)=0;
  BoutReal min3 = min(n,true);
  BoutReal max3 = max(n,true);
  BoutReal min2 = min(n2,true);
  BoutReal max2 = max(n2,true);
  //BoutReal minl2 = min(xlow2,true);
  //BoutReal maxl2 = max(xlow2,true);
  BoutReal minl3 = min(xlow3,true);
  BoutReal maxl3 = max(xlow3,true);

  // Field2D and Field3D should be same
  if ( ( abs(min2-min3)   > 1e-8 )
       || ( abs(max2-max3)   > 1e-8 )
       // || ( abs(minl3-minl2) > 1e-8 )
       // || ( abs(maxl3-maxl2) > 1e-8 )
       ) {
    output.write("%g %g\n%g %g\n",min2,min3,max2,max3);
    //output.write("%g %g\n%g %g\n",minl2,minl3,maxl2,maxl3);
    throw BoutException("Something is wrong");
  }
  // And we also now what these values should be
  int nx = mesh->GlobalNx- mesh->xstart*2;
  BoutReal low=1./nx;
  if ( ( abs(minl3)             > 1e-8 ) ||
       ( abs(maxl3 - 2 + 2*low) > 1e-8 ) ||
       ( abs(min3  - low      ) > 1e-8 ) ||
       ( abs(max3  - 2 + low  ) > 1e-8 ) )
    {
    output.write("%g %g\n%g %g\n",min2,min3,max2,max3);
    output.write("%g %g\n",minl3,maxl3);
    //output.write("%g %g\n%g %g\n",minl2,minl3,maxl2,maxl3);
    throw BoutException("Something other is wrong");
  }
  return 0;
}
