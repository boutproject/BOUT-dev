#include <bout.hxx>
#include <boutmain.hxx>

#include <bout/region.hxx>
#include <bout/assert.hxx>

Field3D n;

int physics_init(bool restarting) {

  Field3D a=1.0, b=1.0, c=2.0;

  Region<ind3D> reg (0,mesh->LocalNx-1,
		    0,mesh->LocalNy-1,
		    0,mesh->LocalNz-1,
		    mesh->LocalNy,
		    mesh->LocalNz);

  BLOCK_REGION_LOOP_SERIAL(reg,i,
			   a[i] = 3.0;
			   b[i] = c[i];
			   );

  //Check expected results
  int nerr=0;
  for (const auto &i: a.region(RGN_ALL)){
    if (a[i] != 3.0) nerr++;
  }
  if(nerr != 0 ) throw BoutException("Unexpected values found in 'a', count %d",nerr);
  for (const auto &i: b.region(RGN_ALL)){
    if (b[i] != c[i]) nerr++;
  }
  if(nerr != 0 ) throw BoutException("Unexpected values found in 'b', count %d",nerr);

  SOLVE_FOR(n);
  return 0;
}

int physics_run(BoutReal t) {
  ddt(n) = 0.;
  return 0;
}
