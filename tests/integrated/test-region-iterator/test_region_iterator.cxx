#include <bout.hxx>
#include <boutmain.hxx>

#include <bout/region.hxx>
#include <bout/assert.hxx>

Field3D n;

int physics_init(bool UNUSED(restarting)) {

  Field3D a=1.0, b=1.0, c=2.0;

  using bout::globals::mesh;
  Region<Ind3D> reg(0, mesh->LocalNx - 1, 0, mesh->LocalNy - 1, 0, mesh->LocalNz - 1,
                    mesh->LocalNy, mesh->LocalNz);

  BOUT_FOR(i, reg) {
    a[i] = 3.0;
    b[i] = c[i];
  }

  //Check expected results
  int nerr=0;
  for (const auto &i : a.getRegion(RGN_ALL)) {
    if (a[i] != 3.0) nerr++;
  }
  if(nerr != 0 ) throw BoutException("Unexpected values found in 'a', count {:d}",nerr);
  nerr=0;
  for (const auto &i : b.getRegion(RGN_ALL)) {
    if (b[i] != c[i]) nerr++;
  }
  if(nerr != 0 ) throw BoutException("Unexpected values found in 'b', count {:d}",nerr);


  Field3D d=1.0, e=1.0, f=2.0;
  BOUT_FOR(i, d.getRegion("RGN_NOBNDRY")) {
    d[i] = 3.0;
    e[i] = f[i];
  }

  //Check expected results
  nerr=0;
  //Expect to find differences just in the boundaries, so work out how many boundary points there area
  const int nerrExpected = (2*mesh->xstart*mesh->LocalNy + 2*mesh->ystart*(mesh->LocalNx-mesh->xstart*2))*mesh->LocalNz;
  for (const auto &i : d.getRegion(RGN_ALL)) {
    if (d[i] != 3.0) nerr++;
  }
  if(nerr != nerrExpected ) throw BoutException("Unexpected values found in 'd', count {:d}",nerr);
  nerr=0;
  for (const auto &i : e.getRegion(RGN_ALL)) {
    if (e[i] != f[i]) nerr++;
  }
  if(nerr != nerrExpected ) throw BoutException("Unexpected values found in 'e', count {:d}",nerr);

  SOLVE_FOR(n);
  return 0;
}

int physics_run(BoutReal UNUSED(t)) {
  ddt(n) = 0.;
  return 0;
}
